"""
Sentinel Service — Job Supervisor
Responsibilities:
  - Poll all job queues (docking, md, qsar, pharmacophore)
  - Detect failures, stale jobs, timeouts
  - Retry with exponential backoff
  - Validate outputs before passing forward
  - Fallback strategies (e.g., Vina failed → switch to GNINA)
  - Escalate to notifications on unrecoverable failures
"""

import os
import logging
import asyncio
import httpx
import redis
import json
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List
from pathlib import Path

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from contextlib import asynccontextmanager

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("sentinel")

REDIS_URL = os.getenv("REDIS_URL", "redis://redis:6379")
DOCKING_SERVICE_URL = os.getenv("DOCKING_SERVICE_URL", "http://docking-service:8002")
MD_SERVICE_URL = os.getenv("MD_SERVICE_URL", "http://md-service:8000")
QSAR_SERVICE_URL = os.getenv("QSAR_SERVICE_URL", "http://qsar-service:8005")
PHARMACOPHORE_SERVICE_URL = os.getenv(
    "PHARMACOPHORE_SERVICE_URL", "http://pharmacophore-service:8004"
)
API_BACKEND_URL = os.getenv("API_BACKEND_URL", "http://api-backend:8000")
NOTIFY_URL = os.getenv("NOTIFY_URL", "http://api-backend:8000/notify")

POLL_INTERVAL = int(os.getenv("SENTINEL_POLL_INTERVAL", "30"))
JOB_TIMEOUT_SECONDS = int(os.getenv("SENTINEL_JOB_TIMEOUT", "3600"))
MAX_RETRIES = int(os.getenv("SENTINEL_MAX_RETRIES", "3"))

try:
    redis_client = redis.from_url(REDIS_URL, decode_responses=True)
    redis_client.ping()
    redis_available = True
except Exception:
    redis_client = None
    redis_available = False

NOTIFICATION_MANAGER_URL = os.getenv(
    "NOTIFICATION_MANAGER_URL", "http://api-backend:8000/notify"
)


def _get_redis():
    return redis_client


class RetryStrategy(BaseModel):
    job_id: str
    service: str
    original_params: Dict[str, Any]
    retry_count: int = 0
    last_error: Optional[str] = None
    strategy: str = "retry"


FALLBACK_STRATEGIES: Dict[str, Dict[str, Any]] = {
    "docking:vina_failed": {
        "strategy": "switch_engine",
        "engine": "gnina",
        "description": "AutoDock Vina failed — switching to GNINA",
    },
    "docking:timeout": {
        "strategy": "reduce_complexity",
        "description": "Docking timed out — reducing exhaustiveness",
        "params": {"exhaustiveness": 4},
    },
    "md:unstable": {
        "strategy": "reduce_sim_time",
        "description": "MD simulation unstable — reducing simulation time",
        "params": {"steps_multiply": 0.5},
    },
    "md:timeout": {
        "strategy": "shorten_sim",
        "description": "MD timed out — shortening simulation",
        "params": {"steps_multiply": 0.25},
    },
    "qsar:training_failed": {
        "strategy": "simpler_model",
        "description": "QSAR training failed — switching to simpler model",
        "params": {"model_type": "Ridge"},
    },
    # max_retries_exceeded keys — match the lookup format used in monitor_job()
    "docking:max_retries_exceeded": {
        "strategy": "switch_engine",
        "engine": "gnina",
        "description": "Docking failed after max retries — switching to GNINA",
    },
    "md:max_retries_exceeded": {
        "strategy": "reduce_sim_time",
        "description": "MD failed after max retries — reducing simulation time",
        "params": {"steps_multiply": 0.25},
    },
    "qsar:max_retries_exceeded": {
        "strategy": "simpler_model",
        "description": "QSAR failed after max retries — switching to Ridge regression",
        "params": {"model_type": "Ridge"},
    },
}


app = FastAPI(
    title="Sentinel Service",
    version="2.0.0",
    description="Job supervisor: monitor, retry, validate",
)


@app.get("/health")
def health():
    return {
        "service": "sentinel",
        "status": "healthy",
        "redis_available": redis_available,
        "poll_interval": POLL_INTERVAL,
        "job_timeout_seconds": JOB_TIMEOUT_SECONDS,
        "timestamp": datetime.now().isoformat(),
    }


@app.get("/")
def root():
    return {"service": "Sentinel (Supervisor)", "version": "2.0.0"}


@app.post("/monitor/job")
def monitor_job(job_id: str, service: str):
    """
    Check a specific job status, detect failures, and decide action.
    Returns: { status, action, details }
    """
    r = _get_redis()
    if not r:
        raise HTTPException(status_code=503, detail="Redis unavailable")

    job_key = f"{service}_job:{job_id}"
    data = r.get(job_key)
    if not data:
        raise HTTPException(
            status_code=404, detail=f"Job {job_id} not found for service {service}"
        )

    job = json.loads(data)
    status = job.get("status", "unknown")
    error = job.get("error")
    retry_count = job.get("retry_count", 0)
    created_at = job.get("created_at", datetime.now().isoformat())

    age_seconds = (datetime.now() - datetime.fromisoformat(created_at)).total_seconds()
    is_stale = age_seconds > JOB_TIMEOUT_SECONDS and status in ("pending", "running")

    if status == "failed":
        if retry_count < MAX_RETRIES:
            return {
                "status": "failed",
                "action": "retry",
                "retry_count": retry_count,
                "message": f"Job failed (attempt {retry_count + 1}/{MAX_RETRIES}): {error}",
            }
        else:
            fallback_key = f"{service}:max_retries_exceeded"
            fallback = FALLBACK_STRATEGIES.get(fallback_key)
            if fallback:
                return {
                    "status": "failed",
                    "action": "fallback",
                    "strategy": fallback,
                    "message": f"Max retries exceeded. Applying fallback: {fallback.get('description')}",
                }
            return {
                "status": "failed",
                "action": "escalate",
                "message": f"Job failed after {MAX_RETRIES} retries: {error}",
            }

    if is_stale:
        if retry_count < MAX_RETRIES:
            return {
                "status": "stale",
                "action": "retry",
                "retry_count": retry_count,
                "message": f"Job stale for {int(age_seconds)}s — requeueing (attempt {retry_count + 1}/{MAX_RETRIES})",
            }
        return {
            "status": "stale",
            "action": "escalate",
            "message": f"Job stale for {int(age_seconds)}s and max retries exceeded",
        }

    if status == "completed":
        result = job.get("result")
        validation = _validate_result(service, result)
        if validation["valid"]:
            return {
                "status": "completed",
                "action": "pass",
                "message": "Job completed and validated",
                "result_summary": _summarize_result(service, result),
            }
        return {
            "status": "completed",
            "action": "revalidate",
            "message": f"Result validation warning: {validation['reason']}",
        }

    return {
        "status": status,
        "action": "monitor",
        "message": f"Job {status}, running for {int(age_seconds)}s",
    }


@app.post("/retry/job")
def retry_job(job_id: str, service: str):
    """Re-queue a failed or stale job for retry with exponential backoff."""
    r = _get_redis()
    if not r:
        raise HTTPException(status_code=503, detail="Redis unavailable")

    job_key = f"{service}_job:{job_id}"
    data = r.get(job_key)
    if not data:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    job = json.loads(data)
    retry_count = job.get("retry_count", 0) + 1
    backoff_seconds = min(2**retry_count * 10, 300)

    job["retry_count"] = retry_count
    job["last_retry_at"] = datetime.now().isoformat()
    job["status"] = "pending"
    job["error"] = None

    r.setex(job_key, 7200, json.dumps(job))

    r.sadd(f"{service}_queue", json.dumps({"job_id": job_id, "retry": retry_count}))

    return {
        "success": True,
        "job_id": job_id,
        "retry_count": retry_count,
        "backoff_seconds": backoff_seconds,
        "message": f"Job re-queued for retry {retry_count}/{MAX_RETRIES} after {backoff_seconds}s backoff",
    }


@app.post("/fallback")
def apply_fallback(job_id: str, service: str, strategy: str):
    """Apply a fallback strategy to a failed job."""
    fallback_key = f"{service}:{strategy}"
    fallback = FALLBACK_STRATEGIES.get(fallback_key)

    if not fallback:
        raise HTTPException(
            status_code=400, detail=f"No fallback strategy for {fallback_key}"
        )

    strategy_type = fallback.get("strategy")

    if strategy_type == "switch_engine":
        return _switch_engine(job_id, service, fallback)
    elif strategy_type == "reduce_complexity":
        return _reduce_complexity(job_id, service, fallback)
    elif strategy_type == "reduce_sim_time":
        return _reduce_sim_time(job_id, service, fallback)
    elif strategy_type == "simpler_model":
        return _simpler_model(job_id, service, fallback)
    else:
        return {"success": False, "message": f"Unknown strategy type: {strategy_type}"}


def _switch_engine(job_id: str, service: str, fallback: Dict) -> Dict:
    if service == "docking":
        new_engine = fallback.get("engine", "gnina")
        return {
            "success": True,
            "action": "switch_engine",
            "new_engine": new_engine,
            "job_id": job_id,
            "message": fallback.get("description"),
            "requeue": True,
        }
    return {"success": False, "message": f"Switch engine not supported for {service}"}


def _reduce_complexity(job_id: str, service: str, fallback: Dict) -> Dict:
    params = fallback.get("params", {})
    return {
        "success": True,
        "action": "reduce_complexity",
        "adjusted_params": params,
        "job_id": job_id,
        "message": fallback.get("description"),
        "requeue": True,
    }


def _reduce_sim_time(job_id: str, service: str, fallback: Dict) -> Dict:
    params = fallback.get("params", {})
    return {
        "success": True,
        "action": "reduce_sim_time",
        "adjusted_params": params,
        "job_id": job_id,
        "message": fallback.get("description"),
        "requeue": True,
    }


def _simpler_model(job_id: str, service: str, fallback: Dict) -> Dict:
    params = fallback.get("params", {})
    return {
        "success": True,
        "action": "simpler_model",
        "adjusted_params": params,
        "job_id": job_id,
        "message": fallback.get("description"),
        "requeue": True,
    }


@app.post("/escalate")
def escalate(job_id: str, service: str, reason: str):
    """Escalate a failed job to notifications and log."""
    r = _get_redis()
    job_data = None
    if r:
        job_key = f"{service}_job:{job_id}"
        data = r.get(job_key)
        if data:
            job_data = json.loads(data)

    title = f"Sentinel Escalation: {service} job failed"
    message = f"Job `{job_id}` ({service}) has failed and cannot be recovered.\nReason: {reason}"
    details = {
        "job_id": job_id,
        "service": service,
        "reason": reason,
        "job_data": job_data,
        "timestamp": datetime.now().isoformat(),
    }

    try:
        with httpx.Client(timeout=10.0) as client:
            client.post(
                f"{API_BACKEND_URL}/notify",
                json={
                    "event": "error",
                    "title": title,
                    "message": message,
                    "details": details,
                },
            )
        notified = True
    except Exception:
        notified = False

    return {
        "success": True,
        "job_id": job_id,
        "service": service,
        "reason": reason,
        "notified": notified,
        "timestamp": datetime.now().isoformat(),
    }


@app.get("/queue/status")
def queue_status():
    """Get status of all job queues."""
    r = _get_redis()
    if not r:
        raise HTTPException(status_code=503, detail="Redis unavailable")

    services = ["docking", "md", "qsar", "pharmacophore", "rdkit"]
    status = {}
    total_pending = 0
    total_running = 0
    total_failed = 0

    for svc in services:
        queue_key = f"{svc}_queue"
        pending = r.scard(queue_key)

        job_keys = list(r.scan_iter(f"{svc}_job:*", count=100))
        running = 0
        failed = 0
        completed = 0
        stale = 0

        for key in job_keys:
            data = r.get(key)
            if data:
                job = json.loads(data)
                s = job.get("status", "unknown")
                created = job.get("created_at", datetime.now().isoformat())
                age = (datetime.now() - datetime.fromisoformat(created)).total_seconds()
                if s == "running":
                    running += 1
                elif s == "failed":
                    failed += 1
                elif s == "completed":
                    completed += 1
                if age > JOB_TIMEOUT_SECONDS and s in ("pending", "running"):
                    stale += 1

        total_pending += pending
        total_running += running
        total_failed += failed

        status[svc] = {
            "pending": pending,
            "running": running,
            "completed": completed,
            "failed": failed,
            "stale": stale,
        }

    return {
        "services": status,
        "totals": {
            "pending": total_pending,
            "running": total_running,
            "failed": total_failed,
        },
        "timestamp": datetime.now().isoformat(),
    }


@app.post("/validate/result")
def validate_result(service: str, result: Dict[str, Any]):
    """Validate a job result before passing it forward."""
    validation = _validate_result(service, result)
    return validation


def _validate_result(service: str, result: Optional[Dict]) -> Dict[str, Any]:
    """Internal: check if a result is valid for the given service type."""
    if not result:
        return {"valid": False, "reason": "No result data"}

    if service == "docking":
        if not isinstance(result, dict):
            return {"valid": False, "reason": "Result is not a dictionary"}
        if "vina_score" not in result and "binding_affinity" not in result:
            return {"valid": False, "reason": "No binding affinity in docking result"}
        if result.get("vina_score", 0) > 0:
            return {"valid": False, "reason": "Docking score out of expected range"}
        return {"valid": True}

    elif service == "md":
        if not isinstance(result, dict):
            return {"valid": False, "reason": "Result is not a dictionary"}
        required = ["sim_time_ns", "n_frames", "avg_energy_kj_mol"]
        missing = [k for k in required if k not in result]
        if missing:
            return {"valid": False, "reason": f"Missing MD fields: {missing}"}
        if result.get("n_frames", 0) == 0:
            return {"valid": False, "reason": "No frames generated"}
        return {"valid": True}

    elif service == "qsar":
        if not isinstance(result, dict):
            return {"valid": False, "reason": "Result is not a dictionary"}
        if "predicted_activity" not in result and "cv_r2" not in result:
            return {"valid": False, "reason": "No prediction or metrics in QSAR result"}
        return {"valid": True}

    elif service == "pharmacophore":
        if not isinstance(result, dict):
            return {"valid": False, "reason": "Result is not a dictionary"}
        if "features" not in result and not result.get("smiles"):
            return {"valid": False, "reason": "No pharmacophore features generated"}
        return {"valid": True}

    return {"valid": True}


def _summarize_result(service: str, result: Optional[Dict]) -> Dict[str, Any]:
    """Extract a summary of a result for logging/notifications."""
    if not result:
        return {}
    if service == "docking":
        return {
            "vina_score": result.get("vina_score"),
            "num_poses": result.get("num_poses", 1),
        }
    elif service == "md":
        return {
            "sim_time_ns": result.get("sim_time_ns"),
            "n_frames": result.get("n_frames"),
            "avg_energy_kj_mol": result.get("avg_energy_kj_mol"),
        }
    elif service == "qsar":
        return {
            "cv_r2": result.get("cv_r2"),
            "n_compounds": result.get("n_compounds"),
        }
    return {}


@app.post("/notify")
def send_notification(
    event: str, title: str, message: str, details: Optional[Dict] = None
):
    """Send a notification through the notification system."""
    try:
        with httpx.Client(timeout=10.0) as client:
            resp = client.post(
                f"{MD_SERVICE_URL}/notify",
                json={
                    "event": event,
                    "title": title,
                    "message": message,
                    "details": details or {},
                },
            )
            resp.raise_for_status()
            return {"success": True}
    except Exception as e:
        logger.warning(f"Notification failed: {e}")
        return {"success": False, "error": str(e)}


async def _poll_jobs():
    """Background task to poll Redis for stale/failed jobs"""
    logger.info("Sentinel watchdog started - polling for stale jobs")

    while True:
        try:
            if redis_client and redis_available:
                r = _get_redis()
                if r:
                    now = datetime.utcnow()
                    job_patterns = [
                        "job:*",
                        "docking_job:*",
                        "md_job:*",
                        "qsar_job:*",
                        "pharmacophore_job:*",
                    ]

                    for pattern in job_patterns:
                        keys = r.keys(pattern)
                        for key in keys:
                            job_data = r.get(key)
                            if not job_data:
                                continue
                            try:
                                job = json.loads(job_data)
                                status = job.get("status", "")
                                updated = job.get("updated_at")

                                if updated and status in ["running", "pending"]:
                                    updated_time = datetime.fromisoformat(
                                        updated.replace("Z", "+00:00")
                                    )
                                    age = (
                                        now - updated_time.replace(tzinfo=None)
                                    ).total_seconds()

                                    if age > JOB_TIMEOUT_SECONDS:
                                        logger.warning(
                                            f"Job {key} timed out (age: {age}s)"
                                        )
                                        job_id = key.split(":")[-1]
                                        await _escalate_job(
                                            job_id, job, f"timeout_after_{int(age)}s"
                                        )

                            except (json.JSONDecodeError, KeyError, ValueError) as e:
                                logger.debug(f"Skipping malformed job {key}: {e}")

        except Exception as e:
            logger.error(f"Polling error: {e}")

        await asyncio.sleep(POLL_INTERVAL)


async def _escalate_job(job_id: str, job: dict, reason: str):
    """Escalate a failed job - mark as failed and notify"""
    try:
        job["status"] = "failed"
        job["failure_reason"] = reason
        job["escalated_at"] = datetime.utcnow().isoformat()

        r = _get_redis()
        if r:
            r.set(f"job:{job_id}", json.dumps(job), ex=86400)

        logger.warning(f"Job {job_id} escalated: {reason}")

        notify_data = {
            "job_id": job_id,
            "status": "failed",
            "reason": reason,
            "service": job.get("service", "unknown"),
        }

        async with httpx.AsyncClient(timeout=10.0) as client:
            try:
                await client.post(NOTIFICATION_MANAGER_URL, json=notify_data)
            except Exception:
                pass

    except Exception as e:
        logger.error(f"Failed to escalate job {job_id}: {e}")


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("=" * 60)
    logger.info("Sentinel Service starting up...")
    logger.info(f"Redis available: {redis_available}")
    logger.info(f"Poll interval: {POLL_INTERVAL}s")
    logger.info(f"Job timeout: {JOB_TIMEOUT_SECONDS}s")
    logger.info(f"Max retries: {MAX_RETRIES}")
    logger.info("Role: Job Supervisor (monitor, retry, validate, escalate)")

    task = asyncio.create_task(_poll_jobs())
    logger.info("Watchdog polling task started")

    yield

    task.cancel()
    try:
        await task
    except asyncio.CancelledError:
        pass
    logger.info("Sentinel Service shutting down...")


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8007)
