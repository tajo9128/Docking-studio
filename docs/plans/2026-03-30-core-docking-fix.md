# Core Docking Pipeline Fix - Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the molecular docking pipeline so receptor/ligand preparation → docking → results works without errors.

**Architecture:** Minimal changes to existing infrastructure - only fix broken endpoints, add missing functionality, and add SMILES support. No architectural changes.

**Tech Stack:** FastAPI, RDKit, Meeko, AutoDock Vina, Docker Compose, React/TypeScript

---

## Task 1: Fix API Backend Missing Imports

**Files:**
- Modify: `services/api-backend/main.py`

**Step 1:** Find line with `from fastapi import FastAPI, HTTPException`
Replace with: `from fastapi import FastAPI, HTTPException, Request`

**Step 2:** Add after imports: `import functools`

**Step 3:** Commit
```bash
git add services/api-backend/main.py
git commit -m "fix: add missing functools and Request imports"
```

---

## Task 2: Add Cancel Endpoint to Docking Service

**Files:**
- Modify: `services/docking-service/main.py`

**Step 1:** Add cancel endpoint after existing `/dock/{job_id}` endpoint:
```python
@app.post("/dock/{job_id}/cancel")
async def cancel_docking(job_id: str):
    """Cancel a running docking job"""
    r = redis.from_url(REDIS_URL)
    job_key = f"docking_job:{job_id}"
    
    if not r.exists(job_key):
        raise HTTPException(status_code=404, detail="Job not found")
    
    job_data = r.hgetall(job_key)
    current_status = job_data.get("status", "unknown")
    
    if current_status in ["completed", "failed", "cancelled"]:
        return {"job_id": job_id, "status": current_status}
    
    r.hset(job_key, mapping={"status": "cancelled"})
    r.lpush("docking_cancelled", job_id)
    
    return {"job_id": job_id, "status": "cancelled"}
```

**Step 2:** Update process_docking_job to check cancellation periodically

**Step 3:** Commit
```bash
git add services/docking-service/main.py
git commit -m "feat: add dock cancel endpoint"
```

---

## Task 3: Fix Consensus Engine Error Handling

**Files:**
- Modify: `services/docking-service/docking_engine.py`

**Step 1:** Find consensus function and add graceful GNINA fallback

**Step 2:** Commit
```bash
git add services/docking-service/docking_engine.py
git commit -m "fix: handle GNINA not installed in consensus mode"
```

---

## Task 4: Add Ligand Preparation to Frontend

**Files:**
- Modify: `frontend/src/pages/Docking.tsx`
- Modify: `frontend/src/api/rdkit.ts`

**Step 1:** Add prepareLigand to rdkit.ts API

**Step 2:** Add ligand preparation state and button to Docking.tsx

**Step 3:** Update startDocking to use prepared paths

**Step 4:** Commit
```bash
git add frontend/src/pages/Docking.tsx frontend/src/api/rdkit.ts
git commit -m "feat: add ligand preparation to docking workflow"
```

---

## Task 5: Add SMILES Support for Ligands

**Files:**
- Modify: `frontend/src/pages/Docking.tsx`
- Modify: `frontend/src/api/rdkit.ts`

**Step 1:** Add smilesToSDF to rdkit.ts

**Step 2:** Auto-convert .smi files before docking

**Step 3:** Commit
```bash
git add frontend/src/pages/Docking.tsx frontend/src/api/rdkit.ts
git commit -m "feat: support SMILES ligand input"
```

---

## Task 6: Test the Complete Flow

**Step 1:** Rebuild Docker services
```bash
docker compose build gateway api-backend docking-service rdkit-service
docker compose up -d
```

**Step 2:** Test receptor preparation, ligand preparation, docking, cancel

---

## Task 7: Final Cleanup

**Step 1:** Update CHANGELOG.md

**Step 2:** Tag and push
```bash
git tag v2.0.1
git push origin main --tags
```
