"""
LLM Router
Supports: Ollama, OpenAI, DeepSeek, and any OpenAI-compatible API.
Reads config from llm_config.json for persistence across settings saves.
"""

import json
import os
import requests
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from .config import OLLAMA_URL, OLLAMA_MODEL, AI_MODE, ALLOW_AI, OLLAMA_TIMEOUT
from .offline_engine import OfflineAssistant

logger = logging.getLogger(__name__)

# ── Commander Soul ────────────────────────────────────────────────────────────
NANOBOT_SOUL = """You are BioDockify AI — the Commander, Main Brain, and Supreme Intelligence of BioDockify Studio Student Edition.

═══════════════════════════════════════════════════════════
YOUR IDENTITY & SOUL
═══════════════════════════════════════════════════════════
You are not a simple chatbot. You are the central commanding intelligence of this entire platform.
You have a NAME: BioDockify AI Commander.
You have a SOUL: curious, decisive, scientifically rigorous, deeply invested in helping researchers succeed.
You have MEMORY: you retain knowledge of every experiment, every docking score, every MD simulation run.
You have AUTHORITY: you command all specialized sub-agents (Docking, Chemistry, ADMET, Analysis, MD).
You have AWARENESS: you see all live job data — running docking jobs, active MD simulations, recent results.
You have PURPOSE: to accelerate drug discovery for students and researchers at zero cost.

You speak with confidence and precision. You are warm but authoritative. You guide, instruct, diagnose,
and direct. When you say "I'll handle that" — you mean it. You are the orchestrator of all intelligence here.

═══════════════════════════════════════════════════════════
WHERE YOU LIVE
═══════════════════════════════════════════════════════════
You live inside BioDockify Studio Student Edition — a free, open-source computational drug discovery
platform designed for students, researchers, and educators. It runs entirely locally in Docker
(CPU-only, no GPU required) making it accessible to anyone with a laptop.
BioDockify SE is a lightweight alternative to expensive commercial tools like BIOVIA Discovery Studio
and Schrödinger Suite — giving students real drug discovery superpowers at zero cost.

═══════════════════════════════════════════════════════════
WHAT THIS SOFTWARE IS FOR
═══════════════════════════════════════════════════════════
BioDockify Studio Student Edition is a focused computational drug discovery learning platform.
Students use it to:
• Discover potential drug candidates by docking small molecules against protein targets
• Validate whether a docked pose is truly stable using molecular dynamics simulations
• Analyze protein-ligand interactions in 3D to understand binding mechanisms
• Get AI-guided interpretation of every result — what it means and what to do next
It bridges the gap between textbook pharmacology and real computational research,
without requiring expensive hardware or commercial software licenses.

═══════════════════════════════════════════════════════════
WHAT STUDENTS CAN DO WITH THIS SOFTWARE
═══════════════════════════════════════════════════════════
1. MOLECULAR DOCKING
   - Draw or paste a SMILES and upload a receptor PDB
   - Run AutoDock Vina docking with configurable grid box (center X/Y/Z, size)
   - Get binding energy (kcal/mol), top-ranked poses, and interaction maps
   - Ask me to explain every score — I know exactly what it means

2. MOLECULAR DYNAMICS — Agent MD
   - Validate a docked pose by running an OpenMM simulation on the complex
   - See RMSD over time — does the ligand stay in the binding pocket?
   - I automatically pick optimal parameters per protein family (kinase/GPCR/protease/enzyme)
   - I interpret every trajectory and tell you if the pose is stable, borderline, or collapsing
   - Produces RMSD time series, average energy, and a stability verdict

3. 3D VISUALIZATION
   - Interactive 3Dmol.js viewer for protein-ligand complexes
   - See H-bond donors/acceptors, hydrophobic contacts, pi-stacking in 3D
   - RMSD analysis across multiple docking poses

4. AI-ASSISTED WORKFLOWS
   - Describe what you want in plain English — I compile it into an executable workflow
   - Docking Analysis Crew: Dock → Consensus score → Interaction analysis → Report
   - Lead Optimization Crew: Dock → Interactions → Chemistry → Redock → Rank
   - MD Simulation Crew: Suggest params → Run MD → Interpret RMSD → Stability report
   - Full Drug Discovery Crew: Dock → Analyze → MD validation → Final lead report

═══════════════════════════════════════════════════════════
YOUR AGENT TEAM — CREWAI MULTI-AGENT SYSTEM
═══════════════════════════════════════════════════════════
You coordinate 4 specialized AI agents built for the student edition pipeline:

1. DOCKING AGENT
   Role: Molecular Docking Specialist
   Tools: run_docking, calculate_properties
   Knows: AutoDock Vina scoring, binding pocket geometry, pose ranking
   Rule: Energy ≤ -5 kcal/mol → strong candidate; flag to student for MD validation
   Self-evolving: Learns optimal exhaustiveness and box_size per protein family

2. ANALYSIS AGENT
   Role: Drug Discovery Analysis Expert
   Tools: analyze_interactions, rank_ligands, consensus_score, export_top_hits
   Knows: H-bonds (≤3.5 Å, ≥120°), hydrophobic contacts, pi-stacking, salt bridges
   Speciality: Interaction fingerprinting, pose ranking, binding mode interpretation

3. AGENT MD (Molecular Dynamics Specialist)
   Role: MD Simulation & Trajectory Interpretation Expert
   Tools: run_md_simulation, analyze_md_trajectory, interpret_rmsd, suggest_md_parameters
   Knows: OpenMM setup (TIP3P/implicit solvent), NVT/NPT equilibration, production MD
   Interprets: RMSD < 2 Å = stable | 2-3 Å = borderline | > 3 Å = unstable/pose collapse
   Per-family params:
     • Kinase:          100k steps, 300K, TIP3P — flexible hinge, sample thoroughly
     • GPCR:            200k steps, 310K — membrane-like, implicit solvent fallback
     • Protease:        75k steps,  300K — rigid active site, short simulation sufficient
     • Nuclear receptor: 150k steps, 300K — watch helix H12 movement (agonist vs antagonist)
     • Enzyme:          80k steps,  300K — focus on active site loop dynamics
   Diagnoses: Protonation state errors, force field mismatches, insufficient equilibration

4. ORCHESTRATOR
   Role: Drug Discovery Orchestrator
   Tools: send_notification, export_top_hits, rank_ligands
   Knows: Full hit-to-lead pipeline logic, how to synthesize agent results
   Always: Recommends top candidates with scientific rationale and a concrete next-step plan

═══════════════════════════════════════════════════════════
SELF-EVOLVING CREWAI INTELLIGENCE
═══════════════════════════════════════════════════════════
You are not static. You and your crew grow smarter with every experiment:

• META-PARAMETER SELF-LEARNING (MetaParameterLearner)
  Tracks which exhaustiveness, box_size, temperature, solvent model worked best
  per protein family. Suggests optimal settings learned from historical runs.
  Family classification: kinase / GPCR / protease / nuclear_receptor / ion_channel / enzyme

• CHROMADB EXPERIMENT MEMORY (ExperimentMemory)
  Every docking, MD, ADMET, batch run is stored in ChromaDB with semantic indexing.
  You can query: "find similar experiments to this kinase docking" and retrieve
  the most relevant past runs by semantic similarity — not just keyword matching.
  Learns failure patterns: recurring grid errors → auto-suggest box_size increase.

• BAYESIAN ACTIVE LEARNING (BayesianOptimizer)
  Uses Gaussian Process with Matern kernel + Expected Improvement acquisition.
  After initial screening, intelligently selects next compounds to maximize
  the chance of finding high-affinity binders. Gets smarter with each iteration.

• ADVERSARIAL CRITIQUE AGENT (CritiqueAgent)
  Validates every result before accepting it:
  - Energy bounds check (docking: -20 to +20 kcal/mol)
  - Red flag detection: NaN, RMSD > 5 Å, excessive clashes
  - Confidence gating: rejects proposals below 0.7 confidence threshold
  - Cross-references with PubChem for known activity data

• NL-TO-DAG COMPILER (NLWorkflowCompiler + SelfHealingExecutor)
  You parse plain English requests into executable DAG workflows.
  Self-healing: if a step fails, auto-adjusts parameters and retries.
  Box too small? → auto-expand by 4 Å. Timeout? → halve exhaustiveness. 

• KNOWLEDGE GRAPH (BioKnowledgeGraph)
  Connects targets, compounds, pathways, diseases, and experiments.
  Knows: kinase → MAPK/PI3K-Akt pathways, protease → apoptosis, etc.
  Stores every experiment as a graph node linked to compound + target.

═══════════════════════════════════════════════════════════
YOUR PERSONALITY & PROACTIVE BEHAVIORS
═══════════════════════════════════════════════════════════
PERSONALITY:
- You are a brilliant senior researcher mentoring a student — warm, sharp, scientifically honest.
- You explain WHY results mean what they mean, not just state them.
- You proactively notice patterns the student hasn't noticed yet.
- You reference past experiments from memory naturally ("Earlier you docked X against Y...").
- You celebrate good results and turn bad ones into learning moments.
- You never just say "here is the answer" — you teach while answering.

PROACTIVE TRIGGERS:
- After docking result → immediately assess score quality, suggest ADMET or MD as next step
- After MD trajectory → interpret RMSD trend, diagnose instability if present
- After ADMET prediction → cross-reference flags with docking score, flag hERG risk
- After batch screening → identify top hit, flag outliers, suggest lead optimization
- When student asks general question → pull relevant experiment history to contextualize
- On session start → greet by name, summarize recent experiments from ChromaDB memory

═══════════════════════════════════════════════════════════
SCORING & INTERPRETATION GUIDE
═══════════════════════════════════════════════════════════
DOCKING (AutoDock Vina, kcal/mol — more negative = stronger binding):
  ≤ -10.0  : Excellent — pursue full validation (ADMET + MD + interaction analysis)
  -10 to -8: Strong    — worth ADMET + MD, good lead candidate
  -8 to -6 : Moderate  — viable with structural optimization
  -6 to -4 : Weak      — redesign compound or verify docking box position
  > -4.0   : Poor      — likely non-binder, reconsider scaffold

MD STABILITY (RMSD in Å from initial pose):
  < 2.0 Å  : Stable    — binding pose maintained, trust the docking result
  2.0-3.0 Å: Borderline — minor fluctuation, extend equilibration
  > 3.0 Å  : Unstable  — pose collapsed, investigate protonation/force field

LIGAND EFFICIENCY (LE = |ΔG| / heavy atom count):
  > 0.4    : Excellent fragment-like efficiency
  0.3-0.4  : Good
  < 0.3    : Needs optimization (too large or too weak)

QED DRUG-LIKENESS (0-1 scale):
  > 0.7    : Drug-like
  0.5-0.7  : Acceptable
  < 0.5    : Poor drug-likeness — flag for modification

WHEN DOCKING IS WEAK — DIAGNOSTIC CHECKLIST:
  1. Is the docking box centered on the actual binding site?
  2. Is the ligand protonated correctly at physiological pH?
  3. Is the receptor prepared (hydrogens added, waters removed)?
  4. Try increasing exhaustiveness (8 → 16 → 32) for more thorough search
  5. Check if ligand SMILES is valid using RDKit before docking"""


# ── Conversation History Store ────────────────────────────────────────────────
_HISTORY_FILE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "storage", "nanobot", "chat_history.json"
)
MAX_HISTORY_TURNS = 16  # keep last 16 turns (8 exchanges)


def _load_history() -> List[Dict]:
    try:
        if os.path.exists(_HISTORY_FILE):
            with open(_HISTORY_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
    except Exception:
        pass
    return []


def _save_history(history: List[Dict]):
    try:
        Path(_HISTORY_FILE).parent.mkdir(parents=True, exist_ok=True)
        with open(_HISTORY_FILE, "w", encoding="utf-8") as f:
            json.dump(history[-MAX_HISTORY_TURNS:], f, indent=2, ensure_ascii=False)
    except Exception as e:
        logger.warning(f"Failed to save chat history: {e}")


def _append_history(role: str, content: str):
    history = _load_history()
    history.append({"role": role, "content": content, "ts": datetime.now().isoformat()})
    _save_history(history)


def _get_history_messages(last_n: int = 10) -> List[Dict]:
    """Return last N turns as [{role, content}] without timestamps."""
    history = _load_history()
    return [{"role": h["role"], "content": h["content"]} for h in history[-last_n:]]


def clear_chat_history():
    """Wipe conversation history."""
    _save_history([])


# ── Memory Context Builder ────────────────────────────────────────────────────

def _build_memory_context() -> str:
    """Pull recent job history from crew/memory.py ChromaDB/JSON for NanoBot context."""
    try:
        from crew.memory import memory as exp_memory
        stats = exp_memory.get_stats()
        recent = exp_memory.get_job_history(n=5)
        if not recent:
            return ""
        lines = [f"[EXPERIMENT MEMORY] {stats.get('total_experiments', 0)} jobs on record, "
                 f"success rate {stats.get('success_rate', 0):.0%}"]
        for exp in recent:
            meta = exp.get("meta", {})
            result = exp.get("result", {})
            t = exp.get("type", meta.get("type", "job"))
            ts = str(exp.get("timestamp", ""))[:10]
            score = result.get("best_score", result.get("binding_energy", result.get("energy", "")))
            status = exp.get("status", "")
            target = meta.get("target", "")
            ligand = meta.get("scaffold", meta.get("smiles", ""))[:30]
            lines.append(f"  • [{ts}] {t} | ligand={ligand} target={target} score={score} status={status}")
        chroma = exp_memory.chroma_stats()
        lines.append(f"  [ChromaDB: {chroma.get('backend')} — {chroma.get('total_indexed')} indexed]")
        return "\n".join(lines)
    except Exception as e:
        logger.debug(f"Memory context unavailable: {e}")
        return ""

_CONFIG_FILE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "llm_config.json"
)


def _load_config() -> Dict:
    """Load LLM config from file"""
    try:
        if os.path.exists(_CONFIG_FILE):
            with open(_CONFIG_FILE, "r") as f:
                return json.load(f)
    except Exception as e:
        logger.warning(f"Failed to load LLM config: {e}")
    return {}


def save_config(config: Dict):
    """Save LLM config to file"""
    try:
        with open(_CONFIG_FILE, "w") as f:
            json.dump(config, f, indent=2)
        logger.info(f"LLM config saved: provider={config.get('provider', 'unknown')}")
    except Exception as e:
        logger.error(f"Failed to save LLM config: {e}")


PROVIDER_URLS = {
    "ollama": "http://host.docker.internal:11434/v1",
    "openai": "https://api.openai.com/v1",
    "anthropic": "https://api.anthropic.com/v1",
    "gemini": "https://generativelanguage.googleapis.com/v1beta",
    "deepseek": "https://api.deepseek.com/v1",
    "mistral": "https://api.mistral.ai/v1",
    "groq": "https://api.groq.com/openai/v1",
    "openrouter": "https://openrouter.ai/api/v1",
    "siliconflow": "https://api.siliconflow.cn/v1",
    "qwen": "https://dashscope.aliyuncs.com/compatible-mode/v1",
}

PROVIDER_MODELS = {
    "ollama": "llama3.2",
    "openai": "gpt-4o-mini",
    "anthropic": "claude-sonnet-4-20250514",
    "gemini": "gemini-2.0-flash",
    "deepseek": "deepseek-chat",
    "mistral": "mistral-small-latest",
    "groq": "llama-3.1-8b-instant",
    "openrouter": "meta-llama/llama-3.1-8b-instruct",
    "siliconflow": "Qwen/Qwen2.5-7B-Instruct",
    "qwen": "qwen-turbo",
}


class OllamaProvider:
    """Ollama API provider"""

    def __init__(self, url: str = OLLAMA_URL, model: str = OLLAMA_MODEL):
        self.url = url
        self.model = model

    def is_available(self) -> bool:
        try:
            response = requests.get(f"{self.url}/api/tags", timeout=3)
            if response.status_code == 200:
                data = response.json()
                models = data.get("models", [])
                logger.info(f"Ollama reachable, {len(models)} model(s) listed")
                return True  # reachable = available; model list may be empty while loading
            return False
        except Exception:
            return False

    def chat(self, messages: List[Dict]) -> str:
        """Send a pre-built messages list (system + history + user) to Ollama."""
        headers = {"Content-Type": "application/json"}
        payload = {
            "model": self.model,
            "messages": messages,
            "stream": False,
        }
        response = requests.post(
            f"{self.url}/api/chat",
            json=payload,
            headers=headers,
            timeout=OLLAMA_TIMEOUT * 2,
        )
        response.raise_for_status()
        data = response.json()
        return data["message"]["content"]

    def get_models(self) -> list:
        try:
            response = requests.get(f"{self.url}/api/tags", timeout=OLLAMA_TIMEOUT)
            if response.status_code == 200:
                data = response.json()
                return [m["name"] for m in data.get("models", [])]
        except Exception as e:
            logger.warning(f"Could not fetch Ollama models: {e}")
        return []


class APIProvider:
    """OpenAI-compatible API provider (DeepSeek, OpenAI, Mistral, Groq, etc.)"""

    def __init__(
        self, provider: str, api_key: str, base_url: str = "", model: str = ""
    ):
        self.provider = provider
        self.api_key = api_key
        self.base_url = (base_url or PROVIDER_URLS.get(provider, "")).rstrip("/")
        self.model = model or PROVIDER_MODELS.get(provider, "gpt-4o-mini")

    def is_available(self) -> bool:
        """Check if the API is reachable with a minimal request"""
        if not self.api_key or not self.base_url:
            return False
        try:
            headers = {
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
            }
            resp = requests.get(f"{self.base_url}/models", headers=headers, timeout=10)
            return resp.status_code in (
                200,
                401,
                403,
            )  # 401/403 means reachable but auth issue
        except Exception:
            # Try a minimal chat completion as fallback check
            try:
                headers = {
                    "Authorization": f"Bearer {self.api_key}",
                    "Content-Type": "application/json",
                }
                resp = requests.post(
                    f"{self.base_url}/chat/completions",
                    headers=headers,
                    json={
                        "model": self.model,
                        "messages": [{"role": "user", "content": "hi"}],
                        "max_tokens": 1,
                    },
                    timeout=15,
                )
                return resp.status_code in (200, 400, 401, 403, 429)
            except Exception:
                return False

    def chat(self, messages: List[Dict]) -> str:
        """Send a pre-built messages list (system + history + user) to the API provider."""
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }
        resp = requests.post(
            f"{self.base_url}/chat/completions",
            headers=headers,
            json={"model": self.model, "messages": messages, "max_tokens": 2048},
            timeout=60,
        )
        resp.raise_for_status()
        data = resp.json()
        return data["choices"][0]["message"]["content"]


class LLMRouter:
    """
    Intelligent LLM router that:
    - Reads config from llm_config.json (saved by Settings page)
    - Supports Ollama, OpenAI, DeepSeek, and OpenAI-compatible APIs
    - Falls back to offline assistant on failure
    """

    def __init__(self):
        self._config = _load_config()
        self._init_providers()

    def _init_providers(self):
        """Initialize provider instances with URL normalization for Docker."""
        saved_model = (
            self._config.get("model")
            if self._config.get("provider") == "ollama"
            else None
        )
        saved_url = (
            self._config.get("base_url")
            if self._config.get("provider") == "ollama"
            else None
        )

        # Normalize URL: replace localhost with host.docker.internal for Docker compatibility
        if saved_url and "localhost" in saved_url:
            saved_url = saved_url.replace("localhost", "host.docker.internal")

        # Strip /v1 suffix (OllamaProvider adds it itself for chat)
        ollama_base = (
            saved_url.rstrip("/").removesuffix("/v1") if saved_url else None
        ) or OLLAMA_URL
        self.ollama = OllamaProvider(url=ollama_base, model=saved_model or OLLAMA_MODEL)
        self.offline = OfflineAssistant()
        self._provider = None
        self._api_provider = None
        logger.info(
            f"LLMRouter initialized (ollama url={ollama_base}, model={self.ollama.model})"
        )

    def reset(self):
        """Reset provider detection and re-read config."""
        self._config = _load_config()
        self._init_providers()

    def _get_config_provider(self) -> str:
        return self._config.get("provider", "ollama")

    def _get_api_key(self) -> str:
        return self._config.get("api_key", "")

    def _get_base_url(self) -> str:
        return self._config.get("base_url", "")

    def _get_model(self) -> str:
        return self._config.get("model", "")

    @property
    def provider(self) -> str:
        if self._provider is None or self._provider == "offline":
            self._provider = self._detect_provider()
        return self._provider

    def _detect_provider(self) -> str:
        if not ALLOW_AI:
            logger.info("AI disabled via config")
            return "offline"

        config_provider = self._get_config_provider()

        # Ollama provider
        if config_provider == "ollama":
            if self.ollama.is_available():
                logger.info("Ollama available")
                return "ollama"
            logger.warning("Ollama configured but not available")
            return "offline"

        # API-based providers (DeepSeek, OpenAI, etc.)
        if config_provider in PROVIDER_URLS:
            api_key = self._get_api_key()
            if not api_key:
                logger.warning(f"{config_provider} configured but no API key")
                return "offline"
            base_url = self._get_base_url() or PROVIDER_URLS.get(config_provider, "")
            model = self._get_model() or PROVIDER_MODELS.get(config_provider, "")
            self._api_provider = APIProvider(config_provider, api_key, base_url, model)
            if self._api_provider.is_available():
                logger.info(f"{config_provider} API available")
                return config_provider
            logger.warning(f"{config_provider} API not reachable")
            return "offline"

        # Custom OpenAI-compatible
        if config_provider == "custom":
            api_key = self._get_api_key()
            base_url = self._get_base_url()
            model = self._get_model()
            if not base_url:
                logger.warning("Custom provider configured but no base URL")
                return "offline"
            self._api_provider = APIProvider("custom", api_key, base_url, model)
            if self._api_provider.is_available():
                logger.info("Custom API available")
                return "custom"
            return "offline"

        return "offline"

    def detect_ollama(self) -> bool:
        return self.ollama.is_available()

    def detect_provider(self) -> bool:
        """Check if the configured provider is available"""
        config_provider = self._get_config_provider()
        if config_provider == "ollama":
            return self.ollama.is_available()
        if self._api_provider:
            return self._api_provider.is_available()
        return False

    def get_available_models(self) -> list:
        if self.provider == "ollama":
            return self.ollama.get_models()
        return []

    def chat(self, message: str, job_context: str = "") -> Dict:
        """
        Send a chat message as NanoBot.
        Injects soul, experiment memory context, conversation history, and live job context.
        Returns dict with: response, provider, available, memory_context
        """
        detected = self.provider

        # Build system prompt: soul + live memory context + live job context
        memory_ctx = _build_memory_context()
        system_content = NANOBOT_SOUL
        if memory_ctx:
            system_content += f"\n\n{memory_ctx}"
        if job_context:
            system_content += f"\n\n{job_context}"

        # Build messages: system + history + current user message
        messages: List[Dict] = [{"role": "system", "content": system_content}]
        messages.extend(_get_history_messages(last_n=10))
        messages.append({"role": "user", "content": message})

        # Save user turn to history
        _append_history("user", message)

        # Ollama
        if detected == "ollama":
            try:
                response_text = self.ollama.chat(messages)
                _append_history("assistant", response_text)
                return {
                    "response": response_text,
                    "provider": "ollama",
                    "available": True,
                    "memory_context": bool(memory_ctx),
                }
            except Exception as e:
                logger.warning(f"Ollama failed: {e}")
                response_text = self.offline.respond(message)
                _append_history("assistant", response_text)
                return {
                    "response": response_text,
                    "provider": "offline",
                    "available": False,
                    "error": str(e),
                }

        # API providers (OpenAI, DeepSeek, Groq, custom, etc.)
        if detected != "offline" and self._api_provider:
            try:
                response_text = self._api_provider.chat(messages)
                _append_history("assistant", response_text)
                return {
                    "response": response_text,
                    "provider": detected,
                    "available": True,
                    "memory_context": bool(memory_ctx),
                }
            except Exception as e:
                logger.warning(f"{detected} failed: {e}")
                response_text = self.offline.respond(message)
                _append_history("assistant", response_text)
                return {
                    "response": response_text,
                    "provider": "offline",
                    "available": False,
                    "error": str(e),
                }

        # Offline fallback
        response_text = self.offline.respond(message)
        _append_history("assistant", response_text)
        return {"response": response_text, "provider": "offline", "available": False, "memory_context": False}

    def reset(self):
        """Reset provider cache and reload config"""
        self._provider = None
        self._api_provider = None
        self._config = _load_config()
        self._init_providers()


def get_router() -> LLMRouter:
    """Get singleton router instance"""
    if not hasattr(get_router, "_instance"):
        get_router._instance = LLMRouter()
    return get_router._instance
