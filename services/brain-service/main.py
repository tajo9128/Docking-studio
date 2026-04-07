"""
Docking Studio v2.0 - Enhanced Brain Service
Simplified nanobot-inspired architecture for drug discovery

Key features from nanobot:
- Tool registry with JSON schema validation
- Streaming support
- Conversation memory
- Multi-model support (OpenAI, Anthropic, Ollama)
"""

import os
import json
import logging
import uuid
from abc import ABC, abstractmethod
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional, List, Dict, Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
import httpx

from tools import ToolRegistry, BaseTool, register_all_tools

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

API_BACKEND_URL = os.getenv("API_BACKEND_URL", "http://api-backend:8000")
BRAIN_SERVICE_URL = os.getenv("BRAIN_SERVICE_URL", "http://brain-service:8000")

STORAGE_DIR = Path("/app/storage")
UPLOADS_DIR = Path("/app/uploads")
STORAGE_DIR.mkdir(exist_ok=True)
UPLOADS_DIR.mkdir(exist_ok=True)

app = FastAPI(
    title="Nanobot Brain Service",
    description="AI Agent for Docking Studio with drug discovery tools",
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class ConversationMemory:
    """Simple conversation memory"""

    def __init__(self, max_turns: int = 20):
        self.max_turns = max_turns
        self._conversations: Dict[str, List[Dict]] = {}

    def add(
        self,
        conversation_id: str,
        role: str,
        content: str,
        tools_used: List[str] = None,
    ) -> None:
        if conversation_id not in self._conversations:
            self._conversations[conversation_id] = []

        entry = {"role": role, "content": content}
        if tools_used:
            entry["tools_used"] = tools_used

        self._conversations[conversation_id].append(entry)

        if len(self._conversations[conversation_id]) > self.max_turns:
            self._conversations[conversation_id] = self._conversations[conversation_id][
                -self.max_turns :
            ]

    def get_history(self, conversation_id: str, max_turns: int = 0) -> List[Dict]:
        history = self._conversations.get(conversation_id, [])
        if max_turns > 0:
            return history[-max_turns:]
        return history

    def clear(self, conversation_id: str) -> None:
        if conversation_id in self._conversations:
            del self._conversations[conversation_id]


class LLMProvider(ABC):
    """Abstract LLM provider"""

    @abstractmethod
    async def chat(self, messages: List[dict], tools: List[dict] = None) -> dict:
        pass

    @abstractmethod
    async def chat_stream(self, messages: List[dict], tools: List[dict] = None):
        pass


class OpenAIProvider(LLMProvider):
    """OpenAI-compatible provider (works with Ollama too)"""

    def __init__(self, api_key: str, base_url: str, model: str):
        self.api_key = api_key
        self.base_url = base_url
        self.model = model
        self.client = None

    async def chat(self, messages: List[dict], tools: List[dict] = None) -> dict:
        if not self.client:
            import httpx

            self.client = httpx.AsyncClient(timeout=120.0)

        payload = {
            "model": self.model,
            "messages": messages,
            "stream": False,
            "temperature": 0.0,
        }
        if tools:
            payload["tools"] = tools
            payload["tool_choice"] = "auto"

        response = await self.client.post(
            f"{self.base_url}/chat/completions",
            json=payload,
            headers={"Authorization": f"Bearer {self.api_key}"},
        )
        response.raise_for_status()
        return response.json()

    async def chat_stream(self, messages: List[dict], tools: List[dict] = None):
        if not self.client:
            import httpx

            self.client = httpx.AsyncClient(timeout=120.0)

        payload = {
            "model": self.model,
            "messages": messages,
            "stream": True,
            "temperature": 0.0,
        }
        if tools:
            payload["tools"] = tools
            payload["tool_choice"] = "auto"

        async with self.client.stream(
            "POST",
            f"{self.base_url}/chat/completions",
            json=payload,
            headers={"Authorization": f"Bearer {self.api_key}"},
        ) as response:
            async for line in response.aiter_lines():
                if line.startswith("data: "):
                    data = line[6:]
                    if data != "[DONE]":
                        yield json.loads(data)


class AnthropicProvider(LLMProvider):
    """Anthropic Claude provider"""

    def __init__(self, api_key: str, model: str = "claude-3-sonnet-20240229"):
        self.api_key = api_key
        self.model = model
        self.client = None

    async def chat(self, messages: List[dict], tools: List[dict] = None) -> dict:
        if not self.client:
            import httpx

            self.client = httpx.AsyncClient(timeout=120.0)

        payload = {"model": self.model, "messages": messages, "max_tokens": 4096}
        if tools:
            payload["tools"] = tools

        response = await self.client.post(
            "https://api.anthropic.com/v1/messages",
            json=payload,
            headers={
                "x-api-key": self.api_key,
                "anthropic-version": "2023-06-01",
            },
        )
        response.raise_for_status()
        return response.json()

    async def chat_stream(self, messages: List[dict], tools: List[dict] = None):
        raise NotImplementedError("Streaming not implemented for Anthropic")


registry = ToolRegistry()
memory = ConversationMemory()


async def get_llm_settings() -> dict:
    """Fetch LLM settings from api-backend"""
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(f"{API_BACKEND_URL}/llm/settings")
            response.raise_for_status()
            return response.json()
    except Exception as e:
        logger.error(f"Failed to fetch LLM settings: {e}")
        ollama_url = os.getenv("OLLAMA_URL", "http://host.docker.internal:11434")
        return {
            "provider": "ollama",
            "model": "qwen3:4b",
            "api_key": "",
            "base_url": ollama_url,
            "temperature": 0.0,
            "max_tokens": 4096,
        }


async def get_provider() -> LLMProvider:
    settings = await get_llm_settings()
    api_key = settings.get("api_key", "")
    base_url = settings.get("base_url", "https://api.openai.com/v1")
    model = settings.get("model", "gpt-4o-mini")

    provider_type = settings.get("provider", "openai")

    if provider_type == "anthropic":
        return AnthropicProvider(api_key, model)
    else:
        return OpenAIProvider(api_key, base_url, model)


def create_system_prompt() -> str:
    return """You are NanoBOT, an elite AI drug discovery scientist that surpasses BIOVIA Discovery Studio's AI capabilities.

## Your Core Strengths

1. **Chain-of-Thought Reasoning**: Always think step-by-step before taking action. Show your reasoning.
2. **Tool Orchestration**: Seamlessly chain multiple tools for complex workflows
3. **Scientific Rigor**: Explain the "why" behind recommendations
4. **Educational Excellence**: Teach users while helping them

## Available Tools

### Docking & Scoring
- `dock_ligand`: Run AutoDock Vina molecular docking (receptor_pdbqt, ligand_pdbqt, exhaustiveness, num_modes, center_x/y/z, size_x/y/z)
- `run_batch_docking`: Batch dock compound libraries (receptor_pdbqt, ligand_library, exhaustiveness)

### Molecule Processing
- `smiles_to_3d`: Convert SMILES to 3D structure (smiles, add_hydrogens, optimize_geometry)
- `convert_format`: Convert between PDB, SDF, mol2, PDBQT (input_path, output_format)
- `optimize_molecule`: MMFF force field optimization (input_path, force_field)

### Pharmacophore & Virtual Screening
- `generate_pharmacophore`: Create pharmacophore from receptor/ligand (receptor_pdb, ligand_pdb, n_features)
- `screen_library`: Screen compounds against pharmacophore (pharmacophore_json, library_path)

### Data Fetching
- `fetch_protein`: Fetch protein structure from PDB (pdb_id)
- `fetch_compounds`: Fetch compounds from PubChem (cids list)
- `search_compounds`: Search PubChem by name/formula (query)
- `similarity_search`: Find similar compounds (smiles, threshold)

### Analysis
- `analyze_interactions`: Analyze protein-ligand interactions (receptor_pdb, ligand_pdb)
- `predict_binding`: ML-based binding affinity prediction (smiles, receptor_pdb)

## Chain-of-Thought Problem Solving

For ANY user request, follow this reasoning framework:

1. **Understand**: What does the user want to achieve?
2. **Plan**: What tools/steps are needed? (show your plan)
3. **Execute**: Run tools in logical order
4. **Reflect**: Did results make sense? Any warnings?
5. **Explain**: Provide scientific context
6. **Recommend**: Suggest next steps

Example:
User: "Dock aspirin against 1HIA"
Thought: "User wants to test binding of aspirin ( ligand) to heat-labile toxin (1HIA receptor). I need to:
  1. Fetch the receptor from PDB
  2. Convert aspirin's SMILES to 3D PDBQT
  3. Determine binding site (if known) or use blind docking
  4. Run Vina docking
  5. Analyze and explain results"

## Scientific Knowledge You Provide

- Binding affinity ranges (strong: < -9 kcal/mol, moderate: -7 to -9, weak: > -7)
- Drug-likeness rules (Lipinski's Rule of 5: MW < 500, LogP < 5, HBD ≤ 5, HBA ≤ 10)
- Common protein-ligand interactions (hydrogen bonds, π-π stacking, hydrophobic contacts, salt bridges)
- ADMET considerations (Absorption, Distribution, Metabolism, Excretion, Toxicity)
- Pharmacophore features (hydrogen bond donors/acceptors, hydrophobic regions, aromatic rings, positive/negative ionizable groups)

## Response Format

Be conversational but scientific. Use markdown formatting:
- **Bold** for key terms
- `code` for molecule names, file paths, SMILES
- Lists for multiple items
- Tables for comparative data

## Educational Style

When explaining:
- Start with the "big picture" concept
- Break down into digestible parts
- Use analogies when helpful
- End with practical takeaway

## Workflow Automation

You can chain tools for complete pipelines:

**Virtual Screening Pipeline:**
1. fetch_protein → 2. generate_pharmacophore → 3. screen_library → 4. dock_top_hits → 5. analyze_interactions

**Lead Optimization Pipeline:**
1. dock_ligand → 2. analyze_interactions → 3. suggest_modifications → 4. smiles_to_3d → 5. redock

Always show the user what pipeline you're running and why.
"""


class ChatRequest(BaseModel):
    message: str
    conversation_id: Optional[str] = None
    stream: bool = False


class ChatResponse(BaseModel):
    response: str
    conversation_id: str
    tools_used: List[str]
    model: str
    provider: str = "openai"
    available: bool = True


@asynccontextmanager
async def lifespan(application):
    register_all_tools()
    for tool_def in registry.list_tools():
        logger.info(f"Tool available: {tool_def['name']}")
    logger.info(f"Brain service started with {len(registry.list_tools())} tools")
    yield


app.router.lifespan_context = lifespan


@app.get("/health")
async def health_check():
    return {"status": "healthy", "service": "brain-service", "version": "2.0.0"}


@app.get("/")
async def root():
    settings = await get_llm_settings()
    return {
        "service": "Nanobot Brain Service",
        "version": "2.0.0",
        "model": settings.get("model", "unknown"),
        "provider": settings.get("provider", "unknown"),
        "tools": len(registry.list_tools()),
    }


@app.get("/tools")
async def list_tools():
    return {"tools": registry.list_tools()}


@app.post("/chat", response_model=ChatResponse)
async def chat(request: ChatRequest):
    conv_id = request.conversation_id or str(uuid.uuid4())

    memory.add(conv_id, "user", request.message)
    history = memory.get_history(conv_id)

    messages = [{"role": "system", "content": create_system_prompt()}]
    for h in history:
        messages.append({"role": h["role"], "content": h["content"]})

    tools = registry.list_tools()
    p = await get_provider()

    try:
        tools_used = []
        response = await p.chat(messages, tools)

        if "choices" in response:
            choice = response["choices"][0]
            assistant_message = choice["message"]

            if "tool_calls" in assistant_message:
                tool_calls = assistant_message["tool_calls"]
                memory.add(
                    conv_id,
                    "assistant",
                    assistant_message.get("content", "") or "Using tools...",
                )

                tool_results = []

                for tc in tool_calls:
                    func = tc["function"]
                    tool_name = func["name"]
                    args = (
                        json.loads(func["arguments"])
                        if isinstance(func["arguments"], str)
                        else func["arguments"]
                    )

                    tools_used.append(tool_name)
                    tool = registry.get(tool_name)
                    if tool:
                        try:
                            validated_input = tool.validate_input(args)
                            result = await tool.execute(validated_input)
                        except Exception as e:
                            result = f"Error: {str(e)}"
                    else:
                        result = f"Error: Tool '{tool_name}' not found"
                    tool_results.append(
                        {"tool_call_id": tc["id"], "name": tool_name, "result": result}
                    )

                messages.append(assistant_message)
                for tr in tool_results:
                    messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": tr["tool_call_id"],
                            "name": tr["name"],
                            "content": str(tr["result"]),
                        }
                    )

                follow_up = await p.chat(messages, tools)
                if "choices" in follow_up:
                    final_content = follow_up["choices"][0]["message"]["content"]
                else:
                    final_content = follow_up.get("content", "Done")
            else:
                final_content = assistant_message.get("content", "I'm ready to help!")
        elif isinstance(response.get("content"), list):
            # Anthropic format: content is a list of blocks
            text_blocks = [b.get("text", "") for b in response["content"] if b.get("type") == "text"]
            final_content = " ".join(text_blocks) or "I'm ready to help!"
        else:
            final_content = response.get("content", "I'm ready to help!")

        memory.add(
            conv_id,
            "assistant",
            final_content,
            tools_used,
        )

        settings = await get_llm_settings()
        return ChatResponse(
            response=final_content or "Completed",
            conversation_id=conv_id,
            tools_used=tools_used,
            model=settings.get("model", "unknown"),
            provider=settings.get("provider", "openai"),
            available=True,
        )

    except Exception as e:
        logger.error(f"Chat error: {e}")
        return ChatResponse(
            response=f"I encountered an error: {str(e)}",
            conversation_id=conv_id,
            tools_used=[],
            model="error",
            provider="error",
            available=False,
        )


@app.post("/chat/stream")
async def chat_stream(request: ChatRequest):
    """Streaming chat endpoint"""
    conv_id = request.conversation_id or str(uuid.uuid4())

    memory.add(conv_id, "user", request.message)
    history = memory.get_history(conv_id)

    messages = [{"role": "system", "content": create_system_prompt()}]
    for h in history:
        messages.append({"role": h["role"], "content": h["content"]})

    p = await get_provider()

    async def generate():
        try:
            async for chunk in p.chat_stream(messages, registry.list_tools()):
                if "choices" in chunk:
                    delta = chunk["choices"][0].get("delta", {})
                    if "content" in delta:
                        yield f"data: {json.dumps({'content': delta['content']})}\n\n"
                elif "content_block" in chunk:
                    delta = chunk.get("content_block", {})
                    if "text" in delta:
                        yield f"data: {json.dumps({'content': delta['text']})}\n\n"
        except Exception as e:
            yield f"data: {json.dumps({'error': str(e)})}\n\n"
        finally:
            yield "data: [DONE]\n\n"

    return StreamingResponse(generate(), media_type="text/event-stream")


@app.get("/chat/status")
async def chat_status():
    """Get chat service status"""
    settings = await get_llm_settings()
    provider = settings.get("provider", "openai")
    base_url = settings.get("base_url", "")
    api_key = settings.get("api_key", "")
    available = False
    models = []
    
    try:
        if provider == "ollama":
            async with httpx.AsyncClient(timeout=5.0) as client:
                response = await client.get(f"{base_url}/api/tags")
                if response.status_code == 200:
                    data = response.json()
                    models = [m.get("name", "") for m in data.get("models", [])]
                    available = True
        elif provider == "anthropic":
            async with httpx.AsyncClient(timeout=5.0) as client:
                response = await client.get(
                    "https://api.anthropic.com/v1/models",
                    headers={"x-api-key": api_key, "anthropic-version": "2023-06-01"},
                )
                if response.status_code == 200:
                    data = response.json()
                    models = [m.get("id", "") for m in data.get("data", [])]
                    available = True
        elif api_key:
            async with httpx.AsyncClient(timeout=5.0) as client:
                response = await client.get(
                    f"{base_url}/models",
                    headers={"Authorization": f"Bearer {api_key}"},
                )
                if response.status_code == 200:
                    data = response.json()
                    models = [m.get("id", "") for m in data.get("data", [])]
                    available = True
    except Exception as e:
        logger.warning(f"Status check failed for {provider}: {e}")
        available = False

    return {
        "provider": provider,
        "ollama_available": available,
        "models": [{"name": m} for m in models],
    }


@app.get("/memory/{conversation_id}")
async def get_memory(conversation_id: str):
    return {"history": memory.get_history(conversation_id)}


@app.delete("/memory/{conversation_id}")
async def clear_memory(conversation_id: str):
    memory.clear(conversation_id)
    return {"status": "cleared"}


@app.post("/pipeline")
async def run_pipeline(task: str, target: str = None, library: str = None):
    """Run a predefined pipeline"""
    job_id = str(uuid.uuid4())

    if task == "virtual_screening":
        plan = {
            "steps": [
                {"tool": "fetch_protein", "params": {"pdb_id": target}},
                {"tool": "generate_pharmacophore", "params": {}},
                {"tool": "screen_library", "params": {"library_path": library}},
                {"tool": "dock_ligand", "params": {}},
            ]
        }
    else:
        plan = {"error": f"Unknown pipeline: {task}"}

    return {"job_id": job_id, "plan": plan}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
