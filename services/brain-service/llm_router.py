from enum import Enum
from typing import Dict, List, Optional, Any, AsyncIterator
import httpx
import os
import json

class LLMProvider(Enum):
    OPENAI = "openai"
    ANTHROPIC = "anthropic"
    GEMINI = "gemini"
    OPENROUTER = "openrouter"
    MISTRAL = "mistral"
    SILICONFLOW = "siliconflow"
    DEEPSEEK = "deepseek"
    QWEN = "qwen"
    OLLAMA = "ollama"

PROVIDER_CONFIG = {
    LLMProvider.OPENAI: {
        "models": ["gpt-4o", "gpt-4o-mini", "gpt-4-turbo"],
        "base_url": "https://api.openai.com/v1",
    },
    LLMProvider.ANTHROPIC: {
        "models": ["claude-3-5-sonnet-20240620", "claude-3-opus-20240229", "claude-3-haiku-20240307"],
        "base_url": "https://api.anthropic.com",
    },
    LLMProvider.GEMINI: {
        "models": ["gemini-1.5-pro", "gemini-1.5-flash", "gemini-pro"],
        "base_url": "https://generativelanguage.googleapis.com/v1beta",
    },
    LLMProvider.OPENROUTER: {
        "models": ["anthropic/claude-3.5-sonnet", "google/gemini-pro", "mistralai/mistral-large"],
        "base_url": "https://openrouter.ai/api/v1",
    },
    LLMProvider.MISTRAL: {
        "models": ["mistral-large-latest", "mistral-7b-instruct-v0.2"],
        "base_url": "https://api.mistral.ai/v1",
    },
    LLMProvider.SILICONFLOW: {
        "models": ["Qwen/Qwen2-72B-Instruct", "deepseek-ai/DeepSeek-V2.5"],
        "base_url": "https://api.siliconflow.cn/v1",
    },
    LLMProvider.DEEPSEEK: {
        "models": ["deepseek-chat", "deepseek-coder"],
        "base_url": "https://api.deepseek.com/v1",
    },
    LLMProvider.QWEN: {
        "models": ["qwen-turbo", "qwen-plus", "qwen-max"],
        "base_url": "https://dashscope.aliyuncs.com/compatible-mode/v1",
    },
    LLMProvider.OLLAMA: {
        "models": ["llama3", "llama3.1", "mistral", "codellama"],
        "base_url": os.getenv("OLLAMA_BASE_URL", "http://localhost:11434/v1"),
    },
}

class LLMClient:
    def __init__(self, provider: LLMProvider, api_key: str, model: str = None):
        self.provider = provider
        self.api_key = api_key
        self.config = PROVIDER_CONFIG[provider]
        self.model = model or self.config["models"][0]
        self.base_url = self.config["base_url"]
        self._client = httpx.AsyncClient(timeout=120.0)
    
    def get_models(self) -> List[str]:
        return self.config["models"]
    
    async def _get_headers(self) -> Dict[str, str]:
        headers = {"Content-Type": "application/json"}
        if self.provider == LLMProvider.ANTHROPIC:
            headers["x-api-key"] = self.api_key
            headers["anthropic-version"] = "2023-06-01"
        elif self.provider == LLMProvider.OLLAMA:
            pass
        else:
            headers["Authorization"] = f"Bearer {self.api_key}"
        return headers
    
    def _format_openai_messages(self, messages: List[Dict]) -> List[Dict]:
        formatted = []
        for msg in messages:
            role = msg.get("role", "user")
            if role == "system":
                role = "system"
            formatted.append({"role": role, "content": msg.get("content", "")})
        return formatted
    
    async def chat(
        self,
        messages: List[Dict],
        temperature: float = 0.7,
        max_tokens: int = 4096,
        stream: bool = False,
        **kwargs
    ) -> Dict | AsyncIterator:
        if self.provider == LLMProvider.ANTHROPIC:
            return await self._chat_anthropic(messages, temperature, max_tokens, stream, **kwargs)
        elif self.provider == LLMProvider.GEMINI:
            return await self._chat_gemini(messages, temperature, max_tokens, stream, **kwargs)
        else:
            return await self._chat_openai_compatible(messages, temperature, max_tokens, stream, **kwargs)
    
    async def _chat_anthropic(
        self,
        messages: List[Dict],
        temperature: float,
        max_tokens: int,
        stream: bool,
        **kwargs
    ) -> Dict | AsyncIterator:
        system_msg = ""
        chat_messages = []
        for msg in messages:
            if msg.get("role") == "system":
                system_msg = msg.get("content", "")
            else:
                chat_messages.append({
                    "role": msg.get("role", "user"),
                    "content": msg.get("content", "")
                })
        
        payload = {
            "model": self.model,
            "messages": chat_messages,
            "temperature": temperature,
            "max_tokens": max_tokens,
            "stream": stream,
        }
        if system_msg:
            payload["system"] = system_msg
        
        headers = await self._get_headers()
        
        if stream:
            async def stream_response():
                async with self._client.stream(
                    "POST",
                    f"{self.base_url}/messages",
                    headers=headers,
                    json=payload,
                ) as response:
                    async for line in response.aiter_lines():
                        if line.startswith("data: "):
                            data = line[6:]
                            if data == "[DONE]":
                                break
                            try:
                                yield json.loads(data)
                            except json.JSONDecodeError:
                                continue
            return stream_response()
        
        async with self._client as client:
            response = await client.post(
                f"{self.base_url}/messages",
                headers=headers,
                json=payload,
            )
            response.raise_for_status()
            result = response.json()
            return {
                "content": result.get("content", [{}])[0].get("text", ""),
                "model": result.get("model", self.model),
                "usage": result.get("usage", {}),
            }
    
    async def _chat_gemini(
        self,
        messages: List[Dict],
        temperature: float,
        max_tokens: int,
        stream: bool,
        **kwargs
    ) -> Dict | AsyncIterator:
        contents = []
        for msg in messages:
            if msg.get("role") == "user":
                contents.append({
                    "role": "user",
                    "parts": [{"text": msg.get("content", "")}]
                })
            elif msg.get("role") == "model":
                contents.append({
                    "role": "model",
                    "parts": [{"text": msg.get("content", "")}]
                })
        
        payload = {
            "contents": contents,
            "generationConfig": {
                "temperature": temperature,
                "maxOutputTokens": max_tokens,
            },
        }
        
        headers = await self._get_headers()
        
        if stream:
            async def stream_response():
                url = f"{self.base_url}/models/{self.model}:streamGenerateContent?key={self.api_key}&alt=sse"
                async with self._client.stream(
                    "POST",
                    url,
                    headers=headers,
                    json=payload,
                ) as response:
                    async for line in response.aiter_lines():
                        if line.startswith("data: "):
                            data = line[6:]
                            try:
                                yield json.loads(data)
                            except json.JSONDecodeError:
                                continue
            return stream_response()
        
        async with self._client as client:
            response = await client.post(
                f"{self.base_url}/models/{self.model}:generateContent?key={self.api_key}",
                headers=headers,
                json=payload,
            )
            response.raise_for_status()
            result = response.json()
            return {
                "content": result.get("candidates", [{}])[0].get("content", {}).get("parts", [{}])[0].get("text", ""),
       
