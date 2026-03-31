from typing import Optional, Dict, List
import httpx
import os

class AIProvider:
    def __init__(self, api_key: Optional[str] = None, base_url: Optional[str] = None, model: str = ""):
        self.api_key = api_key or os.environ.get(f"{self.provider_name.upper()}_API_KEY", "")
        self.base_url = base_url
        self.model = model
    
    async def chat(self, messages: List[dict]) -> str:
        raise NotImplementedError
    
    def _get_headers(self) -> dict:
        return {"Content-Type": "application/json"}


class OpenAIProvider(AIProvider):
    provider_name = "openai"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "https://api.openai.com/v1") + "/chat/completions"
        headers = self._get_headers()
        headers["Authorization"] = f"Bearer {self.api_key}"
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "gpt-3.5-turbo",
                "messages": messages
            })
            response.raise_for_status()
            data = response.json()
            return data["choices"][0]["message"]["content"]


class ClaudeProvider(AIProvider):
    provider_name = "claude"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "https://api.anthropic.com/v1") + "/messages"
        headers = self._get_headers()
        headers["x-api-key"] = self.api_key
        headers["anthropic-version"] = "2023-06-01"
        
        content = "\n".join([f"{m['role']}: {m['content']}" for m in messages])
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "claude-3-sonnet-20240229",
                "messages": [{"role": "user", "content": content}]
            })
            response.raise_for_status()
            data = response.json()
            return data["content"][0]["text"]


class GeminiProvider(AIProvider):
    provider_name = "gemini"
    
    async def chat(self, messages: List[dict]) -> str:
        url = f"https://generativelanguage.googleapis.com/v1beta/models/{self.model or 'gemini-pro'}:generateContent"
        params = {"key": self.api_key}
        
        contents = [{"parts": [{"text": m["content"]}]} for m in messages if m["role"] != "system"]
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, params=params, json={"contents": contents})
            response.raise_for_status()
            data = response.json()
            return data["candidates"][0]["content"]["parts"][0]["text"]


class MistralProvider(AIProvider):
    provider_name = "mistral"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "https://api.mistral.ai/v1") + "/chat/completions"
        headers = self._get_headers()
        headers["Authorization"] = f"Bearer {self.api_key}"
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "mistral-medium",
                "messages": messages
            })
            response.raise_for_status()
            data = response.json()
            return data["choices"][0]["message"]["content"]


class DeepSeekProvider(AIProvider):
    provider_name = "deepseek"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "https://api.deepseek.com/v1") + "/chat/completions"
        headers = self._get_headers()
        headers["Authorization"] = f"Bearer {self.api_key}"
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "deepseek-chat",
                "messages": messages
            })
            response.raise_for_status()
            data = response.json()
            return data["choices"][0]["message"]["content"]


class QwenProvider(AIProvider):
    provider_name = "qwen"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "https://dashscope.aliyuncs.com/compatible-mode/v1") + "/chat/completions"
        headers = self._get_headers()
        headers["Authorization"] = f"Bearer {self.api_key}"
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "qwen-plus",
                "messages": messages
            })
            response.raise_for_status()
            data = response.json()
            return data["choices"][0]["message"]["content"]


class SiliconFlowProvider(AIProvider):
    provider_name = "siliconflow"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "https://api.siliconflow.cn/v1") + "/chat/completions"
        headers = self._get_headers()
        headers["Authorization"] = f"Bearer {self.api_key}"
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "Qwen/Qwen2-7B-Instruct",
                "messages": messages
            })
            response.raise_for_status()
            data = response.json()
            return data["choices"][0]["message"]["content"]


class OpenRouterProvider(AIProvider):
    provider_name = "openrouter"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "https://openrouter.ai/api/v1") + "/chat/completions"
        headers = self._get_headers()
        headers["Authorization"] = f"Bearer {self.api_key}"
        headers["HTTP-Referer"] = "https://biodockify.ai"
        headers["X-Title"] = "BioDockify"
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "anthropic/claude-3-sonnet",
                "messages": messages
            })
            response.raise_for_status()
            data = response.json()
            return data["choices"][0]["message"]["content"]


class OllamaProvider(AIProvider):
    provider_name = "ollama"
    
    async def chat(self, messages: List[dict]) -> str:
        url = (self.base_url or "http://localhost:11434") + "/api/chat"
        headers = self._get_headers()
        
        async with httpx.AsyncClient() as client:
            response = await client.post(url, headers=headers, json={
                "model": self.model or "llama2",
                "messages": messages
            })
            response.raise_for_status()
            data = response.json()
            return data["message"]["content"]


PROVIDERS = {
    "openai": OpenAIProvider,
    "claude": ClaudeProvider,
    "gemini": GeminiProvider,
    "mistral": MistralProvider,
    "deepseek": DeepSeekProvider,
    "qwen": QwenProvider,
    "siliconflow": SiliconFlowProvider,
    "openrouter": OpenRouterProvider,
    "ollama": OllamaProvider
}


async def chat_with_provider(provider: str, messages: List[dict], config: dict) -> str:
    provider_class = PROVIDERS.get(provider)
    if not provider_class:
        return f"Unknown provider: {provider}"
    
    try:
        provider_instance = provider_class(
            api_key=config.get("api_key"),
            base_url=config.get("base_url"),
            model=config.get("model", "")
        )
        return await provider_instance.chat(messages)
    except Exception as e:
        return f"Error: {str(e)}"
