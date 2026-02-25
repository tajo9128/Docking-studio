"""
LLM Router
Auto-detects Ollama and falls back to offline assistant
"""

import requests
import logging
from typing import Dict, Optional

from .config import (
    OLLAMA_URL, 
    OLLAMA_MODEL, 
    AI_MODE, 
    ALLOW_AI,
    OLLAMA_TIMEOUT
)
from .offline_engine import OfflineAssistant

logger = logging.getLogger(__name__)


class OllamaProvider:
    """Ollama API provider"""
    
    def __init__(self, url: str = OLLAMA_URL, model: str = OLLAMA_MODEL):
        self.url = url
        self.model = model
    
    def is_available(self) -> bool:
        """Check if Ollama is running"""
        try:
            response = requests.get(f"{self.url}/api/tags", timeout=OLLAMA_TIMEOUT)
            return response.status_code == 200
        except Exception:
            return False
    
    def chat(self, prompt: str) -> str:
        """Send chat request to Ollama"""
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False
        }
        
        response = requests.post(
            f"{self.url}/api/generate",
            json=payload,
            timeout=OLLAMA_TIMEOUT * 2
        )
        
        response.raise_for_status()
        return response.json()["response"]
    
    def get_models(self) -> list:
        """Get available models"""
        try:
            response = requests.get(f"{self.url}/api/tags", timeout=OLLAMA_TIMEOUT)
            if response.status_code == 200:
                data = response.json()
                return [m["name"] for m in data.get("models", [])]
        except Exception as e:
            logger.warning(f"Could not fetch Ollama models: {e}")
        return []


class LLMRouter:
    """
    Intelligent LLM router that:
    - Auto-detects Ollama
    - Falls back to offline assistant on failure
    - Never crashes
    """
    
    def __init__(self):
        self.ollama = OllamaProvider()
        self.offline = OfflineAssistant()
        self._provider = None
        logger.info("LLMRouter initialized")
    
    @property
    def provider(self) -> str:
        """Get current provider"""
        if self._provider is None:
            self._provider = self._detect_provider()
        return self._provider
    
    def _detect_provider(self) -> str:
        """Auto-detect the best available provider"""
        if not ALLOW_AI:
            logger.info("AI disabled via config")
            return "offline"
        
        if AI_MODE == "offline":
            logger.info("Forced offline mode")
            return "offline"
        
        if AI_MODE == "ollama":
            if self.ollama.is_available():
                logger.info("Ollama available, using Ollama")
                return "ollama"
            logger.warning("Ollama forced but not available")
            return "offline"
        
        if AI_MODE == "auto":
            if self.ollama.is_available():
                logger.info("Auto-detected Ollama")
                return "ollama"
            logger.info("Ollama not available, using offline")
            return "offline"
        
        return "offline"
    
    def detect_ollama(self) -> bool:
        """Check if Ollama is available"""
        return self.ollama.is_available()
    
    def get_available_models(self) -> list:
        """Get available Ollama models"""
        if self.provider == "ollama":
            return self.ollama.get_models()
        return []
    
    def chat(self, message: str) -> Dict:
        """
        Send a chat message.
        
        Returns dict with:
        - response: the response text
        - provider: which provider was used
        - available: whether AI is available
        """
        if self.provider == "ollama":
            try:
                response_text = self.ollama.chat(message)
                return {
                    "response": response_text,
                    "provider": "ollama",
                    "available": True
                }
            except Exception as e:
                logger.warning(f"Ollama failed: {e}, falling back to offline")
                response_text = self.offline.respond(message)
                return {
                    "response": response_text,
                    "provider": "offline",
                    "available": False,
                    "error": str(e)
                }
        
        response_text = self.offline.respond(message)
        return {
            "response": response_text,
            "provider": "offline",
            "available": self.detect_ollama()
        }
    
    def reset(self):
        """Reset provider cache"""
        self._provider = None


def get_router() -> LLMRouter:
    """Get singleton router instance"""
    if not hasattr(get_router, "_instance"):
        get_router._instance = LLMRouter()
    return get_router._instance


if __name__ == "__main__":
    router = LLMRouter()
    
    print(f"Provider: {router.provider}")
    print(f"Ollama available: {router.detect_ollama()}")
    
    test_messages = [
        "What is Vina scoring?",
        "How does GNINA work?",
        "What is consensus scoring?",
        "Tell me about hydrogen bonds"
    ]
    
    for msg in test_messages:
        print(f"\nQ: {msg}")
        result = router.chat(msg)
        print(f"A: {result['response'][:100]}...")
        print(f"Provider: {result['provider']}")
