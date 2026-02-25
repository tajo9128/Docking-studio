"""
AI Module
Simplified to Ollama only with offline fallback
"""

from .config import AI_MODE, OLLAMA_URL, OLLAMA_MODEL, ALLOW_AI
from .offline_engine import OfflineAssistant
from .llm_router import LLMRouter, OllamaProvider, get_router

__all__ = [
    'AI_MODE',
    'OLLAMA_URL', 
    'OLLAMA_MODEL', 
    'ALLOW_AI',
    'OfflineAssistant',
    'LLMRouter',
    'OllamaProvider',
    'get_router'
]
