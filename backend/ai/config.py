"""
AI Configuration
Simplified to Ollama only with offline fallback
"""

import os

AI_MODE = os.getenv("AI_MODE", "auto")
# auto = auto-detect Ollama
# ollama = force Ollama
# offline = force offline fallback

OLLAMA_URL = os.getenv("OLLAMA_URL", "http://ollama:11434")
# Docker: http://ollama:11434
# Host: http://host.docker.internal:11434

OLLAMA_MODEL = os.getenv("OLLAMA_MODEL", "llama3")

ALLOW_AI = os.getenv("ALLOW_AI", "true").lower() == "true"

OLLAMA_TIMEOUT = int(os.getenv("OLLAMA_TIMEOUT", "30"))
