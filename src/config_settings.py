"""
Settings Configuration for Docking Studio
Similar to agent-zero's settings system
"""

import os
import json
from typing import Dict, Any, Optional
from dataclasses import dataclass, field, asdict


@dataclass
class LLMProviderConfig:
    """Configuration for an LLM provider"""
    provider: str = "ollama"
    model: str = "llama3.2"
    api_base: str = "http://localhost:11434"
    api_key: str = ""
    temperature: float = 0.7
    max_tokens: int = 4096


@dataclass
class AppSettings:
    """Application settings"""
    # LLM Settings
    chat_provider: LLMProviderConfig = field(default_factory=LLMProviderConfig)
    
    # Paths
    results_dir: str = "results"
    models_dir: str = "models"
    cache_dir: str = ".cache"
    
    # Docker
    docker_timeout: int = 3600
    gpu_enabled: bool = True
    
    # UI
    theme: str = "light"
    language: str = "en"
    
    # Logging
    log_level: str = "INFO"
    log_file: str = "docking_studio.log"


class Settings:
    """
    Settings manager for Docking Studio.
    Loads from environment variables and config file.
    """
    
    _instance = None
    _settings: Optional[AppSettings] = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        if self._settings is None:
            self._settings = self._load_settings()
    
    def _load_settings(self) -> AppSettings:
        """Load settings from environment and config"""
        settings = AppSettings()
        
        # Load from environment variables
        chat_model = os.environ.get("DOCKING_CHAT_MODEL", "llama3.2")
        chat_provider = os.environ.get("DOCKING_CHAT_PROVIDER", "ollama")
        api_base = os.environ.get("DOCKING_API_BASE", "http://localhost:11434")
        
        if chat_provider == "lm_studio":
            api_base = os.environ.get("DOCKING_API_BASE", "http://localhost:1234/v1")
        
        settings.chat_provider = LLMProviderConfig(
            provider=chat_provider,
            model=chat_model,
            api_base=api_base,
            api_key=os.environ.get("DOCKING_API_KEY", ""),
            temperature=float(os.environ.get("DOCKING_TEMPERATURE", "0.7")),
            max_tokens=int(os.environ.get("DOCKING_MAX_TOKENS", "4096")),
        )
        
        # Other settings
        settings.results_dir = os.environ.get("DOCKING_RESULTS_DIR", "results")
        settings.models_dir = os.environ.get("DOCKING_MODELS_DIR", "models")
        settings.docker_timeout = int(os.environ.get("DOCKING_DOCKER_TIMEOUT", "3600"))
        settings.gpu_enabled = os.environ.get("DOCKING_GPU_ENABLED", "true").lower() == "true"
        settings.theme = os.environ.get("DOCKING_THEME", "light")
        settings.log_level = os.environ.get("DOCKING_LOG_LEVEL", "INFO")
        
        return settings
    
    def get(self) -> AppSettings:
        """Get current settings"""
        return self._settings
    
    def update(self, **kwargs):
        """Update settings"""
        for key, value in kwargs.items():
            if hasattr(self._settings, key):
                setattr(self._settings, key, value)
    
    def save(self, path: str = "settings.json"):
        """Save settings to file"""
        with open(path, 'w') as f:
            json.dump(asdict(self._settings), f, indent=2)
    
    def load(self, path: str = "settings.json"):
        """Load settings from file"""
        if os.path.exists(path):
            with open(path, 'r') as f:
                data = json.load(f)
                for key, value in data.items():
                    if hasattr(self._settings, key):
                        setattr(self._settings, key, value)
    
    @classmethod
    def reset(cls):
        """Reset settings to defaults"""
        cls._settings = None
        cls._instance = None


# Convenience functions
def get_settings() -> AppSettings:
    """Get application settings"""
    return Settings().get()


def update_settings(**kwargs):
    """Update settings"""
    Settings().update(**kwargs)


def save_settings(path: str = "settings.json"):
    """Save settings"""
    Settings().save(path)


def load_settings(path: str = "settings.json"):
    """Load settings"""
    Settings().load(path)


# Default prompts for docking assistant
DOCKING_SYSTEM_PROMPT = """You are an AI assistant for BioDockify Docking Studio, a molecular docking tool.

You can help users with:
- Explaining molecular docking concepts
- Analyzing docking results
- Suggesting binding site optimizations
- Interpreting interaction patterns
- Recommending compound modifications
- Answering questions about drug discovery

Be concise, accurate, and helpful. Use scientific terminology when appropriate."""


DOCKING_ANALYSIS_PROMPT = """You are analyzing molecular docking results.

For the current docking results:
- Consider the binding affinity scores
- Evaluate hydrogen bonds and hydrophobic interactions
- Assess the ligand's fit in the binding pocket
- Note any polar or charged interactions
- Consider drug-likeness properties

Provide insights about:
1. Binding mode quality
2. Key interactions
3. Potential for optimization
4. Recommendations for next steps"""
