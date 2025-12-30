"""
BioDockify Docking Studio - Configuration Management
Handles application settings, preferences, and configuration persistence
"""

import json
import os
from pathlib import Path
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)

class Config:
    """Configuration management for BioDockify Docking Studio"""
    
    def __init__(self):
        """Initialize configuration with default values"""
        self.config_dir = self._get_config_dir()
        self.config_file = self.config_dir / "config.json"
        self.defaults = self._get_defaults()
        self.config = self._load_config()
        
    def _get_config_dir(self) -> Path:
        """Get configuration directory"""
        if os.name == "nt":  # Windows
            config_dir = Path.home() / "AppData" / "Local" / "BioDockify"
        else:  # macOS/Linux
            config_dir = Path.home() / ".config" / "BioDockify"
        
        config_dir.mkdir(parents=True, exist_ok=True)
        return config_dir
    
    def _get_defaults(self) -> Dict[str, Any]:
        """Get default configuration values"""
        return {
            "exhaustiveness": 8,
            "num_modes": 9,
            "box_size": 20.0,
            "flexible_residues": [],
            "docker_memory": 4096,  # 4 GB in MB
            "max_jobs": 1,
            "log_level": "INFO",
            "auto_save": True,
            "checkpoints_enabled": True,
            "agent_zero_enabled": True,
            "ui_theme": "dark"
        }
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file"""
        if not self.config_file.exists():
            return self.defaults.copy()
        
        try:
            with open(self.config_file, 'r') as f:
                config = json.load(f)
            
            # Merge with defaults
            merged_config = self.defaults.copy()
            merged_config.update(config)
            
            return merged_config
            
        except Exception as e:
            logger.error(f"Failed to load config: {e}")
            return self.defaults.copy()
    
    def save(self) -> bool:
        """Save configuration to file"""
        try:
            with open(self.config_file, 'w') as f:
                json.dump(self.config, f, indent=2)
            logger.info(f"Configuration saved to {self.config_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to save config: {e}")
            return False
    
    def get(self, key: str, default: Optional[Any] = None) -> Any:
        """Get configuration value"""
        return self.config.get(key, default)
    
    def set(self, key: str, value: Any) -> None:
        """Set configuration value"""
        self.config[key] = value
        logger.debug(f"Config set: {key} = {value}")
    
    def reset_to_defaults(self) -> None:
        """Reset configuration to defaults"""
        self.config = self.defaults.copy()
        self.save()
        logger.info("Configuration reset to defaults")
