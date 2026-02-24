"""
Settings Dialog for Docking Studio
Allows configuration of Ollama/LM Studio and other settings
"""

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QComboBox, QLineEdit, QGroupBox,
    QFormLayout, QTabWidget, QWidget, QSpinBox,
    QCheckBox, QDialogButtonBox
)
from PyQt6.QtCore import Qt, pyqtSignal
import logging

logger = logging.getLogger(__name__)


class SettingsDialog(QDialog):
    """Settings configuration dialog"""
    
    settings_applied = pyqtSignal(dict)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.setWindowTitle("Settings")
        self.setMinimumSize(500, 400)
        self.setModal(True)
        
        self._setup_ui()
        self._load_current_settings()
    
    def _setup_ui(self):
        """Setup the settings UI"""
        layout = QVBoxLayout(self)
        
        # Tab widget
        tabs = QTabWidget()
        
        # LLM Settings tab
        llm_tab = self._create_llm_tab()
        tabs.addTab(llm_tab, "LLM Provider")
        
        # General Settings tab
        general_tab = self._create_general_tab()
        tabs.addTab(general_tab, "General")
        
        # Docker Settings tab
        docker_tab = self._create_docker_tab()
        tabs.addTab(docker_tab, "Docker")
        
        layout.addWidget(tabs)
        
        # Buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | 
            QDialogButtonBox.StandardButton.Cancel |
            QDialogButtonBox.StandardButton.Apply
        )
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(self._on_apply)
        
        layout.addWidget(buttons)
    
    def _create_llm_tab(self) -> QWidget:
        """Create LLM settings tab"""
        widget = QWidget()
        layout = QFormLayout(widget)
        
        # Provider selection
        self.provider_combo = QComboBox()
        self.provider_combo.addItems(["ollama", "lm_studio", "openai", "deepseek", "groq", "mistral"])
        self.provider_combo.setCurrentText("ollama")
        layout.addRow("Provider:", self.provider_combo)
        
        # Model selection
        self.model_edit = QLineEdit()
        self.model_edit.setPlaceholderText("e.g., llama3.2, codellama, mistral")
        self.model_edit.setText("llama3.2")
        layout.addRow("Model:", self.model_edit)
        
        # API Base URL
        self.api_base_edit = QLineEdit()
        self.api_base_edit.setPlaceholderText("http://localhost:11434")
        self.api_base_edit.setText("http://localhost:11434")
        layout.addRow("API Base URL:", self.api_base_edit)
        
        # API Key (optional)
        self.api_key_edit = QLineEdit()
        self.api_key_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self.api_key_edit.setPlaceholderText("Optional for local models")
        layout.addRow("API Key:", self.api_key_edit)
        
        # Temperature
        self.temperature_spin = QSpinBox()
        self.temperature_spin.setRange(0, 100)
        self.temperature_spin.setValue(70)
        self.temperature_spin.setSuffix(" %")
        layout.addRow("Temperature:", self.temperature_spin)
        
        # Max tokens
        self.max_tokens_spin = QSpinBox()
        self.max_tokens_spin.setRange(256, 16384)
        self.max_tokens_spin.setValue(4096)
        self.max_tokens_spin.setSingleStep(256)
        layout.addRow("Max Tokens:", self.max_tokens_spin)
        
        # Connect provider change to update URL
        self.provider_combo.currentTextChanged.connect(self._on_provider_changed)
        
        return widget
    
    def _create_general_tab(self) -> QWidget:
        """Create general settings tab"""
        widget = QWidget()
        layout = QFormLayout(widget)
        
        # Results directory
        self.results_dir_edit = QLineEdit()
        self.results_dir_edit.setText("results")
        layout.addRow("Results Directory:", self.results_dir_edit)
        
        # Models directory
        self.models_dir_edit = QLineEdit()
        self.models_dir_edit.setText("models")
        layout.addRow("Models Directory:", self.models_dir_edit)
        
        # Cache directory
        self.cache_dir_edit = QLineEdit()
        self.cache_dir_edit.setText(".cache")
        layout.addRow("Cache Directory:", self.cache_dir_edit)
        
        # Theme
        self.theme_combo = QComboBox()
        self.theme_combo.addItems(["light", "dark"])
        layout.addRow("Theme:", self.theme_combo)
        
        # Log level
        self.log_level_combo = QComboBox()
        self.log_level_combo.addItems(["DEBUG", "INFO", "WARNING", "ERROR"])
        self.log_level_combo.setCurrentText("INFO")
        layout.addRow("Log Level:", self.log_level_combo)
        
        return widget
    
    def _create_docker_tab(self) -> QWidget:
        """Create Docker settings tab"""
        widget = QWidget()
        layout = QFormLayout(widget)
        
        # Docker timeout
        self.docker_timeout_spin = QSpinBox()
        self.docker_timeout_spin.setRange(60, 7200)
        self.docker_timeout_spin.setValue(3600)
        self.docker_timeout_spin.setSuffix(" seconds")
        layout.addRow("Timeout:", self.docker_timeout_spin)
        
        # GPU enabled
        self.gpu_enabled_check = QCheckBox()
        self.gpu_enabled_check.setChecked(True)
        layout.addRow("GPU Enabled:", self.gpu_enabled_check)
        
        return widget
    
    def _on_provider_changed(self, provider: str):
        """Handle provider change to set default URL"""
        default_urls = {
            "ollama": "http://localhost:11434",
            "lm_studio": "http://localhost:1234/v1",
            "openai": "https://api.openai.com/v1",
            "deepseek": "https://api.deepseek.com/v1",
            "groq": "https://api.groq.com/openai/v1",
            "mistral": "https://api.mistral.ai/v1"
        }
        self.api_base_edit.setText(default_urls.get(provider, ""))
    
    def _load_current_settings(self):
        """Load current settings into the dialog"""
        try:
            from src.config_settings import Settings
            settings = Settings()
            app_settings = settings.get()
            
            # LLM settings
            self.provider_combo.setCurrentText(app_settings.chat_provider.provider)
            self.model_edit.setText(app_settings.chat_provider.model)
            self.api_base_edit.setText(app_settings.chat_provider.api_base)
            self.api_key_edit.setText(app_settings.chat_provider.api_key)
            self.temperature_spin.setValue(int(app_settings.chat_provider.temperature * 100))
            self.max_tokens_spin.setValue(app_settings.chat_provider.max_tokens)
            
            # General settings
            self.results_dir_edit.setText(app_settings.results_dir)
            self.models_dir_edit.setText(app_settings.models_dir)
            self.cache_dir_edit.setText(app_settings.cache_dir)
            self.theme_combo.setCurrentText(app_settings.theme)
            self.log_level_combo.setCurrentText(app_settings.log_level)
            
            # Docker settings
            self.docker_timeout_spin.setValue(app_settings.docker_timeout)
            self.gpu_enabled_check.setChecked(app_settings.gpu_enabled)
            
        except Exception as e:
            logger.warning(f"Could not load current settings: {e}")
    
    def _get_settings_dict(self) -> dict:
        """Get settings as dictionary"""
        return {
            "chat_provider": {
                "provider": self.provider_combo.currentText(),
                "model": self.model_edit.text(),
                "api_base": self.api_base_edit.text(),
                "api_key": self.api_key_edit.text(),
                "temperature": self.temperature_spin.value() / 100.0,
                "max_tokens": self.max_tokens_spin.value()
            },
            "results_dir": self.results_dir_edit.text(),
            "models_dir": self.models_dir_edit.text(),
            "cache_dir": self.cache_dir_edit.text(),
            "theme": self.theme_combo.currentText(),
            "log_level": self.log_level_combo.currentText(),
            "docker_timeout": self.docker_timeout_spin.value(),
            "gpu_enabled": self.gpu_enabled_check.isChecked()
        }
    
    def _on_apply(self):
        """Apply settings"""
        settings_dict = self._get_settings_dict()
        
        try:
            from src.config_settings import Settings, LLMProviderConfig
            settings = Settings()
            app_settings = settings.get()
            
            # Update LLM settings
            cp = settings_dict["chat_provider"]
            app_settings.chat_provider = LLMProviderConfig(
                provider=cp["provider"],
                model=cp["model"],
                api_base=cp["api_base"],
                api_key=cp["api_key"],
                temperature=cp["temperature"],
                max_tokens=cp["max_tokens"]
            )
            
            # Update other settings
            app_settings.results_dir = settings_dict["results_dir"]
            app_settings.models_dir = settings_dict["models_dir"]
            app_settings.cache_dir = settings_dict["cache_dir"]
            app_settings.theme = settings_dict["theme"]
            app_settings.log_level = settings_dict["log_level"]
            app_settings.docker_timeout = settings_dict["docker_timeout"]
            app_settings.gpu_enabled = settings_dict["gpu_enabled"]
            
            logger.info("Settings applied successfully")
            self.settings_applied.emit(settings_dict)
            
        except Exception as e:
            logger.error(f"Failed to apply settings: {e}")
    
    def _on_accept(self):
        """Handle OK button"""
        self._on_apply()
        self.accept()


def show_settings_dialog(parent=None) -> bool:
    """
    Show the settings dialog.
    
    Returns True if settings were applied, False if cancelled.
    """
    dialog = SettingsDialog(parent)
    return dialog.exec() == QDialog.DialogCode.Accepted
