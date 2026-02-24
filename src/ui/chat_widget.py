"""
Chat Widget for Docking Studio
UI for chatting with Ollama/LM Studio
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QTextEdit, QLineEdit, QLabel, QComboBox, QFrame,
    QScrollArea, QSizePolicy
)
from PyQt6.QtCore import Qt, pyqtSignal, QThread, QTimer
from PyQt6.QtGui import QFont, QTextCursor
import logging

logger = logging.getLogger(__name__)


class ChatWorker(QThread):
    """Worker thread for LLM chat"""
    message_ready = pyqtSignal(str)
    finished = pyqtSignal()
    error = pyqtSignal(str)
    
    def __init__(self, chat_service, message, system_prompt=""):
        super().__init__()
        self.chat_service = chat_service
        self.message = message
        self.system_prompt = system_prompt
    
    def run(self):
        try:
            response = self.chat_service.chat(
                self.message,
                system_prompt=self.system_prompt,
                stream=False
            )
            self.message_ready.emit(response)
        except Exception as e:
            self.error.emit(str(e))
        finally:
            self.finished.emit()


class ChatWidget(QWidget):
    """
    Chat widget for interacting with Ollama/LM Studio.
    """
    
    message_sent = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.chat_service = None
        self.worker = None
        self.is_thinking = False
        
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the chat UI"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Header
        header = self._create_header()
        layout.addWidget(header)
        
        # Chat display area
        self.chat_display = QTextEdit()
        self.chat_display.setReadOnly(True)
        self.chat_display.setFont(QFont("Segoe UI", 10))
        self.chat_display.setStyleSheet("""
            QTextEdit {
                background: #fafafa;
                border: none;
                padding: 12px;
            }
        """)
        layout.addWidget(self.chat_display, 1)
        
        # Input area
        input_frame = self._create_input_area()
        layout.addWidget(input_frame)
    
    def _create_header(self) -> QFrame:
        """Create header with provider selection"""
        header = QFrame()
        header.setStyleSheet("""
            QFrame {
                background: white;
                border-bottom: 1px solid #e5e7eb;
                padding: 12px;
            }
        """)
        layout = QHBoxLayout(header)
        layout.setContentsMargins(12, 8, 12, 8)
        
        title = QLabel("AI Assistant")
        title.setStyleSheet("""
            QLabel {
                font-size: 14px;
                font-weight: 600;
                color: #1f2937;
            }
        """)
        layout.addWidget(title)
        
        layout.addStretch()
        
        # Provider selector
        layout.addWidget(QLabel("Provider:"))
        
        self.provider_combo = QComboBox()
        self.provider_combo.addItems(["Ollama", "LM Studio", "OpenAI"])
        self.provider_combo.setFixedWidth(120)
        self.provider_combo.currentTextChanged.connect(self._on_provider_changed)
        layout.addWidget(self.provider_combo)
        
        # Model selector
        self.model_combo = QComboBox()
        self.model_combo.setFixedWidth(150)
        self.model_combo.addItems(["llama3.2", "qwen2.5", "mistral"])
        layout.addWidget(self.model_combo)
        
        # Connect button
        self.connect_btn = QPushButton("Connect")
        self.connect_btn.setFixedWidth(80)
        self.connect_btn.setStyleSheet("""
            QPushButton {
                background: #2E5AAC;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-weight: 600;
            }
            QPushButton:hover {
                background: #254A8F;
            }
        """)
        self.connect_btn.clicked.connect(self._on_connect_clicked)
        layout.addWidget(self.connect_btn)
        
        # Clear button
        clear_btn = QPushButton("Clear")
        clear_btn.setFixedWidth(60)
        clear_btn.setStyleSheet("""
            QPushButton {
                background: #f3f4f6;
                color: #6b7280;
                border: 1px solid #e5e7eb;
                border-radius: 6px;
                padding: 6px 12px;
            }
            QPushButton:hover {
                background: #e5e7eb;
            }
        """)
        clear_btn.clicked.connect(self._clear_chat)
        layout.addWidget(clear_btn)
        
        return header
    
    def _create_input_area(self) -> QFrame:
        """Create message input area"""
        frame = QFrame()
        frame.setStyleSheet("""
            QFrame {
                background: white;
                border-top: 1px solid #e5e7eb;
                padding: 12px;
            }
        """)
        layout = QHBoxLayout(frame)
        layout.setContentsMargins(12, 8, 12, 8)
        
        # Input field
        self.input_field = QLineEdit()
        self.input_field.setPlaceholderText("Ask about your docking results...")
        self.input_field.setFont(QFont("Segoe UI", 10))
        self.input_field.setStyleSheet("""
            QLineEdit {
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                padding: 10px 14px;
                background: #f9fafb;
            }
            QLineEdit:focus {
                border-color: #2E5AAC;
                background: white;
            }
        """)
        self.input_field.returnPressed.connect(self._send_message)
        layout.addWidget(self.input_field, 1)
        
        # Send button
        self.send_btn = QPushButton("Send")
        self.send_btn.setFixedWidth(80)
        self.send_btn.setStyleSheet("""
            QPushButton {
                background: #2E5AAC;
                color: white;
                border: none;
                border-radius: 8px;
                padding: 10px 16px;
                font-weight: 600;
            }
            QPushButton:hover {
                background: #254A8F;
            }
            QPushButton:disabled {
                background: #9ca3af;
            }
        """)
        self.send_btn.clicked.connect(self._send_message)
        layout.addWidget(self.send_btn)
        
        return frame
    
    def set_chat_service(self, service):
        """Set the chat service"""
        self.chat_service = service
        self._update_status("Ready")
    
    def _on_provider_changed(self, provider: str):
        """Handle provider change"""
        if provider == "Ollama":
            self.model_combo.clear()
            self.model_combo.addItems(["llama3.2", "qwen2.5", "mistral", "phi3", "codellama"])
            if self.chat_service:
                from src.services.llm_chat_service import ModelProvider
                self.chat_service.config.provider = ModelProvider.OLLAMA
                self.chat_service.config.api_base = "http://localhost:11434"
        elif provider == "LM Studio":
            self.model_combo.clear()
            self.model_combo.addItems(["llama-3.2-1b-instruct-q4_k_m", "phi3-mini-4k-instruct-q4", "mistral-7b-instruct-v0.2"])
            if self.chat_service:
                from src.services.llm_chat_service import ModelProvider
                self.chat_service.config.provider = ModelProvider.LM_STUDIO
                self.chat_service.config.api_base = "http://localhost:1234/v1"
        elif provider == "OpenAI":
            self.model_combo.clear()
            self.model_combo.addItems(["gpt-4o", "gpt-4-turbo", "gpt-3.5-turbo"])
            if self.chat_service:
                from src.services.llm_chat_service import ModelProvider
                self.chat_service.config.provider = ModelProvider.OPENAI
                self.chat_service.config.api_base = "https://api.openai.com/v1"
    
    def _on_connect_clicked(self):
        """Handle connect button click"""
        self._update_status("Connecting...")
        
        from src.services.llm_chat_service import (
            get_chat_service, create_ollama_service, create_lmstudio_service,
            ChatConfig, ModelProvider
        )
        
        provider = self.provider_combo.currentText()
        model = self.model_combo.currentText()
        
        if provider == "Ollama":
            self.chat_service = create_ollama_service(model=model)
        elif provider == "LM Studio":
            self.chat_service = create_lmstudio_service(model=model)
        else:
            config = ChatConfig(
                provider=ModelProvider.OPENAI,
                model=model,
                api_base="https://api.openai.com/v1"
            )
            self.chat_service = get_chat_service(config)
        
        # Try to connect
        if self.chat_service.check_connection():
            self._update_status(f"Connected to {provider}")
            self._add_message("system", f"Connected to {provider} with model {model}")
        else:
            self._update_status("Connection failed")
            self._add_message("system", f"Could not connect to {provider}. Make sure it's running.")
    
    def _send_message(self):
        """Send a message"""
        if self.is_thinking or not self.input_field.text().strip():
            return
        
        message = self.input_field.text().strip()
        self.input_field.clear()
        
        # Add user message
        self._add_message("user", message)
        
        # Show thinking indicator
        self._show_thinking(True)
        self.send_btn.setEnabled(False)
        
        # Get system prompt based on context
        system_prompt = self._get_system_prompt()
        
        # Start worker thread
        self.worker = ChatWorker(self.chat_service, message, system_prompt)
        self.worker.message_ready.connect(self._on_response)
        self.worker.error.connect(self._on_error)
        self.worker.finished.connect(self._on_finished)
        self.worker.start()
    
    def _get_system_prompt(self) -> str:
        """Get context-aware system prompt"""
        # Could include recent docking results context
        return """You are an AI assistant for BioDockify Docking Studio.
Be helpful, concise, and accurate. You can analyze molecular docking results,
explain binding interactions, and suggest optimizations."""
    
    def _on_response(self, response: str):
        """Handle response received"""
        self._add_message("assistant", response)
    
    def _on_error(self, error: str):
        """Handle error"""
        self._add_message("system", f"Error: {error}")
    
    def _on_finished(self):
        """Handle worker finished"""
        self._show_thinking(False)
        self.send_btn.setEnabled(True)
        self.worker = None
    
    def _add_message(self, role: str, content: str):
        """Add a message to the chat display"""
        if role == "user":
            color = "#2E5AAC"
            align = "right"
            bg = "#eff6ff"
        elif role == "assistant":
            color = "#1f2937"
            align = "left"
            bg = "#f3f4f6"
        else:  # system
            color = "#6b7280"
            align = "center"
            bg = "#fef3c7"
        
        html = f"""
        <div style="margin: 8px 0; text-align: {align};">
            <div style="display: inline-block; max-width: 80%; padding: 10px 14px; 
                        border-radius: 12px; background: {bg}; color: {color};">
                {content.replace('<', '&lt;').replace('>', '&gt;').replace('\n', '<br>')}
            </div>
        </div>
        """
        
        self.chat_display.append(html)
        self.chat_display.moveCursor(QTextCursor.MoveOperation.End)
    
    def _show_thinking(self, show: bool):
        """Show thinking indicator"""
        self.is_thinking = show
        if show:
            self.chat_display.append("""
                <div style="text-align: left; margin: 8px 0;">
                    <span style="color: #6b7280; font-style: italic;">Thinking...</span>
                </div>
            """)
        else:
            # Remove thinking indicator
            cursor = self.chat_display.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            cursor.select(QTextCursor.SelectionType.LineUnderCursor)
            cursor.removeSelectedText()
            cursor.deleteChar()
    
    def _update_status(self, status: str):
        """Update connection status"""
        self.connect_btn.setText(status)
        logger.info(f"Chat status: {status}")
    
    def _clear_chat(self):
        """Clear chat history"""
        self.chat_display.clear()
        if self.chat_service:
            self.chat_service.clear_history()
    
    def add_context(self, context: str):
        """Add docking context to chat"""
        if self.chat_service:
            prompt = f"Context: {context}\n\nUser question:"
            # This would be used in subsequent messages
            logger.info(f"Added context: {context[:100]}...")
