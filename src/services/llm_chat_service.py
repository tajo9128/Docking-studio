"""
LLM Chat Service
Ollama/LM Studio integration for Docking Studio
Based on agent-zero's LiteLLM implementation
"""

import os
import logging
from typing import List, Dict, Optional, Any, AsyncIterator, Callable
from dataclasses import dataclass, field
from enum import Enum

logger = logging.getLogger(__name__)

# Try to import LiteLLM
LITELLM_AVAILABLE = False
try:
    from litellm import completion, acompletion, embedding
    import litellm
    LITELLM_AVAILABLE = True
except ImportError:
    logger.warning("LiteLLM not available, using fallback")

# Try to import langchain
LANGCHAIN_AVAILABLE = False
try:
    from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
    from langchain_core.chat_models import SimpleChatModel
    LANGCHAIN_AVAILABLE = True
except ImportError:
    logger.warning("LangChain not available")


class ModelProvider(Enum):
    """Supported model providers"""
    OLLAMA = "ollama"
    LM_STUDIO = "lm_studio"
    OPENAI = "openai"
    ANTHROPIC = "anthropic"
    DEEPSEEK = "deepseek"
    GROQ = "groq"
    MISTRAL = "mistral"
    HUGGINGFACE = "huggingface"
    OTHER = "other"


@dataclass
class ChatMessage:
    """A chat message"""
    role: str  # 'user', 'assistant', 'system'
    content: str
    

@dataclass
class ChatConfig:
    """Chat model configuration"""
    provider: ModelProvider = ModelProvider.OLLAMA
    model: str = "llama3.2"
    api_base: str = "http://localhost:11434"
    api_key: str = ""
    temperature: float = 0.7
    max_tokens: int = 4096
    streaming: bool = True
    
    def get_litellm_model(self) -> str:
        """Get model string for LiteLLM"""
        if self.provider == ModelProvider.OLLAMA:
            return f"ollama/{self.model}"
        elif self.provider == ModelProvider.LM_STUDIO:
            return f"lm_studio/{self.model}"
        elif self.provider == ModelProvider.OPENAI:
            return f"openai/{self.model}"
        elif self.provider == ModelProvider.DEEPSEEK:
            return f"deepseek/{self.model}"
        elif self.provider == ModelProvider.GROQ:
            return f"groq/{self.model}"
        elif self.provider == ModelProvider.MISTRAL:
            return f"mistral/{self.model}"
        elif self.provider == ModelProvider.HUGGINGFACE:
            return f"huggingface/{self.model}"
        else:
            return self.model


class ChatService:
    """
    LLM Chat Service supporting Ollama, LM Studio, and other providers.
    """
    
    def __init__(self, config: Optional[ChatConfig] = None):
        """
        Initialize chat service.
        
        Args:
            config: Chat configuration, uses defaults if not provided
        """
        self.config = config or ChatConfig()
        self.history: List[ChatMessage] = []
        
        if not LITELLM_AVAILABLE:
            logger.warning("LiteLLM not available. Install with: pip install litellm")
        
        logger.info(f"ChatService initialized with provider: {self.config.provider.value}")
    
    def configure(self, config: ChatConfig):
        """Update chat configuration"""
        self.config = config
        logger.info(f"ChatService reconfigured: {config.provider.value}/{config.model}")
    
    def add_message(self, role: str, content: str):
        """Add a message to history"""
        self.history.append(ChatMessage(role=role, content=content))
    
    def clear_history(self):
        """Clear chat history"""
        self.history = []
    
    def set_system_prompt(self, prompt: str):
        """Set system prompt"""
        # Remove existing system messages
        self.history = [m for m in self.history if m.role != 'system']
        # Add new system message at the beginning
        self.history.insert(0, ChatMessage(role='system', content=prompt))
    
    def _build_messages(self) -> List[Dict[str, str]]:
        """Build messages list for LiteLLM"""
        return [{"role": m.role, "content": m.content} for m in self.history]
    
    def chat(
        self,
        message: str,
        system_prompt: Optional[str] = None,
        stream: bool = False,
        callback: Optional[Callable[[str], None]] = None
    ) -> str:
        """
        Send a chat message and get response.
        
        Args:
            message: User message
            system_prompt: Optional system prompt
            stream: Enable streaming
            callback: Optional callback for streaming
            
        Returns:
            Assistant response
        """
        # Add user message to history
        self.add_message('user', message)
        
        # Build messages
        messages = self._build_messages()
        
        if system_prompt and not any(m['role'] == 'system' for m in messages):
            messages.insert(0, {"role": "system", "content": system_prompt})
        
        if not LITELLM_AVAILABLE:
            return "LiteLLM not available. Please install: pip install litellm"
        
        try:
            model = self.config.get_litellm_model()
            
            kwargs = {
                "model": model,
                "messages": messages,
                "temperature": self.config.temperature,
                "max_tokens": self.config.max_tokens,
            }
            
            # Add API base for local models
            if self.config.provider in [ModelProvider.OLLAMA, ModelProvider.LM_STUDIO]:
                kwargs["api_base"] = self.config.api_base
            
            if self.config.api_key:
                kwargs["api_key"] = self.config.api_key
            
            if stream and callback:
                # Streaming mode
                response = ""
                for chunk in completion(**kwargs):
                    if chunk.choices and len(chunk.choices) > 0:
                        delta = chunk.choices[0].delta.content or ""
                        response += delta
                        callback(delta)
                self.add_message('assistant', response)
                return response
            else:
                # Non-streaming mode
                response = completion(**kwargs)
                if response.choices and len(response.choices) > 0:
                    content = response.choices[0].message.content
                    self.add_message('assistant', content)
                    return content
                return ""
                
        except Exception as e:
            logger.error(f"Chat error: {e}")
            return f"Error: {str(e)}"
    
    async def achat(
        self,
        message: str,
        system_prompt: Optional[str] = None,
        callback: Optional[Callable[[str], None]] = None
    ) -> str:
        """Async version of chat"""
        # Add user message to history
        self.add_message('user', message)
        
        # Build messages
        messages = self._build_messages()
        
        if system_prompt and not any(m['role'] == 'system' for m in messages):
            messages.insert(0, {"role": "system", "content": system_prompt})
        
        if not LITELLM_AVAILABLE:
            return "LiteLLM not available. Please install: pip install litellm"
        
        try:
            model = self.config.get_litellm_model()
            
            kwargs = {
                "model": model,
                "messages": messages,
                "temperature": self.config.temperature,
                "max_tokens": self.config.max_tokens,
            }
            
            if self.config.provider in [ModelProvider.OLLAMA, ModelProvider.LM_STUDIO]:
                kwargs["api_base"] = self.config.api_base
            
            if self.config.api_key:
                kwargs["api_key"] = self.config.api_key
            
            response_text = ""
            
            async for chunk in await acompletion(**kwargs):
                if chunk.choices and len(chunk.choices) > 0:
                    delta = chunk.choices[0].delta.content or ""
                    response_text += delta
                    if callback:
                        callback(delta)
            
            self.add_message('assistant', response_text)
            return response_text
            
        except Exception as e:
            logger.error(f"Async chat error: {e}")
            return f"Error: {str(e)}"
    
    def get_available_models(self) -> List[str]:
        """Get available models from the provider"""
        if self.config.provider == ModelProvider.OLLAMA:
            return self._get_ollama_models()
        elif self.config.provider == ModelProvider.LM_STUDIO:
            return self._get_lmstudio_models()
        return []
    
    def _get_ollama_models(self) -> List[str]:
        """Get models from local Ollama instance"""
        try:
            import requests
            response = requests.get(f"{self.config.api_base}/api/tags", timeout=5)
            if response.status_code == 200:
                data = response.json()
                return [m['name'] for m in data.get('models', [])]
        except Exception as e:
            logger.warning(f"Failed to get Ollama models: {e}")
        return []
    
    def _get_lmstudio_models(self) -> List[str]:
        """Get models from LM Studio"""
        try:
            import requests
            response = requests.get(f"{self.config.api_base}/api/v0/models", timeout=5)
            if response.status_code == 200:
                data = response.json()
                return [m['id'] for m in data.get('data', [])]
        except Exception as e:
            logger.warning(f"Failed to get LM Studio models: {e}")
        return []
    
    def check_connection(self) -> bool:
        """Check if the LLM provider is available"""
        try:
            if self.config.provider == ModelProvider.OLLAMA:
                import requests
                response = requests.get(f"{self.config.api_base}/api/tags", timeout=5)
                return response.status_code == 200
            elif self.config.provider == ModelProvider.LM_STUDIO:
                import requests
                response = requests.get(f"{self.config.api_base}/api/v0/models", timeout=5)
                return response.status_code == 200
            else:
                # For cloud providers, try a simple completion
                return bool(self.chat("test", system_prompt="Reply with 'ok'"))
        except Exception as e:
            logger.warning(f"Connection check failed: {e}")
            return False


class EmbeddingService:
    """Embedding service for Ollama/LM Studio"""
    
    def __init__(
        self,
        provider: ModelProvider = ModelProvider.OLLAMA,
        model: str = "nomic-embed-text",
        api_base: str = "http://localhost:11434"
    ):
        self.provider = provider
        self.model = model
        self.api_base = api_base
    
    def embed_text(self, text: str) -> List[float]:
        """Get embedding for text"""
        if not LITELLM_AVAILABLE:
            logger.warning("LiteLLM not available for embeddings")
            return []
        
        try:
            if self.provider == ModelProvider.OLLAMA:
                response = embedding(
                    model=f"ollama/{self.model}",
                    input=[text],
                    api_base=self.api_base
                )
            elif self.provider == ModelProvider.LM_STUDIO:
                response = embedding(
                    model=f"lm_studio/{self.model}",
                    input=[text],
                    api_base=self.api_base
                )
            else:
                response = embedding(model=self.model, input=[text])
            
            if response.data and len(response.data) > 0:
                return response.data[0].embedding
        except Exception as e:
            logger.error(f"Embedding error: {e}")
        
        return []
    
    def embed_batch(self, texts: List[str]) -> List[List[float]]:
        """Get embeddings for multiple texts"""
        if not LITELLM_AVAILABLE:
            return [[] for _ in texts]
        
        try:
            if self.provider == ModelProvider.OLLAMA:
                response = embedding(
                    model=f"ollama/{self.model}",
                    input=texts,
                    api_base=self.api_base
                )
            elif self.provider == ModelProvider.LM_STUDIO:
                response = embedding(
                    model=f"lm_studio/{self.model}",
                    input=texts,
                    api_base=self.api_base
                )
            else:
                response = embedding(model=self.model, input=texts)
            
            return [item.embedding for item in response.data]
        except Exception as e:
            logger.error(f"Batch embedding error: {e}")
            return [[] for _ in texts]


# Default chat service instance
_default_chat_service: Optional[ChatService] = None


def get_chat_service(config: Optional[ChatConfig] = None) -> ChatService:
    """Get or create default chat service"""
    global _default_chat_service
    if config:
        _default_chat_service = ChatService(config)
    elif _default_chat_service is None:
        _default_chat_service = ChatService()
    return _default_chat_service


def create_ollama_service(
    model: str = "llama3.2",
    api_base: str = "http://localhost:11434"
) -> ChatService:
    """Create chat service for Ollama"""
    config = ChatConfig(
        provider=ModelProvider.OLLAMA,
        model=model,
        api_base=api_base
    )
    return ChatService(config)


def create_lmstudio_service(
    model: str = "llama-3.2-1b-instruct-q4_k_m",
    api_base: str = "http://localhost:1234/v1"
) -> ChatService:
    """Create chat service for LM Studio"""
    config = ChatConfig(
        provider=ModelProvider.LM_STUDIO,
        model=model,
        api_base=api_base,
        api_key="not-needed"
    )
    return ChatService(config)
