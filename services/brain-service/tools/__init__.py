"""
Docking Studio v2.0 - Nanobot Brain Tool System
Plugin-based architecture for drug discovery tools
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field
import logging

logger = logging.getLogger(__name__)


class ToolInput(BaseModel):
    """Base input schema for all tools"""

    class Config:
        extra = "allow"


class ToolOutput(BaseModel):
    """Base output schema for all tools"""

    success: bool = True
    data: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class BaseTool(ABC):
    """Abstract base class for all Nanobot tools"""

    name: str = ""
    description: str = ""
    category: str = "general"
    input_schema: Dict[str, Any] = {}
    output_schema: Dict[str, Any] = {}

    def __init__(self):
        self.logger = logging.getLogger(f"tools.{self.name}")

    @abstractmethod
    async def execute(self, input_data: ToolInput) -> ToolOutput:
        """Execute the tool with given input"""
        pass

    def validate_input(self, input_data: Dict[str, Any]) -> ToolInput:
        """Validate and convert input to ToolInput"""
        return ToolInput(**input_data)

    def to_definition(self) -> Dict[str, Any]:
        """Return tool definition in OpenAI function-calling format"""
        return {
            "type": "function",
            "function": {
                "name": self.name,
                "description": self.description,
                "parameters": self.input_schema if self.input_schema else {
                    "type": "object",
                    "properties": {},
                },
            },
        }


class ToolRegistry:
    """Dynamic registry for Nanobot tools"""

    _instance = None
    _tools: Dict[str, BaseTool] = {}

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._tools = {}
            cls._instance.logger = logging.getLogger("tool_registry")
        return cls._instance

    def register(self, tool: BaseTool) -> None:
        """Register a tool"""
        self._tools[tool.name] = tool
        self.logger.info(f"Registered tool: {tool.name}")

    def unregister(self, name: str) -> None:
        """Unregister a tool"""
        if name in self._tools:
            del self._tools[name]
            self.logger.info(f"Unregistered tool: {name}")

    def get(self, name: str) -> Optional[BaseTool]:
        """Get a tool by name"""
        return self._tools.get(name)

    def list_tools(self) -> List[Dict[str, Any]]:
        """List all registered tools"""
        return [tool.to_definition() for tool in self._tools.values()]

    def list_by_category(self, category: str) -> List[BaseTool]:
        """List tools by category"""
        return [t for t in self._tools.values() if t.category == category]

    async def execute(self, tool_name: str, input_data: Dict[str, Any]) -> ToolOutput:
        """Execute a tool by name"""
        tool = self.get(tool_name)
        if not tool:
            return ToolOutput(success=False, error=f"Tool not found: {tool_name}")

        try:
            validated_input = tool.validate_input(input_data)
            return await tool.execute(validated_input)
        except Exception as e:
            self.logger.error(f"Tool execution failed: {tool_name} - {e}")
            return ToolOutput(success=False, error=str(e))


def register_all_tools():
    """Register all available tools"""
    from tools.docking import DockingTool, BatchDockingTool
    from tools.rdkit import (
        SmilesTo3DTool,
        ConvertFormatTool,
        OptimizeMoleculeTool,
        CalculatePropertiesTool,
        GenerateMoleculeVariantsTool,
    )
    from tools.pharmacophore import GeneratePharmacophoreTool, ScreenLibraryTool
    from integrations.pubchem import FetchCompoundsTool, SearchCompoundsTool
    from integrations.pdb import FetchProteinTool
    from tools.analysis import (
        AnalyzeInteractionsTool,
        PredictBindingTool,
    )
    from tools.rdkit import PredictADMETTool, SuggestOptimizationTool

    registry = ToolRegistry()

    # Docking tools
    registry.register(DockingTool())
    registry.register(BatchDockingTool())

    # RDKit tools
    registry.register(SmilesTo3DTool())
    registry.register(ConvertFormatTool())
    registry.register(OptimizeMoleculeTool())
    registry.register(CalculatePropertiesTool())
    registry.register(GenerateMoleculeVariantsTool())

    # Pharmacophore tools
    registry.register(GeneratePharmacophoreTool())
    registry.register(ScreenLibraryTool())

    # Integration tools
    registry.register(FetchCompoundsTool())
    registry.register(SearchCompoundsTool())
    registry.register(FetchProteinTool())

    # Analysis tools
    registry.register(AnalyzeInteractionsTool())
    registry.register(PredictBindingTool())
    registry.register(PredictADMETTool())
    registry.register(SuggestOptimizationTool())

    return registry
