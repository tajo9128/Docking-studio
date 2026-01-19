"""Custom Exceptions for Error Propagation"""

from typing import Optional

class BioDockVizException(Exception):
    """Base exception for BioDockViz"""
    
    def __init__(
        self,
        message: str,
        code: str = "UNKNOWN_ERROR",
        details: Optional[str] = None,
        correlation_id: Optional[str] = None,
    ):
        self.message = message
        self.code = code
        self.details = details
        self.correlation_id = correlation_id
        super().__init__(self.message)
    
    def to_dict(self) -> dict:
        """Convert exception to dictionary"""
        return {
            'type': self.__class__.__name__,
            'code': self.code,
            'message': self.message,
            'details': self.details,
            'correlation_id': self.correlation_id,
        }

class UploadException(BioDockVizException):
    """Upload exception"""
    pass

class ParseException(BioDockVizException):
    """Parsing exception"""
    pass

class ValidationException(BioDockVizException):
    """Validation exception"""
    pass

class ProcessingException(BioDockVizException):
    """Processing exception"""
    pass

class AnalysisException(BioDockVizException):
    """Analysis exception"""
    pass

class VisualizationException(BioDockVizException):
    """Visualization exception"""
    pass

class ExportException(BioDockVizException):
    """Export exception"""
    pass

class DatabaseException(BioDockVizException):
    """Database exception"""
    pass

class MolecularException(BioDockVizException):
    """Molecular engine exception"""
    pass

class SystemException(BioDockVizException):
    """System exception"""
    pass

class RateLimitException(BioDockVizException):
    """Rate limit exception"""
    pass

class AuthenticationException(BioDockVizException):
    """Authentication exception"""
    pass

class ForbiddenException(BioDockVizException):
    """Forbidden exception"""
    pass

class NotFoundException(BioDockVizException):
    """Not found exception"""
    pass
