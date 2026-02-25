"""
Security Scanner Engine
Provides unified interface for security scanning tools
"""

import subprocess
import json
import os
import logging
from typing import Dict, Optional, Tuple, List
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class Severity(Enum):
    """Security severity levels"""
    SAFE = "SAFE"
    LOW = "LOW"
    MEDIUM = "MEDIUM"
    HIGH = "HIGH"
    CRITICAL = "CRITICAL"


@dataclass
class SecurityScanResult:
    """Result of a security scan"""
    scan_type: str
    severity: Severity
    issues_count: int
    issues: List[Dict]
    raw_output: str
    success: bool
    error_message: Optional[str] = None


def run_command(cmd: List[str], timeout: int = 300) -> Tuple[str, bool]:
    """
    Run a shell command and return output.
    
    Args:
        cmd: Command and arguments as list
        timeout: Timeout in seconds
    
    Returns:
        Tuple of (output, success)
    """
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        )
        return result.stdout + result.stderr, result.returncode == 0
    except subprocess.TimeoutExpired:
        return f"Command timed out after {timeout}s", False
    except FileNotFoundError:
        return f"Command not found: {cmd[0]}", False
    except Exception as e:
        return str(e), False


def parse_trivy_output(output: str) -> Tuple[Severity, int, List[Dict]]:
    """Parse Trivy output to extract issues"""
    issues = []
    severity = Severity.SAFE
    
    lines = output.split('\n')
    for line in lines:
        if 'CRITICAL' in line.upper():
            severity = Severity.CRITICAL
            issues.append({'severity': 'CRITICAL', 'message': line.strip()})
        elif 'HIGH' in line.upper() and severity != Severity.CRITICAL:
            severity = Severity.HIGH
            issues.append({'severity': 'HIGH', 'message': line.strip()})
        elif 'MEDIUM' in line.upper() and severity not in [Severity.CRITICAL, Severity.HIGH]:
            severity = Severity.MEDIUM
            issues.append({'severity': 'MEDIUM', 'message': line.strip()})
    
    return severity, len(issues), issues


def parse_bandit_output(output: str) -> Tuple[Severity, int, List[Dict]]:
    """Parse Bandit output"""
    issues = []
    severity = Severity.SAFE
    
    try:
        if 'Issue' in output or 'Vulnerability' in output:
            severity = Severity.HIGH
            for line in output.split('\n'):
                if 'Issue' in line:
                    issues.append({'severity': 'HIGH', 'message': line.strip()})
    except Exception as e:
        logger.warning(f"Failed to parse bandit output: {e}")
    
    return severity, len(issues), issues


def parse_safety_output(output: str) -> Tuple[Severity, int, List[Dict]]:
    """Parse Safety output"""
    issues = []
    severity = Severity.SAFE
    
    if 'vulnerable' in output.lower() or 'insecure' in output.lower():
        severity = Severity.CRITICAL
        for line in output.split('\n'):
            if 'vulnerable' in line.lower() or 'insecure' in line.lower():
                issues.append({'severity': 'CRITICAL', 'message': line.strip()})
    
    return severity, len(issues), issues


def parse_gitleaks_output(output: str) -> Tuple[Severity, int, List[Dict]]:
    """Parse Gitleaks output"""
    issues = []
    severity = Severity.SAFE
    
    if 'Findings:' in output or 'secrets found' in output.lower():
        severity = Severity.CRITICAL
        for line in output.split('\n'):
            if ' Findings' in line or 'secret' in line.lower():
                issues.append({'severity': 'CRITICAL', 'message': line.strip()})
    
    return severity, len(issues), issues


class SecurityScanner:
    """Unified security scanner interface"""
    
    def __init__(self):
        self.base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    def scan_trivy_filesystem(self) -> SecurityScanResult:
        """Scan filesystem with Trivy"""
        output, success = run_command([
            'trivy', 'fs', '.',
            '--severity', 'HIGH,CRITICAL',
            '--format', 'json'
        ])
        
        severity, count, issues = parse_trivy_output(output)
        
        return SecurityScanResult(
            scan_type='trivy_fs',
            severity=severity,
            issues_count=count,
            issues=issues,
            raw_output=output,
            success=success
        )
    
    def scan_trivy_image(self, image: str = "docking-studio-backend:latest") -> SecurityScanResult:
        """Scan Docker image with Trivy"""
        output, success = run_command([
            'trivy', 'image', image,
            '--severity', 'HIGH,CRITICAL',
            '--format', 'json'
        ], timeout=600)
        
        severity, count, issues = parse_trivy_output(output)
        
        return SecurityScanResult(
            scan_type='trivy_image',
            severity=severity,
            issues_count=count,
            issues=issues,
            raw_output=output,
            success=success
        )
    
    def scan_bandit(self) -> SecurityScanResult:
        """Scan Python code with Bandit"""
        output, success = run_command([
            'bandit', '-r', 'backend/',
            '-f', 'json'
        ])
        
        severity, count, issues = parse_bandit_output(output)
        
        return SecurityScanResult(
            scan_type='bandit',
            severity=severity,
            issues_count=count,
            issues=issues,
            raw_output=output,
            success=success
        )
    
    def scan_safety(self) -> SecurityScanResult:
        """Scan dependencies with Safety"""
        output, success = run_command([
            'safety', 'check',
            '--json'
        ])
        
        severity, count, issues = parse_safety_output(output)
        
        return SecurityScanResult(
            scan_type='safety',
            severity=severity,
            issues_count=count,
            issues=issues,
            raw_output=output,
            success=success
        )
    
    def scan_gitleaks(self) -> SecurityScanResult:
        """Scan for secrets with Gitleaks"""
        output, success = run_command([
            'gitleaks', 'detect',
            '--no-git',
            '--report-format', 'json'
        ])
        
        severity, count, issues = parse_gitleaks_output(output)
        
        return SecurityScanResult(
            scan_type='gitleaks',
            severity=severity,
            issues_count=count,
            issues=issues,
            raw_output=output,
            success=success
        )
    
    def scan_all(self) -> Dict[str, SecurityScanResult]:
        """Run all security scans"""
        results = {}
        
        scanners = [
            ('trivy_fs', self.scan_trivy_filesystem),
            ('bandit', self.scan_bandit),
            ('safety', self.scan_safety),
            ('gitleaks', self.scan_gitleaks),
        ]
        
        for name, scanner in scanners:
            try:
                result = scanner()
                results[name] = result
                logger.info(f"{name}: {result.severity.value} ({result.issues_count} issues)")
            except Exception as e:
                logger.error(f"{name} scan failed: {e}")
                results[name] = SecurityScanResult(
                    scan_type=name,
                    severity=Severity.UNKNOWN if hasattr(Severity, 'UNKNOWN') else Severity.SAFE,
                    issues_count=0,
                    issues=[],
                    raw_output=str(e),
                    success=False,
                    error_message=str(e)
                )
        
        return results
    
    def is_secure(self, results: Dict[str, SecurityScanResult]) -> bool:
        """Check if all scans passed (no CRITICAL or HIGH issues)"""
        for result in results.values():
            if result.severity in [Severity.CRITICAL, Severity.HIGH]:
                return False
        return True
    
    def get_summary(self, results: Dict[str, SecurityScanResult]) -> Dict:
        """Get summary of all scans"""
        total_issues = sum(r.issues_count for r in results.values())
        worst = max((r.severity for r in results.values()), 
                   key=lambda s: [Severity.SAFE, Severity.LOW, Severity.MEDIUM, Severity.HIGH, Severity.CRITICAL].index(s))
        
        return {
            'total_issues': total_issues,
            'worst_severity': worst.value,
            'is_secure': self.is_secure(results),
            'scans': {
                name: {
                    'severity': r.severity.value,
                    'issues': r.issues_count,
                    'success': r.success
                }
                for name, r in results.items()
            }
        }


if __name__ == "__main__":
    scanner = SecurityScanner()
    results = scanner.scan_all()
    summary = scanner.get_summary(results)
    print(json.dumps(summary, indent=2))
