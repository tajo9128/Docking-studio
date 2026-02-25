"""
Security Monitor Module
Agent Zero security brain for runtime security monitoring
"""

import sqlite3
import json
import logging
from datetime import datetime
from typing import Dict, List, Optional
from dataclasses import asdict

from security.scanner import (
    SecurityScanner, 
    SecurityScanResult, 
    Severity
)

logger = logging.getLogger(__name__)

DB_PATH = "storage/jobs.db"


class SecurityMonitor:
    """
    Agent Zero Security Monitor
    Runs security scans and manages security state
    """
    
    def __init__(self, db_path: str = DB_PATH):
        self.db_path = db_path
        self.scanner = SecurityScanner()
        self._init_db()
    
    def _init_db(self):
        """Initialize security tables"""
        import os
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        
        cur.execute("""
            CREATE TABLE IF NOT EXISTS security_reports (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                scan_type TEXT NOT NULL,
                severity TEXT NOT NULL,
                issues_count INTEGER DEFAULT 0,
                issues_json TEXT,
                raw_output TEXT,
                is_secure INTEGER DEFAULT 1,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        cur.execute("""
            CREATE TABLE IF NOT EXISTS security_status (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                last_scan_at DATETIME,
                overall_severity TEXT,
                is_secure INTEGER DEFAULT 1,
                total_issues INTEGER DEFAULT 0,
                scan_results_json TEXT
            )
        """)
        
        conn.commit()
        conn.close()
    
    def store_report(self, result: SecurityScanResult) -> int:
        """Store a scan result in the database"""
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        
        cur.execute("""
            INSERT INTO security_reports 
            (scan_type, severity, issues_count, issues_json, raw_output, is_secure)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (
            result.scan_type,
            result.severity.value,
            result.issues_count,
            json.dumps(result.issues),
            result.raw_output[:10000],
            1 if result.severity in [Severity.SAFE, Severity.LOW] else 0
        ))
        
        report_id = cur.lastrowid
        conn.commit()
        conn.close()
        
        return report_id
    
    def update_overall_status(self, results: Dict[str, SecurityScanResult]):
        """Update overall security status"""
        summary = self.scanner.get_summary(results)
        
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        
        cur.execute("""
            INSERT INTO security_status 
            (last_scan_at, overall_severity, is_secure, total_issues, scan_results_json)
            VALUES (?, ?, ?, ?, ?)
        """, (
            datetime.now().isoformat(),
            summary['worst_severity'],
            1 if summary['is_secure'] else 0,
            summary['total_issues'],
            json.dumps(summary)
        ))
        
        conn.commit()
        conn.close()
    
    def run_full_scan(self) -> Dict:
        """
        Run all security scans and store results.
        
        Returns:
            Summary dictionary
        """
        logger.info("Starting full security scan...")
        
        results = self.scanner.scan_all()
        
        for result in results.values():
            self.store_report(result)
        
        summary = self.scanner.get_summary(results)
        self.update_overall_status(results)
        
        logger.info(f"Security scan complete: {summary['worst_severity']} ({summary['total_issues']} issues)")
        
        return summary
    
    def get_latest_status(self) -> Optional[Dict]:
        """Get the latest security status"""
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        
        cur.execute("""
            SELECT last_scan_at, overall_severity, is_secure, total_issues, scan_results_json
            FROM security_status
            ORDER BY id DESC
            LIMIT 1
        """)
        
        row = cur.fetchone()
        conn.close()
        
        if row:
            return {
                'last_scan_at': row[0],
                'overall_severity': row[1],
                'is_secure': bool(row[2]),
                'total_issues': row[3],
                'scan_results': json.loads(row[4]) if row[4] else {}
            }
        return None
    
    def get_recent_reports(self, limit: int = 10) -> List[Dict]:
        """Get recent security reports"""
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        
        cur.execute("""
            SELECT scan_type, severity, issues_count, created_at
            FROM security_reports
            ORDER BY created_at DESC
            LIMIT ?
        """, (limit,))
        
        rows = cur.fetchall()
        conn.close()
        
        return [
            {
                'scan_type': r[0],
                'severity': r[1],
                'issues_count': r[2],
                'created_at': r[3]
            }
            for r in rows
        ]
    
    def block_if_critical(self) -> bool:
        """
        Check if there are critical issues and raise exception if so.
        
        Returns:
            True if secure, raises RuntimeError if not secure
        """
        status = self.get_latest_status()
        
        if status and not status['is_secure']:
            if status['overall_severity'] in ['CRITICAL', 'HIGH']:
                raise RuntimeError(
                    f"Security blocked: {status['total_issues']} issues found "
                    f"(worst: {status['overall_severity']})"
                )
        
        return True
    
    def is_docking_allowed(self) -> bool:
        """Check if docking is allowed based on security status"""
        status = self.get_latest_status()
        
        if status is None:
            logger.warning("No security scan performed yet")
            return True
        
        return status['is_secure']
    
    def get_security_issues(self, scan_type: Optional[str] = None) -> List[Dict]:
        """Get detailed security issues"""
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        
        if scan_type:
            cur.execute("""
                SELECT scan_type, severity, issues_json, raw_output, created_at
                FROM security_reports
                WHERE scan_type = ?
                ORDER BY created_at DESC
                LIMIT 1
            """, (scan_type,))
        else:
            cur.execute("""
                SELECT scan_type, severity, issues_json, raw_output, created_at
                FROM security_reports
                ORDER BY created_at DESC
                LIMIT 10
            """)
        
        rows = cur.fetchall()
        conn.close()
        
        return [
            {
                'scan_type': r[0],
                'severity': r[1],
                'issues': json.loads(r[2]) if r[2] else [],
                'raw_output': r[3],
                'created_at': r[4]
            }
            for r in rows
        ]
    
    def clear_old_reports(self, days: int = 30):
        """Clear security reports older than specified days"""
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        
        cur.execute("""
            DELETE FROM security_reports
            WHERE created_at < datetime('now', '-' || ? || ' days')
        """, (days,))
        
        deleted = cur.rowcount
        conn.commit()
        conn.close()
        
        logger.info(f"Cleared {deleted} old security reports")
        return deleted


class SecurityPolicy:
    """
    Security policy engine for automated decision making
    """
    
    ALLOW_INSECURE = False
    CRITICAL_THRESHOLD = 0
    HIGH_THRESHOLD = 5
    
    @classmethod
    def evaluate(cls, summary: Dict) -> tuple[bool, str]:
        """
        Evaluate security summary against policy.
        
        Returns:
            (allowed, message)
        """
        worst = summary.get('worst_severity', 'SAFE')
        total = summary.get('total_issues', 0)
        
        if worst == 'CRITICAL':
            return False, f"CRITICAL security issues detected ({total} total)"
        
        if worst == 'HIGH' and total > cls.HIGH_THRESHOLD:
            return False, f"Too many HIGH severity issues ({total})"
        
        if worst == 'HIGH':
            if cls.ALLOW_INSECURE:
                return True, f"HIGH severity issues found ({total}), but ALLOW_INSECURE is True"
            return False, f"HIGH severity issues require attention ({total})"
        
        return True, "Security check passed"


if __name__ == "__main__":
    monitor = SecurityMonitor()
    summary = monitor.run_full_scan()
    print(json.dumps(summary, indent=2))
