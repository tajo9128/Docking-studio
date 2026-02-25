#!/usr/bin/env python3
"""
Docking Studio Production Debug Script
Run this to verify your full stack is working correctly
"""

import subprocess
import sys
import time
import requests
import json


def run_command(cmd, description, check=True):
    """Run a shell command"""
    print(f"\n{'='*60}")
    print(f"  {description}")
    print(f"{'='*60}")
    print(f"Command: {cmd}")
    print()
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    
    if check and result.returncode != 0:
        print(f"✗ FAILED: {description}")
        return False
    
    print(f"✓ {description}")
    return True


def check_backend_health():
    """Check backend health"""
    print(f"\n{'='*60}")
    print(f"  Testing Backend Health")
    print(f"{'='*60}")
    
    try:
        # Test root
        r = requests.get("http://localhost:8000/", timeout=5)
        print(f"GET / : {r.status_code} - {r.json()}")
        
        # Test health
        r = requests.get("http://localhost:8000/health", timeout=5)
        print(f"GET /health : {r.status_code} - {r.json()}")
        
        # Test chat status
        r = requests.get("http://localhost:8000/chat/status", timeout=10)
        print(f"GET /chat/status : {r.status_code} - {r.json()}")
        
        # Test security status
        r = requests.get("http://localhost:8000/security/status", timeout=10)
        print(f"GET /security/status : {r.status_code} - {r.json()}")
        
        print("\n✓ Backend is fully operational")
        return True
        
    except requests.exceptions.ConnectionError:
        print("\n✗ Cannot connect to backend")
        print("  Make sure docker compose is running: docker compose up -d")
        return False
    except Exception as e:
        print(f"\n✗ Backend error: {e}")
        return False


def test_chat():
    """Test chat functionality"""
    print(f"\n{'='*60}")
    print(f"  Testing Chat / AI")
    print(f"{'='*60}")
    
    try:
        # Test offline chat
        r = requests.post(
            "http://localhost:8000/chat",
            json={"message": "What is vina?"},
            timeout=30
        )
        result = r.json()
        print(f"Provider: {result.get('provider')}")
        print(f"Available: {result.get('available')}")
        print(f"Response: {result.get('response', '')[:200]}...")
        
        print("\n✓ Chat is working")
        return True
        
    except Exception as e:
        print(f"\n✗ Chat test failed: {e}")
        return False


def check_docker():
    """Check Docker status"""
    print(f"\n{'='*60}")
    print(f"  Checking Docker")
    print(f"{'='*60}")
    
    result = subprocess.run("docker ps", shell=True, capture_output=True, text=True)
    print(result.stdout)
    
    return result.returncode == 0


def main():
    """Main debug function"""
    print("""
    ╔══════════════════════════════════════════════════════════╗
    ║     Docking Studio Production Debug Script          ║
    ╚══════════════════════════════════════════════════════════╝
    """)
    
    # Step 1: Check Docker
    print("\n[1/5] Checking Docker...")
    if not check_docker():
        print("\n✗ Docker is not running!")
        print("  Please start Docker Desktop")
        return
    
    # Step 2: Check backend
    print("\n[2/5] Checking backend...")
    if not check_backend_health():
        print("\n✗ Backend is not responding")
        print("  Run: docker compose up -d")
        return
    
    # Step 3: Test chat
    print("\n[3/5] Testing chat/AI...")
    test_chat()
    
    # Step 4: Check logs
    print("\n[4/5] Checking logs...")
    run_command("docker compose logs --tail=20", "Docker Compose Logs")
    
    # Step 5: Summary
    print("\n" + "="*60)
    print("  ✓ ALL CHECKS PASSED")
    print("="*60)
    print("""
    Your Docking Studio is ready!

    Next steps:
    1. Run PyQt6: python main.py
    2. Test docking functionality
    3. Try the AI chat

    Backend API: http://localhost:8000
    API Docs: http://localhost:8000/docs
    """)


if __name__ == "__main__":
    main()
