import subprocess
import sys
from typing import Tuple, Optional

def detect_gpu_platform() -> Tuple[str, bool, bool]:
    cuda_available = False
    opencl_available = False
    platform = "CPU"
    
    try:
        result = subprocess.run(
            ['nvidia-smi'],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            cuda_available = True
            platform = "CUDA"
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    if not cuda_available:
        try:
            import numpy as np
            # Check for OpenCL using a simple test
            # In production, would check OpenCL devices directly
            opencl_available = False
            if opencl_available:
                platform = "OpenCL"
        except ImportError:
            pass
    
    return platform, cuda_available, opencl_available


def get_openmm_platform() -> str:
    try:
        import openmm as mm
        from openmm import Platform
        
        platform, cuda, opencl = detect_gpu_platform()
        
        if cuda:
            cuda_platform = Platform.getPlatformByName('CUDA')
            if cuda_platform.getSpeed() > 0:
                return 'CUDA'
        
        if opencl:
            opencl_platform = Platform.getPlatformByName('OpenCL')
            if opencl_platform.getSpeed() > 0:
                return 'OpenCL'
        
        return 'CPU'
    except Exception as e:
        print(f"Error detecting OpenMM platform: {e}")
        return 'CPU'


def get_gpu_info() -> dict:
    platform, cuda, opencl = detect_gpu_platform()
    
    gpu_info = {
        'platform': platform,
        'cuda_available': cuda,
        'opencl_available': opencl,
        'gpus': []
    }
    
    if cuda:
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=name,memory.total,driver_version', '--format=csv,noheader'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        parts = line.split(',')
                        if len(parts) >= 2:
                            gpu_info['gpus'].append({
                                'name': parts[0].strip(),
                                'memory': parts[1].strip(),
                                'driver': parts[2].strip() if len(parts) > 2 else 'Unknown'
                            })
        except Exception:
            pass
    
    return gpu_info
