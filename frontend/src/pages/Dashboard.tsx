import { useState, useEffect } from 'react'

interface Stats {
  total_jobs: number;
  completed_jobs: number;
  active_jobs: number;
}

interface GPUInfo {
  cuda_available: boolean;
  opencl_available: boolean;
  platform: string;
}

export default function Dashboard() {
  const [stats, setStats] = useState<Stats>({ total_jobs: 0, completed_jobs: 0, active_jobs: 0 })
  const [gpuInfo, setGpuInfo] = useState<GPUInfo | null>(null)

  useEffect(() => {
    fetch('/api/stats')
      .then(res => res.json())
      .then(setStats)
      .catch(console.error)

    fetch('/api/md/gpu-info')
      .then(res => res.json())
      .then(setGpuInfo)
      .catch(console.error)
  }, [])

  return (
    <div className="page">
      <div className="page-header">
        <h2>Dashboard</h2>
        <p>Welcome to BioDockify - Molecular Docking Studio</p>
      </div>

      <div className="stats-grid">
        <div className="stat-card">
          <h3>Total Jobs</h3>
          <div className="value">{stats.total_jobs}</div>
        </div>
        <div className="stat-card">
          <h3>Completed</h3>
          <div className="value">{stats.completed_jobs}</div>
        </div>
        <div className="stat-card">
          <h3>Active</h3>
          <div className="value">{stats.active_jobs}</div>
        </div>
        <div className="stat-card">
          <h3>GPU Platform</h3>
          <div className="value" style={{ fontSize: '1rem' }}>{gpuInfo?.platform || 'Loading...'}</div>
        </div>
      </div>

      <div className="grid-2">
        <div className="card">
          <h3>Quick Actions</h3>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '0.5rem', marginTop: '1rem' }}>
            <button className="btn btn-primary" onClick={() => window.location.hash = '#/docking'}>
              New Docking Job
            </button>
            <button className="btn btn-secondary" onClick={() => window.location.hash = '#/md'}>
              Start MD Simulation
            </button>
            <button className="btn btn-secondary" onClick={() => window.location.hash = '#/ai'}>
              Ask BioDockify AI
            </button>
          </div>
        </div>

        <div className="card">
          <h3>System Status</h3>
          <div style={{ marginTop: '1rem' }}>
            <div className="status-item" style={{ marginBottom: '0.5rem' }}>
              <span className="status-dot"></span>
              <span>Backend API</span>
            </div>
            <div className="status-item" style={{ marginBottom: '0.5rem' }}>
              <span className="status-dot"></span>
              <span>File System</span>
            </div>
            <div className="status-item">
              <span className={`status-dot ${gpuInfo?.cuda_available ? '' : 'warning'}`}></span>
              <span>GPU Acceleration: {gpuInfo?.cuda_available ? 'CUDA' : gpuInfo?.opencl_available ? 'OpenCL' : 'CPU'}</span>
            </div>
          </div>
        </div>
      </div>

      <div className="card" style={{ marginTop: '1rem' }}>
        <h3>Recent Activity</h3>
        <p style={{ color: 'var(--text-secondary)', marginTop: '1rem' }}>No recent activity. Start a new docking job to get started.</p>
      </div>
    </div>
  )
}
