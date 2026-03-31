import { useState, useEffect } from 'react'

interface MDJob {
  job_id: string;
  status: string;
  message?: string;
}

export default function MDSimulation() {
  const [pdbFile, setPdbFile] = useState<File | null>(null)
  const [forcefield, setForcefield] = useState('amber14-all.xml')
  const [duration, setDuration] = useState(10)
  const [temperature, setTemperature] = useState(300)
  const [jobs, setJobs] = useState<MDJob[]>([])
  const [isSubmitting, setIsSubmitting] = useState(false)

  const handleSubmit = async () => {
    if (!pdbFile) {
      alert('Please upload a PDB file')
      return
    }

    setIsSubmitting(true)
    try {
      const response = await fetch('/api/md/jobs', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          pdb_file: pdbFile.name,
          forcefield,
          duration_ns: duration,
          temperature_k: temperature
        })
      })

      const job = await response.json()
      setJobs([job, ...jobs])
    } catch (error) {
      console.error('Failed to submit job:', error)
    }
    setIsSubmitting(false)
  }

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault()
    const file = e.dataTransfer.files[0]
    if (file) setPdbFile(file)
  }

  return (
    <div className="page">
      <div className="page-header">
        <h2>Molecular Dynamics</h2>
        <p>Run OpenMM MD simulations with GPU acceleration</p>
      </div>

      <div className="grid-2">
        <div>
          <div className="card">
            <h3>Structure File</h3>
            <div
              className="drop-zone"
              onDragOver={e => e.preventDefault()}
              onDrop={handleDrop}
              onClick={() => document.getElementById('pdb-input')?.click()}
            >
              <input
                id="pdb-input"
                type="file"
                accept=".pdb"
                style={{ display: 'none' }}
                onChange={e => e.target.files?.[0] && setPdbFile(e.target.files[0])}
              />
              <p>{pdbFile ? pdbFile.name : 'Drop PDB file or click to upload'}</p>
            </div>
          </div>

          <div className="card">
            <h3>Parameters</h3>
            <div className="form-group">
              <label>Force Field</label>
              <select value={forcefield} onChange={e => setForcefield(e.target.value)}>
                <option value="amber14-all.xml">AMBER14 All</option>
                <option value="amber14/tip3p.xml">AMBER14 TIP3P</option>
                <option value="charmm36.xml">CHARMM36</option>
              </select>
            </div>
            <div className="form-group">
              <label>Duration (ns): {duration}</label>
              <input
                type="range"
                min="1"
                max="100"
                value={duration}
                onChange={e => setDuration(Number(e.target.value))}
                style={{ width: '100%' }}
              />
            </div>
            <div className="form-group">
              <label>Temperature (K): {temperature}</label>
              <input
                type="range"
                min="200"
                max="400"
                value={temperature}
                onChange={e => setTemperature(Number(e.target.value))}
                style={{ width: '100%' }}
              />
            </div>
            <button
              className="btn btn-primary"
              style={{ width: '100%' }}
              onClick={handleSubmit}
              disabled={isSubmitting}
            >
              {isSubmitting ? 'Starting...' : 'Start MD Simulation'}
            </button>
          </div>
        </div>

        <div>
          <div className="card">
            <h3>Active Simulations</h3>
            <div className="job-list">
              {jobs.length === 0 ? (
                <p style={{ color: 'var(--text-secondary)' }}>No active simulations</p>
              ) : (
                jobs.map(job => (
                  <div key={job.job_id} className="job-item">
                    <span className="job-id">Job {job.job_id}</span>
                    <span className={`job-status ${job.status}`}>{job.status}</span>
                  </div>
                ))
              )}
            </div>
          </div>

          <div className="card">
            <h3>GPU Acceleration</h3>
            <p style={{ color: 'var(--text-secondary)', marginTop: '0.5rem' }}>
              OpenMM will automatically detect and use the best available platform:
            </p>
            <ul style={{ marginTop: '1rem', paddingLeft: '1.5rem', color: 'var(--text-secondary)' }}>
              <li>CUDA (NVIDIA GPUs)</li>
              <li>OpenCL (AMD/Intel GPUs)</li>
              <li>CPU (fallback)</li>
            </ul>
          </div>
        </div>
      </div>
    </div>
  )
}
