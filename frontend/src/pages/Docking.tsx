import { useState } from 'react'

interface DockingJob {
  job_id: string;
  status: string;
  message?: string;
}

export default function Docking() {
  const [receptorFile, setReceptorFile] = useState<File | null>(null)
  const [ligandFile, setLigandFile] = useState<File | null>(null)
  const [ligandSmiles, setLigandSmiles] = useState('')
  const [exhaustiveness, setExhaustiveness] = useState(32)
  const [nPoses, setNPoses] = useState(10)
  const [jobs, setJobs] = useState<DockingJob[]>([])
  const [selectedJob, setSelectedJob] = useState<string | null>(null)
  const [isSubmitting, setIsSubmitting] = useState(false)

  const handleSubmit = async () => {
    if (!receptorFile) {
      alert('Please upload a receptor file')
      return
    }
    if (!ligandFile && !ligandSmiles) {
      alert('Please upload a ligand file or enter a SMILES string')
      return
    }

    setIsSubmitting(true)
    try {
      const formData = new FormData()
      formData.append('receptor', receptorFile)
      if (ligandFile) formData.append('ligand', ligandFile)
      formData.append('smiles', ligandSmiles)
      formData.append('exhaustiveness', exhaustiveness.toString())
      formData.append('n_poses', nPoses.toString())

      const response = await fetch('/api/docking/jobs', {
        method: 'POST',
        body: formData
      })

      const job = await response.json()
      setJobs([job, ...jobs])
      setSelectedJob(job.job_id)
    } catch (error) {
      console.error('Failed to submit job:', error)
    }
    setIsSubmitting(false)
  }

  const handleDrop = (e: React.DragEvent, setter: (f: File) => void) => {
    e.preventDefault()
    const file = e.dataTransfer.files[0]
    if (file) setter(file)
  }

  return (
    <div className="page">
      <div className="page-header">
        <h2>Molecular Docking</h2>
        <p>Run AutoDock Vina docking simulations</p>
      </div>

      <div className="grid-2">
        <div>
          <div className="card">
            <h3>Receptor</h3>
            <div
              className="drop-zone"
              onDragOver={e => e.preventDefault()}
              onDrop={e => handleDrop(e, setReceptorFile)}
              onClick={() => document.getElementById('receptor-input')?.click()}
            >
              <input
                id="receptor-input"
                type="file"
                accept=".pdb,.pdbqt"
                style={{ display: 'none' }}
                onChange={e => e.target.files?.[0] && setReceptorFile(e.target.files[0])}
              />
              <p>{receptorFile ? receptorFile.name : 'Drop PDB/PDBQT file or click to upload'}</p>
            </div>
          </div>

          <div className="card">
            <h3>Ligand</h3>
            <div
              className="drop-zone"
              onDragOver={e => e.preventDefault()}
              onDrop={e => handleDrop(e, setLigandFile)}
              onClick={() => document.getElementById('ligand-input')?.click()}
            >
              <input
                id="ligand-input"
                type="file"
                accept=".pdb,.pdbqt,.sdf,.mol"
                style={{ display: 'none' }}
                onChange={e => e.target.files?.[0] && setLigandFile(e.target.files[0])}
              />
              <p>{ligandFile ? ligandFile.name : 'Drop PDB/PDBQT/SDF file or click to upload'}</p>
            </div>

            <div style={{ marginTop: '1rem', textAlign: 'center', color: 'var(--text-secondary)' }}>
              - OR -
            </div>

            <div className="form-group" style={{ marginTop: '1rem' }}>
              <label>SMILES String</label>
              <input
                type="text"
                placeholder="Enter SMILES (e.g., CC(=O)Oc1ccccc1C(=O)O)"
                value={ligandSmiles}
                onChange={e => setLigandSmiles(e.target.value)}
              />
            </div>
          </div>

          <div className="card">
            <h3>Parameters</h3>
            <div className="form-group">
              <label>Exhaustiveness: {exhaustiveness}</label>
              <input
                type="range"
                min="1"
                max="64"
                value={exhaustiveness}
                onChange={e => setExhaustiveness(Number(e.target.value))}
                style={{ width: '100%' }}
              />
            </div>
            <div className="form-group">
              <label>Number of Poses: {nPoses}</label>
              <input
                type="range"
                min="1"
                max="20"
                value={nPoses}
                onChange={e => setNPoses(Number(e.target.value))}
                style={{ width: '100%' }}
              />
            </div>
            <button
              className="btn btn-primary"
              style={{ width: '100%' }}
              onClick={handleSubmit}
              disabled={isSubmitting}
            >
              {isSubmitting ? 'Submitting...' : 'Start Docking'}
            </button>
          </div>
        </div>

        <div>
          <div className="card">
            <h3>Job Queue</h3>
            <div className="job-list">
              {jobs.length === 0 ? (
                <p style={{ color: 'var(--text-secondary)' }}>No jobs yet</p>
              ) : (
                jobs.map(job => (
                  <div
                    key={job.job_id}
                    className={`job-item ${selectedJob === job.job_id ? 'active' : ''}`}
                    onClick={() => setSelectedJob(job.job_id)}
                  >
                    <span className="job-id">Job {job.job_id}</span>
                    <span className={`job-status ${job.status}`}>{job.status}</span>
                  </div>
                ))
              )}
            </div>
          </div>

          {selectedJob && (
            <div className="card">
              <h3>Results - Job {selectedJob}</h3>
              <div className="results-viewer">
                <p style={{ color: 'var(--text-secondary)' }}>
                  Results will appear here when the job completes.
                </p>
                <div className="pose-grid">
                  {[1, 2, 3, 4].map(rank => (
                    <div key={rank} className="pose-card">
                      <div className="pose-rank">Rank #{rank}</div>
                      <div className="pose-score">-8.5</div>
                      <div className="pose-energy">kcal/mol</div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  )
}
