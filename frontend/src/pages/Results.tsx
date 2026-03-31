import { useState, useEffect } from 'react'

interface Job {
  job_id: string;
  type: string;
  status: string;
  created_at: string;
}

export default function Results() {
  const [jobs, setJobs] = useState<Job[]>([])
  const [selectedJob, setSelectedJob] = useState<Job | null>(null)

  useEffect(() => {
    // Load jobs from API
    fetch('/api/docking/jobs')
      .then(res => res.json())
      .then(data => setJobs(data.jobs || []))
      .catch(console.error)
  }, [])

  return (
    <div className="page">
      <div className="page-header">
        <h2>Results</h2>
        <p>View completed docking and simulation results</p>
      </div>

      <div className="grid-2">
        <div className="card">
          <h3>All Jobs</h3>
          <div className="job-list">
            {jobs.length === 0 ? (
              <p style={{ color: 'var(--text-secondary)' }}>No jobs found. Run a docking simulation to see results here.</p>
            ) : (
              jobs.map(job => (
                <div
                  key={job.job_id}
                  className={`job-item ${selectedJob?.job_id === job.job_id ? 'active' : ''}`}
                  onClick={() => setSelectedJob(job)}
                >
                  <div>
                    <div className="job-id">{job.job_id}</div>
                    <div style={{ fontSize: '0.75rem', color: 'var(--text-secondary)' }}>
                      {job.type} - {job.created_at}
                    </div>
                  </div>
                  <span className={`job-status ${job.status}`}>{job.status}</span>
                </div>
              ))
            )}
          </div>
        </div>

        <div className="card">
          <h3>Job Details</h3>
          {selectedJob ? (
            <div className="results-viewer">
              <p><strong>Job ID:</strong> {selectedJob.job_id}</p>
              <p><strong>Type:</strong> {selectedJob.type}</p>
              <p><strong>Status:</strong> {selectedJob.status}</p>
              <p><strong>Created:</strong> {selectedJob.created_at}</p>
              <div style={{ marginTop: '1rem' }}>
                <h4>Output Files</h4>
                <p style={{ color: 'var(--text-secondary)' }}>Download links will appear here when available.</p>
              </div>
            </div>
          ) : (
            <p style={{ color: 'var(--text-secondary)' }}>Select a job to view details</p>
          )}
        </div>
      </div>
    </div>
  )
}
