import { useState, useEffect, useRef } from 'react'
import { useTheme } from '@/contexts/ThemeContext'

interface DockingResult {
  id: number
  job_uuid: string
  pose_id: number
  ligand_name: string
  vina_score: number | null
  gnina_score: number | null
  rf_score: number | null
  consensus: number | null
  pdb_data: string | null
  hydrophobic_term?: number
  rotatable_penalty?: number
  lipo_contact?: number
  final_score?: number
  composite_score?: number
  constraint_penalty?: number
}

interface Job {
  id: number
  job_uuid: string
  job_name: string
  receptor_file: string | null
  ligand_file: string | null
  status: string
  created_at: string
  completed_at: string | null
  binding_energy: number | null
  confidence_score: number | null
  engine: string | null
}

export function Results() {
  const { theme } = useTheme()
  const isDark = theme === 'dark'
  const viewer3dRef = useRef<HTMLDivElement>(null)
  const viewerRef = useRef<any>(null)
  const [jobs, setJobs] = useState<Job[]>([])
  const [selectedJob, setSelectedJob] = useState<Job | null>(null)
  const [selectedResults, setSelectedResults] = useState<DockingResult[]>([])
  const [selectedPose, setSelectedPose] = useState<DockingResult | null>(null)
  const [loading, setLoading] = useState(true)
  const [activeTab, setActiveTab] = useState<'all' | 'completed' | 'running' | 'failed'>('all')
  const [jobFiles, setJobFiles] = useState<Record<string, { filename: string; url: string | null; size_bytes?: number; exists: boolean }>>({}) 
  const [logText, setLogText] = useState<string>('')
  const [aiExplanation, setAiExplanation] = useState<string>('')
  const [aiLoading, setAiLoading] = useState(false)
  const [aiConvId, setAiConvId] = useState<string | null>(null)

  const bgClass = isDark ? 'bg-gray-900' : 'bg-gray-50'
  const cardBg = isDark ? 'bg-gray-800' : 'bg-white'
  const borderClass = isDark ? 'border-gray-700' : 'border-gray-200'
  const textClass = isDark ? 'text-white' : 'text-gray-900'
  const subtextClass = isDark ? 'text-gray-400' : 'text-gray-500'
  const hoverClass = isDark ? 'hover:bg-gray-700' : 'hover:bg-gray-50'

  useEffect(() => {
    fetchJobs()
  }, [])

  useEffect(() => {
    setAiExplanation('')
    setAiConvId(null)
  }, [selectedPose?.id])

  useEffect(() => {
    if (!viewer3dRef.current || viewerRef.current) return

    const initViewer = () => {
      if (!(window as any).$3Dmol || !viewer3dRef.current) return false
      const viewer = (window as any).$3Dmol.createViewer(viewer3dRef.current, {
        backgroundColor: '#1a1a2e',
      })
      viewerRef.current = viewer
      viewer.addModel('', 'pdb')
      viewer.zoomTo()
      viewer.render()
      return true
    }

    if (!initViewer()) {
      const interval = setInterval(() => {
        if (initViewer()) clearInterval(interval)
      }, 200)
      return () => clearInterval(interval)
    }
  }, [selectedJob])

  useEffect(() => {
    if (!selectedPose || !viewerRef.current) return
    const viewer = viewerRef.current
    viewer.clear()
    if (selectedPose.pdb_data) {
      const text = selectedPose.pdb_data.trim()
      const isPDBQT = text.includes('ROOT') || text.includes('BRANCH')
      viewer.addModel(text, isPDBQT ? 'pdbqt' : 'pdb')
    }
    viewer.setStyle({}, { stick: {} })
    viewer.zoomTo()
    viewer.render()
  }, [selectedPose])

  const fetchJobs = async () => {
    setLoading(true)
    try {
      const res = await fetch('/jobs')
      const data = await res.json()
      const jobList = data.jobs || []
      setJobs(jobList)
      if (jobList.length > 0 && !selectedJob) {
        selectJob(jobList[0])
      }
    } catch (err) {
      console.error('Failed to fetch jobs:', err)
    } finally {
      setLoading(false)
    }
  }

  const selectJob = async (job: Job) => {
    setSelectedJob(job)
    setSelectedPose(null)
    setJobFiles({})
    setLogText('')
    try {
      const [resResults, resFiles] = await Promise.all([
        fetch(`/jobs/${job.job_uuid}/results`),
        fetch(`/jobs/${job.job_uuid}/files`),
      ])
      const dataResults = await resResults.json()
      const results = dataResults.results || []
      setSelectedResults(results)
      if (results.length > 0) setSelectedPose(results[0])

      if (resFiles.ok) {
        const dataFiles = await resFiles.json()
        setJobFiles(dataFiles.files || {})
        setLogText(dataFiles.log_text || '')
      }
    } catch (err) {
      console.error('Failed to fetch results:', err)
      setSelectedResults([])
    }
  }

  const deleteJob = async (jobUuid: string) => {
    if (!confirm('Delete this job and all its results?')) return
    try {
      await fetch(`/jobs/${jobUuid}`, { method: 'DELETE' })
      setJobs(prev => prev.filter(j => j.job_uuid !== jobUuid))
      if (selectedJob?.job_uuid === jobUuid) {
        setSelectedJob(null)
        setSelectedResults([])
        setSelectedPose(null)
      }
    } catch (err) {
      console.error('Failed to delete job:', err)
    }
  }

  const handleDownloadPDB = (pose: DockingResult) => {
    if (!pose.pdb_data) return
    const blob = new Blob([pose.pdb_data], { type: 'chemical/x-pdb' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `${pose.ligand_name || `pose_${pose.pose_id}`}.pdb`
    a.click()
    URL.revokeObjectURL(url)
  }

  const handleSnapshot = () => {
    if (!viewerRef.current) return
    try {
      const uri = viewerRef.current.pngURI()
      const a = document.createElement('a')
      a.href = uri
      a.download = `snapshot_${selectedJob?.job_uuid || 'pose'}.png`
      a.click()
    } catch (e) {
      console.error('Snapshot failed', e)
    }
  }

  const askAI = async () => {
    if (!selectedPose) return
    setAiLoading(true)
    setAiExplanation('')
    const scoreLines = [
      `Ligand: ${selectedPose.ligand_name || `Pose ${selectedPose.pose_id}`}`,
      selectedPose.vina_score != null ? `Vina Score: ${selectedPose.vina_score.toFixed(2)} kcal/mol` : null,
      selectedPose.gnina_score != null ? `GNINA CNN Score: ${selectedPose.gnina_score.toFixed(2)}` : null,
      selectedPose.rf_score != null ? `RF Score: ${selectedPose.rf_score.toFixed(2)}` : null,
      `Composite Score: ${(selectedPose.composite_score ?? selectedPose.final_score ?? selectedPose.vina_score)?.toFixed(2) ?? 'N/A'} kcal/mol`,
    ].filter(Boolean).join('\n')
    try {
      const res = await fetch('/brain/chat', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          message: `You are a drug discovery assistant helping a student. Explain these molecular docking results in clear, educational terms:\n\n${scoreLines}\n\nCover: what these scores mean, whether this is a strong/moderate/weak binder, and what the student's next steps should be.`,
          conversation_id: aiConvId || undefined,
        }),
      })
      const data = await res.json()
      setAiExplanation(data.response || 'No explanation returned.')
      if (data.conversation_id) setAiConvId(data.conversation_id)
    } catch {
      setAiExplanation('AI explanation unavailable. Check that the brain service is running.')
    }
    setAiLoading(false)
  }

  const handleDownloadFile = (url: string, filename: string) => {
    const a = document.createElement('a')
    a.href = url
    a.download = filename
    a.target = '_blank'
    a.click()
  }

  const filteredJobs = jobs.filter(job => {
    if (activeTab === 'all') return true
    if (activeTab === 'completed') return job.status === 'completed'
    if (activeTab === 'running') return job.status === 'running' || job.status === 'pending'
    if (activeTab === 'failed') return job.status === 'failed' || job.status === 'cancelled'
    return true
  })

  const getStatusBadge = (status: string) => {
    const baseClass = 'px-2 py-0.5 rounded-full text-xs font-medium'
    switch (status) {
      case 'completed':
        return `${baseClass} ${isDark ? 'bg-green-900 text-green-300' : 'bg-green-100 text-green-700'}`
      case 'running':
      case 'pending':
        return `${baseClass} ${isDark ? 'bg-blue-900 text-blue-300' : 'bg-blue-100 text-blue-700'}`
      case 'failed':
      case 'cancelled':
        return `${baseClass} ${isDark ? 'bg-red-900 text-red-300' : 'bg-red-100 text-red-700'}`
      default:
        return `${baseClass} ${isDark ? 'bg-gray-700 text-gray-300' : 'bg-gray-100 text-gray-600'}`
    }
  }

  const getStatusDot = (status: string) => {
    const baseClass = 'w-2 h-2 rounded-full'
    switch (status) {
      case 'completed':
        return `${baseClass} bg-green-500`
      case 'running':
      case 'pending':
        return `${baseClass} bg-blue-500 animate-pulse`
      case 'failed':
      case 'cancelled':
        return `${baseClass} bg-red-500`
      default:
        return `${baseClass} bg-gray-400`
    }
  }

  return (
    <div className={`h-full flex ${bgClass}`}>
      <div className={`w-80 ${cardBg} ${borderClass} border-r overflow-hidden flex flex-col`}>
        <div className={`p-4 ${borderClass} border-b`}>
          <div className="flex items-center justify-between mb-3">
            <h2 className={`font-semibold ${textClass}`}>Job History</h2>
            <button onClick={fetchJobs} className={`p-1.5 rounded-lg ${isDark ? 'hover:bg-gray-700' : 'hover:bg-gray-100'}`} title="Refresh">
              <svg className={`w-4 h-4 ${subtextClass}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
              </svg>
            </button>
          </div>
          <div className="flex gap-1 flex-wrap">
            {(['all', 'completed', 'running', 'failed'] as const).map(tab => (
              <button
                key={tab}
                onClick={() => setActiveTab(tab)}
                className={`px-2.5 py-1 text-xs font-medium rounded-lg transition-colors ${
                  activeTab === tab
                    ? isDark ? 'bg-blue-600 text-white' : 'bg-blue-500 text-white'
                    : isDark ? 'text-gray-400 hover:bg-gray-700' : 'text-gray-600 hover:bg-gray-100'
                }`}
              >
                {tab.charAt(0).toUpperCase() + tab.slice(1)}
              </button>
            ))}
          </div>
        </div>

        <div className="flex-1 overflow-y-auto">
          {loading ? (
            <div className="flex items-center justify-center h-32">
              <div className="animate-spin text-2xl">⟳</div>
            </div>
          ) : filteredJobs.length === 0 ? (
            <div className={`text-center py-12 ${subtextClass}`}>
              <div className="text-4xl mb-3">📋</div>
              <p className="text-sm">No jobs found</p>
              <p className="text-xs mt-1">Run a docking simulation to see history</p>
            </div>
          ) : (
            <div className="divide-y divide-gray-200 dark:divide-gray-700">
              {filteredJobs.map(job => (
                <button
                  key={job.job_uuid}
                  onClick={() => selectJob(job)}
                  className={`w-full p-4 text-left transition-colors ${
                    selectedJob?.job_uuid === job.job_uuid
                      ? isDark ? 'bg-blue-900/30 border-l-4 border-blue-500' : 'bg-blue-50 border-l-4 border-blue-500'
                      : `${hoverClass}`
                  }`}
                >
                  <div className="flex items-start justify-between">
                    <div className="flex-1 min-w-0">
                      <p className={`font-medium text-sm truncate ${textClass}`}>{job.job_name}</p>
                      <p className={`text-xs ${subtextClass} mt-1`}>
                        {new Date(job.created_at).toLocaleDateString()} {new Date(job.created_at).toLocaleTimeString()}
                      </p>
                    </div>
                    <div className="flex items-center gap-2 ml-2">
                      <span className={getStatusDot(job.status)} />
                    </div>
                  </div>
                  <div className="flex items-center gap-2 mt-2">
                    {getStatusBadge(job.status)}
                    {job.binding_energy && (
                      <span className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                        {job.binding_energy.toFixed(2)} kcal/mol
                      </span>
                    )}
                  </div>
                  {job.engine && (
                    <p className={`text-xs ${subtextClass} mt-1`}>Engine: {job.engine.toUpperCase()}</p>
                  )}
                </button>
              ))}
            </div>
          )}
        </div>
        <div className={`p-3 ${borderClass} border-t ${subtextClass} text-xs text-center`}>
          {jobs.length} total | {jobs.filter(j => j.status === 'completed').length} completed
        </div>
      </div>

      <div className="flex-1 overflow-y-auto">
        {selectedJob ? (
          <div className="p-6">
            <div className={`${cardBg} rounded-xl ${borderClass} border mb-6`}>
              <div className={`p-6 ${borderClass} border-b`}>
                <div className="flex items-start justify-between">
                  <div>
                    <h1 className={`text-2xl font-bold ${textClass}`}>{selectedJob.job_name}</h1>
                    <p className={`mt-1 ${subtextClass}`}>
                      Job ID: <code className={`${isDark ? 'bg-gray-700' : 'bg-gray-100'} px-1 rounded`}>{selectedJob.job_uuid}</code>
                    </p>
                  </div>
                  <div className="flex items-center gap-3">
                    {getStatusBadge(selectedJob.status)}
                    <button onClick={() => deleteJob(selectedJob.job_uuid)} className={`p-2 rounded-lg ${isDark ? 'hover:bg-red-900/50' : 'hover:bg-red-100'}`} title="Delete Job">
                      <svg className={`w-5 h-5 ${isDark ? 'text-red-400' : 'text-red-500'}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 7l-.867 12.142A2 2 0 0116.138 21H7.862a2 2 0 01-1.995-1.858L5 7m5 4v6m4-6v6m1-10V4a1 1 0 00-1-1h-4a1 1 0 00-1 1v3M4 7h16" />
                      </svg>
                    </button>
                  </div>
                </div>
              </div>
              <div className="p-6 grid grid-cols-2 md:grid-cols-4 gap-4">
                <div>
                  <p className={`text-xs ${subtextClass} uppercase tracking-wide`}>Created</p>
                  <p className={`font-medium ${textClass} mt-1`}>{new Date(selectedJob.created_at).toLocaleDateString()}</p>
                  <p className={`text-sm ${subtextClass}`}>{new Date(selectedJob.created_at).toLocaleTimeString()}</p>
                </div>
                {selectedJob.completed_at && (
                  <div>
                    <p className={`text-xs ${subtextClass} uppercase tracking-wide`}>Completed</p>
                    <p className={`font-medium ${textClass} mt-1`}>{new Date(selectedJob.completed_at).toLocaleDateString()}</p>
                    <p className={`text-sm ${subtextClass}`}>{new Date(selectedJob.completed_at).toLocaleTimeString()}</p>
                  </div>
                )}
                <div>
                  <p className={`text-xs ${subtextClass} uppercase tracking-wide`}>Engine</p>
                  <p className={`font-medium ${textClass} mt-1`}>{selectedJob.engine?.toUpperCase() || 'VINA'}</p>
                </div>
                <div>
                  <p className={`text-xs ${subtextClass} uppercase tracking-wide`}>Best Score</p>
                  <p className={`font-bold text-lg mt-1 ${selectedJob.binding_energy ? 'text-green-600' : subtextClass}`}>
                    {selectedJob.binding_energy ? `${selectedJob.binding_energy.toFixed(2)} kcal/mol` : '-'}
                  </p>
                </div>
              </div>
            </div>

            {selectedResults.length > 0 && (
              <>
                <div className={`${cardBg} rounded-xl ${borderClass} border mb-6`}>
                  <div className={`p-4 ${borderClass} border-b`}>
                    <h2 className={`font-semibold ${textClass}`}>Docking Poses — {selectedResults.length} pose(s)</h2>
                  </div>
                  <div className="overflow-x-auto">
                    <table className="w-full">
                      <thead className={isDark ? 'bg-gray-700/50' : 'bg-gray-50'}>
                        <tr>
                          <th className={`px-6 py-3 text-left text-xs font-medium ${subtextClass} uppercase`}>#</th>
                          <th className={`px-6 py-3 text-left text-xs font-medium ${subtextClass} uppercase`}>Pose</th>
                          <th className={`px-6 py-3 text-left text-xs font-medium ${subtextClass} uppercase`}>Vina Score</th>
                          <th className={`px-6 py-3 text-left text-xs font-medium ${subtextClass} uppercase`}>CNN Score</th>
                          <th className={`px-6 py-3 text-left text-xs font-medium ${subtextClass} uppercase`}>RF Score</th>
                          <th className={`px-6 py-3 text-left text-xs font-medium ${subtextClass} uppercase`}>Consensus</th>
                          <th className={`px-6 py-3 text-left text-xs font-medium ${subtextClass} uppercase`}>Composite</th>
                        </tr>
                      </thead>
                      <tbody className={`divide-y ${borderClass}`}>
                        {selectedResults.map((result, i) => (
                          <tr
                            key={result.id}
                            className={`transition-colors cursor-pointer ${selectedPose?.id === result.id ? (isDark ? 'bg-blue-900/30' : 'bg-blue-50') : ''} ${i === 0 ? (isDark ? 'bg-green-900/20' : 'bg-green-50') : ''}`}
                            onClick={() => setSelectedPose(result)}
                          >
                            <td className="px-6 py-4">
                              <span className={`inline-flex items-center justify-center w-6 h-6 rounded-full text-xs font-bold ${i === 0 ? 'bg-green-500 text-white' : isDark ? 'bg-gray-600 text-white' : 'bg-gray-200 text-gray-700'}`}>
                                {i + 1}
                              </span>
                            </td>
                            <td className={`px-6 py-4 text-sm ${textClass}`}>{result.ligand_name || `Pose ${result.pose_id}`}</td>
                            <td className="px-6 py-4">
                              <span className={`font-bold ${(result.vina_score || 0) < -8 ? 'text-green-600' : (result.vina_score || 0) < -6 ? 'text-yellow-600' : isDark ? 'text-gray-300' : 'text-gray-700'}`}>
                                {result.vina_score?.toFixed(2) || '-'}
                              </span>
                            </td>
                            <td className={`px-6 py-4 text-sm ${result.gnina_score ? 'text-purple-600' : subtextClass}`}>
                              {result.gnina_score?.toFixed(2) || '-'}
                            </td>
                            <td className={`px-6 py-4 text-sm ${result.rf_score ? 'text-blue-600' : subtextClass}`}>
                              {result.rf_score?.toFixed(2) || '-'}
                            </td>
                            <td className={`px-6 py-4 text-sm ${result.consensus ? 'text-cyan-600' : subtextClass}`}>
                              {result.consensus?.toFixed(2) || '-'}
                            </td>
                            <td className="px-6 py-4">
                              <span className={`font-bold text-sm ${(result.composite_score ?? result.final_score ?? 0) < -8 ? 'text-green-600' : (result.composite_score ?? result.final_score ?? 0) < -6 ? 'text-yellow-600' : isDark ? 'text-gray-300' : 'text-gray-700'}`}>
                                {(result.composite_score ?? result.final_score ?? result.vina_score)?.toFixed(2) || '-'}
                              </span>
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>

                {Object.keys(jobFiles).length > 0 && (
                  <div className={`${cardBg} rounded-xl ${borderClass} border mb-6`}>
                    <div className={`p-4 ${borderClass} border-b flex items-center justify-between`}>
                      <h2 className={`font-semibold ${textClass}`}>📦 Download All Files</h2>
                      {logText && (
                        <button
                          onClick={() => {
                            const blob = new Blob([logText], { type: 'text/plain' })
                            handleDownloadFile(URL.createObjectURL(blob), `log_${selectedJob?.job_uuid}.txt`)
                          }}
                          className={`text-xs px-2.5 py-1 rounded font-medium ${isDark ? 'bg-gray-700 hover:bg-gray-600 text-white' : 'bg-gray-200 hover:bg-gray-300 text-gray-800'}`}
                        >
                          ⬇ Log Text
                        </button>
                      )}
                    </div>
                    <div className="p-4 flex flex-wrap gap-2">
                      {(Object.entries(jobFiles) as [string, { filename: string; url: string | null; size_bytes?: number; exists: boolean }][]).map(([key, info]) =>
                        info.exists && info.url ? (
                          <a
                            key={key}
                            href={info.url}
                            download={info.filename}
                            className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-medium border transition-colors ${
                              isDark
                                ? 'bg-gray-700 hover:bg-gray-600 border-gray-600 text-gray-200'
                                : 'bg-gray-50 hover:bg-gray-100 border-gray-200 text-gray-700'
                            }`}
                          >
                            <span>{key === 'docking' ? '🧬' : key === 'log' ? '📄' : key === 'receptor' ? '🔵' : key === 'ligand' ? '🟢' : '📁'}</span>
                            <span className="capitalize">{key}</span>
                            <span className={`${subtextClass}`}>({info.filename.split('.').pop()?.toUpperCase()})</span>
                            {info.size_bytes && <span className={subtextClass}>· {(info.size_bytes / 1024).toFixed(1)}KB</span>}
                          </a>
                        ) : null
                      )}
                    </div>
                    {logText && (
                      <details className={`border-t ${borderClass}`}>
                        <summary className={`px-4 py-2 cursor-pointer text-xs font-medium ${subtextClass} hover:opacity-80`}>View Docking Log</summary>
                        <pre className={`px-4 pb-4 text-xs overflow-x-auto whitespace-pre-wrap ${isDark ? 'text-green-300' : 'text-gray-600'}`}>{logText}</pre>
                      </details>
                    )}
                  </div>
                )}

                {selectedPose && (
                  <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                    <div className={`${cardBg} rounded-xl ${borderClass} border`}>
                      <div className={`p-4 ${borderClass} border-b flex items-center justify-between`}>
                        <h2 className={`font-semibold ${textClass}`}>3D Viewer — {selectedPose.ligand_name || `Pose ${selectedPose.pose_id}`}</h2>
                        <div className="flex gap-2 flex-wrap">
                          {selectedPose.pdb_data && (
                            <button
                              onClick={() => handleDownloadPDB(selectedPose)}
                              className={`text-xs px-2.5 py-1 rounded font-medium ${isDark ? 'bg-blue-700 hover:bg-blue-600 text-white' : 'bg-blue-500 hover:bg-blue-600 text-white'}`}
                            >
                              ⬇ PDB
                            </button>
                          )}
                          <button
                            onClick={handleSnapshot}
                            className={`text-xs px-2.5 py-1 rounded font-medium ${isDark ? 'bg-purple-700 hover:bg-purple-600 text-white' : 'bg-purple-500 hover:bg-purple-600 text-white'}`}
                          >
                            📷 Snapshot
                          </button>
                        </div>
                      </div>
                      <div
                        ref={viewer3dRef}
                        style={{ height: '400px', position: 'relative', background: '#1a1a2e' }}
                      />
                      {!selectedPose.pdb_data && (
                        <div className={`p-6 text-center text-sm ${subtextClass}`}>
                          No PDBQT data available for this pose
                        </div>
                      )}
                    </div>

                    <div className="space-y-6">
                      <div className={`${cardBg} rounded-xl ${borderClass} border`}>
                        <div className={`p-4 ${borderClass} border-b`}>
                          <h2 className={`font-semibold ${textClass}`}>Scoring Breakdown</h2>
                          <span className={`text-sm ${subtextClass}`}>Vina Score</span>
                          <span className={`font-mono font-bold ${(selectedPose.vina_score || 0) < -8 ? 'text-green-600' : (selectedPose.vina_score || 0) < -6 ? 'text-yellow-600' : isDark ? 'text-gray-300' : 'text-gray-700'}`}>
                            {selectedPose.vina_score?.toFixed(2) || '-'} kcal/mol
                          </span>
                        </div>
                        {selectedPose.gnina_score != null && (
                          <div className="flex justify-between items-center">
                            <span className={`text-sm ${subtextClass}`}>GNINA CNN</span>
                            <span className="font-mono text-purple-600 font-bold">{selectedPose.gnina_score.toFixed(2)}</span>
                          </div>
                        )}
                        {selectedPose.rf_score != null && (
                          <div className="flex justify-between items-center">
                            <span className={`text-sm ${subtextClass}`}>RF Score</span>
                            <span className="font-mono text-blue-600 font-bold">{selectedPose.rf_score.toFixed(2)}</span>
                          </div>
                        )}
                        {selectedPose.hydrophobic_term != null && (
                          <div className="flex justify-between items-center">
                            <span className={`text-sm ${subtextClass}`}>Hydrophobic</span>
                            <span className="font-mono font-bold">{selectedPose.hydrophobic_term.toFixed(3)}</span>
                          </div>
                        )}
                        {selectedPose.rotatable_penalty != null && (
                          <div className="flex justify-between items-center">
                            <span className={`text-sm ${subtextClass}`}>Rotatable Penalty</span>
                            <span className="font-mono font-bold">{selectedPose.rotatable_penalty.toFixed(3)}</span>
                          </div>
                        )}
                        {selectedPose.constraint_penalty != null && (
                          <div className="flex justify-between items-center">
                            <span className={`text-sm ${subtextClass}`}>Constraint Penalty</span>
                            <span className="font-mono font-bold text-orange-500">{selectedPose.constraint_penalty.toFixed(3)}</span>
                          </div>
                        )}
                        <div className={`border-t pt-2 flex justify-between items-center ${isDark ? 'border-gray-700' : 'border-gray-200'}`}>
                          <span className={`font-semibold ${textClass}`}>Composite Score</span>
                          <span className={`font-mono font-bold text-lg ${(selectedPose.composite_score ?? selectedPose.final_score ?? 0) < -8 ? 'text-green-600' : (selectedPose.composite_score ?? selectedPose.final_score ?? 0) < -6 ? 'text-yellow-600' : isDark ? 'text-gray-300' : 'text-gray-700'}`}>
                            {(selectedPose.composite_score ?? selectedPose.final_score ?? selectedPose.vina_score)?.toFixed(2) || '-'} kcal/mol
                          </span>
                        </div>
                        <div className={`mt-3 rounded-lg p-3 ${isDark ? 'bg-gray-700/50' : 'bg-gray-50'}`}>
                          <div className="space-y-2">
                            <div className="flex items-center gap-2">
                              <span className={`w-3 h-3 rounded-full ${(selectedPose.vina_score || 0) < -8 ? 'bg-green-500' : (selectedPose.vina_score || 0) < -6 ? 'bg-yellow-500' : 'bg-red-500'}`} />
                              <span className={`text-sm ${textClass}`}>
                                {(selectedPose.vina_score || 0) < -8 ? 'Strong binding affinity' : (selectedPose.vina_score || 0) < -6 ? 'Moderate binding affinity' : 'Weak binding affinity'}
                              </span>
                            </div>
                            {selectedPose.gnina_score != null && (
                              <div className="flex items-center gap-2">
                                <span className={`w-3 h-3 rounded-full ${selectedPose.gnina_score > 0.5 ? 'bg-green-500' : selectedPose.gnina_score > 0 ? 'bg-yellow-500' : 'bg-red-500'}`} />
                                <span className={`text-sm ${textClass}`}>
                                  CNN quality: {selectedPose.gnina_score > 0.5 ? 'High' : selectedPose.gnina_score > 0 ? 'Moderate' : 'Low'}
                                </span>
                              </div>
                            )}
                          </div>
                        </div>
                      </div>

                    {/* AI Explanation Panel */}
                    <div className={`${cardBg} rounded-xl ${borderClass} border`}>
                      <div className={`p-4 ${borderClass} border-b flex items-center justify-between`}>
                        <div className="flex items-center gap-2">
                          <span className="text-lg">🤖</span>
                          <h2 className={`font-semibold ${textClass}`}>AI Result Explanation</h2>
                        </div>
                        <button
                          onClick={askAI}
                          disabled={aiLoading}
                          className={`flex items-center gap-1.5 px-3 py-1.5 rounded-lg text-xs font-medium transition-colors disabled:opacity-50 ${
                            isDark ? 'bg-purple-700 hover:bg-purple-600 text-white' : 'bg-purple-600 hover:bg-purple-700 text-white'
                          }`}
                        >
                          {aiLoading ? (
                            <>
                              <svg className="animate-spin w-3 h-3" fill="none" viewBox="0 0 24 24">
                                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z" />
                              </svg>
                              Thinking…
                            </>
                          ) : (
                            <>
                              <svg className="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
                              </svg>
                              {aiExplanation ? 'Re-explain' : 'Explain with AI'}
                            </>
                          )}
                        </button>
                      </div>
                      <div className="p-4">
                        {!aiExplanation && !aiLoading && (
                          <p className={`text-sm ${subtextClass}`}>
                            Click <strong>Explain with AI</strong> to get an educational breakdown of these docking scores.
                          </p>
                        )}
                        {aiLoading && (
                          <div className="flex items-center gap-2">
                            <div className="animate-pulse w-2 h-2 rounded-full bg-purple-500" />
                            <span className={`text-sm ${subtextClass}`}>Generating explanation…</span>
                          </div>
                        )}
                        {aiExplanation && (
                          <div className={`text-sm leading-relaxed whitespace-pre-wrap ${textClass}`}>
                            {aiExplanation}
                          </div>
                        )}
                      </div>
                    </div>
                    </div>
                  </div>
                )}
              </>
            )}

            {selectedResults.length === 0 && selectedJob.status === 'completed' && (
              <div className={`${cardBg} rounded-xl ${borderClass} border p-12 text-center`}>
                <div className="text-4xl mb-3">📊</div>
                <p className={subtextClass}>No results available for this job</p>
              </div>
            )}
          </div>
        ) : (
          <div className="flex items-center justify-center h-full">
            <div className="text-center">
              <svg className={`w-20 h-20 mx-auto ${isDark ? 'text-gray-600' : 'text-gray-300'}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1} d="M9 17v-2m3 2v-4m3 4v-6m2 10H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
              </svg>
              <p className={`mt-4 ${subtextClass}`}>Select a job from the history to view results</p>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
