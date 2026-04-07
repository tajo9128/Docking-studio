import { useState, useEffect } from 'react'
import { useTheme } from '@/contexts/ThemeContext'
import { Link } from 'react-router-dom'

const PC = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'

interface PcResult {
  cid: number
  name: string
  formula: string
  mw: string
  smiles: string
}

export function Docking() {
  const { theme } = useTheme()
  const isDark = theme === 'dark'
  const [receptorFile, setReceptorFile] = useState<File | null>(null)
  const [ligandFile, setLigandFile] = useState<File | null>(null)
  const [center, setCenter] = useState({ x: 0, y: 0, z: 0 })
  const [size, setSize] = useState({ x: 20, y: 20, z: 20 })
  const [exhaustiveness, setExhaustiveness] = useState(8)
  const [running, setRunning] = useState(false)
  const [result, setResult] = useState<any>(null)
  const [error, setError] = useState('')
  const [jobs, setJobs] = useState<any[]>([])
  const [jobFiles, setJobFiles] = useState<Record<string, { filename: string; url: string | null; size_bytes?: number; exists: boolean }>>({})

  // Progress tracking state
  const [progress, setProgress] = useState(0)
  const [progressMessage, setProgressMessage] = useState('')
  const [currentStep, setCurrentStep] = useState(0)

  // PubChem lookup
  const [pcOpen, setPcOpen]     = useState(false)
  const [pcQuery, setPcQuery]   = useState('')
  const [pcBusy, setPcBusy]     = useState(false)
  const [pcResult, setPcResult] = useState<PcResult | null>(null)
  const [pcErr, setPcErr]       = useState<string | null>(null)
  const [copied, setCopied]     = useState(false)

  const searchPubChem = async () => {
    const q = pcQuery.trim()
    if (!q) return
    setPcBusy(true); setPcErr(null); setPcResult(null)
    try {
      const cidRes = await fetch(`${PC}/compound/name/${encodeURIComponent(q)}/cids/JSON`)
      if (!cidRes.ok) throw new Error('Compound not found')
      const cid: number = (await cidRes.json()).IdentifierList.CID[0]
      const propRes = await fetch(
        `${PC}/compound/cid/${cid}/property/IUPACName,MolecularFormula,MolecularWeight,IsomericSMILES/JSON`
      )
      const p = (await propRes.json()).PropertyTable.Properties[0]
      setPcResult({ cid, name: p.IUPACName ?? q, formula: p.MolecularFormula ?? '', mw: String(p.MolecularWeight ?? ''), smiles: p.IsomericSMILES ?? '' })
    } catch (e: any) { setPcErr(e.message ?? 'Not found') }
    finally { setPcBusy(false) }
  }

  const copySMILES = () => {
    if (!pcResult) return
    navigator.clipboard.writeText(pcResult.smiles)
    setCopied(true); setTimeout(() => setCopied(false), 2000)
  }

  useEffect(() => { fetchJobs() }, [])

  const fetchJobs = async () => {
    try {
      const res = await fetch('/jobs')
      const data = await res.json()
      setJobs((data.jobs || []).filter((j: any) => j.status === 'completed').slice(0, 10))
    } catch { /* no-op */ }
  }

  const runDocking = async () => {
    if (!receptorFile || !ligandFile) {
      setError('Upload both receptor and ligand files')
      return
    }
    const recExt = receptorFile.name.split('.').pop()?.toLowerCase()
    const ligExt = ligandFile.name.split('.').pop()?.toLowerCase()
    const recOk = ['pdb', 'pdbqt', 'ent'].includes(recExt ?? '')
    const ligOk = ['pdb', 'pdbqt', 'sdf', 'mol2', 'ent'].includes(ligExt ?? '')
    if (!recOk) { setError(`Receptor format ".${recExt}" not supported. Use .pdb, .pdbqt, or .ent`); return }
    if (!ligOk) { setError(`Ligand format ".${ligExt}" not supported. Use .pdb, .pdbqt, .sdf, .mol2, or .ent`); return }
    setRunning(true)
    setError('')
    setResult(null)
    
    // Reset progress state
    setProgress(0)
    setProgressMessage('Initialising...')
    setCurrentStep(0)
    
    try {
      // Read file contents
      const receptorContent = await readFileContent(receptorFile)
      const ligandContent = await readFileContent(ligandFile)
      
      // Start docking via API
      const response = await fetch('/api/docking/run', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          receptor_content: receptorContent,
          ligand_content: ligandContent,
          receptor_filename: receptorFile.name,
          ligand_filename: ligandFile.name,
          center_x: center.x,
          center_y: center.y,
          center_z: center.z,
          size_x: size.x,
          size_y: size.y,
          size_z: size.z,
          exhaustiveness,
          scoring: 'vina'
        })
      })
      
      const accepted = await response.json()
      if (accepted.error) throw new Error(accepted.error)
      
      const jobId = accepted.job_id
      
      // Connect to SSE stream for real-time progress
      const eventSource = new EventSource(`/dock/${jobId}/stream`)
      
      await new Promise<void>((resolve, reject) => {
        eventSource.onmessage = (event) => {
          try {
            const data = JSON.parse(event.data)
            setProgress(data.progress)
            setProgressMessage(data.message)
            setCurrentStep(getStepFromProgress(data.progress))
            
            if (data.status === 'completed') {
              eventSource.close()
              resolve()
            } else if (data.status === 'failed') {
              eventSource.close()
              reject(new Error(data.message || 'Docking failed'))
            }
          } catch (e) {
            // Ignore parse errors
          }
        }
        
        eventSource.onerror = () => {
          eventSource.close()
          // Fall back to polling if SSE fails
          pollUntilComplete(jobId).then(resolve).catch(reject)
        }
      })
      
      // Fetch final results
      const resultRes = await fetch(`/api/docking/result/${jobId}`)
      const data = await resultRes.json()
      
      if (data.error) setError(data.error)
      else { 
        setResult(data)
        setJobFiles({})
        // Fetch downloadable files
        try {
          const filesRes = await fetch(`/jobs/${jobId}/files`)
          if (filesRes.ok) {
            const filesData = await filesRes.json()
            setJobFiles(filesData.files || {})
          }
        } catch { /* files may not exist yet */ }
        setJobs(prev => [{ 
          job_uuid: jobId, 
          job_name: `Docking ${new Date().toLocaleTimeString()}`,
          status: 'completed',
          binding_energy: data.best_score,
          engine: data.engine
        }, ...prev])
      }
    } catch (e: any) {
      setError(e.message || 'Docking failed')
    }
    
    setRunning(false)
  }
  
  const pollUntilComplete = async (jobId: string): Promise<void> => {
    for (let attempt = 0; attempt < 300; attempt++) {
      await new Promise(r => setTimeout(r, 2000))
      const res = await fetch(`/dock/${jobId}/status`)
      const status = await res.json()
      
      setProgress(status.progress)
      setProgressMessage(status.message)
      setCurrentStep(getStepFromProgress(status.progress))
      
      if (status.status === 'completed') return
      if (status.status === 'failed') throw new Error(status.message)
    }
    throw new Error('Docking timed out')
  }
  
  const getStepFromProgress = (pct: number): number => {
    if (pct < 5) return 0
    if (pct < 25) return 1  // Protein prep
    if (pct < 45) return 2  // Ligand prep
    if (pct < 50) return 3  // Grid config
    if (pct < 85) return 4  // Docking
    if (pct < 100) return 5 // Processing
    return 6 // Complete
  }
  
  const readFileContent = (file: File): Promise<string> => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = e => resolve(e.target?.result as string)
      reader.onerror = reject
      reader.readAsText(file)
    })
  }

  return (
    <div className={`h-full flex ${isDark ? 'bg-gray-900 text-white' : 'bg-gray-50 text-gray-900'}`}>
      {/* Left - Form */}
      <div className="w-96 p-6 border-r overflow-y-auto" style={{ borderColor: isDark ? '#374151' : '#e5e7eb' }}>
        <h1 className="text-xl font-bold mb-1">AutoDock Vina</h1>
        <p className={`text-xs mb-4 ${isDark ? 'text-gray-500' : 'text-gray-400'}`}>
          Upload prepared PDBQT files. Use AutoDock Tools or other software to prepare files.
        </p>

        {/* ── PubChem Lookup ── */}
        <div className={`mb-4 rounded-lg border ${isDark ? 'border-blue-800 bg-blue-950/40' : 'border-blue-200 bg-blue-50'}`}>
          <button
            onClick={() => setPcOpen(o => !o)}
            className={`w-full flex items-center justify-between px-3 py-2.5 text-sm font-medium ${
              isDark ? 'text-blue-300' : 'text-blue-700'
            }`}
          >
            <span>🔬 Find ligand on PubChem by name</span>
            <span>{pcOpen ? '▲' : '▼'}</span>
          </button>

          {pcOpen && (
            <div className="px-3 pb-3 space-y-2">
              <div className="flex gap-2">
                <input
                  className={`flex-1 px-2.5 py-1.5 rounded border text-sm focus:outline-none focus:ring-2 focus:ring-blue-500 ${
                    isDark ? 'bg-gray-800 border-gray-600 text-white placeholder-gray-500' : 'bg-white border-gray-300 text-gray-900 placeholder-gray-400'
                  }`}
                  placeholder="e.g. Ibuprofen, Aspirin, Caffeine…"
                  value={pcQuery}
                  onChange={e => setPcQuery(e.target.value)}
                  onKeyDown={e => e.key === 'Enter' && searchPubChem()}
                />
                <button
                  onClick={searchPubChem} disabled={pcBusy}
                  className="px-3 py-1.5 bg-blue-600 hover:bg-blue-700 text-white text-sm rounded disabled:opacity-50"
                >
                  {pcBusy ? '…' : 'Search'}
                </button>
              </div>

              {pcErr && <p className="text-xs text-red-500">{pcErr}</p>}

              {pcResult && (
                <div className={`rounded-lg p-2.5 flex gap-3 ${
                  isDark ? 'bg-gray-800' : 'bg-white border border-gray-200'
                }`}>
                  <img
                    src={`${PC}/compound/cid/${pcResult.cid}/PNG?image_size=small`}
                    alt={pcResult.name}
                    className="w-20 h-20 object-contain rounded bg-white shrink-0"
                  />
                  <div className="flex-1 min-w-0 space-y-1">
                    <p className={`text-xs font-semibold truncate ${
                      isDark ? 'text-gray-200' : 'text-gray-800'
                    }`}>{pcResult.name}</p>
                    <p className={`text-xs ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>
                      {pcResult.formula} · {pcResult.mw} g/mol
                    </p>
                    <p className={`text-xs font-mono break-all leading-tight ${
                      isDark ? 'text-green-400' : 'text-green-700'
                    }`}>{pcResult.smiles}</p>
                    <div className="flex gap-2 pt-1">
                      <button
                        onClick={copySMILES}
                        className={`text-xs px-2 py-1 rounded font-medium ${
                          copied
                            ? 'bg-green-600 text-white'
                            : isDark ? 'bg-gray-700 text-gray-200 hover:bg-gray-600' : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
                        }`}
                      >
                        {copied ? '✓ Copied!' : 'Copy SMILES'}
                      </button>
                      <a
                        href={`https://pubchem.ncbi.nlm.nih.gov/compound/${pcResult.cid}#section=3D-Conformer`}
                        target="_blank"
                        rel="noreferrer"
                        className="text-xs px-2 py-1 rounded font-medium bg-blue-600 text-white hover:bg-blue-700"
                      >
                        View on PubChem ↗
                      </a>
                    </div>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>

        <div className="space-y-4">
          <div>
            <label className={`text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Receptor (.pdb, .pdbqt, .ent)</label>
            <input type="file" accept=".pdb,.pdbqt,.ent" onChange={e => setReceptorFile(e.target.files?.[0] || null)}
              className={`w-full mt-1 p-2 rounded border text-sm ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-300'}`} />
            {receptorFile && <p className="text-xs text-green-500 mt-1">✓ {receptorFile.name}</p>}
          </div>

          <div>
            <label className={`text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Ligand (.pdb, .pdbqt, .sdf, .mol2, .ent)</label>
            <input type="file" accept=".pdb,.pdbqt,.sdf,.mol2,.ent" onChange={e => setLigandFile(e.target.files?.[0] || null)}
              className={`w-full mt-1 p-2 rounded border text-sm ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-300'}`} />
            {ligandFile && <p className="text-xs text-green-500 mt-1">✓ {ligandFile.name}</p>}
          </div>

          <div>
            <label className={`text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Grid Center (X, Y, Z)</label>
            <div className="grid grid-cols-3 gap-2 mt-1">
              {(['x', 'y', 'z'] as const).map(axis => (
                <input key={axis} type="number" value={center[axis]} onChange={e => setCenter(p => ({ ...p, [axis]: Number(e.target.value) }))}
                  className={`p-2 rounded border text-sm text-center ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-300'}`} />
              ))}
            </div>
          </div>

          <div>
            <label className={`text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Grid Size (X, Y, Z)</label>
            <div className="grid grid-cols-3 gap-2 mt-1">
              {(['x', 'y', 'z'] as const).map(axis => (
                <input key={axis} type="number" value={size[axis]} onChange={e => setSize(p => ({ ...p, [axis]: Number(e.target.value) }))}
                  className={`p-2 rounded border text-sm text-center ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-300'}`} />
              ))}
            </div>
          </div>

          <div>
            <label className={`text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Exhaustiveness: {exhaustiveness}</label>
            <input type="range" min="1" max="64" value={exhaustiveness} onChange={e => setExhaustiveness(Number(e.target.value))}
              className="w-full mt-1" />
          </div>

          <button onClick={runDocking} disabled={running || !receptorFile || !ligandFile}
            className="w-full py-2.5 bg-blue-600 text-white rounded-lg font-medium hover:bg-blue-700 disabled:opacity-50 transition-colors">
            {running ? 'Docking…' : 'Run Docking'}
          </button>
        </div>

        {error && <p className="mt-4 text-sm text-red-500">{error}</p>}

        {/* Recent jobs */}
        {jobs.length > 0 && (
          <div className="mt-6">
            <h3 className={`text-xs font-semibold mb-2 ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Recent Jobs</h3>
            <div className="space-y-1">
              {jobs.map((j: any) => (
                <div key={j.job_uuid} className={`text-xs p-2 rounded ${isDark ? 'bg-gray-800' : 'bg-white border'}`}>
                  <span className="font-medium">{j.job_name}</span>
                  <span className={`ml-2 ${j.binding_energy < -7 ? 'text-green-500' : 'text-gray-400'}`}>
                    {j.binding_energy?.toFixed(2)} kcal/mol
                  </span>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>

      {/* Right - Results */}
      <div className="flex-1 p-6 overflow-y-auto">
        {!result && !running && (
          <div className="flex items-center justify-center h-full">
            <div className="text-center">
              <div className="text-4xl mb-3">⬡</div>
              <p className={isDark ? 'text-gray-500' : 'text-gray-400'}>Upload PDBQT files and run docking</p>
            </div>
          </div>
        )}
        {running && (
          <div className="flex flex-col items-center justify-center h-full p-8">
            <div className="w-full max-w-md">
              {/* Progress header */}
              <div className="flex items-center justify-between mb-2">
                <span className={`text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
                  {currentStep === 0 && 'Initialising...'}
                  {currentStep === 1 && 'Preparing Protein'}
                  {currentStep === 2 && 'Preparing Ligand'}
                  {currentStep === 3 && 'Configuring Grid'}
                  {currentStep === 4 && 'Running Docking'}
                  {currentStep === 5 && 'Processing Results'}
                  {currentStep === 6 && 'Complete'}
                </span>
                <span className={`text-sm font-bold ${isDark ? 'text-blue-400' : 'text-blue-600'}`}>
                  {progress}%
                </span>
              </div>
              
              {/* Progress bar */}
              <div className={`w-full h-3 rounded-full overflow-hidden ${isDark ? 'bg-gray-700' : 'bg-gray-200'}`}>
                <div 
                  className="h-full bg-gradient-to-r from-blue-500 to-blue-600 transition-all duration-300 ease-out rounded-full"
                  style={{ width: `${progress}%` }}
                />
              </div>
              
              {/* Current message from backend */}
              <div className={`mt-4 p-3 rounded-lg border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'}`}>
                <div className="flex items-start gap-2">
                  <div className="animate-pulse w-2 h-2 rounded-full bg-blue-500 mt-1.5" />
                  <p className={`text-sm ${isDark ? 'text-gray-300' : 'text-gray-600'}`}>
                    {progressMessage || 'Initialising docking pipeline...'}
                  </p>
                </div>
              </div>
              
              {/* Steps indicator */}
              <div className="flex justify-between mt-4 px-1">
                {['Init', 'Protein', 'Ligand', 'Grid', 'Dock', 'Save'].map((label, idx) => (
                  <div key={label} className="flex flex-col items-center">
                    <div className={`w-2 h-2 rounded-full transition-colors ${
                      idx < currentStep ? 'bg-green-500' : 
                      idx === currentStep ? 'bg-blue-500 animate-pulse' : 
                      isDark ? 'bg-gray-600' : 'bg-gray-300'
                    }`} />
                    <span className={`text-xs mt-1 ${
                      idx <= currentStep ? (isDark ? 'text-gray-400' : 'text-gray-600') : 
                      (isDark ? 'text-gray-600' : 'text-gray-400')
                    }`}>{label}</span>
                  </div>
                ))}
              </div>
            </div>
          </div>
        )}
        {result && !running && (
          <div className="max-w-2xl">
            <div className="flex items-center justify-between mb-4">
              <h2 className="text-lg font-bold">Docking Complete</h2>
              <Link
                to="/results"
                className="flex items-center gap-1.5 px-3 py-1.5 bg-blue-600 hover:bg-blue-700 text-white text-sm rounded-lg font-medium transition-colors"
              >
                <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 12a3 3 0 11-6 0 3 3 0 016 0zM2.458 12C3.732 7.943 7.523 5 12 5c4.478 0 8.268 2.943 9.542 7-1.274 4.057-5.064 7-9.542 7-4.477 0-8.268-2.943-9.542-7z" />
                </svg>
                View Full Results & 3D Viewer
              </Link>
            </div>

            {/* Best score banner */}
            {result.best_score && (
              <div className={`mb-4 p-4 rounded-xl border-2 ${
                result.best_score < -8
                  ? isDark ? 'bg-green-900/30 border-green-600' : 'bg-green-50 border-green-400'
                  : result.best_score < -6
                  ? isDark ? 'bg-yellow-900/30 border-yellow-600' : 'bg-yellow-50 border-yellow-400'
                  : isDark ? 'bg-gray-800 border-gray-600' : 'bg-gray-50 border-gray-300'
              }`}>
                <p className={`text-xs uppercase tracking-wide mb-1 ${
                  isDark ? 'text-gray-400' : 'text-gray-500'
                }`}>Best Binding Affinity</p>
                <p className={`text-3xl font-bold font-mono ${
                  result.best_score < -8 ? 'text-green-500'
                  : result.best_score < -6 ? 'text-yellow-500'
                  : isDark ? 'text-gray-300' : 'text-gray-700'
                }`}>{result.best_score?.toFixed(2)} kcal/mol</p>
                <p className={`text-xs mt-1 ${
                  isDark ? 'text-gray-400' : 'text-gray-500'
                }`}>
                  {result.best_score < -8 ? '🟢 Strong binder — worth further validation'
                  : result.best_score < -6 ? '🟡 Moderate binder — consider optimisation'
                  : '🔴 Weak binder — structural modifications recommended'}
                </p>
              </div>
            )}

            {/* Poses table */}
            {result.results && result.results.length > 0 && (
              <div className={`mb-4 rounded-xl border overflow-hidden ${
                isDark ? 'border-gray-700' : 'border-gray-200'
              }`}>
                <div className={`px-4 py-2.5 text-xs font-semibold uppercase tracking-wide ${
                  isDark ? 'bg-gray-700 text-gray-300' : 'bg-gray-50 text-gray-500'
                }`}>
                  Docking Poses ({result.results.length})
                </div>
                <table className="w-full text-sm">
                  <thead className={isDark ? 'bg-gray-800/50' : 'bg-white'}>
                    <tr>
                      <th className={`px-4 py-2 text-left text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Pose</th>
                      <th className={`px-4 py-2 text-left text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Vina Score (kcal/mol)</th>
                      <th className={`px-4 py-2 text-left text-xs font-medium ${isDark ? 'text-gray-400' : 'text-gray-500'}`}>Rating</th>
                    </tr>
                  </thead>
                  <tbody className={`divide-y ${isDark ? 'divide-gray-700' : 'divide-gray-200'}`}>
                    {result.results.map((r: any, i: number) => (
                      <tr key={i} className={i === 0 ? isDark ? 'bg-green-900/20' : 'bg-green-50' : ''}>
                        <td className={`px-4 py-2.5 font-medium ${isDark ? 'text-gray-200' : 'text-gray-700'}`}>
                          <span className={`inline-flex items-center justify-center w-5 h-5 rounded-full text-xs font-bold mr-2 ${
                            i === 0 ? 'bg-green-500 text-white' : isDark ? 'bg-gray-600 text-white' : 'bg-gray-200 text-gray-700'
                          }`}>{i + 1}</span>
                          {r.mode ?? `Pose ${i + 1}`}
                        </td>
                        <td className={`px-4 py-2.5 font-mono font-bold ${
                          r.vina_score < -8 ? 'text-green-500' : r.vina_score < -6 ? 'text-yellow-500' : isDark ? 'text-gray-300' : 'text-gray-700'
                        }`}>{r.vina_score?.toFixed(2)}</td>
                        <td className={`px-4 py-2.5 text-xs ${
                          isDark ? 'text-gray-400' : 'text-gray-500'
                        }`}>
                          {r.vina_score < -8 ? '🟢 Strong' : r.vina_score < -6 ? '🟡 Moderate' : '🔴 Weak'}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}

            {/* Download section */}
            {Object.keys(jobFiles).length > 0 && (
              <div className={`rounded-xl border ${
                isDark ? 'border-gray-700 bg-gray-800' : 'border-gray-200 bg-white'
              }`}>
                <div className={`px-4 py-2.5 border-b text-sm font-semibold flex items-center gap-2 ${
                  isDark ? 'border-gray-700 text-white' : 'border-gray-200 text-gray-800'
                }`}>
                  <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
                  </svg>
                  Download Result Files
                </div>
                <div className="p-3 flex flex-wrap gap-2">
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
                        <span className={isDark ? 'text-gray-500' : 'text-gray-400'}>
                          ({info.filename.split('.').pop()?.toUpperCase()})
                        </span>
                        {info.size_bytes && (
                          <span className={isDark ? 'text-gray-500' : 'text-gray-400'}>
                            · {(info.size_bytes / 1024).toFixed(1)}KB
                          </span>
                        )}
                      </a>
                    ) : null
                  )}
                </div>
              </div>
            )}

            {!result.results?.length && (
              <p className={isDark ? 'text-gray-400' : 'text-gray-500'}>No pose data available</p>
            )}
          </div>
        )}
      </div>
    </div>
  )
}
