import { useState, useRef } from 'react'
import { useTheme } from '@/contexts/ThemeContext'

interface Config {
  center_x: number
  center_y: number
  center_z: number
  size_x: number
  size_y: number
  size_z: number
  exhaustiveness: number
  num_modes: number
}

interface DockingResult {
  mode: number
  vina_score: number
  gnina_score?: number
  rf_score?: number
  rmsd_lb?: number
  pdb_path?: string
}

interface DownloadUrls {
  log_file?: string
  docking_file?: string
  grid_file?: string
  vina_log?: string
  vina_docking?: string
  gnina_log?: string
  gnina_docking?: string
}

type StageStatus = 'pending' | 'running' | 'completed' | 'failed'

interface DockingStage {
  id: string
  label: string
  status: StageStatus
  progress: number
  message?: string
}

const INITIAL_STAGES: DockingStage[] = [
  { id: 'protein_prep', label: 'Protein Preparation', status: 'pending', progress: 0 },
  { id: 'ligand_prep', label: 'Ligand Preparation', status: 'pending', progress: 0 },
  { id: 'grid_generation', label: 'Grid Box Generation', status: 'pending', progress: 0 },
  { id: 'docking', label: 'Molecular Docking', status: 'pending', progress: 0 },
  { id: 'analysis', label: 'Results Analysis', status: 'pending', progress: 0 },
]

export function Docking() {
  const { theme } = useTheme()
  const isDark = theme === 'dark'
  const eventSourceRef = useRef<EventSource | null>(null)
  
  const [receptor, setReceptor] = useState<File | null>(null)
  const [receptorContent, setReceptorContent] = useState<string | null>(null)
  const [ligand, setLigand] = useState<File | null>(null)
  const [ligandContent, setLigandContent] = useState<string | null>(null)
  const [smiles, setSmiles] = useState('')
  const [config, setConfig] = useState<Config>({
    center_x: 0, center_y: 0, center_z: 0,
    size_x: 20, size_y: 20, size_z: 20,
    exhaustiveness: 32, num_modes: 10
  })
  const [scoringFunction, setScoringFunction] = useState<'vina' | 'gnina' | 'rf' | 'vinardo'>('vina')
  const [stages, setStages] = useState<DockingStage[]>(INITIAL_STAGES)
  const [overallProgress, setOverallProgress] = useState(0)
  const [isRunning, setIsRunning] = useState(false)
  const [error, setError] = useState('')
  const [results, setResults] = useState<DockingResult[]>([])
  const [_jobId, setJobId] = useState<string | null>(null)
  const [downloadUrls, setDownloadUrls] = useState<DownloadUrls | null>(null)
  const [routingDecision, setRoutingDecision] = useState<string>('')

  const resetStages = () => {
    setStages(INITIAL_STAGES.map(s => ({ ...s })))
    setOverallProgress(0)
    setResults([])
    setError('')
    setJobId(null)
    setDownloadUrls(null)
    setRoutingDecision('')
  }

  const updateStage = (id: string, updates: Partial<DockingStage>) => {
    setStages(prev => prev.map(s => s.id === id ? { ...s, ...updates } : s))
    recalculateOverall()
  }

  const recalculateOverall = () => {
    setStages(prev => {
      const total = prev.reduce((sum, s) => sum + s.progress, 0)
      setOverallProgress(Math.round(total / prev.length))
      return prev
    })
  }

  const completeStage = (id: string, message?: string) => {
    updateStage(id, { status: 'completed', progress: 100, message })
  }

  const failStage = (id: string, message?: string) => {
    updateStage(id, { status: 'failed', message })
    setError(message || 'Stage failed')
  }

  const runStage = (id: string, message?: string) => {
    updateStage(id, { status: 'running', message })
  }

  const readFileContent = (file: File): Promise<string> => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = e => resolve(e.target?.result as string)
      reader.onerror = reject
      reader.readAsText(file)
    })
  }

  const handleReceptorUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    
    try {
      const formData = new FormData()
      formData.append('file', file)
      
      const res = await fetch('/upload', {
        method: 'POST',
        body: formData,
      })
      const data = await res.json()
      
      if (data.path) {
        const content = await readFileContent(file)
        setReceptor(file)
        setReceptorContent(content)
        setError('')
      }
    } catch {
      const content = await readFileContent(file)
      setReceptor(file)
      setReceptorContent(content)
      setError('')
    }
  }

  const handleLigandUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    
    try {
      const formData = new FormData()
      formData.append('file', file)
      
      const res = await fetch('/upload', {
        method: 'POST',
        body: formData,
      })
      const data = await res.json()
      
      if (data.path) {
        const content = await readFileContent(file)
        setLigand(file)
        setLigandContent(content)
        setError('')
      }
    } catch {
      const content = await readFileContent(file)
      setLigand(file)
      setLigandContent(content)
      setError('')
    }
  }

  const handleDocking = async () => {
    resetStages()
    setIsRunning(true)
    setError('')

    try {
      if (smiles) {
        runStage('ligand_prep', 'Converting SMILES to 3D structure...')
        updateStage('ligand_prep', { progress: 25 })
        
        completeStage('protein_prep', 'Using default receptor')
        updateStage('ligand_prep', { progress: 50 })
        
        runStage('grid_generation', 'Setting up grid box...')
        updateStage('grid_generation', { progress: 50 })
        completeStage('grid_generation', `Grid: ${config.size_x}x${config.size_y}x${config.size_z} Å`)
        updateStage('ligand_prep', { progress: 75 })
        
        runStage('docking', 'Running AutoDock Vina...')
        
        const res = await fetch('/api/docking/run', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            smiles: smiles,
            receptor_content: receptorContent || '',
            ligand_content: '',
            center_x: config.center_x,
            center_y: config.center_y,
            center_z: config.center_z,
            size_x: config.size_x,
            size_y: config.size_y,
            size_z: config.size_z,
            exhaustiveness: config.exhaustiveness,
            num_modes: config.num_modes,
            scoring: scoringFunction,
          })
        })
        
        updateStage('ligand_prep', { progress: 100 })
        updateStage('docking', { progress: 50 })
        
        const data = await res.json()
        
        if (data.error) {
          failStage('docking', data.error)
          setError(data.error)
          setIsRunning(false)
          return
        }
        
        updateStage('docking', { progress: 100 })
        completeStage('docking', `Best score: ${data.best_score?.toFixed(2) || 'N/A'} kcal/mol`)
        completeStage('ligand_prep', 'SMILES converted to PDBQT')
        
        runStage('analysis', 'Analyzing binding modes...')
        updateStage('analysis', { progress: 50 })
        
        if (data.results && data.results.length > 0) {
          setResults(data.results)
        }
        
        if (data.download_urls) {
          setDownloadUrls(data.download_urls)
        }
        if (data.routing_decision) {
          setRoutingDecision(data.routing_decision)
        }
        
        updateStage('analysis', { progress: 100 })
        completeStage('analysis', `${data.results?.length || 0} poses generated`)
        
      } else if (receptorContent && ligandContent) {
        runStage('protein_prep', 'Preparing protein with RDKit...')
        updateStage('protein_prep', { progress: 25 })
        
        runStage('ligand_prep', 'Preparing ligand with RDKit...')
        updateStage('ligand_prep', { progress: 25 })
        
        runStage('grid_generation', 'Setting up grid box...')
        updateStage('grid_generation', { progress: 25 })
        
        completeStage('grid_generation', `Grid: ${config.size_x}x${config.size_y}x${config.size_z} Å`)
        updateStage('protein_prep', { progress: 50 })
        updateStage('ligand_prep', { progress: 50 })
        
        runStage('docking', 'Running AutoDock Vina...')
        
        const res = await fetch('/api/docking/run', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            smiles: '',
            receptor_content: receptorContent,
            ligand_content: ligandContent,
            center_x: config.center_x,
            center_y: config.center_y,
            center_z: config.center_z,
            size_x: config.size_x,
            size_y: config.size_y,
            size_z: config.size_z,
            exhaustiveness: config.exhaustiveness,
            num_modes: config.num_modes,
            scoring: scoringFunction,
          })
        })
        
        updateStage('protein_prep', { progress: 75 })
        updateStage('ligand_prep', { progress: 75 })
        updateStage('docking', { progress: 50 })
        
        const data = await res.json()
        
        if (data.error) {
          failStage('docking', data.error)
          setError(data.error)
          setIsRunning(false)
          return
        }
        
        updateStage('protein_prep', { progress: 100 })
        updateStage('ligand_prep', { progress: 100 })
        completeStage('protein_prep', 'Protein prepared')
        completeStage('ligand_prep', 'Ligand prepared to PDBQT')
        updateStage('docking', { progress: 100 })
        completeStage('docking', `Best: ${data.best_score?.toFixed(2) || 'N/A'} kcal/mol`)
        
        runStage('analysis', 'Analyzing binding modes...')
        updateStage('analysis', { progress: 50 })
        
        if (data.results && data.results.length > 0) {
          setResults(data.results)
        }
        
        if (data.download_urls) {
          setDownloadUrls(data.download_urls)
        }
        if (data.routing_decision) {
          setRoutingDecision(data.routing_decision)
        }
        
        updateStage('analysis', { progress: 100 })
        completeStage('analysis', `${data.results?.length || 0} poses generated`)
        
      } else {
        setError('Please upload receptor (PDB) + ligand (SDF/MOL2), OR enter a SMILES string')
        setIsRunning(false)
        return
      }
      
    } catch (err: any) {
      failStage('docking', 'Error: ' + err.message)
      setError(err.message)
    }

    setIsRunning(false)
  }

  const cancelDocking = () => {
    if (eventSourceRef.current) {
      eventSourceRef.current.close()
    }
    setIsRunning(false)
    setStages(prev => prev.map(s => s.status === 'running' ? { ...s, status: 'failed', message: 'Cancelled by user' } : s))
  }

  const handleConfigChange = (key: keyof Config, value: number) => {
    setConfig(prev => ({ ...prev, [key]: parseFloat(String(value)) }))
  }

  const getStageIcon = (status: StageStatus) => {
    switch (status) {
      case 'completed': return '✓'
      case 'running': return '⟳'
      case 'failed': return '✗'
      default: return '○'
    }
  }

  const getStageColor = (status: StageStatus, isDark: boolean) => {
    switch (status) {
      case 'completed': return isDark ? 'text-green-400' : 'text-green-600'
      case 'running': return isDark ? 'text-blue-400' : 'text-blue-600'
      case 'failed': return isDark ? 'text-red-400' : 'text-red-600'
      default: return isDark ? 'text-gray-500' : 'text-gray-400'
    }
  }

  const bgClass = isDark ? 'bg-gray-900' : 'bg-gray-50'
  const cardBg = isDark ? 'bg-gray-800' : 'bg-white'
  const borderClass = isDark ? 'border-gray-700' : 'border-gray-200'
  const textClass = isDark ? 'text-white' : 'text-gray-900'
  const subtextClass = isDark ? 'text-gray-400' : 'text-gray-500'
  const inputClass = isDark 
    ? 'bg-gray-700 border-gray-600 text-white focus:ring-blue-500'
    : 'bg-white border-gray-300 text-gray-900 focus:ring-blue-500'

  return (
    <div className={`h-full flex ${bgClass}`}>
      <div className={`w-80 ${cardBg} ${borderClass} border-r overflow-y-auto`}>
        <div className={`p-4 ${borderClass} border-b sticky top-0 ${cardBg} z-10`}>
          <h2 className={`font-semibold ${textClass}`}>Docking Parameters</h2>
        </div>

        <div className="p-4 space-y-5">
          <div>
            <label className={`block text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-1`}>SMILES (Ligand)</label>
            <input
              type="text"
              value={smiles}
              onChange={e => { setSmiles(e.target.value); setLigand(null); setLigandContent(null); setError(''); }}
              placeholder="CC(=O)Oc1ccccc1C(=O)O"
              disabled={isRunning}
              className={`w-full px-3 py-2 text-xs font-mono border rounded-lg focus:ring-2 focus:border-transparent ${inputClass}`}
            />
            <p className={`text-xs ${subtextClass} mt-1`}>Or upload files below</p>
          </div>

          <div>
            <label className={`block text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-2`}>Receptor (PDB)</label>
            <div className={`border-2 border-dashed ${borderClass} rounded-lg p-3 text-center ${isRunning ? 'opacity-50' : 'hover:border-blue-500'} transition-colors cursor-pointer`}>
              <input type="file" accept=".pdb" className="hidden" id="receptor" onChange={handleReceptorUpload} disabled={isRunning} />
              <label htmlFor="receptor" className="cursor-pointer text-sm">
                {receptor ? (
                  <span className="text-green-600">✓ {receptor.name}</span>
                ) : (
                  <span className={subtextClass}>Click to upload receptor PDB</span>
                )}
              </label>
            </div>
          </div>

          <div>
            <label className={`block text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-2`}>Ligand (SDF/MOL2)</label>
            <div className={`border-2 border-dashed ${borderClass} rounded-lg p-3 text-center ${isRunning ? 'opacity-50' : 'hover:border-blue-500'} transition-colors cursor-pointer`}>
              <input type="file" accept=".sdf,.mol2" className="hidden" id="ligand" onChange={handleLigandUpload} disabled={isRunning} />
              <label htmlFor="ligand" className="cursor-pointer text-sm">
                {ligand ? (
                  <span className="text-green-600">✓ {ligand.name}</span>
                ) : (
                  <span className={subtextClass}>Click to upload ligand SDF/MOL2</span>
                )}
              </label>
            </div>
          </div>

          <div>
            <label className={`block text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-2`}>Scoring Function</label>
            <select
              value={scoringFunction}
              onChange={e => setScoringFunction(e.target.value as any)}
              disabled={isRunning}
              className={`w-full px-3 py-2 border rounded-lg text-sm focus:ring-2 focus:border-transparent ${inputClass}`}
            >
              <option value="vina">AutoDock Vina (default)</option>
              <option value="gnina">GNINA (CNN scoring)</option>
              <option value="rf">RF-Score (Random Forest)</option>
              <option value="vinardo">Vinardo</option>
            </select>
          </div>

          <div>
            <h3 className={`text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-2`}>Grid Box Center</h3>
            <div className="grid grid-cols-3 gap-2">
              {(['x', 'y', 'z'] as const).map(axis => (
                <div key={axis}>
                  <label className={`text-xs ${subtextClass} uppercase`}>{axis}</label>
                  <input
                    type="number"
                    value={config[`center_${axis}`]}
                    onChange={e => handleConfigChange(`center_${axis}`, parseFloat(e.target.value))}
                    disabled={isRunning}
                    className={`w-full px-2 py-1 border rounded text-sm ${inputClass}`}
                  />
                </div>
              ))}
            </div>
          </div>

          <div>
            <h3 className={`text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-2`}>Grid Box Size (Å)</h3>
            <div className="grid grid-cols-3 gap-2">
              {(['x', 'y', 'z'] as const).map(axis => (
                <div key={axis}>
                  <label className={`text-xs ${subtextClass} uppercase`}>{axis}</label>
                  <input
                    type="number"
                    value={config[`size_${axis}`]}
                    onChange={e => handleConfigChange(`size_${axis}`, parseFloat(e.target.value))}
                    disabled={isRunning}
                    className={`w-full px-2 py-1 border rounded text-sm ${inputClass}`}
                  />
                </div>
              ))}
            </div>
          </div>

          <div>
            <label className={`block text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-1`}>
              Exhaustiveness: <span className="font-mono">{config.exhaustiveness}</span>
            </label>
            <input
              type="range" min="1" max="64" value={config.exhaustiveness}
              onChange={e => handleConfigChange('exhaustiveness', parseInt(e.target.value))}
              disabled={isRunning}
              className="w-full"
            />
          </div>

          <div>
            <label className={`block text-sm font-medium ${isDark ? 'text-gray-300' : 'text-gray-700'} mb-1`}>
              Number of Modes: <span className="font-mono">{config.num_modes}</span>
            </label>
            <input
              type="range" min="1" max="50" value={config.num_modes}
              onChange={e => handleConfigChange('num_modes', parseInt(e.target.value))}
              disabled={isRunning}
              className="w-full"
            />
          </div>

          {isRunning ? (
            <button
              onClick={cancelDocking}
              className="w-full bg-red-600 text-white py-2.5 px-4 rounded-lg font-medium hover:bg-red-700 transition-colors"
            >
              Cancel
            </button>
          ) : (
            <button
              onClick={handleDocking}
              className="w-full bg-blue-600 text-white py-2.5 px-4 rounded-lg font-medium hover:bg-blue-700 transition-colors disabled:opacity-50"
            >
              Start Docking
            </button>
          )}

          {error && (
            <div className={`p-3 rounded-lg bg-red-100 border border-red-300 text-red-700 text-sm`}>
              {error}
            </div>
          )}
        </div>
      </div>

      <div className="flex-1 flex flex-col">
        {isRunning && (
          <div className={`${cardBg} ${borderClass} border-b p-4`}>
            <div className="flex items-center justify-between mb-3">
              <h3 className={`font-semibold ${textClass}`}>Docking Progress</h3>
              <span className={`text-sm ${subtextClass}`}>{overallProgress}%</span>
            </div>
            
            <div className={`w-full h-2 ${isDark ? 'bg-gray-700' : 'bg-gray-200'} rounded-full overflow-hidden mb-4`}>
              <div 
                className="h-full bg-gradient-to-r from-blue-500 to-blue-600 transition-all duration-300 ease-out"
                style={{ width: `${overallProgress}%` }}
              />
            </div>

            <div className="grid grid-cols-1 md:grid-cols-5 gap-2">
              {stages.map(stage => (
                <div 
                  key={stage.id}
                  className={`p-2 rounded-lg ${isDark ? 'bg-gray-700' : 'bg-gray-100'}`}
                >
                  <div className="flex items-center gap-2 mb-1">
                    <span className={`text-sm ${getStageColor(stage.status, isDark)}`}>
                      {stage.status === 'running' ? (
                        <span className="animate-spin">⟳</span>
                      ) : (
                        getStageIcon(stage.status)
                      )}
                    </span>
                    <span className={`text-xs font-medium ${textClass}`}>{stage.label}</span>
                  </div>
                  
                  <div className={`w-full h-1 ${isDark ? 'bg-gray-600' : 'bg-gray-200'} rounded-full overflow-hidden`}>
                    <div 
                      className={`h-full transition-all duration-300 ${
                        stage.status === 'completed' ? 'bg-green-500' :
                        stage.status === 'running' ? 'bg-blue-500 animate-pulse' :
                        stage.status === 'failed' ? 'bg-red-500' :
                        isDark ? 'bg-gray-600' : 'bg-gray-200'
                      }`}
                      style={{ width: `${stage.progress}%` }}
                    />
                  </div>
                  
                  {stage.message && (
                    <p className={`text-xs mt-1 ${subtextClass}`}>{stage.message}</p>
                  )}
                </div>
              ))}
            </div>
          </div>
        )}

        <div className="flex-1 bg-gray-200 dark:bg-gray-900 flex items-center justify-center">
          <div className="text-center">
            <div className="text-6xl mb-4">🧬</div>
            <p className={`${subtextClass}`}>
              {isRunning ? 'Docking in progress...' : receptor || ligand || smiles ? 'Ready to dock' : 'Upload files or enter SMILES to begin'}
            </p>
          </div>
        </div>
      </div>

      <div className={`w-80 ${cardBg} ${borderClass} border-l overflow-y-auto`}>
        <div className={`p-4 ${borderClass} border-b sticky top-0 ${cardBg} z-10`}>
          <div className="flex items-center justify-between">
            <h2 className={`font-semibold ${textClass}`}>Docking Results</h2>
            {results.length > 0 && (
              <span className={`px-2 py-0.5 rounded text-xs ${isDark ? 'bg-blue-900 text-blue-300' : 'bg-blue-100 text-blue-700'}`}>
                {results.length} poses
              </span>
            )}
          </div>
        </div>

        <div className="p-4">
          {results.length === 0 ? (
            <div className={`text-center py-12 ${subtextClass}`}>
              <div className="text-4xl mb-3">📊</div>
              <p className="text-sm">No results yet</p>
              <p className="text-xs mt-1">Run a docking simulation to see results</p>
            </div>
          ) : (
            <div className="space-y-3">
              
              {downloadUrls && (
                <div className={`p-3 rounded-lg ${isDark ? 'bg-gray-700' : 'bg-gray-50'} mb-4`}>
                  <div className={`text-xs font-medium ${textClass} mb-2`}>
                    📥 Download Files
                  </div>
                  <div className="flex flex-wrap gap-2">
                    {downloadUrls.log_file && (
                      <a href={downloadUrls.log_file} download className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-blue-800 text-blue-300 hover:bg-blue-700' : 'bg-blue-100 text-blue-700 hover:bg-blue-200'} flex items-center gap-1`}>
                        📄 Log
                      </a>
                    )}
                    {downloadUrls.docking_file && (
                      <a href={downloadUrls.docking_file} download className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-green-800 text-green-300 hover:bg-green-700' : 'bg-green-100 text-green-700 hover:bg-green-200'} flex items-center gap-1`}>
                        🧬 Docking
                      </a>
                    )}
                    {downloadUrls.grid_file && (
                      <a href={downloadUrls.grid_file} download className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-purple-800 text-purple-300 hover:bg-purple-700' : 'bg-purple-100 text-purple-700 hover:bg-purple-200'} flex items-center gap-1`}>
                        📐 Grid
                      </a>
                    )}
                    {downloadUrls.vina_log && downloadUrls.vina_log !== downloadUrls.log_file && (
                      <a href={downloadUrls.vina_log} download className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-blue-800 text-blue-300 hover:bg-blue-700' : 'bg-blue-100 text-blue-700 hover:bg-blue-200'} flex items-center gap-1`}>
                        📄 Vina Log
                      </a>
                    )}
                    {downloadUrls.vina_docking && downloadUrls.vina_docking !== downloadUrls.docking_file && (
                      <a href={downloadUrls.vina_docking} download className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-green-800 text-green-300 hover:bg-green-700' : 'bg-green-100 text-green-700 hover:bg-green-200'} flex items-center gap-1`}>
                        🧬 Vina Docking
                      </a>
                    )}
                    {downloadUrls.gnina_log && (
                      <a href={downloadUrls.gnina_log} download className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-orange-800 text-orange-300 hover:bg-orange-700' : 'bg-orange-100 text-orange-700 hover:bg-orange-200'} flex items-center gap-1`}>
                        📄 GNINA Log
                      </a>
                    )}
                    {downloadUrls.gnina_docking && (
                      <a href={downloadUrls.gnina_docking} download className={`text-xs px-2 py-1 rounded ${isDark ? 'bg-red-800 text-red-300 hover:bg-red-700' : 'bg-red-100 text-red-700 hover:bg-red-200'} flex items-center gap-1`}>
                        🧬 GNINA Docking
                      </a>
                    )}
                  </div>
                </div>
              )}
              
              {routingDecision && (
                <div className={`text-xs ${subtextClass} mb-2 px-2 py-1 rounded ${isDark ? 'bg-gray-700' : 'bg-gray-100'}`}>
                  🔀 {routingDecision}
                </div>
              )}
              
              <div className={`text-xs ${subtextClass} mb-2`}>
                Sorted by {scoringFunction.toUpperCase()} score (lowest = best)
              </div>
              
              {results.map((result, i) => (
                <div 
                  key={i} 
                  className={`${isDark ? 'bg-gray-700' : 'bg-gray-50'} rounded-lg p-3 ${i === 0 ? 'ring-2 ring-green-500' : ''}`}
                >
                  <div className="flex justify-between items-center">
                    <div className="flex items-center gap-2">
                      {i === 0 && <span className="text-xs bg-green-500 text-white px-1.5 py-0.5 rounded">Best</span>}
                      <span className={`font-medium text-sm ${textClass}`}>Pose {i + 1}</span>
                    </div>
                    <span className="font-bold text-blue-600">{result.vina_score.toFixed(2)}</span>
                  </div>
                  <div className={`text-xs ${subtextClass} mt-1`}>kcal/mol</div>
                  
                  <div className="flex gap-2 mt-2">
                    {result.gnina_score !== undefined && (
                      <span className={`text-xs px-2 py-0.5 rounded ${isDark ? 'bg-purple-900 text-purple-300' : 'bg-purple-100 text-purple-700'}`}>
                        CNN: {result.gnina_score.toFixed(2)}
                      </span>
                    )}
                    {result.rf_score !== undefined && (
                      <span className={`text-xs px-2 py-0.5 rounded ${isDark ? 'bg-green-900 text-green-300' : 'bg-green-100 text-green-700'}`}>
                        RF: {result.rf_score.toFixed(2)}
                      </span>
                    )}
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  )
}
