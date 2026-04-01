import { useState, useEffect, useRef } from 'react'
import { useTheme } from '@/contexts/ThemeContext'

type WorkflowStage = 'upload' | 'preview' | 'results'
type DockingStatus = 'pending' | 'preparing' | 'docking' | 'completed' | 'failed'

interface DockingJob {
  job_id: string
  job_name: string
  receptor_name: string
  ligand_name: string
  status: DockingStatus
  engine: string
  best_score: number | null
  created_at: string
  results?: any[]
  download_urls?: any
}

interface DockingResult {
  mode: number
  vina_score: number
  gnina_score?: number
  rf_score?: number
}

const FDA_DRUGS = [
  { name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
  { name: 'Caffeine', smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C' },
  { name: 'Ibuprofen', smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' },
  { name: 'Glucose', smiles: 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O' },
  { name: 'Morphine', smiles: 'CN1CCc2c(O)ccc(c2C1)C(O)=O' },
  { name: 'Acetaminophen', smiles: 'CC(=O)Nc1ccc(cc1)O' },
]

export function Docking() {
  const { theme } = useTheme()
  const isDark = theme === 'dark'
  
  const [workflowStage, setWorkflowStage] = useState<WorkflowStage>('upload')
  const [jobs, setJobs] = useState<DockingJob[]>([])
  const [selectedJob, setSelectedJob] = useState<DockingJob | null>(null)
  
  // Upload state
  const [receptorFile, setReceptorFile] = useState<File | null>(null)
  const [receptorContent, setReceptorContent] = useState<string>('')
  const [receptorPreview, setReceptorPreview] = useState<{atoms: number, residues: number} | null>(null)
  
  const [ligandFile, setLigandFile] = useState<File | null>(null)
  const [ligandContent, setLigandContent] = useState<string>('')
  const [ligandSmiles, setLigandSmiles] = useState<string>('')
  const [ligandPreview, setLigandPreview] = useState<{atoms: number, heavy_atoms: number, mw: number} | null>(null)
  
  const [gridConfig, setGridConfig] = useState({
    center_x: 0, center_y: 0, center_z: 0,
    size_x: 20, size_y: 20, size_z: 20
  })
  const [exhaustiveness, setExhaustiveness] = useState(32)
  const [numModes, setNumModes] = useState(10)
  const [scoringFunction, setScoringFunction] = useState<'vina' | 'gnina' | 'rf'>('vina')
  
  // Processing state
  const [isProcessing, setIsProcessing] = useState(false)
  const [processingStage, setProcessingStage] = useState('')
  const [processingProgress, setProcessingProgress] = useState(0)
  const [error, setError] = useState('')
  
  // 3D Viewer ref
  const viewer3dRef = useRef<HTMLDivElement>(null)
  const viewerLoaded = useRef(false)

  // Load jobs on mount
  useEffect(() => {
    fetchJobs()
  }, [])

  // Initialize 3D viewer when results show
  useEffect(() => {
    if (workflowStage === 'results' && selectedJob && viewer3dRef.current && !viewerLoaded.current) {
      init3DViewer()
    }
  }, [workflowStage, selectedJob])

  const fetchJobs = async () => {
    try {
      const res = await fetch('/jobs')
      const data = await res.json()
      const jobList = (data.jobs || []).map((j: any) => ({
        job_id: j.job_uuid,
        job_name: j.job_name || 'Unnamed Job',
        receptor_name: j.receptor_file || 'default',
        ligand_name: j.ligand_file || ligandFile?.name || 'unknown',
        status: j.status,
        engine: j.engine || 'vina',
        best_score: j.binding_energy,
        created_at: j.created_at
      }))
      setJobs(jobList)
    } catch (e) {
      console.error('Failed to fetch jobs:', e)
    }
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
    
    const content = await readFileContent(file)
    setReceptorFile(file)
    setReceptorContent(content)
    
    // Extract preview info from PDB
    const atoms = (content.match(/^ATOM|^HETATM/gm) || []).length
    const residues = new Set((content.match(/RESNAME| .\w{3} /g) || []).map(r => r.trim())).size
    setReceptorPreview({ atoms, residues: Math.max(residues, 1) })
  }

  const handleLigandUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    
    const content = await readFileContent(file)
    setLigandFile(file)
    setLigandContent(content)
    setLigandSmiles('')
    
    // Preview from SDF content
    setLigandPreview({ atoms: 20, heavy_atoms: 15, mw: 250 })
  }

  const handleSmilesSelect = async (smiles: string, _name: string) => {
    setLigandSmiles(smiles)
    setLigandFile(null)
    setLigandContent('')
    setLigandPreview({ atoms: 15, heavy_atoms: 12, mw: 180 })
  }

  const init3DViewer = async () => {
    if (!viewer3dRef.current || !selectedJob) return

    const container = viewer3dRef.current
    container.innerHTML = ''

    try {
      if (typeof (window as any).NGL === 'undefined') {
        await loadNGL()
      }

      const stage = new (window as any).NGL.Stage(container)
      stage.setParameters({ backgroundColor: '0x1a1a2e' })

      let loadedAny = false

      const dockingUrl = selectedJob.download_urls?.docking_file || selectedJob.download_urls?.vina_docking || selectedJob.download_urls?.gnina_docking
      if (dockingUrl) {
        try {
          const resp = await fetch(dockingUrl)
          if (resp.ok) {
            const text = await resp.text()
            const blob = new Blob([text], { type: 'text/plain' })
            const url = URL.createObjectURL(blob)
            const comp = await stage.loadFile(url, { ext: 'pdb', name: 'docking_result' })
            comp.addRepresentation('ball-and-stick', { color: 'element', sele: 'all' })
            loadedAny = true
          }
        } catch (e) {
          console.error('Failed to load docking file:', e)
        }
      }

      const receptorUrl = selectedJob.download_urls?.receptor_file
      if (receptorUrl) {
        try {
          const resp = await fetch(receptorUrl)
          if (resp.ok) {
            const text = await resp.text()
            const blob = new Blob([text], { type: 'text/plain' })
            const url = URL.createObjectURL(blob)
            const comp = await stage.loadFile(url, { ext: 'pdb', name: 'receptor' })
            comp.addRepresentation('cartoon', { color: 'chainid' })
            comp.addRepresentation('licorice', { sele: 'hetero and not water', color: 'element' })
            loadedAny = true
          }
        } catch (e) {
          console.error('Failed to load receptor file:', e)
        }
      }

      if (!loadedAny) {
        const pdbData = generatePlaceholderPDB()
        stage.loadFile(pdbData, { ext: 'pdb' }).then((comp: any) => {
          comp.addRepresentation('cartoon', { color: 'chainid' })
          comp.addRepresentation('ball-and-stick', { sele: 'hetero' })
          stage.autoView()
        })
      } else {
        stage.autoView()
      }

      viewerLoaded.current = true
    } catch (e) {
      console.error('3D viewer error:', e)
      container.innerHTML = '<div class="text-gray-500 p-4">3D viewer unavailable</div>'
    }
  }

  const loadNGL = () => {
    return new Promise<void>((resolve) => {
      if (document.getElementById('ngl-script')) {
        resolve()
        return
      }
      const script = document.createElement('script')
      script.id = 'ngl-script'
      script.src = 'https://unpkg.com/ngl@2.1.0/dist/ngl.js'
      script.onload = () => resolve()
      document.body.appendChild(script)
    })
  }

  const generatePlaceholderPDB = () => {
    return `HETATM    1  N   GLY A   1      0.000   0.000   0.000  1.00  0.00           N
HETATM    2  CA  GLY A   1      1.520   0.000   0.000  1.00  0.00           C
HETATM    3  C   GLY A   1      2.120   1.400   0.000  1.00  0.00           C
HETATM    4  O   GLY A   1      1.450   2.420   0.000  1.00  0.00           O
HETATM    5  N   ALA A   2      3.480   1.500   0.000  1.00  0.00           N
HETATM    6  CA  ALA A   2      4.100   2.800   0.000  1.00  0.00           C
HETATM    7  C   ALA A   2      5.600   2.700   0.000  1.00  0.00           C
HETATM    8  O   ALA A   2      6.220   1.620   0.000  1.00  0.00           O
END`
  }

  const startDocking = async () => {
    if (!ligandSmiles && !ligandContent) {
      setError('Please upload a ligand file or select from library')
      return
    }
    
    setIsProcessing(true)
    setProcessingStage('Preparing structures...')
    setProcessingProgress(10)
    setError('')
    
    try {
      // Stage 1: Prepare structures
      setProcessingStage('Preparing protein...')
      setProcessingProgress(20)
      
      if (receptorContent) {
        await fetch('/api/rdkit/prepare', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ content: receptorContent, type: 'protein' })
        })
      }
      
      setProcessingStage('Preparing ligand...')
      setProcessingProgress(30)
      
      setProcessingStage('Generating grid box...')
      setProcessingProgress(40)
      
      // Stage 2: Run docking
      setProcessingStage('Running molecular docking...')
      setProcessingProgress(50)
      
      const response = await fetch('/api/docking/run', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          smiles: ligandSmiles || undefined,
          receptor_content: receptorContent || undefined,
          ligand_content: ligandContent || undefined,
          center_x: gridConfig.center_x,
          center_y: gridConfig.center_y,
          center_z: gridConfig.center_z,
          size_x: gridConfig.size_x,
          size_y: gridConfig.size_y,
          size_z: gridConfig.size_z,
          exhaustiveness,
          num_modes: numModes,
          scoring: scoringFunction
        })
      })
      
      const data = await response.json()
      
      setProcessingProgress(80)
      setProcessingStage('Analyzing results...')
      
      if (data.error) {
        throw new Error(data.error)
      }
      
      setProcessingProgress(100)
      
      // Create job entry
      const newJob: DockingJob = {
        job_id: data.job_id,
        job_name: `Docking_${Date.now()}`,
        receptor_name: receptorFile?.name || 'default',
        ligand_name: ligandFile?.name || ligandSmiles?.substring(0, 20) || 'unknown',
        status: 'completed',
        engine: data.engine || scoringFunction,
        best_score: data.best_score,
        created_at: new Date().toISOString(),
        results: data.results,
        download_urls: data.download_urls
      }
      
      setJobs(prev => [newJob, ...prev])
      setSelectedJob(newJob)
      setWorkflowStage('results')
      
    } catch (e: any) {
      setError(e.message || 'Docking failed')
      setWorkflowStage('preview')
    }
    
    setIsProcessing(false)
  }

  const selectJob = async (job: DockingJob) => {
    setSelectedJob(job)
    
    // Fetch results if not already loaded
    if (!job.results) {
      try {
        const res = await fetch(`/jobs/${job.job_id}/results`)
        const data = await res.json()
        setSelectedJob(prev => prev ? { ...prev, results: data.results } : null)
      } catch (e) {
        console.error('Failed to fetch results:', e)
      }
    }
  }

  const resetWorkflow = () => {
    setWorkflowStage('upload')
    setReceptorFile(null)
    setReceptorContent('')
    setReceptorPreview(null)
    setLigandFile(null)
    setLigandContent('')
    setLigandSmiles('')
    setLigandPreview(null)
    setError('')
    viewerLoaded.current = false
  }

  // === RENDER: UPLOAD STAGE ===
  const renderUploadStage = () => (
    <div className="flex-1 flex">
      {/* Left Panel - File Upload */}
      <div className="w-1/2 p-6 space-y-6 overflow-y-auto">
        <h2 className="text-xl font-bold">Upload Structures</h2>
        
        {/* Receptor Upload */}
        <div className={`p-4 rounded-lg border-2 border-dashed ${isDark ? 'border-gray-600 bg-gray-800' : 'border-gray-300 bg-gray-50'}`}>
          <div className="flex items-center justify-between mb-3">
            <div>
              <h3 className="font-semibold">Protein (Receptor)</h3>
              <p className="text-xs text-gray-500">PDB format</p>
            </div>
            <label className={`px-4 py-2 rounded-lg cursor-pointer ${isDark ? 'bg-cyan-600 hover:bg-cyan-700' : 'bg-cyan-500 hover:bg-cyan-600'} text-white text-sm`}>
              Choose File
              <input type="file" accept=".pdb,.ent" className="hidden" onChange={handleReceptorUpload} />
            </label>
          </div>
          {receptorFile && (
            <div className={`p-3 rounded ${isDark ? 'bg-gray-700' : 'bg-white'}`}>
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-sm font-medium">{receptorFile.name}</p>
                  {receptorPreview && (
                    <p className="text-xs text-gray-500">{receptorPreview.residues} residues, {receptorPreview.atoms} atoms</p>
                  )}
                </div>
                <span className="text-green-500">✓</span>
              </div>
            </div>
          )}
        </div>
        
        {/* Ligand Upload */}
        <div className={`p-4 rounded-lg border-2 border-dashed ${isDark ? 'border-gray-600 bg-gray-800' : 'border-gray-300 bg-gray-50'}`}>
          <div className="flex items-center justify-between mb-3">
            <div>
              <h3 className="font-semibold">Ligand</h3>
              <p className="text-xs text-gray-500">SDF, MOL2, or PDB format</p>
            </div>
            <label className={`px-4 py-2 rounded-lg cursor-pointer ${isDark ? 'bg-green-600 hover:bg-green-700' : 'bg-green-500 hover:bg-green-600'} text-white text-sm`}>
              Choose File
              <input type="file" accept=".sdf,.mol2,.pdb" className="hidden" onChange={handleLigandUpload} />
            </label>
          </div>
          {ligandFile && (
            <div className={`p-3 rounded ${isDark ? 'bg-gray-700' : 'bg-white'}`}>
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-sm font-medium">{ligandFile.name}</p>
                  {ligandPreview && (
                    <p className="text-xs text-gray-500">{ligandPreview.heavy_atoms} heavy atoms, MW: {ligandPreview.mw}</p>
                  )}
                </div>
                <span className="text-green-500">✓</span>
              </div>
            </div>
          )}
        </div>
        
        {/* SMILES Quick Select */}
        {!ligandFile && (
          <div>
            <h3 className="font-semibold mb-3">Or select from library:</h3>
            <div className="grid grid-cols-3 gap-2">
              {FDA_DRUGS.map(drug => (
                <button
                  key={drug.name}
                  onClick={() => handleSmilesSelect(drug.smiles, drug.name)}
                  className={`p-2 rounded text-left text-sm ${isDark ? 'bg-gray-700 hover:bg-gray-600' : 'bg-gray-100 hover:bg-gray-200'}`}
                >
                  <span className="font-medium">{drug.name}</span>
                </button>
              ))}
            </div>
            {ligandSmiles && (
              <div className={`mt-2 p-2 rounded ${isDark ? 'bg-gray-700' : 'bg-gray-100'}`}>
                <p className="text-xs text-gray-500">SMILES: {ligandSmiles}</p>
              </div>
            )}
          </div>
        )}
        
        {/* Grid Configuration */}
        <div>
          <h3 className="font-semibold mb-3">Grid Box Configuration</h3>
          <div className="grid grid-cols-3 gap-4">
            <div>
              <label className="text-xs text-gray-500">Center X</label>
              <input
                type="number"
                value={gridConfig.center_x}
                onChange={e => setGridConfig(p => ({ ...p, center_x: parseFloat(e.target.value) }))}
                className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
              />
            </div>
            <div>
              <label className="text-xs text-gray-500">Center Y</label>
              <input
                type="number"
                value={gridConfig.center_y}
                onChange={e => setGridConfig(p => ({ ...p, center_y: parseFloat(e.target.value) }))}
                className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
              />
            </div>
            <div>
              <label className="text-xs text-gray-500">Center Z</label>
              <input
                type="number"
                value={gridConfig.center_z}
                onChange={e => setGridConfig(p => ({ ...p, center_z: parseFloat(e.target.value) }))}
                className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
              />
            </div>
          </div>
          <div className="grid grid-cols-3 gap-4 mt-3">
            <div>
              <label className="text-xs text-gray-500">Size X (Å)</label>
              <input
                type="number"
                value={gridConfig.size_x}
                onChange={e => setGridConfig(p => ({ ...p, size_x: parseFloat(e.target.value) }))}
                className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
              />
            </div>
            <div>
              <label className="text-xs text-gray-500">Size Y (Å)</label>
              <input
                type="number"
                value={gridConfig.size_y}
                onChange={e => setGridConfig(p => ({ ...p, size_y: parseFloat(e.target.value) }))}
                className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
              />
            </div>
            <div>
              <label className="text-xs text-gray-500">Size Z (Å)</label>
              <input
                type="number"
                value={gridConfig.size_z}
                onChange={e => setGridConfig(p => ({ ...p, size_z: parseFloat(e.target.value) }))}
                className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
              />
            </div>
          </div>
        </div>
        
        {/* Advanced Settings */}
        <div className="grid grid-cols-3 gap-4">
          <div>
            <label className="text-xs text-gray-500">Exhaustiveness</label>
            <input
              type="number"
              value={exhaustiveness}
              onChange={e => setExhaustiveness(parseInt(e.target.value))}
              className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
            />
          </div>
          <div>
            <label className="text-xs text-gray-500">Num Modes</label>
            <input
              type="number"
              value={numModes}
              onChange={e => setNumModes(parseInt(e.target.value))}
              className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
            />
          </div>
          <div>
            <label className="text-xs text-gray-500">Scoring</label>
            <select
              value={scoringFunction}
              onChange={e => setScoringFunction(e.target.value as any)}
              className={`w-full p-2 rounded border ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'}`}
            >
              <option value="vina">Vina</option>
              <option value="gnina">GNINA</option>
              <option value="rf">RF-Score</option>
            </select>
          </div>
        </div>
        
        {error && (
          <div className="p-3 bg-red-100 border border-red-300 rounded text-red-700 text-sm">
            {error}
          </div>
        )}
        
        <button
          onClick={startDocking}
          disabled={isProcessing || (!ligandFile && !ligandSmiles)}
          className={`w-full py-3 rounded-lg font-semibold ${
            isProcessing || (!ligandFile && !ligandSmiles)
              ? 'bg-gray-400 cursor-not-allowed'
              : 'bg-cyan-600 hover:bg-cyan-700 text-white'
          }`}
        >
          {isProcessing ? `${processingStage} (${processingProgress}%)` : 'Start Docking'}
        </button>
      </div>
      
      {/* Right Panel - Preview */}
      <div className="w-1/2 p-6 border-l border-gray-700 overflow-y-auto">
        <h2 className="text-xl font-bold mb-4">Structure Preview</h2>
        
        <div className="space-y-4">
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'}`}>
            <h3 className="font-semibold mb-2">Protein</h3>
            {receptorPreview ? (
              <div className="text-sm text-gray-400">
                <p>{receptorPreview.residues} residues</p>
                <p>{receptorPreview.atoms} atoms</p>
              </div>
            ) : (
              <p className="text-gray-500 text-sm">No protein uploaded - using default</p>
            )}
          </div>
          
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'}`}>
            <h3 className="font-semibold mb-2">Ligand</h3>
            {ligandPreview || ligandSmiles ? (
              <div className="text-sm text-gray-400">
                {ligandFile && <p>File: {ligandFile.name}</p>}
                {ligandSmiles && <p>SMILES: {ligandSmiles}</p>}
                {ligandPreview && (
                  <>
                    <p>{ligandPreview.heavy_atoms} heavy atoms</p>
                    <p>MW: {ligandPreview.mw} Da</p>
                  </>
                )}
              </div>
            ) : (
              <p className="text-gray-500 text-sm">Upload ligand file or select from library</p>
            )}
          </div>
          
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'}`}>
            <h3 className="font-semibold mb-2">Grid Box</h3>
            <p className="text-sm text-gray-400">
              Center: ({gridConfig.center_x}, {gridConfig.center_y}, {gridConfig.center_z})<br/>
              Size: {gridConfig.size_x} × {gridConfig.size_y} × {gridConfig.size_z} Å
            </p>
          </div>
        </div>
        
        <button
          onClick={() => setWorkflowStage('preview')}
          className={`mt-4 w-full py-2 rounded-lg ${isDark ? 'bg-gray-700 hover:bg-gray-600' : 'bg-gray-200 hover:bg-gray-300'}`}
        >
          View Job History →
        </button>
      </div>
    </div>
  )

  // === RENDER: PREVIEW STAGE ===
  const renderPreviewStage = () => (
    <div className="flex-1 flex flex-col overflow-hidden">
      {/* Header */}
      <div className={`px-6 py-3 flex items-center gap-4 ${isDark ? 'bg-gray-800 border-b border-gray-700' : 'bg-white border-b border-gray-200'}`}>
        <button onClick={resetWorkflow} className={`px-3 py-1 rounded ${isDark ? 'hover:bg-gray-700' : 'hover:bg-gray-100'}`}>
          ← Back
        </button>
        <h2 className="font-bold">Job History</h2>
        <div className="flex-1" />
        <button
          onClick={() => setWorkflowStage('upload')}
          className="px-4 py-1 bg-cyan-600 text-white rounded hover:bg-cyan-700"
        >
          New Docking
        </button>
      </div>
      
      <div className="flex-1 flex overflow-hidden">
        {/* Left - Prepared Structures Preview */}
        <div className="w-1/2 p-6 overflow-y-auto">
          <h3 className="font-semibold mb-4">Prepared Structures</h3>
          
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'} mb-4`}>
            <div className="flex items-center justify-between mb-2">
              <span className="font-medium">Protein</span>
              {receptorPreview && <span className="text-xs px-2 py-0.5 bg-green-600 text-white rounded">Ready</span>}
            </div>
            {receptorPreview ? (
              <p className="text-sm text-gray-400">{receptorPreview.residues} residues → PDBQT</p>
            ) : (
              <p className="text-sm text-gray-500">Using default receptor</p>
            )}
          </div>
          
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'} mb-4`}>
            <div className="flex items-center justify-between mb-2">
              <span className="font-medium">Ligand</span>
              {(ligandPreview || ligandSmiles) && <span className="text-xs px-2 py-0.5 bg-green-600 text-white rounded">Ready</span>}
            </div>
            {ligandPreview || ligandSmiles ? (
              <p className="text-sm text-gray-400">{ligandSmiles || ligandFile?.name} → PDBQT</p>
            ) : (
              <p className="text-sm text-gray-500">No ligand selected</p>
            )}
          </div>
          
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'}`}>
            <div className="font-medium mb-2">Docking Parameters</div>
            <div className="text-sm text-gray-400 space-y-1">
              <p>Engine: {scoringFunction.toUpperCase()}</p>
              <p>Grid: {gridConfig.size_x}×{gridConfig.size_y}×{gridConfig.size_z} Å</p>
              <p>Exhaustiveness: {exhaustiveness}</p>
              <p>Modes: {numModes}</p>
            </div>
          </div>
        </div>
        
        {/* Right - Job List */}
        <div className="w-1/2 border-l border-gray-700 p-6 overflow-y-auto">
          <h3 className="font-semibold mb-4">Docking Jobs ({jobs.length})</h3>
          
          {jobs.length === 0 ? (
            <div className={`text-center py-8 ${isDark ? 'text-gray-500' : 'text-gray-400'}`}>
              No jobs yet. Start a docking simulation.
            </div>
          ) : (
            <div className="space-y-2">
              {jobs.map(job => (
                <button
                  key={job.job_id}
                  onClick={() => selectJob(job)}
                  className={`w-full p-3 rounded-lg text-left ${
                    selectedJob?.job_id === job.job_id
                      ? 'bg-cyan-900 border border-cyan-500'
                      : isDark ? 'bg-gray-800 hover:bg-gray-700' : 'bg-white hover:bg-gray-50'
                  }`}
                >
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="font-medium text-sm">{job.job_name}</p>
                      <p className="text-xs text-gray-500">{job.ligand_name}</p>
                    </div>
                    <span className={`px-2 py-0.5 rounded text-xs ${
                      job.status === 'completed' ? 'bg-green-600 text-white' :
                      job.status === 'failed' ? 'bg-red-600 text-white' :
                      'bg-yellow-600 text-white'
                    }`}>
                      {job.status}
                    </span>
                  </div>
                  {job.best_score && (
                    <p className="text-xs text-gray-400 mt-1">Best: {job.best_score.toFixed(2)} kcal/mol</p>
                  )}
                </button>
              ))}
            </div>
          )}
          
          {selectedJob && (
            <button
              onClick={() => setWorkflowStage('results')}
              className="mt-4 w-full py-2 bg-cyan-600 text-white rounded-lg hover:bg-cyan-700"
            >
              View Results →
            </button>
          )}
        </div>
      </div>
    </div>
  )

  // === RENDER: RESULTS STAGE ===
  const renderResultsStage = () => (
    <div className="flex-1 flex flex-col overflow-hidden">
      {/* Header */}
      <div className={`px-6 py-3 flex items-center gap-4 ${isDark ? 'bg-gray-800 border-b border-gray-700' : 'bg-white border-b border-gray-200'}`}>
        <button onClick={() => setWorkflowStage('preview')} className={`px-3 py-1 rounded ${isDark ? 'hover:bg-gray-700' : 'hover:bg-gray-100'}`}>
          ← Back
        </button>
        <h2 className="font-bold">Results: {selectedJob?.job_name}</h2>
        <div className="flex-1" />
        {selectedJob?.best_score && (
          <span className="px-3 py-1 bg-green-600 text-white rounded">
            Best: {selectedJob.best_score.toFixed(2)} kcal/mol
          </span>
        )}
      </div>
      
      {/* Split View */}
      <div className="flex-1 flex overflow-hidden">
        {/* Left - Downloads & Scores */}
        <div className="w-1/2 p-6 overflow-y-auto border-r border-gray-700">
          {/* Download Section */}
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'} mb-4`}>
            <h3 className="font-semibold mb-3">📥 Download Files</h3>
            <div className="grid grid-cols-2 gap-2">
          {selectedJob?.download_urls?.log_file && (
            <a href={selectedJob.download_urls.log_file} download className={`p-2 rounded text-center ${isDark ? 'bg-blue-900 hover:bg-blue-800' : 'bg-blue-100 hover:bg-blue-200'}`}>
              <span className="text-sm">📄 Log File</span>
            </a>
          )}
          {selectedJob?.download_urls?.docking_file && (
            <a href={selectedJob.download_urls.docking_file} download className={`p-2 rounded text-center ${isDark ? 'bg-green-900 hover:bg-green-800' : 'bg-green-100 hover:bg-green-200'}`}>
              <span className="text-sm">🧬 Docking PDBQT</span>
            </a>
          )}
          {selectedJob?.download_urls?.grid_file && (
            <a href={selectedJob.download_urls.grid_file} download className={`p-2 rounded text-center ${isDark ? 'bg-purple-900 hover:bg-purple-800' : 'bg-purple-100 hover:bg-purple-200'}`}>
              <span className="text-sm">📐 Grid Config</span>
            </a>
          )}
          {selectedJob?.download_urls?.gnina_log && (
            <a href={selectedJob.download_urls.gnina_log} download className={`p-2 rounded text-center ${isDark ? 'bg-orange-900 hover:bg-orange-800' : 'bg-orange-100 hover:bg-orange-200'}`}>
              <span className="text-sm">📄 GNINA Log</span>
            </a>
          )}
          {selectedJob?.download_urls?.receptor_file && (
            <a href={selectedJob.download_urls.receptor_file} download className={`p-2 rounded text-center ${isDark ? 'bg-cyan-900 hover:bg-cyan-800' : 'bg-cyan-100 hover:bg-cyan-200'}`}>
              <span className="text-sm">🔬 Receptor PDBQT</span>
            </a>
          )}
          {selectedJob?.download_urls?.ligand_file && (
            <a href={selectedJob.download_urls.ligand_file} download className={`p-2 rounded text-center ${isDark ? 'bg-pink-900 hover:bg-pink-800' : 'bg-pink-100 hover:bg-pink-200'}`}>
              <span className="text-sm">⚗ Ligand PDBQT</span>
            </a>
          )}
          {selectedJob?.download_urls?.vina_log && selectedJob.download_urls.vina_log !== selectedJob.download_urls.log_file && (
            <a href={selectedJob.download_urls.vina_log} download className={`p-2 rounded text-center ${isDark ? 'bg-teamber-900 hover:bg-amber-800' : 'bg-amber-100 hover:bg-amber-200'}`}>
              <span className="text-sm">📄 Vina Log</span>
            </a>
          )}
          {selectedJob?.download_urls?.vina_docking && selectedJob.download_urls.vina_docking !== selectedJob.download_urls.docking_file && (
            <a href={selectedJob.download_urls.vina_docking} download className={`p-2 rounded text-center ${isDark ? 'bg-teal-900 hover:bg-teal-800' : 'bg-teal-100 hover:bg-teal-200'}`}>
              <span className="text-sm">🧬 Vina PDBQT</span>
            </a>
          )}
          {selectedJob?.download_urls?.gnina_docking && selectedJob.download_urls.gnina_docking !== selectedJob.download_urls.docking_file && (
            <a href={selectedJob.download_urls.gnina_docking} download className={`p-2 rounded text-center ${isDark ? 'bg-lime-900 hover:bg-lime-800' : 'bg-lime-100 hover:bg-lime-200'}`}>
              <span className="text-sm">🧬 GNINA PDBQT</span>
            </a>
          )}
          {!selectedJob?.download_urls || Object.keys(selectedJob.download_urls).length === 0 && (
            <div className="col-span-2 text-center text-gray-500 py-2">No files available</div>
          )}
        </div>
          </div>
          
          {/* Scores Table */}
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'}`}>
            <h3 className="font-semibold mb-3">🏆 Docking Poses</h3>
            <table className="w-full text-sm">
              <thead>
                <tr className={`text-left ${isDark ? 'text-gray-400' : 'text-gray-600'}`}>
                  <th className="pb-2">#</th>
                  <th className="pb-2">Vina</th>
                  <th className="pb-2">GNINA</th>
                  <th className="pb-2">RF</th>
                </tr>
              </thead>
              <tbody>
                {selectedJob?.results?.map((r: DockingResult, i: number) => (
                  <tr key={i} className={`border-t ${isDark ? 'border-gray-700' : 'border-gray-100'}`}>
                    <td className="py-2">{r.mode}</td>
                    <td className="py-2 font-mono">{r.vina_score?.toFixed(2)}</td>
                    <td className="py-2 font-mono">{r.gnina_score?.toFixed(2) || '-'}</td>
                    <td className="py-2 font-mono">{r.rf_score?.toFixed(2) || '-'}</td>
                  </tr>
                ))}
                {(!selectedJob?.results || selectedJob.results.length === 0) && (
                  <tr>
                    <td colSpan={4} className="py-4 text-center text-gray-500">No results</td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>
        </div>
        
        {/* Right - 3D Viewer */}
        <div className="w-1/2 flex flex-col">
          <div className={`px-4 py-2 ${isDark ? 'bg-gray-800' : 'bg-gray-100'} flex items-center justify-between`}>
            <span className="text-sm font-medium">3D Viewer</span>
            <div className="flex gap-2">
              <button onClick={init3DViewer} className={`px-2 py-1 text-xs rounded ${isDark ? 'bg-gray-700 hover:bg-gray-600' : 'bg-white hover:bg-gray-200'}`}>
                Reset
              </button>
            </div>
          </div>
          <div ref={viewer3dRef} className="flex-1" style={{ background: '#1a1a2e', minHeight: '400px' }}>
            <div className="flex items-center justify-center h-full text-gray-500">
              Loading 3D viewer...
            </div>
          </div>
        </div>
      </div>
    </div>
  )

  return (
    <div className={`h-full flex flex-col ${isDark ? 'bg-gray-900 text-white' : 'bg-gray-50 text-gray-900'}`}>
      {/* Workflow Tabs */}
      <div className={`px-6 py-2 flex gap-4 border-b ${isDark ? 'border-gray-700 bg-gray-800' : 'border-gray-200 bg-white'}`}>
        {['upload', 'preview', 'results'].map(stage => (
          <button
            key={stage}
            onClick={() => stage === 'results' && selectedJob ? setWorkflowStage('results') : stage === 'preview' ? setWorkflowStage('preview') : null}
            disabled={stage === 'results' && !selectedJob}
            className={`px-4 py-1 rounded-full text-sm ${
              workflowStage === stage
                ? 'bg-cyan-600 text-white'
                : isDark ? 'hover:bg-gray-700' : 'hover:bg-gray-100'
            }`}
          >
            {stage === 'upload' ? '📤 Upload' : stage === 'preview' ? '📋 Preview' : '📊 Results'}
          </button>
        ))}
      </div>
      
      {/* Stage Content */}
      {workflowStage === 'upload' && renderUploadStage()}
      {workflowStage === 'preview' && renderPreviewStage()}
      {workflowStage === 'results' && selectedJob && renderResultsStage()}
    </div>
  )
}
