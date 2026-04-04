import { useState, useEffect, useRef } from 'react'
import { useTheme } from '@/contexts/ThemeContext'
import { InteractionPanel } from '@/components/InteractionPanel'
import { ExportPanel } from '@/components/ExportPanel'

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
  completed_at?: string
  results?: any[]
  download_urls?: any
  receptor_content?: string
  ligand_pdb?: string
  log_text?: string
}

interface DockingResult {
  mode: number
  vina_score: number
  gnina_score?: number
  rf_score?: number
  hydrophobic_term?: number
  rotatable_penalty?: number
  lipo_contact?: number
  final_score?: number
  composite_score?: number
  constraint_penalty?: number
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
  const [receptorFormat, setReceptorFormat] = useState<string>('')
  const [receptorFilename, setReceptorFilename] = useState<string>('')
  const [receptorPreview, setReceptorPreview] = useState<{atoms: number, residues: number} | null>(null)
  
  const [ligandFile, setLigandFile] = useState<File | null>(null)
  const [ligandContent, setLigandContent] = useState<string>('')
  const [ligandFormat, setLigandFormat] = useState<string>('')
  const [ligandFilename, setLigandFilename] = useState<string>('')
  const [ligandSmiles, setLigandSmiles] = useState<string>('')
  const [ligandPreview, setLigandPreview] = useState<{atoms: number, heavy_atoms: number, mw: number} | null>(null)
  
  const [gridConfig, setGridConfig] = useState({
    center_x: 0, center_y: 0, center_z: 0,
    size_x: 20, size_y: 20, size_z: 20
  })
  const [exhaustiveness, setExhaustiveness] = useState(32)
  const [numModes, setNumModes] = useState(10)
  const [scoringFunction, setScoringFunction] = useState<'vina' | 'gnina' | 'rf'>('vina')
  const [enableFlexibility, setEnableFlexibility] = useState(false)
  const [constraints, setConstraints] = useState<any[]>([])
  const [showCartoon, setShowCartoon] = useState(true)
  const [showSurface, setShowSurface] = useState(false)
  const [selectedPoseIdx, setSelectedPoseIdx] = useState(0)
  const nglStageRef = useRef<any>(null)

  // AI Job Assistant state
  const [aiQuestion, setAiQuestion] = useState('')
  const [aiAnswer, setAiAnswer] = useState('')
  const [aiLoading, setAiLoading] = useState(false)
  const [aiJobId, setAiJobId] = useState('')
  const [showAIPanel, setShowAIPanel] = useState(false)
  
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

  // Re-render viewer when toggles change
  useEffect(() => {
    if (nglStageRef.current) {
      nglStageRef.current.eachComponent((comp: any) => {
        comp.eachRepresentation((rep: any) => {
          if (rep.name === 'cartoon') rep.setVisibility(showCartoon)
          if (rep.name === 'surface') rep.setVisibility(showSurface)
        })
      })
    }
  }, [showCartoon, showSurface])

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

  const detectFormat = (filename: string, content: string, isLigand: boolean): string => {
    const ext = filename.split('.').pop()?.toLowerCase() || ''
    const extMap: Record<string, string> = {
      pdb: 'pdb', pdbqt: 'pdbqt', mol: 'mol', mol2: 'mol2',
      sdf: 'sdf', smi: 'smiles', smiles: 'smiles', inchi: 'inchi'
    }
    if (extMap[ext]) return extMap[ext]
    const trimmed = content.trim()
    if (isLigand && trimmed.length < 500 && !trimmed.includes('\n')) {
      if (/^[A-Za-z0-9@+\-\[\]\(\)\\%=#]+$/.test(trimmed)) return 'smiles'
    }
    if (trimmed.startsWith('InChI=')) return 'inchi'
    if (/^ATOM |^HETATM /m.test(trimmed)) {
      const lines = trimmed.split('\n').filter(l => l.startsWith('ATOM') || l.startsWith('HETATM'))
      if (lines.length > 0 && lines[0].length >= 78) {
        const type = lines[0].substring(76, 78).trim()
        if (['C','A','N','NA','OA','S','SA','P','F','CL','BR','I','HD'].includes(type)) return 'pdbqt'
      }
      return 'pdb'
    }
    if (trimmed.includes('V2000') || trimmed.includes('V3000')) return trimmed.includes('$$$$') ? 'sdf' : 'mol'
    if (trimmed.includes('@<TRIPOS>')) return 'mol2'
    return 'unknown'
  }

  const handleReceptorUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    
    const content = await readFileContent(file)
    const format = detectFormat(file.name, content, false)
    setReceptorFile(file)
    setReceptorContent(content)
    setReceptorFormat(format)
    setReceptorFilename(file.name)
    
    const atoms = (content.match(/^ATOM|^HETATM/gm) || []).length
    const residues = new Set((content.match(/ .\w{3} /g) || []).map(r => r.trim())).size
    setReceptorPreview({ atoms, residues: Math.max(residues, 1) })
  }

  const fetchLigandProps = async (smiles: string) => {
    try {
      const res = await fetch('/api/chem/properties', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const d = await res.json()
      if (d && !d.error) {
        setLigandPreview({ atoms: d.num_atoms || 0, heavy_atoms: d.heavy_atoms || 0, mw: Math.round(d.mw || 0) })
      }
    } catch { /* keep existing preview */ }
  }

  const handleLigandUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    
    const content = await readFileContent(file)
    const format = detectFormat(file.name, content, true)
    setLigandFile(file)
    setLigandContent(content)
    setLigandFormat(format)
    setLigandFilename(file.name)
    setLigandSmiles('')
    setLigandPreview({ atoms: 0, heavy_atoms: 0, mw: 0 })

    const sdfSmiles = content.match(/^> +<SMILES>\s*\n([^\n]+)/m)?.[1]
    if (sdfSmiles) fetchLigandProps(sdfSmiles.trim())
  }

  const handleSmilesSelect = async (smiles: string, _name: string) => {
    setLigandSmiles(smiles)
    setLigandFile(null)
    setLigandContent('')
    setLigandFormat('smiles')
    setLigandFilename('')
    setLigandPreview({ atoms: 0, heavy_atoms: 0, mw: 0 })
    fetchLigandProps(smiles)
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
      nglStageRef.current = stage

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
    setProcessingStage('Submitting job...')
    setProcessingProgress(5)
    setError('')

    try {
      const response = await fetch('/api/docking/run', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          smiles: ligandSmiles || undefined,
          receptor_content: receptorContent || undefined,
          receptor_format: receptorFormat || undefined,
          receptor_filename: receptorFilename || undefined,
          ligand_content: ligandContent || undefined,
          ligand_format: ligandFormat || (ligandSmiles ? 'smiles' : undefined),
          ligand_filename: ligandFilename || undefined,
          center_x: gridConfig.center_x,
          center_y: gridConfig.center_y,
          center_z: gridConfig.center_z,
          size_x: gridConfig.size_x,
          size_y: gridConfig.size_y,
          size_z: gridConfig.size_z,
          exhaustiveness,
          num_modes: numModes,
          scoring: scoringFunction,
          enable_flexibility: enableFlexibility,
          constraints: constraints.length > 0 ? constraints : undefined
        })
      })

      const accepted = await response.json()
      if (accepted.error) throw new Error(accepted.error)

      const jobId = accepted.job_id
      setProcessingStage('Docking in progress...')

      // Poll /dock/{id}/status until done
      let data: any = null
      for (let attempt = 0; attempt < 300; attempt++) {
        await new Promise(r => setTimeout(r, 2000))
        const statusRes = await fetch(`/dock/${jobId}/status`)
        const status = await statusRes.json()
        const pct = Math.min(10 + Math.round((status.progress || 0) * 0.85), 90)
        setProcessingProgress(pct)
        setProcessingStage(status.message || 'Running...')

        if (status.status === 'completed') {
          const resultRes = await fetch(`/api/docking/result/${jobId}`)
          data = await resultRes.json()
          break
        }
        if (status.status === 'failed') {
          throw new Error(status.message || 'Docking failed')
        }
      }

      if (!data) throw new Error('Docking timed out')

      setProcessingProgress(100)

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
        download_urls: data.download_urls,
        receptor_content: receptorContent || undefined,
        ligand_pdb: data.files?.docking || undefined,
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
    setAiJobId(job.job_id)
    setAiAnswer('')

    if (job.status === 'completed' && !job.results?.length) {
      try {
        const res = await fetch(`/api/jobs/${job.job_id}/full`)
        if (res.ok) {
          const data = await res.json()
          setSelectedJob({
            job_id: data.job_id,
            job_name: data.job_name || job.job_name,
            receptor_name: data.receptor_name || job.receptor_name,
            ligand_name: data.ligand_name || job.ligand_name,
            status: 'completed',
            engine: data.engine || job.engine,
            best_score: data.best_score ?? job.best_score,
            created_at: data.created_at || job.created_at,
            completed_at: data.completed_at,
            results: data.results,
            download_urls: data.download_urls,
            log_text: data.log_text,
          })
        }
      } catch (e) {
        console.error('Failed to restore job from history:', e)
      }
    }
  }

  const askAI = async () => {
    if (!aiJobId || !aiQuestion.trim()) return
    setAiLoading(true)
    setAiAnswer('')
    try {
      const res = await fetch('/api/ai/job-explain', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ job_id: aiJobId, question: aiQuestion })
      })
      const data = await res.json()
      setAiAnswer(data.answer || data.error || 'No response')
    } catch (e: any) {
      setAiAnswer('Error: ' + e.message)
    }
    setAiLoading(false)
  }

  const resetWorkflow = () => {
    setWorkflowStage('upload')
    setReceptorFile(null)
    setReceptorContent('')
    setReceptorFormat('')
    setReceptorFilename('')
    setReceptorPreview(null)
    setLigandFile(null)
    setLigandContent('')
    setLigandFormat('')
    setLigandFilename('')
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
              <p className="text-xs text-gray-500">PDB, PDBQT, MOL2, SDF, PDB ID</p>
            </div>
            <label className={`px-4 py-2 rounded-lg cursor-pointer ${isDark ? 'bg-cyan-600 hover:bg-cyan-700' : 'bg-cyan-500 hover:bg-cyan-600'} text-white text-sm`}>
              Choose File
              <input type="file" accept=".pdb,.pdbqt,.mol2,.sdf,.mol" className="hidden" onChange={handleReceptorUpload} />
            </label>
          </div>
          {receptorFile && (
            <div className={`p-3 rounded ${isDark ? 'bg-gray-700' : 'bg-white'}`}>
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-sm font-medium">{receptorFile.name}</p>
                  <p className="text-xs text-gray-500">
                    {receptorFormat && <span className="inline-block px-1.5 py-0.5 bg-gray-100 rounded text-xs mr-1">{receptorFormat.toUpperCase()}</span>}
                    {receptorPreview ? `${receptorPreview.residues} residues, ${receptorPreview.atoms} atoms` : ''}
                  </p>
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
              <p className="text-xs text-gray-500">SMILES, SDF, MOL, MOL2, PDB, PDBQT</p>
            </div>
            <label className={`px-4 py-2 rounded-lg cursor-pointer ${isDark ? 'bg-green-600 hover:bg-green-700' : 'bg-green-500 hover:bg-green-600'} text-white text-sm`}>
              Choose File
              <input type="file" accept=".sdf,.mol2,.pdb,.pdbqt,.mol,.smi" className="hidden" onChange={handleLigandUpload} />
            </label>
          </div>
          {ligandFile && (
            <div className={`p-3 rounded ${isDark ? 'bg-gray-700' : 'bg-white'}`}>
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-sm font-medium">{ligandFile.name}</p>
                  <p className="text-xs text-gray-500">
                    {ligandFormat && <span className="inline-block px-1.5 py-0.5 bg-gray-100 rounded text-xs mr-1">{ligandFormat.toUpperCase()}</span>}
                    {ligandPreview ? `${ligandPreview.heavy_atoms} heavy atoms, MW: ${ligandPreview.mw}` : ''}
                  </p>
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
          <div className="flex items-center justify-between mb-3">
            <h3 className="font-semibold">Grid Box Configuration</h3>
            {receptorContent && (
              <button
                onClick={async () => {
                  try {
                    const res = await fetch('/api/docking/binding-site', {
                      method: 'POST',
                      headers: { 'Content-Type': 'application/json' },
                      body: JSON.stringify({ receptor_content: receptorContent })
                    })
                    const d = await res.json()
                    if (d.center_x != null) setGridConfig(d)
                  } catch { /* no-op */ }
                }}
                className={`text-xs px-3 py-1 rounded ${isDark ? 'bg-cyan-900 hover:bg-cyan-800 text-cyan-300' : 'bg-cyan-50 hover:bg-cyan-100 text-cyan-700'} border border-cyan-500`}
              >
                ⚡ Auto-detect
              </button>
            )}
          </div>
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
        
        <div className={`p-3 rounded-lg border ${isDark ? 'border-gray-600 bg-gray-800' : 'border-gray-200 bg-gray-50'}`}>
          <label className="flex items-center gap-3 cursor-pointer">
            <input
              type="checkbox"
              checked={enableFlexibility}
              onChange={e => setEnableFlexibility(e.target.checked)}
              className="w-5 h-5 text-cyan-600 rounded"
            />
            <div>
              <p className="font-medium text-sm">Enable Side-Chain Flexibility</p>
              <p className="text-xs text-gray-500">
                Generates multiple receptor conformations for better binding poses
              </p>
            </div>
          </label>
          {enableFlexibility && (
            <div className="mt-2 p-2 bg-yellow-50 border border-yellow-200 rounded text-xs text-yellow-700">
              ⚠️ Runtime will increase 3-27x. Estimated: ~{Math.ceil(exhaustiveness * numModes * 3 / 60)}-{Math.ceil(exhaustiveness * numModes * 27 / 60)} min
            </div>
          )}
        </div>
        
        <div className={`p-3 rounded-lg border ${isDark ? 'border-gray-600 bg-gray-800' : 'border-gray-200 bg-gray-50'}`}>
          <h4 className="font-medium text-sm mb-2">Docking Constraints</h4>
          <button
            onClick={() => setConstraints([...constraints, { type: 'hydrogen_bond', weight: 2.0 }])}
            className={`text-xs px-3 py-1 rounded ${isDark ? 'bg-cyan-900 hover:bg-cyan-800' : 'bg-cyan-100 hover:bg-cyan-200'}`}
          >
            + Add Constraint
          </button>
          <div className="space-y-2 mt-2">
            {constraints.map((c, idx) => (
              <div key={idx} className={`p-2 rounded border ${isDark ? 'border-gray-600' : 'border-gray-300'}`}>
                <div className="flex justify-between items-center mb-1">
                  <select
                    value={c.type}
                    onChange={e => {
                      const nc = [...constraints];
                      nc[idx] = { ...nc[idx], type: e.target.value };
                      setConstraints(nc);
                    }}
                    className={`text-xs p-1 rounded ${isDark ? 'bg-gray-700' : 'bg-white'}`}
                  >
                    <option value="hydrogen_bond">H-Bond</option>
                    <option value="positional">Positional</option>
                    <option value="metal_coordination">Metal Coordination</option>
                  </select>
                  <button
                    onClick={() => setConstraints(constraints.filter((_, i) => i !== idx))}
                    className="text-red-500 hover:text-red-400 text-xs"
                  >
                    ✕
                  </button>
                </div>
                {c.type === 'hydrogen_bond' && (
                  <>
                    <input type="number" placeholder="Ligand atom idx" value={c.ligand_atom_idx || ''}
                      onChange={e => { const nc = [...constraints]; nc[idx] = { ...nc[idx], ligand_atom_idx: parseInt(e.target.value) }; setConstraints(nc); }}
                      className={`w-full text-xs p-1 rounded mb-1 ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'} border`} />
                    <input type="text" placeholder="Receptor pattern (e.g. ASP:OD1)" value={c.receptor_pattern || ''}
                      onChange={e => { const nc = [...constraints]; nc[idx] = { ...nc[idx], receptor_pattern: e.target.value }; setConstraints(nc); }}
                      className={`w-full text-xs p-1 rounded ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'} border`} />
                  </>
                )}
                {c.type === 'positional' && (
                  <>
                    <input type="text" placeholder="Atom indices (comma-sep)" value={c.atom_indices?.join(',') || ''}
                      onChange={e => { const nc = [...constraints]; nc[idx] = { ...nc[idx], atom_indices: e.target.value.split(',').map(Number).filter((n: number) => !isNaN(n)) }; setConstraints(nc); }}
                      className={`w-full text-xs p-1 rounded mb-1 ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'} border`} />
                    <div className="grid grid-cols-3 gap-1 mb-1">
                      {['x', 'y', 'z'].map((axis, i) => (
                        <input key={axis} type="number" placeholder={`Center ${axis}`} value={c.center?.[i] || ''}
                          onChange={e => { const nc = [...constraints]; const center = [...(nc[idx].center || [0,0,0])]; center[i] = parseFloat(e.target.value) || 0; nc[idx] = { ...nc[idx], center }; setConstraints(nc); }}
                          className={`text-xs p-1 rounded ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'} border`} />
                      ))}
                    </div>
                    <input type="number" placeholder="Radius (Å)" value={c.radius || ''}
                      onChange={e => { const nc = [...constraints]; nc[idx] = { ...nc[idx], radius: parseFloat(e.target.value) }; setConstraints(nc); }}
                      className={`w-full text-xs p-1 rounded ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'} border`} />
                  </>
                )}
                {c.type === 'metal_coordination' && (
                  <>
                    <div className="grid grid-cols-3 gap-1 mb-1">
                      {['x', 'y', 'z'].map((axis, i) => (
                        <input key={axis} type="number" placeholder={`Metal ${axis}`} value={c.metal_coords?.[i] || ''}
                          onChange={e => { const nc = [...constraints]; const mc = [...(nc[idx].metal_coords || [0,0,0])]; mc[i] = parseFloat(e.target.value) || 0; nc[idx] = { ...nc[idx], metal_coords: mc }; setConstraints(nc); }}
                          className={`text-xs p-1 rounded ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'} border`} />
                      ))}
                    </div>
                    <input type="text" placeholder="Donor indices (comma-sep)" value={c.donor_indices?.join(',') || ''}
                      onChange={e => { const nc = [...constraints]; nc[idx] = { ...nc[idx], donor_indices: e.target.value.split(',').map(Number).filter((n: number) => !isNaN(n)) }; setConstraints(nc); }}
                      className={`w-full text-xs p-1 rounded ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-white border-gray-300'} border`} />
                  </>
                )}
                <div className="flex items-center gap-2 mt-1">
                  <span className="text-xs text-gray-500">Weight: {c.weight}</span>
                  <input type="range" min="0.5" max="5.0" step="0.5" value={c.weight || 2.0}
                    onChange={e => { const nc = [...constraints]; nc[idx] = { ...nc[idx], weight: parseFloat(e.target.value) }; setConstraints(nc); }}
                    className="flex-1" />
                </div>
              </div>
            ))}
          </div>
        </div>
        
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
          {(!selectedJob?.download_urls || Object.keys(selectedJob.download_urls).length === 0) && (
            <div className="col-span-2 text-center text-gray-500 py-2">No files available</div>
          )}
        </div>
          </div>
          
          {/* Simulated results banner (d10) */}
          {selectedJob?.engine === 'simulated' && (
            <div className="mb-4 px-3 py-2 rounded bg-yellow-900 border border-yellow-600 text-yellow-300 text-xs">
              ⚠ Simulated results — Vina/GNINA unavailable. Scores are illustrative only.
            </div>
          )}

          {/* Scores Table */}
          <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'}`}>
            <h3 className="font-semibold mb-3">🏆 Docking Poses <span className="text-xs font-normal text-gray-500">(click row to view pose)</span></h3>
            <table className="w-full text-sm">
              <thead>
                <tr className={`text-left ${isDark ? 'text-gray-400' : 'text-gray-600'}`}>
                  <th className="pb-2">#</th>
                  <th className="pb-2">Vina</th>
                  <th className="pb-2">GNINA</th>
                  <th className="pb-2">RF</th>
                  <th className="pb-2">Hydrophobic</th>
                  <th className="pb-2">Rotatable</th>
                  <th className="pb-2">Lipo</th>
                  <th className="pb-2">Composite</th>
                  {selectedJob?.results?.some((r: DockingResult) => r.constraint_penalty) && (
                    <th className="pb-2">Constraint</th>
                  )}
                </tr>
              </thead>
              <tbody>
                {selectedJob?.results?.map((r: DockingResult, i: number) => (
                  <tr
                    key={i}
                    onClick={() => {
                      setSelectedPoseIdx(i)
                      if (nglStageRef.current) {
                        nglStageRef.current.eachComponent((comp: any) => {
                          if (comp.name === 'docking_result') {
                            comp.setFrame(i)
                          }
                        })
                      }
                    }}
                    className={`border-t cursor-pointer ${
                      i === selectedPoseIdx
                        ? isDark ? 'bg-cyan-900' : 'bg-cyan-50'
                        : isDark ? 'border-gray-700 hover:bg-gray-700' : 'border-gray-100 hover:bg-gray-50'
                    }`}
                  >
                    <td className="py-2">{r.mode}</td>
                    <td className="py-2 font-mono">{r.vina_score?.toFixed(2)}</td>
                    <td className="py-2 font-mono">{r.gnina_score?.toFixed(2) || '-'}</td>
                    <td className="py-2 font-mono">{r.rf_score?.toFixed(2) || '-'}</td>
                    <td className="py-2 font-mono">{r.hydrophobic_term?.toFixed(3) || '-'}</td>
                    <td className="py-2 font-mono">{r.rotatable_penalty?.toFixed(3) || '-'}</td>
                    <td className="py-2 font-mono">{r.lipo_contact?.toFixed(3) || '-'}</td>
                    <td className="py-2 font-mono font-bold">{(r.final_score ?? r.composite_score ?? r.vina_score)?.toFixed(3)}</td>
                    {selectedJob?.results?.some((r: DockingResult) => r.constraint_penalty) && (
                      <td className="py-2 font-mono">{r.constraint_penalty?.toFixed(3) || '-'}</td>
                    )}
                  </tr>
                ))}
                {(!selectedJob?.results || selectedJob.results.length === 0) && (
                  <tr>
                    <td colSpan={9} className="py-4 text-center text-gray-500">No results</td>
                  </tr>
                )}
              </tbody>
            </table>
          </div>

          {/* Interaction Panel */}
          {selectedJob?.receptor_content && selectedJob?.ligand_pdb && (
            <InteractionPanel
              ligandPdb={selectedJob.ligand_pdb}
              receptorPdb={selectedJob.receptor_content}
              isDark={isDark}
            />
          )}

          {/* Download Files */}
          {selectedJob?.download_urls && Object.keys(selectedJob.download_urls).length > 0 && (
            <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-white'}`}>
              <h3 className="font-semibold mb-3">📥 Download Files</h3>
              <div className="flex flex-wrap gap-2">
                {Object.entries(selectedJob.download_urls).map(([key, url]: [string, any]) => (
                  <a
                    key={key}
                    href={url}
                    download
                    className="px-3 py-1.5 text-xs rounded bg-cyan-700 hover:bg-cyan-600 text-white font-mono"
                  >
                    ↓ {key}
                  </a>
                ))}
              </div>
            </div>
          )}

          {/* AI Job Assistant */}
          <div className={`rounded-lg border ${isDark ? 'bg-gray-800 border-gray-700' : 'bg-white border-gray-200'}`}>
            <button
              onClick={() => setShowAIPanel(v => !v)}
              className="w-full px-4 py-3 flex items-center justify-between text-sm font-semibold"
            >
              <span>🤖 Ask BioDockify AI about this job</span>
              <span className="text-gray-400">{showAIPanel ? '▲' : '▼'}</span>
            </button>
            {showAIPanel && (
              <div className="px-4 pb-4 space-y-3">
                <div className="flex gap-2">
                  <input
                    type="text"
                    value={aiJobId}
                    onChange={e => setAiJobId(e.target.value)}
                    placeholder="Job ID"
                    className={`w-40 px-2 py-1 text-xs rounded border font-mono ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-gray-50 border-gray-300'}`}
                  />
                  <span className="text-xs text-gray-500 self-center">(pre-filled)</span>
                </div>
                <textarea
                  value={aiQuestion}
                  onChange={e => setAiQuestion(e.target.value)}
                  onKeyDown={e => { if (e.key === 'Enter' && e.ctrlKey) askAI() }}
                  placeholder="Ask anything: explain my binding scores, why did it fail, is -8.2 kcal/mol a good result, what does GNINA score mean…"
                  rows={3}
                  className={`w-full px-3 py-2 text-xs rounded border resize-none ${isDark ? 'bg-gray-700 border-gray-600' : 'bg-gray-50 border-gray-300'}`}
                />
                <button
                  onClick={askAI}
                  disabled={aiLoading || !aiQuestion.trim()}
                  className="px-4 py-1.5 text-xs rounded bg-cyan-600 hover:bg-cyan-700 text-white disabled:opacity-50"
                >
                  {aiLoading ? '⏳ Thinking…' : 'Ask AI (Ctrl+Enter)'}
                </button>
                {aiAnswer && (
                  <div className={`p-3 rounded text-xs leading-relaxed whitespace-pre-wrap border ${isDark ? 'bg-gray-900 border-gray-700 text-gray-200' : 'bg-gray-50 border-gray-200 text-gray-800'}`}>
                    {aiAnswer}
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
        
        {/* Right - 3D Viewer */}
        <div className="w-1/2 flex flex-col">
          <div className={`px-4 py-2 ${isDark ? 'bg-gray-800' : 'bg-gray-100'} flex items-center justify-between`}>
            <span className="text-sm font-medium">3D Viewer</span>
            <div className="flex gap-2">
              <button onClick={() => setShowCartoon(v => !v)} className={`px-2 py-1 text-xs rounded ${showCartoon ? 'bg-cyan-600 text-white' : isDark ? 'bg-gray-700' : 'bg-white'}`}>
                Cartoon
              </button>
              <button onClick={() => setShowSurface(v => !v)} className={`px-2 py-1 text-xs rounded ${showSurface ? 'bg-cyan-600 text-white' : isDark ? 'bg-gray-700' : 'bg-white'}`}>
                Surface
              </button>
              <button onClick={() => { nglStageRef.current = null; viewerLoaded.current = false; init3DViewer() }} className={`px-2 py-1 text-xs rounded ${isDark ? 'bg-gray-700 hover:bg-gray-600' : 'bg-white hover:bg-gray-200'}`}>
                Reset
              </button>
            </div>
          </div>
          <div ref={viewer3dRef} className="flex-1" style={{ background: '#1a1a2e', minHeight: '400px' }}>
            <div className="flex items-center justify-center h-full text-gray-500">
              Loading 3D viewer...
            </div>
          </div>
          <ExportPanel viewerRef={viewer3dRef} isDark={isDark} />
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
