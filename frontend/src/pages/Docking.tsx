import { useState, useCallback } from 'react'
import { useNavigate } from 'react-router-dom'
import { Card, Button, Input, ProgressBar, Badge, Tabs, TabPanel } from '@/components/ui'
import { useDockingStream } from '@/hooks'
import { startDocking, cancelDocking } from '@/api/docking'
import { uploadFile, downloadFile } from '@/api/upload'
import { prepareProtein, prepareReceptorPDBQT, prepareLigand, smilesToSDF } from '@/api/rdkit'
import type { DockingConfig } from '@/lib/types'

const SAMPLE_RECEPTOR_PDB = `HEADER    IMMUNE SYSTEM                           15-APR-98   1HIA
TITLE     CRYSTAL STRUCTURE OF HUMAN HEAT-LABILE ENTEROTOXIN
EXPDTA    X-RAY DIFFRACTION
HETATM    1  N   MET A   1       9.287  11.614  -1.060  1.00 28.12           N
HETATM    2  CA  MET A   1      10.286  11.038  -1.989  1.00 27.95           C
HETATM    3  C   MET A   1      11.669  11.451  -1.661  1.00 28.26           C
HETATM    4  O   MET A   1      12.039  12.619  -1.589  1.00 30.27           O
HETATM    5  N   ILE A   2      12.516  10.421  -1.472  1.00 27.21           N
HETATM    6  CA  ILE A   2      13.922  10.668  -1.148  1.00 27.26           C
HETATM    7  C   ILE A   2      14.802  10.296  -2.322  1.00 27.85           C
HETATM    8  O   ILE A   2      15.995  10.583  -2.399  1.00 29.59           O
HETATM    9  N   GLN A   3      14.261   9.687  -3.311  1.00 26.84           N
HETATM   10  CA  GLN A   3      14.996   9.257  -4.500  1.00 27.01           C
HETATM   11  C   GLN A   3      14.103   8.443  -5.420  1.00 28.07           C
HETATM   12  O   GLN A   3      13.020   7.954  -5.086  1.00 28.77           O
HETATM   13  N   VAL A   4      14.578   8.279  -6.640  1.00 26.58           N
HETATM   14  CA  VAL A   4      13.816   7.510  -7.636  1.00 26.74           C
HETATM   15  C   VAL A   4      12.339   7.896  -7.742  1.00 27.03           C
HETATM   16  O   VAL A   4      11.513   7.663  -6.843  1.00 27.72           O
HETATM   17  N   GLU A   5      12.000   8.508  -8.879  1.00 26.07           N
HETATM   18  CA  GLU A   5      10.630   8.932  -9.112  1.00 26.23           C
HETATM   19  C   GLU A   5      10.523  10.453  -9.246  1.00 27.04           C
HETATM   20  O   GLU A   5       9.424  11.008  -9.378  1.00 28.03           O
HETATM   21  N   TYR A   6       9.876  11.657 -10.220  1.00 25.97           N
HETATM   22  CA  TYR A   6       9.847  13.112 -10.446  1.00 26.46           C
HETATM   23  C   TYR A   6       8.484  13.715 -10.165  1.00 27.33           C
HETATM   24  O   TYR A   6       7.453  13.062  -9.997  1.00 27.68           O
HETATM   25  N   CYS A   7       8.484  15.049 -10.107  1.00 26.46           N
HETATM   26  CA  CYS A   7       7.236  15.757  -9.845  1.00 27.05           C
HETATM   27  C   CYS A   7       7.390  17.269  -9.673  1.00 27.96           C
HETATM   28  O   CYS A   7       8.511  17.810  -9.730  1.00 29.03           O
CONECT    1    2
CONECT    2    1    3
CONECT    3    2    4    5
END
`

const SAMPLE_LIGAND_SMI = 'CC(=O)Oc1ccccc1C(=O)O'

interface TooltipProps {
  text: string
  children: React.ReactNode
}

function Tooltip({ text, children }: TooltipProps) {
  return (
    <div className="relative inline-block group">
      {children}
      <div className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 bg-gray-900 text-white text-xs rounded-lg w-64 opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-50 shadow-lg">
        {text}
        <div className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-900"></div>
      </div>
    </div>
  )
}

export function Docking() {
  const navigate = useNavigate()
  const [activeTab, setActiveTab] = useState('input')
  const [receptorFile, setReceptorFile] = useState<File | null>(null)
  const [ligandFiles, setLigandFiles] = useState<File[]>([])
  const [ligandFile, setLigandFile] = useState<File | null>(null)
  const [receptorName, setReceptorName] = useState('')
  const [ligandName, setLigandName] = useState('')
  const [uploadError, setUploadError] = useState<string | null>(null)
  const [currentJobId, setCurrentJobId] = useState<string | null>(null)
  const [preparingReceptor, setPreparingReceptor] = useState(false)
  const [preparedReceptorPath, setPreparedReceptorPath] = useState<string | null>(null)
  const [receptorPrepInfo, setReceptorPrepInfo] = useState<{watersRemoved: number; hydrogensAdded: number} | null>(null)
  const [preparingReceptorPDBQT, setPreparingReceptorPDBQT] = useState(false)
  const [preparedReceptorPDBQTPath, setPreparedReceptorPDBQTPath] = useState<string | null>(null)
  const [preparedLigandPath, setPreparedLigandPath] = useState<string | null>(null)
  const [preparingLigand, setPreparingLigand] = useState(false)
  const [showSampleDataHint, setShowSampleDataHint] = useState(true)
  const [config, setConfig] = useState<DockingConfig>({
    center_x: 0,
    center_y: 0,
    center_z: 0,
    size_x: 20,
    size_y: 20,
    size_z: 20,
    exhaustiveness: 8,
    num_modes: 9,
    engine: 'vina',
    batch_size: 5,
  })

  const { progress, error } = useDockingStream(currentJobId)

  const handleLoadSampleData = useCallback(() => {
    const receptorBlob = new Blob([SAMPLE_RECEPTOR_PDB], { type: 'text/plain' })
    const receptorFakeFile = new File([receptorBlob], '1HIA_protein.pdb', { type: 'text/plain' })
    setReceptorFile(receptorFakeFile)
    setReceptorName('1HIA_protein.pdb (Sample)')
    
    const ligandBlob = new Blob([SAMPLE_LIGAND_SMI], { type: 'text/plain' })
    const ligandFakeFile = new File([ligandBlob], 'aspirin.smi', { type: 'text/plain' })
    setLigandFiles([ligandFakeFile])
    setLigandName('Aspirin (Sample)')
    setShowSampleDataHint(false)
    setUploadError(null)
  }, [])

  const handleReceptorUpload = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    setReceptorFile(file)
    setReceptorName(file.name)
    setPreparedReceptorPath(null)
    setReceptorPrepInfo(null)
    setPreparedReceptorPDBQTPath(null)
    setUploadError(null)
    setShowSampleDataHint(false)
  }, [])

  const handleLigandUpload = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(e.target.files || [])
    if (!files.length) return
    setLigandFiles(files)
    setLigandFile(files[0])
    setLigandName(`${files.length} files selected`)
    setPreparedLigandPath(null)
    setUploadError(null)
    setShowSampleDataHint(false)
  }, [])

  const handlePrepareReceptor = async () => {
    if (!receptorFile) return
    setPreparingReceptor(true)
    setUploadError(null)
    try {
      const uploadResult = await uploadFile(receptorFile)
      const receptorPath = uploadResult.path
      const fileContent = await downloadFile(receptorPath)
      const pdbContent = typeof fileContent === 'string' ? fileContent : fileContent.content || fileContent
      const prepResult = await prepareProtein(
        pdbContent as string,
        receptorName.replace(/\.[^.]+$/, ''),
        true,
        true
      )
      if (prepResult.success) {
        setPreparedReceptorPath(prepResult.pdb_path)
        setReceptorPrepInfo({
          watersRemoved: prepResult.original_atoms - prepResult.final_atoms,
          hydrogensAdded: prepResult.final_atoms - prepResult.original_atoms,
        })
      }
    } catch (err) {
      console.error('Failed to prepare receptor:', err)
      setUploadError(err instanceof Error ? err.message : 'Failed to prepare receptor')
    } finally {
      setPreparingReceptor(false)
    }
  }

  const handlePrepareReceptorPDBQT = async () => {
    if (!receptorFile) return
    setPreparingReceptorPDBQT(true)
    setUploadError(null)
    try {
      const uploadResult = await uploadFile(receptorFile)
      const receptorPath = uploadResult.path
      const fileContent = await downloadFile(receptorPath)
      const pdbContent = typeof fileContent === 'string' ? fileContent : fileContent.content || fileContent
      const prepResult = await prepareReceptorPDBQT(
        pdbContent as string,
        receptorName.replace(/\.[^.]+$/, ''),
        true
      )
      if (prepResult.success) {
        setPreparedReceptorPDBQTPath(prepResult.pdbqt_path)
      } else {
        throw new Error('Receptor PDBQT preparation failed')
      }
    } catch (err) {
      console.error('Failed to prepare receptor PDBQT:', err)
      setUploadError(err instanceof Error ? err.message : 'Failed to prepare receptor PDBQT')
    } finally {
      setPreparingReceptorPDBQT(false)
    }
  }

  const handlePrepareLigand = async () => {
    if (ligandFiles.length === 0) return
    setPreparingLigand(true)
    setUploadError(null)
    try {
      const currentLigandFile = ligandFile || ligandFiles[0]
      
      if (currentLigandFile.name.toLowerCase().endsWith('.smi')) {
        const content = await downloadFile(currentLigandFile)
        const smiles = content.trim().split('\n')[0]
        
        const sdfResult = await smilesToSDF(smiles, currentLigandFile.name)
        
        if (sdfResult.sdf_content) {
          const sdfFileName = currentLigandFile.name.replace(/\.smi$/i, '.sdf')
          const prepResult = await prepareLigand(sdfResult.sdf_content, sdfFileName)
          if (prepResult.success && prepResult.pdbqt_path) {
            setPreparedLigandPath(prepResult.pdbqt_path)
          } else {
            throw new Error(prepResult.message || 'Ligand preparation failed')
          }
        } else {
          throw new Error('Failed to convert SMILES to SDF format')
        }
        setPreparingLigand(false)
        return
      }
      
      setLigandFile(currentLigandFile)
      const content = await downloadFile(currentLigandFile)
      const pdbContent = typeof content === 'string' ? content : content.content || content
      const prepResult = await prepareLigand(pdbContent as string, ligandFiles[0].name)
      if (prepResult.success && prepResult.pdbqt_path) {
        setPreparedLigandPath(prepResult.pdbqt_path)
      } else {
        throw new Error(prepResult.message || 'Ligand preparation failed')
      }
    } catch (err) {
      console.error('Failed to prepare ligand:', err)
      setUploadError(err instanceof Error ? err.message : 'Failed to prepare ligand')
    } finally {
      setPreparingLigand(false)
    }
  }

  const handleStartDocking = async () => {
    if (!receptorFile) {
      setUploadError('Please upload a receptor file first')
      return
    }

    const ligandToUse = preparedLigandPath || (ligandFiles.length > 0 ? ligandFiles[0] : null)
    if (!ligandToUse) {
      setUploadError('Please upload a ligand file')
      return
    }

    if (!preparedReceptorPDBQTPath && !preparedReceptorPath) {
      setUploadError('Please prepare the receptor first')
      return
    }

    setUploadError(null)
    try {
      let receptorPath: string
      if (preparedReceptorPDBQTPath) {
        receptorPath = preparedReceptorPDBQTPath
      } else if (preparedReceptorPath) {
        receptorPath = preparedReceptorPath
      } else {
        const uploadResult = await uploadFile(receptorFile)
        receptorPath = uploadResult.path
      }
      
      const ligandUploadResult = await uploadFile(ligandToUse)
      const ligandPath = ligandUploadResult.path

      const result = await startDocking(receptorPath, ligandPath, {
        center_x: config.center_x,
        center_y: config.center_y,
        center_z: config.center_z,
        size_x: config.size_x,
        size_y: config.size_y,
        size_z: config.size_z,
        exhaustiveness: config.exhaustiveness,
        num_modes: config.num_modes,
        engine: config.engine
      })
      setCurrentJobId(result.job_id)
      setActiveTab('progress')
    } catch (err) {
      console.error('Failed to start docking:', err)
      setUploadError(err instanceof Error ? err.message : 'Failed to start docking')
    }
  }

  const handleCancel = async () => {
    if (!currentJobId) return
    try {
      await cancelDocking(currentJobId)
      setCurrentJobId(null)
      setActiveTab('input')
    } catch (err) {
      console.error('Failed to cancel:', err)
    }
  }

  const tabs = [
    { id: 'input', label: '📁 Input Files' },
    { id: 'settings', label: '⚙ Docking Settings' },
    { id: 'hardware', label: '🖥 Hardware' },
    { id: 'progress', label: '📊 Progress' },
  ]

  return (
    <div className="p-6">
      {/* Header */}
      <div className="mb-6">
        <h1 className="text-2xl font-bold text-text-primary">New Docking Experiment</h1>
        <p className="text-text-secondary mt-1">Configure and run molecular docking simulation</p>
      </div>

      {showSampleDataHint && (
        <Card padding="lg" className="mb-6 bg-gradient-to-r from-blue-50 to-purple-50 border-blue-200">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-4">
              <span className="text-4xl">🚀</span>
              <div>
                <h3 className="font-bold text-text-primary">New to Molecular Docking?</h3>
                <p className="text-sm text-text-secondary mt-1">
                  Start quickly with sample data — protein + ligand pre-loaded for you!
                </p>
              </div>
            </div>
            <Button variant="primary" onClick={handleLoadSampleData}>
              Load Sample Data
            </Button>
          </div>
        </Card>
      )}

      <Tabs tabs={tabs} activeTab={activeTab} onChange={setActiveTab} />

      <TabPanel>
        {activeTab === 'input' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mt-6">
            {/* Receptor Card */}
            <Card padding="lg">
              <div className="flex items-center gap-3 mb-4">
                <span className="w-8 h-8 rounded-full bg-primary-50 text-primary flex items-center justify-center font-bold">1</span>
                <div className="flex-1">
                  <div className="flex items-center gap-2">
                    <h3 className="font-bold text-text-primary">Target Receptor</h3>
                    <Tooltip text="The receptor is the protein target (usually from PDB). We prepare it by removing water molecules and adding hydrogen atoms before docking.">
                      <span className="text-gray-400 cursor-help text-xs">ℹ️</span>
                    </Tooltip>
                  </div>
                  <p className="text-xs text-text-tertiary">Protein structure (.pdb, .pdbqt)</p>
                </div>
              </div>

              <label
                className="border-2 border-dashed border-border-light rounded-xl p-8 text-center cursor-pointer hover:border-primary hover:bg-primary-50/30 transition-all block"
                style={{ borderColor: '#e2e8f0' }}
              >
                <input
                  type="file"
                  accept=".pdb,.pdbqt,.mol2,.mmcif"
                  className="hidden"
                  onChange={handleReceptorUpload}
                />
                <div className="text-3xl mb-2">📁</div>
                <p className="font-semibold text-text-primary">Upload Receptor</p>
                <p className="text-xs text-text-tertiary mt-1">.pdb, .mmcif, .mol2</p>
              </label>

              {receptorName && (
                <div className="mt-4 p-3 bg-surface-secondary rounded-lg flex items-center gap-3">
                  <span>📄</span>
                  <span className="text-sm font-medium text-text-primary flex-1 truncate">{receptorName}</span>
                  <Badge variant="success">Ready</Badge>
                </div>
              )}

              {receptorName && !preparedReceptorPath && !preparedReceptorPDBQTPath && (
                <Button
                  variant="primary"
                  className="w-full mt-4"
                  onClick={handlePrepareReceptor}
                  disabled={preparingReceptor}
                >
                  {preparingReceptor ? '⏳ Preparing...' : '🧪 Prepare Receptor'}
                </Button>
              )}

              {receptorName && !preparedReceptorPDBQTPath && (
                <Button
                  variant="primary"
                  className="w-full mt-2"
                  onClick={handlePrepareReceptorPDBQT}
                  disabled={preparingReceptorPDBQT}
                >
                  {preparingReceptorPDBQT ? '⏳ Preparing AMBER...' : '⚡ Prepare Receptor (AMBER)'}
                </Button>
              )}

              {preparedReceptorPDBQTPath && (
                <div className="mt-4 p-3 bg-green-50 border border-green-200 rounded-lg">
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="success">✓ AMBER Prepared</Badge>
                  </div>
                  <p className="text-xs text-green-700">
                    Receptor ready for Vina/GNINA with AMBER charges
                  </p>
                </div>
              )}

              {preparedReceptorPath && (
                <div className="mt-4 p-3 bg-green-50 border border-green-200 rounded-lg">
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="success">✓ Prepared</Badge>
                  </div>
                  {receptorPrepInfo && (
                    <p className="text-xs text-green-700">
                      Waters removed: {receptorPrepInfo.watersRemoved > 0 ? receptorPrepInfo.watersRemoved : 0}, 
                      Hydrogens added: {receptorPrepInfo.hydrogensAdded}
                    </p>
                  )}
                </div>
              )}

              <Button
                variant="outline"
                className="w-full mt-4"
                onClick={() => document.querySelector<HTMLInputElement>('input[type=file]')?.click()}
              >
                Select Receptor File
              </Button>
            </Card>

            {/* Ligand Card */}
            <Card padding="lg">
              <div className="flex items-center gap-3 mb-4">
                <span className="w-8 h-8 rounded-full bg-success-bg text-success flex items-center justify-center font-bold">2</span>
                <div className="flex-1">
                  <div className="flex items-center gap-2">
                    <h3 className="font-bold text-text-primary">Ligand Molecule</h3>
                    <Tooltip text="The ligand is the small molecule drug or compound you want to test against the receptor. Can be a single molecule or a library.">
                      <span className="text-gray-400 cursor-help text-xs">ℹ️</span>
                    </Tooltip>
                  </div>
                  <p className="text-xs text-text-tertiary">Small molecules, drugs (.sdf, .mol2, .smi)</p>
                </div>
              </div>

              <label
                className="border-2 border-dashed border-border-light rounded-xl p-8 text-center cursor-pointer hover:border-success hover:bg-success-bg/30 transition-all block"
                style={{ borderColor: '#e2e8f0' }}
              >
                <input
                  type="file"
                  accept=".sdf,.mol2,.pdbqt,.smi,.smiles"
                  multiple
                  className="hidden"
                  onChange={handleLigandUpload}
                />
                <div className="text-3xl mb-2">🧪</div>
                <p className="font-semibold text-text-primary">Upload Ligands</p>
                <p className="text-xs text-text-tertiary mt-1">.sdf, .mol2, .smi</p>
              </label>

              {ligandName && (
                <div className="mt-4 p-3 bg-surface-secondary rounded-lg flex items-center gap-3">
                  <span>📄</span>
                  <span className="text-sm font-medium text-text-primary flex-1 truncate">{ligandName}</span>
                  <Badge variant="success">Ready</Badge>
                </div>
              )}

              {ligandName && !preparedLigandPath && (
                <Button
                  variant="primary"
                  className="w-full mt-4"
                  onClick={handlePrepareLigand}
                  disabled={preparingLigand}
                >
                  {preparingLigand ? '⏳ Preparing Ligand...' : '⚡ Prepare Ligand (PDBQT)'}
                </Button>
              )}

              {preparedLigandPath && (
                <div className="mt-4 p-3 bg-green-50 border border-green-200 rounded-lg">
                  <div className="flex items-center gap-2 mb-2">
                    <Badge variant="success">✓ Prepared</Badge>
                  </div>
                  <p className="text-xs text-green-700">Ligand ready for docking</p>
                </div>
              )}

              <Button
                variant="outline"
                className="w-full mt-4"
                onClick={() => document.querySelectorAll<HTMLInputElement>('input[type=file]')[1]?.click()}
              >
                Select Ligand Files
              </Button>
            </Card>
          </div>
        )}

        {activeTab === 'settings' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mt-6">
            {/* Engine Selection */}
            <Card padding="lg">
              <div className="flex items-center gap-2 mb-2">
                <h3 className="font-bold text-text-primary">🔬 Docking Engine</h3>
                <Tooltip text="AutoDock Vina uses physics-based scoring. GNINA adds deep learning CNN scores. Consensus combines both for more reliable predictions.">
                  <span className="text-gray-400 cursor-help text-xs">ℹ️</span>
                </Tooltip>
              </div>
              <p className="text-xs text-text-secondary mb-4">Choose how to score binding poses</p>

              <div className="space-y-3">
                {[
                  { id: 'vina', label: 'AutoDock Vina', desc: 'Physics-based scoring', detail: 'Fast, well-validated', icon: '⚡' },
                  { id: 'gnina', label: 'GNINA', desc: 'Deep learning CNN', detail: 'Neural network accuracy', icon: '🧠' },
                  { id: 'consensus', label: 'Consensus (Recommended)', desc: 'Vina + GNINA combined', detail: 'Best reliability', icon: '🎯' },
                ].map((engine) => (
                  <label
                    key={engine.id}
                    className="flex items-center gap-3 p-3 rounded-lg border cursor-pointer transition-all"
                    style={{
                      borderColor: config.engine === engine.id ? '#2e5aac' : '#e2e8f0',
                      background: config.engine === engine.id ? '#f0f4ff' : 'white',
                    }}
                  >
                    <input
                      type="radio"
                      name="engine"
                      value={engine.id}
                      checked={config.engine === engine.id}
                      onChange={(e) => setConfig({ ...config, engine: e.target.value as any })}
                      className="sr-only"
                    />
                    <div className="text-2xl">{engine.icon}</div>
                    <div className="flex-1">
                      <p className="font-semibold text-text-primary">{engine.label}</p>
                      <p className="text-xs text-text-primary">{engine.desc}</p>
                      <p className="text-xs text-text-tertiary">{engine.detail}</p>
                    </div>
                  </label>
                ))}
              </div>
            </Card>

            {/* Grid Configuration */}
            <Card padding="lg">
              <div className="flex items-center gap-2 mb-2">
                <h3 className="font-bold text-text-primary">📦 Search Space (Grid Box)</h3>
                <Tooltip text="The grid box defines where Vina will search for ligand binding positions. The box should cover the active site or region of interest.">
                  <span className="text-gray-400 cursor-help text-xs">ℹ️</span>
                </Tooltip>
              </div>
              <p className="text-xs text-text-secondary mb-4">Define the region to search for binding</p>

              <div className="space-y-4">
                <div>
                  <p className="text-sm font-semibold text-text-secondary mb-2">Center (Å)</p>
                  <div className="grid grid-cols-3 gap-3">
                    {['center_x', 'center_y', 'center_z'].map((key, i) => (
                      <Input
                        key={key}
                        label={`${['X', 'Y', 'Z'][i]}`}
                        type="number"
                        value={config[key as keyof typeof config]}
                        onChange={(e) => setConfig({ ...config, [key]: parseFloat(e.target.value) || 0 })}
                      />
                    ))}
                  </div>
                </div>

                <div>
                  <p className="text-sm font-semibold text-text-secondary mb-2">Size (Å)</p>
                  <div className="grid grid-cols-3 gap-3">
                    {['size_x', 'size_y', 'size_z'].map((key, i) => (
                      <Input
                        key={key}
                        label={`${['X', 'Y', 'Z'][i]}`}
                        type="number"
                        value={config[key as keyof typeof config]}
                        onChange={(e) => setConfig({ ...config, [key]: parseInt(e.target.value) || 20 })}
                      />
                    ))}
                  </div>
                </div>

                <div>
                  <Input
                    label="Exhaustiveness"
                    type="number"
                    value={config.exhaustiveness}
                    min={1}
                    max={32}
                    onChange={(e) => setConfig({ ...config, exhaustiveness: parseInt(e.target.value) || 8 })}
                    hint="Sampling exhaustiveness (1-32)"
                  />
                </div>

                <div>
                  <Input
                    label="Batch Size"
                    type="number"
                    value={config.batch_size}
                    min={1}
                    max={10}
                    onChange={(e) => setConfig({ ...config, batch_size: parseInt(e.target.value) || 5 })}
                    hint="Ligands per batch (1-10)"
                  />
                </div>
              </div>
            </Card>
          </div>
        )}

        {activeTab === 'hardware' && (
          <div className="mt-6">
            <Card padding="lg">
              <div className="flex items-center gap-2 mb-4">
                <h3 className="font-bold text-text-primary">🖥 Hardware Acceleration</h3>
                <Tooltip text="GPU acceleration significantly speeds up docking calculations using NVIDIA CUDA.">
                  <span className="text-gray-400 cursor-help text-xs">ℹ️</span>
                </Tooltip>
              </div>
              <div className="flex items-center gap-4 p-4 bg-surface-secondary rounded-lg mb-4">
                <span className="text-3xl">🚀</span>
                <div>
                  <p className="font-semibold text-text-primary">GPU Acceleration</p>
                  <p className="text-sm text-text-secondary">Powered by NVIDIA CUDA</p>
                </div>
              </div>
              <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
                <h4 className="font-semibold text-blue-800 mb-2">How Docking Works</h4>
                <p className="text-xs text-blue-700 mb-2">
                  Molecular docking predicts how a small molecule (ligand) binds to a protein (receptor) by searching through possible orientations and scoring them.
                </p>
                <p className="text-xs text-blue-600">
                  <strong>AutoDock Vina</strong> uses a smart search algorithm. <strong>GNINA</strong> adds deep learning to improve accuracy.
                </p>
              </div>
            </Card>
          </div>
        )}

        {activeTab === 'progress' && (
          <div className="mt-6">
            <Card padding="lg">
              <h3 className="font-bold text-text-primary mb-4">📊 Docking Progress</h3>

              {progress ? (
                <div className="space-y-4">
                  <ProgressBar value={progress.progress} size="lg" />
                  <div className="flex items-center justify-between">
                    <Badge
                      variant={
                        progress.status === 'completed' ? 'success' :
                        progress.status === 'failed' ? 'error' :
                        progress.status === 'cancelled' ? 'warning' : 'info'
                      }
                    >
                      {progress.status.toUpperCase()}
                    </Badge>
                    <span className="text-sm text-text-secondary">{progress.message}</span>
                  </div>

                  {progress.status === 'completed' && (
                    <div className="flex gap-3 mt-4">
                      <Button onClick={() => navigate(`/results?job=${currentJobId}`)}>View Results</Button>
                      <Button variant="outline" onClick={() => setActiveTab('input')}>New Experiment</Button>
                    </div>
                  )}

                  {progress.status === 'running' && (
                    <Button variant="danger" onClick={handleCancel} className="mt-4">
                      Cancel Experiment
                    </Button>
                  )}
                </div>
              ) : (
                <div className="text-center py-8 text-text-tertiary">
                  <p>No active docking job</p>
                  <p className="text-xs mt-1">Start a new experiment from the Input Files tab</p>
                </div>
              )}

              {error && (
                <div className="mt-4 p-3 bg-red-50 border border-red-200 rounded-lg text-red-700 text-sm">
                  {error}
                </div>
              )}

              {uploadError && (
                <div className="mt-4 p-3 bg-red-50 border border-red-200 rounded-lg text-red-700 text-sm">
                  {uploadError}
                </div>
              )}
            </Card>
          </div>
        )}

        {/* Start Button - always visible */}
        {activeTab !== 'progress' && (
          <div className="mt-6 flex items-center justify-between p-4 bg-surface-secondary rounded-xl">
            <div className="text-sm text-text-secondary">
              {receptorFile && ligandFiles.length > 0 ? (
                <span className="text-success">✓ Ready to start docking</span>
              ) : (
                <span>Upload receptor and ligand files to begin</span>
              )}
            </div>
            <Button
              onClick={handleStartDocking}
              disabled={!receptorFile || ligandFiles.length === 0 || !!uploadError}
              className="bg-gradient-to-r from-primary to-secondary"
            >
              {uploadError ? 'Error' : '▶ Start Docking'}
            </Button>
          </div>
        )}
      </TabPanel>
    </div>
  )
}
