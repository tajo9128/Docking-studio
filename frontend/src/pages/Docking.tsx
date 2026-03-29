import { useState, useCallback } from 'react'
import { useNavigate } from 'react-router-dom'
import { Card, Button, Input, ProgressBar, Badge, Tabs, TabPanel } from '@/components/ui'
import { useDockingStream } from '@/hooks'
import { startDocking, cancelDocking } from '@/api/docking'
import { uploadFile, downloadFile } from '@/api/upload'
import { prepareProtein } from '@/api/rdkit'
import type { DockingConfig } from '@/lib/types'

export function Docking() {
  const navigate = useNavigate()
  const [activeTab, setActiveTab] = useState('input')
  const [receptorFile, setReceptorFile] = useState<File | null>(null)
  const [ligandFiles, setLigandFiles] = useState<File[]>([])
  const [receptorName, setReceptorName] = useState('')
  const [ligandName, setLigandName] = useState('')
  const [uploadError, setUploadError] = useState<string | null>(null)
  const [currentJobId, setCurrentJobId] = useState<string | null>(null)
  const [preparingReceptor, setPreparingReceptor] = useState(false)
  const [preparedReceptorPath, setPreparedReceptorPath] = useState<string | null>(null)
  const [receptorPrepInfo, setReceptorPrepInfo] = useState<{watersRemoved: number; hydrogensAdded: number} | null>(null)
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

  const handleReceptorUpload = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file) return
    setReceptorFile(file)
    setReceptorName(file.name)
    setPreparedReceptorPath(null)
    setReceptorPrepInfo(null)
    setUploadError(null)
  }, [])

  const handleLigandUpload = useCallback(async (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(e.target.files || [])
    if (!files.length) return
    setLigandFiles(files)
    setLigandName(`${files.length} files selected`)
    setUploadError(null)
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

  const handleStartDocking = async () => {
    if (!receptorFile || ligandFiles.length === 0) return

    setUploadError(null)
    try {
      let receptorPath: string
      if (preparedReceptorPath) {
        receptorPath = preparedReceptorPath
      } else {
        const uploadResult = await uploadFile(receptorFile)
        receptorPath = uploadResult.path
      }
      
      const ligandUploadResult = await uploadFile(ligandFiles[0])
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

      <Tabs tabs={tabs} activeTab={activeTab} onChange={setActiveTab} />

      <TabPanel>
        {activeTab === 'input' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mt-6">
            {/* Receptor Card */}
            <Card padding="lg">
              <div className="flex items-center gap-3 mb-4">
                <span className="w-8 h-8 rounded-full bg-primary-50 text-primary flex items-center justify-center font-bold">1</span>
                <div>
                  <h3 className="font-bold text-text-primary">Target Receptor</h3>
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

              {receptorName && !preparedReceptorPath && (
                <Button
                  variant="primary"
                  className="w-full mt-4"
                  onClick={handlePrepareReceptor}
                  disabled={preparingReceptor}
                >
                  {preparingReceptor ? '⏳ Preparing...' : '🧪 Prepare Receptor'}
                </Button>
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
                <div>
                  <h3 className="font-bold text-text-primary">Ligand Library</h3>
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
              <h3 className="font-bold text-text-primary mb-4">🔬 BioDockify Tri-Score Protocol</h3>
              <p className="text-xs text-text-tertiary mb-4">International Standard for molecular docking</p>

              <div className="space-y-3">
                {[
                  { id: 'vina', label: 'Vina', desc: 'Physics-based scoring', checked: true },
                  { id: 'gnina', label: 'GNINA', desc: 'Deep learning CNN', checked: false },
                  { id: 'rfscore', label: 'RF-Score', desc: 'Random Forest ML', checked: false },
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
                    <div
                      className="w-4 h-4 rounded-full border-2 flex items-center justify-center"
                      style={{
                        borderColor: config.engine === engine.id ? '#2e5aac' : '#e2e8f0',
                      }}
                    >
                      {config.engine === engine.id && (
                        <div className="w-2 h-2 rounded-full bg-primary" />
                      )}
                    </div>
                    <div className="flex-1">
                      <p className="font-semibold text-text-primary">{engine.label}</p>
                      <p className="text-xs text-text-tertiary">{engine.desc}</p>
                    </div>
                  </label>
                ))}
              </div>

              <p className="text-xs text-text-tertiary mt-4 italic">
                Powered by Vina (Physics) + GNINA (Deep Learning) + RF-Score (ML)
              </p>
            </Card>

            {/* Grid Configuration */}
            <Card padding="lg">
              <h3 className="font-bold text-text-primary mb-4">📦 Grid Box Configuration</h3>

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
              <h3 className="font-bold text-text-primary mb-4">🖥 Hardware Status</h3>
              <div className="flex items-center gap-4 p-4 bg-surface-secondary rounded-lg mb-4">
                <span className="text-3xl">🖥️</span>
                <div>
                  <p className="font-semibold text-text-primary">GPU Acceleration</p>
                  <p className="text-sm text-text-secondary">Powered by NVIDIA CUDA</p>
                </div>
              </div>
              <p className="text-sm text-text-secondary">
                The system will automatically use GPU acceleration when available for faster docking calculations.
              </p>
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
