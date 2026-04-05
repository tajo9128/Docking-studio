import { useState, useEffect, useRef } from 'react'
import { Card, Button, Tabs, TabPanel } from '@/components/ui'
import Plot from 'react-plotly.js'
import {
  getDescriptorGroups,
  uploadDataset,
  startTraining,
  getTrainingStatus,
  predictSingle,
  predictBatch,
  listModels,
  deleteModel,
  type DatasetUploadResult,
  type TrainResult,
  type SavedModel,
  type PredictionSingle,
} from '@/api/qsar'

const MODEL_TYPES = ['RandomForest', 'GradientBoosting', 'SVR', 'PLS', 'Ridge', 'Lasso']

const SAMPLE_DATASET_CSV = `smiles,activity
CC(=O)OC1=CC=CC=C1C(=O)O,4.8239
Cn1cnc2c1c(=O)n(c(=O)n2C)C,5.3979
CC(C)Cc1ccc(cc1)C(C)C(=O)O,3.9707
CCO,-0.3
CO,-0.77
CC(=O)C,-0.24
c1ccc(cc1)O,2.1666
Nc1ccccc1,2.0622
OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O,-3.0
c1ccccc1,1.0
CC(=O)Oc1ccccc1OC(=O)C,3.5
c1ccc2c(c1)ccc3c2cccc3,3.2
CC1=CC=C(C)C=1C(=O)O,4.0
Nc1ccc(cc1)S(=O)(=O)N,3.8
CC1=CC(=CC=C1C)C(C)C(=O)O,4.2
`

const TABS = [
  { id: 'dataset', label: '📁 Dataset' },
  { id: 'train', label: '🧠 Train Model' },
  { id: 'predict', label: '🔮 Predict' },
  { id: 'models', label: '💾 Saved Models' },
]

export function QSARModeling() {
  const [activeTab, setActiveTab] = useState('dataset')
  const [descriptorGroups, setDescriptorGroups] = useState<string[]>([])
  const [selectedGroups, setSelectedGroups] = useState<string[]>(['all'])
  const [dataset, setDataset] = useState<DatasetUploadResult | null>(null)
  const [csvFile, setCsvFile] = useState<File | null>(null)
  const [smilesCol, setSmilesCol] = useState('smiles')
  const [activityCol, setActivityCol] = useState('activity')
  const [uploadLoading, setUploadLoading] = useState(false)
  const [uploadError, setUploadError] = useState<string | null>(null)
  const [modelName, setModelName] = useState('My QSAR Model')
  const [modelType, setModelType] = useState('RandomForest')
  const [cvFolds, setCvFolds] = useState(5)
  const [trainJobId, setTrainJobId] = useState<string | null>(null)
  const [trainResult, setTrainResult] = useState<TrainResult | null>(null)
  const [trainError, setTrainError] = useState<string | null>(null)
  const [trainStatus, setTrainStatus] = useState<string>('idle')
  const [savedModels, setSavedModels] = useState<SavedModel[]>([])
  const [modelsLoading, setModelsLoading] = useState(false)
  const [predictSmiles, setPredictSmiles] = useState('')
  const [predictBulk, setPredictBulk] = useState('')
  const [predictMode, setPredictMode] = useState<'single' | 'bulk'>('single')
  const [selectedModelId, setSelectedModelId] = useState<string>('')
  const [prediction, setPrediction] = useState<PredictionSingle | null>(null)
  const [bulkPredictions, setBulkPredictions] = useState<any[] | null>(null)
  const [predictLoading, setPredictLoading] = useState(false)
  const [predictError, setPredictError] = useState<string | null>(null)
  const pollingRef = useRef<ReturnType<typeof setInterval> | null>(null)

  useEffect(() => {
    loadDescriptorGroups()
    loadModels()
  }, [])

  useEffect(() => {
    if (trainJobId && trainStatus === 'running') {
      pollingRef.current = setInterval(checkTrainingStatus, 2000)
    }
    return () => {
      if (pollingRef.current) clearInterval(pollingRef.current)
    }
  }, [trainJobId, trainStatus])

  async function loadDescriptorGroups() {
    try {
      const data = await getDescriptorGroups()
      setDescriptorGroups(data.groups)
    } catch (e) {
      console.warn('Could not load descriptor groups')
    }
  }

  async function loadModels() {
    setModelsLoading(true)
    try {
      const data = await listModels()
      setSavedModels(data.models)
      if (data.models.length > 0 && !selectedModelId) {
        setSelectedModelId(data.models[0].model_id)
      }
    } catch (e) {
      console.warn('Could not load models')
    } finally {
      setModelsLoading(false)
    }
  }

  async function handleUpload() {
    if (!csvFile) {
      const blob = new Blob([SAMPLE_DATASET_CSV], { type: 'text/csv' })
      const file = new File([blob], 'sample_qsar.csv', { type: 'text/csv' })
      setCsvFile(file)
      await uploadCsvFile(file)
      return
    }
    await uploadCsvFile(csvFile)
  }

  async function uploadCsvFile(file: File) {
    setUploadLoading(true)
    setUploadError(null)
    try {
      const groupsStr = selectedGroups.join(',')
      const result = await uploadDataset(file, smilesCol, activityCol, groupsStr)
      setDataset(result)
    } catch (err: any) {
      const msg = err?.response?.data?.detail || err?.message || 'Upload failed'
      setUploadError(msg)
    } finally {
      setUploadLoading(false)
    }
  }

  async function handleTrain() {
    if (!dataset) return
    setTrainError(null)
    setTrainResult(null)
    setTrainStatus('pending')
    try {
      const response = await startTraining(
        dataset.X,
        dataset.y,
        dataset.feature_names,
        modelType,
        modelName,
        activityCol,
        selectedGroups,
        cvFolds
      )
      setTrainJobId(response.job_id)
      setTrainStatus('running')
    } catch (err: any) {
      const msg = err?.response?.data?.detail || err?.message || 'Training failed'
      setTrainError(msg)
      setTrainStatus('idle')
    }
  }

  async function checkTrainingStatus() {
    if (!trainJobId) return
    try {
      const status = await getTrainingStatus(trainJobId)
      setTrainStatus(status.status)
      if (status.status === 'completed' && status.result) {
        setTrainResult(status.result)
        if (pollingRef.current) clearInterval(pollingRef.current)
        loadModels()
      } else if (status.status === 'failed') {
        setTrainError(status.error || 'Training failed')
        if (pollingRef.current) clearInterval(pollingRef.current)
      }
    } catch (e) {
      console.warn('Status poll failed')
    }
  }

  async function handlePredict() {
    if (!selectedModelId) return
    setPredictLoading(true)
    setPredictError(null)
    setPrediction(null)
    setBulkPredictions(null)
    try {
      if (predictMode === 'single' && predictSmiles.trim()) {
        const result = await predictSingle(selectedModelId, predictSmiles.trim())
        setPrediction(result)
      } else if (predictMode === 'bulk' && predictBulk.trim()) {
        const lines = predictBulk.split('\n').filter(l => l.trim())
        const result = await predictBatch(selectedModelId, lines)
        setBulkPredictions(result.predictions)
      }
    } catch (err: any) {
      const msg = err?.response?.data?.detail || err?.message || 'Prediction failed'
      setPredictError(msg)
    } finally {
      setPredictLoading(false)
    }
  }

  async function handleDeleteModel(modelId: string) {
    try {
      await deleteModel(modelId)
      loadModels()
    } catch (e) {
      console.warn('Delete failed')
    }
  }

  function downloadPredictions() {
    if (!bulkPredictions) return
    const csv = ['smiles,predicted_activity,ad_status'] +
      bulkPredictions.map((p: any) => `${p.smiles},${p.predicted_activity ?? ''},${p.ad_status}`).join('\n')
    const blob = new Blob([csv], { type: 'text/csv' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = 'qsar_predictions.csv'
    a.click()
    URL.revokeObjectURL(url)
  }

  function getAdBadge(status: string) {
    switch (status) {
      case 'in_domain': return <span className="px-2 py-0.5 text-xs rounded bg-green-900 text-green-300">In Domain</span>
      case 'warning': return <span className="px-2 py-0.5 text-xs rounded bg-yellow-900 text-yellow-300">Warning</span>
      case 'out_of_domain': return <span className="px-2 py-0.5 text-xs rounded bg-red-900 text-red-300">Out of Domain</span>
      default: return <span className="px-2 py-0.5 text-xs rounded bg-gray-700 text-gray-300">Unknown</span>
    }
  }

  function renderPlotly(data: any, layout: any) {
    try {
      return <Plot data={data} layout={{ ...layout, paper_bgcolor: 'transparent', plot_bgcolor: 'transparent' }} config={{ displayModeBar: false }} style={{ width: '100%', height: 300 }} />
    } catch (e) {
      return <div className="h-[300px] bg-gray-800 rounded flex items-center justify-center text-gray-500">Chart unavailable</div>
    }
  }

  return (
    <div className="p-6">
      <div className="mb-6">
        <h1 className="text-2xl font-bold text-text-primary">QSAR Modeling</h1>
        <p className="text-text-secondary mt-1">Train ML models to predict molecular activity</p>
      </div>

      <Tabs tabs={TABS} activeTab={activeTab} onChange={setActiveTab} />

      <TabPanel>
        {activeTab === 'dataset' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mt-4">
            <Card>
              <h3 className="font-bold text-text-primary mb-4">Upload Dataset</h3>
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">SMILES Column</label>
                  <input
                    type="text"
                    value={smilesCol}
                    onChange={(e) => setSmilesCol(e.target.value)}
                    className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm"
                  />
                </div>
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">Activity Column</label>
                  <input
                    type="text"
                    value={activityCol}
                    onChange={(e) => setActivityCol(e.target.value)}
                    className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm"
                  />
                </div>
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">Descriptor Groups</label>
                  <div className="flex flex-wrap gap-2">
                    {descriptorGroups.map((g: string) => (
                      <label key={g} className="flex items-center gap-1 text-sm text-gray-300 cursor-pointer">
                        <input
                          type="checkbox"
                          checked={selectedGroups.includes(g)}
                          onChange={(e) => {
                            if (e.target.checked) {
                              setSelectedGroups(prev => [...prev, g])
                            } else {
                              setSelectedGroups(prev => prev.filter(x => x !== g))
                            }
                          }}
                          className="rounded"
                        />
                        {g}
                      </label>
                    ))}
                  </div>
                </div>
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">CSV File</label>
                  <input
                    type="file"
                    accept=".csv"
                    onChange={(e) => setCsvFile(e.target.files?.[0] || null)}
                    className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm"
                  />
                </div>
                <div className="flex gap-3">
                  <Button onClick={handleUpload} disabled={uploadLoading}>
                    {uploadLoading ? 'Uploading...' : dataset ? 'Re-upload' : 'Upload Dataset'}
                  </Button>
                  {!csvFile && (
                    <Button variant="secondary" onClick={handleUpload} disabled={uploadLoading}>
                      Load Sample Data
                    </Button>
                  )}
                </div>
                {uploadError && (
                  <div className="p-3 rounded bg-red-900/50 border border-red-700 text-red-300 text-sm">
                    {uploadError}
                  </div>
                )}
              </div>
            </Card>

            {dataset ? (
              <Card>
                <h3 className="font-bold text-text-primary mb-4">Dataset Summary</h3>
                <div className="space-y-3">
                  <div className="grid grid-cols-2 gap-4">
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">Compounds</p>
                      <p className="text-xl font-bold text-white">{dataset.n_compounds}</p>
                    </div>
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">Features</p>
                      <p className="text-xl font-bold text-white">{dataset.n_features}</p>
                    </div>
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">Activity Mean</p>
                      <p className="text-xl font-bold text-white">{dataset.activity_mean.toFixed(3)}</p>
                    </div>
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">Activity Std</p>
                      <p className="text-xl font-bold text-white">{dataset.activity_std.toFixed(3)}</p>
                    </div>
                  </div>
                  <div className="grid grid-cols-2 gap-4">
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">Min Activity</p>
                      <p className="text-lg font-bold text-white">{dataset.activity_min.toFixed(3)}</p>
                    </div>
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">Max Activity</p>
                      <p className="text-lg font-bold text-white">{dataset.activity_max.toFixed(3)}</p>
                    </div>
                  </div>
                  {dataset.failed_count > 0 && (
                    <div className="p-3 rounded bg-yellow-900/30 border border-yellow-700 text-sm text-yellow-300">
                      {dataset.failed_count} molecules failed to parse
                    </div>
                  )}
                  <div>
                    <p className="text-xs text-gray-400 mb-1">Feature Names</p>
                    <div className="flex flex-wrap gap-1 max-h-24 overflow-y-auto">
                      {dataset.feature_names.slice(0, 20).map((f: string) => (
                        <span key={f} className="px-2 py-0.5 text-xs bg-gray-700 rounded text-gray-300">{f}</span>
                      ))}
                      {dataset.feature_names.length > 20 && (
                        <span className="text-xs text-gray-500">+{dataset.feature_names.length - 20} more</span>
                      )}
                    </div>
                  </div>
                </div>
              </Card>
            ) : (
              <Card>
                <h3 className="font-bold text-text-primary mb-4">Dataset Preview</h3>
                <p className="text-gray-500 text-sm">Upload a CSV with SMILES and activity columns to begin</p>
                <div className="mt-4 p-4 bg-gray-800 rounded text-xs text-gray-400 font-mono">
                  <p>smiles,activity</p>
                  <p>CC(=O)OC1=CC=CC=C1C(=O)O,4.8239</p>
                  <p>Cn1cnc2c1c(=O)n(c(=O)n2C)C,5.3979</p>
                  <p>CC(C)Cc1ccc(cc1)C(C)C(=O)O,3.9707</p>
                  <p>...</p>
                </div>
              </Card>
            )}
          </div>
        )}

        {activeTab === 'train' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mt-4">
            <Card>
              <h3 className="font-bold text-text-primary mb-4">Training Configuration</h3>
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">Model Name</label>
                  <input
                    type="text"
                    value={modelName}
                    onChange={(e) => setModelName(e.target.value)}
                    className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm"
                  />
                </div>
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">Algorithm</label>
                  <select
                    value={modelType}
                    onChange={(e) => setModelType(e.target.value)}
                    className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm"
                  >
                    {MODEL_TYPES.map(m => <option key={m} value={m}>{m}</option>)}
                  </select>
                </div>
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">
                    Cross-Validation Folds: {cvFolds}
                  </label>
                  <input
                    type="range"
                    min={3}
                    max={10}
                    value={cvFolds}
                    onChange={(e) => setCvFolds(Number(e.target.value))}
                    className="w-full"
                  />
                </div>
                <Button
                  onClick={handleTrain}
                  disabled={!dataset || trainStatus === 'running'}
                  className="w-full"
                >
                  {trainStatus === 'running' ? 'Training...' : 'Train Model'}
                </Button>
                {trainError && (
                  <div className="p-3 rounded bg-red-900/50 border border-red-700 text-red-300 text-sm">
                    {trainError}
                  </div>
                )}
              </div>
            </Card>

            <Card>
              <h3 className="font-bold text-text-primary mb-4">Training Results</h3>
              {trainResult ? (
                <div className="space-y-4">
                  <div className="grid grid-cols-2 gap-3">
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">CV R²</p>
                      <p className="text-2xl font-bold text-green-400">{trainResult.metrics.cv_r2.toFixed(4)}</p>
                    </div>
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">CV RMSE</p>
                      <p className="text-2xl font-bold text-blue-400">{trainResult.metrics.cv_rmse.toFixed(4)}</p>
                    </div>
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">CV MAE</p>
                      <p className="text-2xl font-bold text-purple-400">{trainResult.metrics.cv_mae.toFixed(4)}</p>
                    </div>
                    <div className="bg-gray-800 rounded p-3">
                      <p className="text-xs text-gray-400">Train R²</p>
                      <p className="text-2xl font-bold text-white">{trainResult.metrics.train_r2.toFixed(4)}</p>
                    </div>
                  </div>
                  {trainResult.scatter_plot && renderPlotly(trainResult.scatter_plot.data, trainResult.scatter_plot.layout)}
                  <div className="flex justify-between items-center text-xs text-gray-400">
                    <span>Model ID: {trainResult.model_id}</span>
                    <span>Saved to Models tab</span>
                  </div>
                </div>
              ) : trainStatus === 'running' ? (
                <div className="flex items-center justify-center h-40">
                  <div className="text-center">
                    <div className="animate-spin text-3xl mb-2">⏳</div>
                    <p className="text-gray-400 text-sm">Training in progress...</p>
                    <p className="text-gray-500 text-xs mt-1">This may take a few minutes</p>
                  </div>
                </div>
              ) : (
                <div className="flex items-center justify-center h-40">
                  <p className="text-gray-500 text-sm">Upload a dataset and configure training to see results</p>
                </div>
              )}
            </Card>
          </div>
        )}

        {activeTab === 'predict' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mt-4">
            <Card>
              <h3 className="font-bold text-text-primary mb-4">Predict Activity</h3>
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-text-secondary mb-1">Select Model</label>
                  <select
                    value={selectedModelId}
                    onChange={(e) => setSelectedModelId(e.target.value)}
                    className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm"
                  >
                    {savedModels.length === 0 && <option value="">No models saved</option>}
                    {savedModels.map((m: SavedModel) => (
                      <option key={m.model_id} value={m.model_id}>
                        {m.name} ({m.model_type}) - R²={m.metrics.cv_r2?.toFixed(3) ?? 'N/A'}
                      </option>
                    ))}
                  </select>
                </div>
                <div className="flex gap-2 border-b pb-2">
                  <button
                    onClick={() => setPredictMode('single')}
                    className={`px-3 py-1 text-sm rounded ${predictMode === 'single' ? 'bg-primary text-white' : 'text-gray-400'}`}
                  >
                    Single SMILES
                  </button>
                  <button
                    onClick={() => setPredictMode('bulk')}
                    className={`px-3 py-1 text-sm rounded ${predictMode === 'bulk' ? 'bg-primary text-white' : 'text-gray-400'}`}
                  >
                    Bulk
                  </button>
                </div>

                {predictMode === 'single' ? (
                  <div>
                    <div className="flex gap-2 mb-2">
                      <input
                        type="text"
                        value={predictSmiles}
                        onChange={(e) => setPredictSmiles(e.target.value)}
                        placeholder="Enter SMILES, e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
                        className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm"
                      />
                    </div>
                  </div>
                ) : (
                  <div>
                    <label className="block text-sm font-medium text-text-secondary mb-1">
                      One SMILES per line
                    </label>
                    <textarea
                      value={predictBulk}
                      onChange={(e) => setPredictBulk(e.target.value)}
                      placeholder="CC(=O)OC1=CC=CC=C1C(=O)O&#10;Cn1cnc2c1c(=O)n(c(=O)n2C)C&#10;CC(C)Cc1ccc(cc1)C(C)C(=O)O"
                      rows={8}
                      className="w-full bg-gray-800 border border-gray-600 rounded px-3 py-2 text-white text-sm font-mono"
                    />
                  </div>
                )}

                <Button
                  onClick={handlePredict}
                  disabled={predictLoading || !selectedModelId || (predictMode === 'single' && !predictSmiles.trim())}
                  className="w-full"
                >
                  {predictLoading ? 'Predicting...' : 'Predict'}
                </Button>
                {predictError && (
                  <div className="p-3 rounded bg-red-900/50 border border-red-700 text-red-300 text-sm">
                    {predictError}
                  </div>
                )}
              </div>
            </Card>

            <Card>
              <h3 className="font-bold text-text-primary mb-4">Prediction Results</h3>
              {prediction ? (
                <div className="space-y-3">
                  <div className="bg-gray-800 rounded p-4">
                    <p className="text-xs text-gray-400">Input SMILES</p>
                    <p className="text-sm text-white font-mono break-all">{prediction.smiles}</p>
                  </div>
                  <div className="bg-gray-800 rounded p-4">
                    <p className="text-xs text-gray-400">Predicted Activity</p>
                    <p className="text-3xl font-bold text-green-400">{prediction.predicted_activity.toFixed(4)}</p>
                  </div>
                  <div className="flex items-center gap-3">
                    <p className="text-sm text-gray-400">Applicability Domain:</p>
                    {getAdBadge(prediction.ad_status)}
                  </div>
                  {prediction.ad_leverage !== undefined && (
                    <div className="text-xs text-gray-500">
                      Leverage: {prediction.ad_leverage.toFixed(4)}
                      {prediction.ad_warning_threshold && ` (warn: ${prediction.ad_warning_threshold.toFixed(4)})`}
                    </div>
                  )}
                </div>
              ) : bulkPredictions ? (
                <div className="space-y-3">
                  <div className="flex justify-between items-center">
                    <p className="text-sm text-gray-400">
                      {bulkPredictions.filter((p: any) => p.predicted_activity !== null).length} / {bulkPredictions.length} predicted
                    </p>
                    <Button variant="secondary" onClick={downloadPredictions} className="text-xs">
                      Download CSV
                    </Button>
                  </div>
                  <div className="max-h-80 overflow-y-auto space-y-1">
                    {bulkPredictions.map((p: any, i: number) => (
                      <div key={i} className="flex justify-between items-center bg-gray-800 rounded px-3 py-2 text-sm">
                        <span className="text-gray-300 font-mono truncate flex-1 mr-3">{p.smiles}</span>
                        <span className="text-white font-mono mr-3">
                          {p.predicted_activity !== null ? p.predicted_activity.toFixed(4) : '—'}
                        </span>
                        {getAdBadge(p.ad_status)}
                      </div>
                    ))}
                  </div>
                </div>
              ) : (
                <div className="flex items-center justify-center h-40">
                  <p className="text-gray-500 text-sm">Select a model and enter SMILES to predict</p>
                </div>
              )}
            </Card>
          </div>
        )}

        {activeTab === 'models' && (
          <div className="mt-4">
            <Card>
              <h3 className="font-bold text-text-primary mb-4">Saved Models</h3>
              {modelsLoading ? (
                <div className="flex justify-center py-8"><span className="animate-spin text-2xl">⏳</span></div>
              ) : savedModels.length === 0 ? (
                <div className="text-center py-8 text-gray-500">
                  <p className="mb-2">No models saved yet</p>
                  <p className="text-sm">Train a model and save it to see it here</p>
                </div>
              ) : (
                <div className="overflow-x-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="border-b border-gray-700 text-left text-gray-400">
                        <th className="pb-2 font-medium">Name</th>
                        <th className="pb-2 font-medium">Algorithm</th>
                        <th className="pb-2 font-medium">CV R²</th>
                        <th className="pb-2 font-medium">CV RMSE</th>
                        <th className="pb-2 font-medium">Dataset Size</th>
                        <th className="pb-2 font-medium">Created</th>
                        <th className="pb-2 font-medium">Actions</th>
                      </tr>
                    </thead>
                    <tbody>
                      {savedModels.map((m: SavedModel) => (
                        <tr key={m.model_id} className="border-b border-gray-800">
                          <td className="py-3 text-white font-medium">{m.name}</td>
                          <td className="py-3 text-gray-300">{m.model_type}</td>
                          <td className="py-3 text-green-400">{m.metrics.cv_r2?.toFixed(4) ?? '—'}</td>
                          <td className="py-3 text-blue-400">{m.metrics.cv_rmse?.toFixed(4) ?? '—'}</td>
                          <td className="py-3 text-gray-300">{m.n_features} features</td>
                          <td className="py-3 text-gray-500">{new Date(m.created_at).toLocaleDateString()}</td>
                          <td className="py-3">
                            <button
                              onClick={() => handleDeleteModel(m.model_id)}
                              className="text-red-400 hover:text-red-300 text-xs"
                            >
                              Delete
                            </button>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </Card>
          </div>
        )}
      </TabPanel>

    </div>
  )
}
