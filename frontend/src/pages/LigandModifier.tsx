import { useState } from 'react'
import { ligandModifierAPI, type LigandModifierRequest, type ModificationResult } from '@/api/ligandModifier'

const MODES = [
  { key: 'similarity_search', label: 'Similarity Search', icon: '🔍', desc: 'Find similar compounds in PubChem' },
  { key: 'prompt_based', label: 'Prompt-Based', icon: '✏️', desc: 'Describe modifications in natural language' },
  { key: 'autonomous', label: 'Autonomous', icon: '🤖', desc: 'AI-driven iterative optimization' },
] as const

const SAMPLE_PROMPTS = [
  'Add polar groups but keep MW < 450',
  'Increase lipophilicity',
  'Reduce molecular weight',
  'Add H-bond donors for better binding',
  'Improve drug-likeness (Lipinski rules)',
]

export function LigandModifier() {
  const [parentSmiles, setParentSmiles] = useState('')
  const [receptorPdb, setReceptorPdb] = useState('')
  const [mode, setMode] = useState<'similarity_search' | 'prompt_based' | 'autonomous'>('prompt_based')
  const [prompt, setPrompt] = useState('')
  const [similarityThreshold, setSimilarityThreshold] = useState(0.85)
  const [maxVariants, setMaxVariants] = useState(50)

  const [jobStatus, setJobStatus] = useState<string | null>(null)
  const [progress, setProgress] = useState(0)
  const [results, setResults] = useState<ModificationResult[]>([])
  const [error, setError] = useState<string | null>(null)
  const [isRunning, setIsRunning] = useState(false)
  const [selectedSmiles, setSelectedSmiles] = useState<string | null>(null)

  async function handleStart() {
    if (!parentSmiles.trim()) {
      setError('Please enter a parent ligand SMILES')
      return
    }

    try {
      setIsRunning(true)
      setError(null)
      setResults([])
      setProgress(0)
      setJobStatus(null)

      const req: LigandModifierRequest = {
        parent_smiles: parentSmiles.trim(),
        receptor_pdb: receptorPdb.trim() || 'DUMMY',
        mode,
        prompt: mode !== 'similarity_search' ? prompt : undefined,
        database: 'pubchem',
        similarity_threshold: similarityThreshold,
        max_variants: maxVariants,
        docking_exhaustiveness: 8,
      }

      const { job_id } = await ligandModifierAPI.optimize(req)

      const finalStatus = await ligandModifierAPI.pollUntilComplete(
        job_id,
        (status) => {
          setJobStatus(status.status)
          setProgress(status.progress)
        }
      )

      if (finalStatus.status === 'completed') {
        setResults(finalStatus.results)
      } else if (finalStatus.status === 'failed') {
        setError(finalStatus.error || 'Optimization failed')
      }
    } catch (e: any) {
      setError(e.message || 'Failed to start optimization')
    } finally {
      setIsRunning(false)
      setJobStatus(null)
      setProgress(0)
    }
  }

  function handleCancel() {
    setIsRunning(false)
    setJobStatus('cancelled')
  }

  function copySmiles(smi: string) {
    navigator.clipboard.writeText(smi)
  }

  const statusColors: Record<string, string> = {
    queued: 'bg-gray-500',
    parsing: 'bg-yellow-500',
    fetching: 'bg-blue-500',
    transforming: 'bg-purple-500',
    validating: 'bg-orange-500',
    docking: 'bg-cyan-500',
    completed: 'bg-green-500',
    failed: 'bg-red-500',
    cancelled: 'bg-gray-400',
  }

  const sourceColors: Record<string, string> = {
    database: 'bg-blue-100 text-blue-800',
    transformation: 'bg-purple-100 text-purple-800',
    autonomous: 'bg-green-100 text-green-800',
  }

  return (
    <div className="p-6 max-w-7xl mx-auto space-y-6">
      {/* Header */}
      <div className="mb-2">
        <h1 className="text-2xl font-bold text-gray-900">Ligand Modifier</h1>
        <p className="text-gray-500 mt-1">Find similar compounds and apply intelligent modifications</p>
      </div>

      {/* Input Section */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4 bg-white rounded-lg shadow-sm p-4">
        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">Parent Ligand (SMILES)</label>
          <textarea
            value={parentSmiles}
            onChange={e => setParentSmiles(e.target.value)}
            placeholder="CCO or paste SMILES..."
            className="w-full p-2 border border-gray-300 rounded-lg text-sm font-mono h-20 resize-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
          />
        </div>
        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">Receptor PDB (optional)</label>
          <textarea
            value={receptorPdb}
            onChange={e => setReceptorPdb(e.target.value)}
            placeholder="Paste PDB content for docking-based ranking..."
            className="w-full p-2 border border-gray-300 rounded-lg text-sm font-mono h-20 resize-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
          />
        </div>
      </div>

      {/* Mode Selection */}
      <div className="bg-white rounded-lg shadow-sm p-4 space-y-4">
        <h3 className="font-semibold text-gray-900">Modification Mode</h3>
        <div className="grid grid-cols-3 gap-3">
          {MODES.map(m => (
            <button
              key={m.key}
              onClick={() => setMode(m.key as typeof mode)}
              className={`p-3 rounded-lg border-2 text-left transition-all ${
                mode === m.key
                  ? 'border-blue-500 bg-blue-50'
                  : 'border-gray-200 hover:border-gray-300'
              }`}
            >
              <div className="text-lg mb-1">{m.icon}</div>
              <div className="font-medium text-sm text-gray-900">{m.label}</div>
              <div className="text-xs text-gray-500 mt-0.5">{m.desc}</div>
            </button>
          ))}
        </div>

        {mode === 'prompt_based' && (
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">Modification Prompt</label>
            <textarea
              value={prompt}
              onChange={e => setPrompt(e.target.value)}
              placeholder="e.g., 'Add polar groups but keep MW < 450'"
              className="w-full p-2 border border-gray-300 rounded-lg text-sm h-16 resize-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
            />
            <div className="flex flex-wrap gap-1.5 mt-2">
              {SAMPLE_PROMPTS.map(p => (
                <button
                  key={p}
                  onClick={() => setPrompt(p)}
                  className="text-xs px-2 py-1 bg-gray-100 hover:bg-gray-200 rounded-full text-gray-600 transition-colors"
                >
                  {p}
                </button>
              ))}
            </div>
          </div>
        )}

        {mode === 'similarity_search' && (
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Similarity Threshold: {(similarityThreshold * 100).toFixed(0)}%
            </label>
            <input
              type="range"
              min="0.8"
              max="0.95"
              step="0.05"
              value={similarityThreshold}
              onChange={e => setSimilarityThreshold(parseFloat(e.target.value))}
              className="w-full"
            />
            <div className="flex justify-between text-xs text-gray-400">
              <span>80% (more results)</span>
              <span>95% (stricter match)</span>
            </div>
          </div>
        )}

        <div className="flex items-center gap-4">
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">Max Variants</label>
            <input
              type="number"
              min="10"
              max="100"
              value={maxVariants}
              onChange={e => setMaxVariants(parseInt(e.target.value) || 50)}
              className="w-20 p-2 border border-gray-300 rounded-lg text-sm"
            />
          </div>
        </div>
      </div>

      {/* Action Buttons */}
      <div className="flex gap-3">
        <button
          onClick={handleStart}
          disabled={isRunning || !parentSmiles.trim()}
          className="px-6 py-2.5 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-400 disabled:cursor-not-allowed text-white rounded-lg font-medium transition-colors"
        >
          {isRunning ? '⏳ Running...' : '▶ Start Optimization'}
        </button>
        {isRunning && (
          <button
            onClick={handleCancel}
            className="px-4 py-2.5 bg-red-600 hover:bg-red-700 text-white rounded-lg transition-colors"
          >
            ⏹ Cancel
          </button>
        )}
      </div>

      {/* Progress Bar */}
      {jobStatus && (
        <div className="bg-white rounded-lg shadow-sm p-4">
          <div className="flex items-center justify-between mb-2">
            <span className="text-sm font-medium text-gray-700 capitalize">{jobStatus.replace('_', ' ')}</span>
            <span className="text-sm text-gray-500">{(progress * 100).toFixed(0)}%</span>
          </div>
          <div className="w-full bg-gray-200 rounded-full h-2">
            <div
              className={`h-2 rounded-full transition-all duration-500 ${statusColors[jobStatus] || 'bg-gray-500'}`}
              style={{ width: `${progress * 100}%` }}
            />
          </div>
        </div>
      )}

      {/* Error Display */}
      {error && (
        <div className="p-4 bg-red-50 border border-red-200 rounded-lg text-red-700 text-sm">
          ⚠️ {error}
        </div>
      )}

      {/* Results Table */}
      {results.length > 0 && (
        <div className="bg-white rounded-lg shadow-sm overflow-hidden">
          <div className="px-4 py-3 border-b border-gray-200">
            <h3 className="font-semibold text-gray-900">
              Results ({results.length} variants)
            </h3>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="bg-gray-50 text-left text-xs text-gray-500 uppercase">
                  <th className="px-4 py-2">#</th>
                  <th className="px-4 py-2">Modified SMILES</th>
                  <th className="px-4 py-2">Transform</th>
                  <th className="px-4 py-2">MW</th>
                  <th className="px-4 py-2">LogP</th>
                  <th className="px-4 py-2">HBD</th>
                  <th className="px-4 py-2">HBA</th>
                  <th className="px-4 py-2">Source</th>
                  <th className="px-4 py-2">Actions</th>
                </tr>
              </thead>
              <tbody>
                {results.map((r, i) => (
                  <tr
                    key={i}
                    onClick={() => setSelectedSmiles(selectedSmiles === r.modified_smiles ? null : r.modified_smiles)}
                    className={`border-t border-gray-100 cursor-pointer transition-colors ${
                      selectedSmiles === r.modified_smiles ? 'bg-blue-50' : 'hover:bg-gray-50'
                    }`}
                  >
                    <td className="px-4 py-2 text-gray-500">{i + 1}</td>
                    <td className="px-4 py-2 font-mono text-xs max-w-xs truncate" title={r.modified_smiles}>
                      {r.modified_smiles}
                    </td>
                    <td className="px-4 py-2 text-gray-600">{r.applied_transform}</td>
                    <td className="px-4 py-2">{r.properties?.mw?.toFixed(1) ?? '-'}</td>
                    <td className="px-4 py-2">{r.properties?.logp?.toFixed(2) ?? '-'}</td>
                    <td className="px-4 py-2">{r.properties?.hbd ?? '-'}</td>
                    <td className="px-4 py-2">{r.properties?.hba ?? '-'}</td>
                    <td className="px-4 py-2">
                      <span className={`px-2 py-0.5 rounded-full text-xs font-medium ${sourceColors[r.source] || 'bg-gray-100 text-gray-700'}`}>
                        {r.source}
                      </span>
                    </td>
                    <td className="px-4 py-2">
                      <button
                        onClick={e => { e.stopPropagation(); copySmiles(r.modified_smiles) }}
                        className="text-blue-600 hover:text-blue-800 text-xs font-medium"
                      >
                        Copy
                      </button>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Empty State */}
      {!isRunning && results.length === 0 && !error && (
        <div className="text-center py-12 text-gray-400">
          <div className="text-4xl mb-2">🧬</div>
          <p className="text-sm">Enter a SMILES string and start optimization to see results</p>
        </div>
      )}
    </div>
  )
}
