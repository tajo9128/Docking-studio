import { useState, useEffect, useRef } from 'react'
import { Button } from '@/components/ui'

const FDA_DRUGS = [
  { name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O', use: 'Anti-inflammatory' },
  { name: 'Caffeine', smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C', use: 'Stimulant' },
  { name: 'Glucose', smiles: 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', use: 'Energy metabolism' },
  { name: 'Ibuprofen', smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', use: 'Analgesic' },
  { name: 'Morphine', smiles: 'CN1CCc2c(O)ccc(c2C1)C(O)=O', use: 'Analgesic' },
  { name: 'Benzene', smiles: 'c1ccccc1', use: 'Aromatic scaffold' },
  { name: 'Acetaminophen', smiles: 'CC(=O)Nc1ccc(cc1)O', use: 'Analgesic' },
  { name: 'Lisinopril', smiles: 'c1ccc2c(c1)CCCC2C(=O)NCCCC(N)C(=O)O', use: 'ACE inhibitor' },
  { name: 'Metformin', smiles: 'CN(C)N=C(N)N', use: 'Diabetes' },
  { name: 'Warfarin', smiles: 'CC(=O)OC(Cc1c(O)c2ccccc2oc1=O)C(c1ccccc1)=O', use: 'Anticoagulant' },
  { name: 'Tamoxifen', smiles: 'CC(C)=c1ccccc1C(OCCCN(C)C)c1ccccc1', use: 'Anticancer' },
  { name: 'Sildenafil', smiles: 'CCCC1=C2N(C(=O)N1CCC)CCCC2c3ccc(cc3)S(=O)(=O)N', use: 'PDE5 inhibitor' },
]

const MUTATION_STRATEGIES = [
  { key: 'add_halogen', label: 'Add Halogen (F, Cl, Br)', icon: '+' },
  { key: 'add_oh', label: 'Add OH Group', icon: '+' },
  { key: 'add_nh2', label: 'Add NH2 Group', icon: '+' },
  { key: 'add_aromatic', label: 'Add Aromatic Ring', icon: '+' },
  { key: 'bioisostere', label: 'Bioisosteric Replace', icon: '~' },
  { key: 'reduce_flex', label: 'Reduce Flexibility', icon: '-' },
]

interface Properties {
  mw: number | null
  logp: number | null
  hbd: number | null
  hba: number | null
  tpsa: number | null
  rotatable: number | null
  formula: string | null
  valid: boolean
}

interface Suggestion {
  text: string
  type: 'good' | 'warning' | 'error' | 'info'
}

export function ChemDraw() {
  const [smiles, setSmiles] = useState('CC(=O)Oc1ccccc1C(=O)O')
  const [molName] = useState('Aspirin')
  const [properties, setProperties] = useState<Properties>({
    mw: 180.16, logp: 1.19, hbd: 1, hba: 3, tpsa: 63.60,
    rotatable: 4, formula: 'C9H8O4', valid: true
  })
  const [suggestions, setSuggestions] = useState<Suggestion[]>([])
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState<'2d' | '3d'>('2d')
  const canvas2dRef = useRef<HTMLCanvasElement>(null)
  const viewer3dRef = useRef<HTMLDivElement>(null)

  useEffect(() => {
    if (activeTab === '2d' && smiles) {
      draw2D()
    }
  }, [smiles, activeTab])

  useEffect(() => {
    if (activeTab === '3d') {
      init3DViewer()
    }
  }, [activeTab, smiles])

  const draw2D = async () => {
    const canvas = canvas2dRef.current
    if (!canvas || !smiles) return

    try {
      if (typeof (window as any).SmilesDrawer === 'undefined') {
        await loadSmilesDrawer()
      }
      const drawer = new (window as any).SmilesDrawer.Drawer({ width: 400, height: 300, bondThickness: 2 })
      ;(window as any).SmilesDrawer.parse(smiles, (tree: any) => {
        drawer.draw(tree, canvas, 'light', false)
      })
    } catch (e) {
      console.error('2D draw error:', e)
    }
  }

  const loadSmilesDrawer = () => {
    return new Promise<void>((resolve) => {
      if (document.getElementById('smiles-drawer-script')) {
        resolve()
        return
      }
      const css = document.createElement('link')
      css.rel = 'stylesheet'
      css.href = 'https://unpkg.com/smiles-drawer@2.2.1/dist/smiles-drawer.min.css'
      document.head.appendChild(css)

      const script = document.createElement('script')
      script.id = 'smiles-drawer-script'
      script.src = 'https://unpkg.com/smiles-drawer@2.2.1/dist/smiles-drawer.min.js'
      script.onload = () => resolve()
      document.body.appendChild(script)
    })
  }

  const init3DViewer = async () => {
    const container = viewer3dRef.current
    if (!container) return

    try {
      if (typeof (window as any).NGL === 'undefined') {
        await loadNGL()
      }
      const stage = new (window as any).NGL.Stage(container)
      stage.setParameters({ backgroundColor: '0x1a1a2e' })

      const res = await fetch(`/api/chem/dock`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles, receptor_id: 'default' })
      })
      const data = await res.json()

      if (data.job_id) {
        const pdbRes = await fetch(`/api/chem/3d/${data.job_id}`)
        const pdbData = await pdbRes.json()
        if (pdbData.pdb) {
          stage.removeAllComponents()
          stage.loadFile(data.job_id + '.pdb', { asTrajectory: false }).then((comp: any) => {
            comp.addRepresentation('ball-and-stick')
            stage.autoView()
          })
        }
      }
    } catch (e) {
      console.error('3D viewer error:', e)
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

  const analyzeMolecule = async () => {
    if (!smiles) return
    setLoading(true)
    try {
      const res = await fetch('/api/chem/properties', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const data = await res.json()
      if (data && !data.error) {
        setProperties({
          mw: data.mw,
          logp: data.logp,
          hbd: data.hbd,
          hba: data.hba,
          tpsa: data.tpsa,
          rotatable: data.rotatable_bonds,
          formula: data.formula,
          valid: data.valid !== false
        })
      }
    } catch (e) {
      console.error('Properties error:', e)
    }

    try {
      const res = await fetch('/api/chem/suggestions', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const data = await res.json()
      if (data.suggestions) {
        setSuggestions(data.suggestions.map((s: string) => {
          let type: Suggestion['type'] = 'info'
          if (s.includes('good') || s.includes('Good') || s.includes('drug-like') || s.includes('Ready')) type = 'good'
          if (s.includes('high') || s.includes('exceeds') || s.includes('low') || s.includes('poor') || s.includes('Reduce')) type = 'warning'
          return { text: s, type }
        }))
      }
    } catch (e) {
      console.error('Suggestions error:', e)
    }
    setLoading(false)
  }

  const mutateMolecule = (strategy: string) => {
    let newSmiles = smiles
    switch (strategy) {
      case 'add_halogen':
        newSmiles = smiles.replace(/([Cc])\)$/, '$1F)')
        if (newSmiles === smiles) newSmiles = smiles + 'F'
        break
      case 'add_oh':
        newSmiles = smiles.replace(/([Cc])\)$/, '$1O)')
        if (newSmiles === smiles) newSmiles = smiles + 'O'
        break
      case 'add_nh2':
        newSmiles = smiles.replace(/([Cc])\)$/, '$1N)')
        if (newSmiles === smiles) newSmiles = smiles + 'N'
        break
      case 'add_aromatic':
        newSmiles = smiles + 'c1ccccc1'
        break
      case 'bioisostere':
        newSmiles = smiles.replace(/CO/, 'CF').replace(/NH2/, 'OH')
        if (newSmiles === smiles) newSmiles = smiles.replace(/C=O/, 'C=S')
        break
      case 'reduce_flex':
        newSmiles = smiles.replace(/\([A-Za-z0-9]+\)/g, '')
        break
    }
    setSmiles(newSmiles)
  }

  const optimizeMolecule = async () => {
    const strategies = ['add_halogen', 'bioisostere', 'reduce_flex']
    const best = strategies[Math.floor(Math.random() * strategies.length)]
    mutateMolecule(best)
    setSuggestions(prev => [
      { text: `AI applied: ${MUTATION_STRATEGIES.find(s => s.key === best)?.label}`, type: 'good' },
      ...prev.slice(0, 4)
    ])
  }

  const copySmiles = () => {
    navigator.clipboard.writeText(smiles)
  }

  const clearEditor = () => {
    setSmiles('')
    setSuggestions([])
    setProperties({ mw: null, logp: null, hbd: null, hba: null, tpsa: null, rotatable: null, formula: null, valid: false })
  }

  const loadExample = (_name: string, smilesStr: string) => {
    setSmiles(smilesStr)
  }

  const dockMolecule = async () => {
    if (!smiles) return
    try {
      await fetch('/api/chem/dock', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles, receptor_id: 'default' })
      })
      window.location.href = '/docking'
    } catch (e) {
      console.error('Dock error:', e)
    }
  }

  const RuleCheck = ({ pass, label, value }: { pass: boolean; label: string; value?: string }) => (
    <div className={`flex justify-between items-center text-xs py-1 px-2 rounded ${pass ? 'text-green-600' : 'text-red-500'}`}>
      <span>{pass ? '✓' : '✗'} {label}</span>
      {value && <span className="font-mono">{value}</span>}
    </div>
  )

  return (
    <div className="h-full flex flex-col bg-gray-100">
      <div className="bg-slate-900 text-white px-6 py-3 flex items-center gap-4">
        <h1 className="text-lg font-bold text-cyan-400">Biodockify ChemDraw</h1>
        <span className="text-xs text-gray-400">Design, analyze & dock molecules</span>
        <div className="flex-1" />
        <Button variant="outline" size="sm" onClick={clearEditor}>Clear</Button>
        <Button variant="outline" size="sm" onClick={copySmiles}>Copy SMILES</Button>
        <Button variant="outline" size="sm" onClick={analyzeMolecule} disabled={loading}>
          {loading ? 'Analyzing...' : 'Analyze'}
        </Button>
        <Button variant="primary" size="sm" onClick={dockMolecule}>Dock</Button>
        <Button variant="secondary" size="sm" onClick={optimizeMolecule}>Optimize</Button>
      </div>

      <div className="flex-1 flex overflow-hidden">
        <div className="w-72 bg-white border-r border-gray-200 flex flex-col overflow-y-auto">
          <div className="p-3 border-b border-gray-200">
            <label className="text-xs font-medium text-gray-600 mb-1 block">SMILES Input</label>
            <textarea
              value={smiles}
              onChange={e => setSmiles(e.target.value)}
              className="w-full px-3 py-2 text-xs font-mono border border-gray-300 rounded-lg resize-none bg-gray-50 focus:ring-2 focus:ring-cyan-400 focus:border-cyan-400"
              rows={3}
              placeholder="Enter SMILES..."
            />
            <Button variant="primary" size="sm" className="w-full mt-2" onClick={analyzeMolecule} disabled={loading}>
              Analyze Structure
            </Button>
          </div>

          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700 mb-2">Molecule Library</div>
            <div className="grid grid-cols-2 gap-1">
              {FDA_DRUGS.map(drug => (
                <button
                  key={drug.name}
                  onClick={() => loadExample(drug.name, drug.smiles)}
                  className="text-xs px-2 py-1.5 bg-gray-50 hover:bg-cyan-50 hover:text-cyan-700 border border-gray-200 rounded text-left transition-colors"
                  title={drug.use}
                >
                  {drug.name}
                </button>
              ))}
            </div>
          </div>

          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700 mb-2">Properties (Lipinski)</div>
            <div className="space-y-0.5">
              <div className="flex justify-between text-xs py-1 px-2 bg-gray-50 rounded">
                <span className="text-gray-600">MW</span>
                <span className="font-mono font-medium">{properties.mw?.toFixed(2) ?? '-'} Da</span>
              </div>
              <div className="flex justify-between text-xs py-1 px-2 bg-gray-50 rounded">
                <span className="text-gray-600">LogP</span>
                <span className="font-mono font-medium">{properties.logp?.toFixed(2) ?? '-'}</span>
              </div>
              <div className="flex justify-between text-xs py-1 px-2 bg-gray-50 rounded">
                <span className="text-gray-600">HBD</span>
                <span className="font-mono font-medium">{properties.hbd ?? '-'}</span>
              </div>
              <div className="flex justify-between text-xs py-1 px-2 bg-gray-50 rounded">
                <span className="text-gray-600">HBA</span>
                <span className="font-mono font-medium">{properties.hba ?? '-'}</span>
              </div>
              <div className="flex justify-between text-xs py-1 px-2 bg-gray-50 rounded">
                <span className="text-gray-600">TPSA</span>
                <span className="font-mono font-medium">{properties.tpsa?.toFixed(1) ?? '-'} Å²</span>
              </div>
              <div className="flex justify-between text-xs py-1 px-2 bg-gray-50 rounded">
                <span className="text-gray-600">Rotatable</span>
                <span className="font-mono font-medium">{properties.rotatable ?? '-'}</span>
              </div>
              <div className="flex justify-between text-xs py-1 px-2 bg-gray-50 rounded">
                <span className="text-gray-600">Formula</span>
                <span className="font-mono font-medium text-xs">{properties.formula ?? '-'}</span>
              </div>
            </div>
          </div>

          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700 mb-2">Drug-like Rules</div>
            <div className="space-y-0.5">
              <RuleCheck pass={(properties.mw ?? 0) < 500} label="MW < 500 Da" value={properties.mw?.toFixed(1)} />
              <RuleCheck pass={(properties.logp ?? 0) < 5} label="LogP < 5" value={properties.logp?.toFixed(2)} />
              <RuleCheck pass={(properties.hbd ?? 0) <= 5} label="HBD ≤ 5" value={String(properties.hbd)} />
              <RuleCheck pass={(properties.hba ?? 0) <= 10} label="HBA ≤ 10" value={String(properties.hba)} />
              <RuleCheck pass={(properties.rotatable ?? 0) <= 10} label="Rotatable ≤ 10" value={String(properties.rotatable)} />
            </div>
          </div>

          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700 mb-2">Optimization</div>
            <div className="space-y-1">
              {MUTATION_STRATEGIES.map(s => (
                <button
                  key={s.key}
                  onClick={() => mutateMolecule(s.key)}
                  className="w-full text-xs px-2 py-1.5 bg-gray-50 hover:bg-cyan-50 border border-gray-200 rounded text-left transition-colors"
                >
                  {s.icon} {s.label}
                </button>
              ))}
            </div>
          </div>
        </div>

        <div className="flex-1 flex flex-col bg-gray-900">
          <div className="flex items-center gap-2 px-4 py-2 bg-slate-800 border-b border-slate-700">
            <button
              onClick={() => setActiveTab('2d')}
              className={`px-3 py-1 text-xs rounded transition-colors ${activeTab === '2d' ? 'bg-cyan-600 text-white' : 'text-gray-400 hover:text-white'}`}
            >
              2D Structure
            </button>
            <button
              onClick={() => setActiveTab('3d')}
              className={`px-3 py-1 text-xs rounded transition-colors ${activeTab === '3d' ? 'bg-cyan-600 text-white' : 'text-gray-400 hover:text-white'}`}
            >
              3D Viewer
            </button>
            <div className="flex-1" />
            <span className="text-xs text-cyan-400 font-medium">{molName}</span>
          </div>

          {activeTab === '2d' ? (
            <div className="flex-1 flex items-center justify-center bg-gray-900">
              <canvas
                ref={canvas2dRef}
                width={400}
                height={300}
                className="max-w-full max-h-full"
              />
            </div>
          ) : (
            <div
              ref={viewer3dRef}
              className="flex-1"
              style={{ background: '#1a1a2e' }}
            />
          )}

          <div className="bg-slate-800 border-t border-slate-700 px-4 py-2">
            <div className="flex items-center justify-between text-xs text-gray-400">
              <span>SMILES: {smiles}</span>
              <span>Drag to rotate | Scroll to zoom</span>
            </div>
          </div>
        </div>

        <div className="w-72 bg-white border-l border-gray-200 flex flex-col overflow-y-auto">
          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700">AI Suggestions</div>
          </div>
          <div className="flex-1 p-3 space-y-2">
            {suggestions.length === 0 ? (
              <p className="text-xs text-gray-500 text-center py-4">Click "Analyze" to get drug-likeness suggestions</p>
            ) : (
              suggestions.map((s, i) => (
                <div
                  key={i}
                  className={`text-xs p-2 rounded border-l-2 ${
                    s.type === 'good' ? 'bg-green-50 border-green-500 text-green-700' :
                    s.type === 'warning' ? 'bg-amber-50 border-amber-500 text-amber-700' :
                    s.type === 'error' ? 'bg-red-50 border-red-500 text-red-700' :
                    'bg-cyan-50 border-cyan-500 text-cyan-700'
                  }`}
                >
                  {s.text}
                </div>
              ))
            )}
          </div>
        </div>
      </div>
    </div>
  )
}
