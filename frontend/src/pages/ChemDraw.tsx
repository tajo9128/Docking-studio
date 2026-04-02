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

const HETEROCYCLES = [
  { name: 'Pyridine', smiles: 'c1ccncc1' },
  { name: 'Pyrimidine', smiles: 'c1cncnc1' },
  { name: 'Piperidine', smiles: 'C1CCNCC1' },
  { name: 'Pyrrole', smiles: 'c1cc[nH]c1' },
  { name: 'Imidazole', smiles: 'c1cnc[nH]1' },
  { name: 'Indole', smiles: 'c1ccc2[nH]ccc2c1' },
  { name: 'Quinoline', smiles: 'c1ccc2ncccc2c1' },
  { name: 'Isoquinoline', smiles: 'c1ccc2cnccc2c1' },
  { name: 'Thiazole', smiles: 'c1cncs1' },
  { name: 'Furan', smiles: 'c1ccoc1' },
  { name: 'Thiophene', smiles: 'c1ccsc1' },
  { name: 'Oxazole', smiles: 'c1cnco1' },
  { name: 'Pyrazine', smiles: 'c1cncnc1' },
  { name: 'Triazole', smiles: 'c1cnnc[nH]1' },
  { name: 'Purine', smiles: 'c1[nH]cnc2c1ncn2' },
  { name: 'Morpholine', smiles: 'C1COCCN1' },
  { name: 'Tetrahydropyran', smiles: 'C1CCOCC1' },
  { name: 'Piperazine', smiles: 'C1CNCCN1' },
  { name: 'Oxetane', smiles: 'C1COC1' },
  { name: 'Azetidine', smiles: 'C1CNC1' },
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
  const [molName, setMolName] = useState('Aspirin')
  const [properties, setProperties] = useState<Properties>({
    mw: 180.16, logp: 1.19, hbd: 1, hba: 3, tpsa: 63.60,
    rotatable: 4, formula: 'C9H8O4', valid: true
  })
  const [suggestions, setSuggestions] = useState<Suggestion[]>([])
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState<'2d' | '3d'>('2d')
  const [loadedFile, setLoadedFile] = useState<string | null>(null)
  const [iupacName, setIupacName] = useState<string>('')
  const [inchi, setInchi] = useState<string>('')
  const [inchiKey, setInchiKey] = useState<string>('')
  const [molBlock, setMolBlock] = useState<string>('')
  const [pdb3d, setPdb3d] = useState<string>('')
  const [conformers, setConformers] = useState<Array<{ energy: number; idx: number }>>([])
  const [selectedConformer, setSelectedConformer] = useState<number>(0)
  const [reactionMode, setReactionMode] = useState(false)
  const [reactionSmiles, setReactionSmiles] = useState('')
  const [showTemplates, setShowTemplates] = useState(false)
  const [templateFilter, setTemplateFilter] = useState('')
  const [showExport, setShowExport] = useState(false)
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

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.target instanceof HTMLInputElement || e.target instanceof HTMLTextAreaElement) return
      if (e.key === 'q' || e.key === 'Q') {
        e.preventDefault()
        cleanupStructure()
      }
      if (e.key === 'a' || e.key === 'A') {
        e.preventDefault()
        analyzeMolecule()
      }
    }
    window.addEventListener('keydown', handleKeyDown)
    return () => window.removeEventListener('keydown', handleKeyDown)
  }, [smiles])

  const draw2D = async () => {
    const canvas = canvas2dRef.current
    if (!canvas || !smiles) return

    try {
      if (typeof (window as any).SmilesDrawer === 'undefined') {
        await loadSmilesDrawer()
      }
      const drawer = new (window as any).SmilesDrawer.Drawer({ width: 600, height: 400, bondThickness: 2, compactDrawing: false })
      ;(window as any).SmilesDrawer.parse(smiles, (tree: any) => {
        drawer.draw(tree, canvas, 'dark', false)
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

      const res = await fetch('/api/chem/to-3d', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const data = await res.json()

      if (data.pdb) {
        stage.removeAllComponents()
        const blob = new Blob([data.pdb], { type: 'text/plain' })
        stage.loadFile(URL.createObjectURL(blob), { ext: 'pdb' }).then((comp: any) => {
          comp.addRepresentation('ball+stick', { color: 'element', radius: 0.3 })
          stage.autoView()
        })
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

  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>, fileType: string) => {
    const file = e.target.files?.[0]
    if (!file) return
    
    setLoading(true)
    try {
      const content = await file.text()
      
      if (fileType === 'sdf' || fileType === 'mol2' || fileType === 'pdb') {
        const res = await fetch('/api/chem/extract-smiles', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ content, format: fileType })
        })
        const data = await res.json()
        
        if (data.smiles) {
          setSmiles(data.smiles)
          setMolName(file.name.replace(/\.[^/.]+$/, ''))
          setLoadedFile(file.name)
        } else {
          const analyzeRes = await fetch('/api/chem/properties', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ content })
          })
          const analyzeData = await analyzeRes.json()
          if (analyzeData.mw) {
            setProperties({
              mw: analyzeData.mw,
              logp: analyzeData.logp,
              hbd: analyzeData.hbd,
              hba: analyzeData.hba,
              tpsa: analyzeData.tpsa,
              rotatable: analyzeData.rotatable_bonds,
              formula: analyzeData.formula,
              valid: true
            })
            setMolName(file.name.replace(/\.[^/.]+$/, ''))
            setLoadedFile(file.name)
          }
        }
      }
    } catch (err) {
      console.error('File upload error:', err)
    }
    setLoading(false)
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

  const cleanupStructure = async () => {
    if (!smiles) return
    try {
      const res = await fetch('/api/chem/cleanup', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const data = await res.json()
      if (data.cleaned) {
        setSmiles(data.cleaned)
        setSuggestions(prev => [{ text: '✓ Structure cleaned and canonicalized', type: 'good' }, ...prev.slice(0, 4)])
      }
    } catch (e) {
      console.error('Cleanup error:', e)
    }
  }

  const fetchIUPAC = async () => {
    if (!smiles) return
    try {
      const res = await fetch('/api/chem/iupac', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const data = await res.json()
      setIupacName(data.iupac || 'Not found')
    } catch (e) {
      setIupacName('OPSIN unavailable')
    }
  }

  const fetchInChI = async () => {
    if (!smiles) return
    try {
      const res = await fetch('/api/chem/inchi', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const data = await res.json()
      setInchi(data.inchi || '')
      setInchiKey(data.inchi_key || '')
    } catch (e) {
      console.error('InChI error:', e)
    }
  }

  const generateMolBlock = async () => {
    if (!smiles) return
    try {
      const res = await fetch('/api/chem/to-mol', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles })
      })
      const data = await res.json()
      setMolBlock(data.mol_block || '')
    } catch (e) {
      console.error('MOL error:', e)
    }
  }

  const generateConformers = async () => {
    if (!smiles) return
    setLoading(true)
    try {
      const res = await fetch('/api/chem/conformers', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles, n_conformers: 5 })
      })
      const data = await res.json()
      if (data.n_conformers) {
        const confList = data.energies.map((e: number, i: number) => ({ energy: e, idx: i }))
        confList.sort((a: any, b: any) => a.energy - b.energy)
        setConformers(confList)
        setSelectedConformer(confList[0].idx)
        setPdb3d(data.pdb || '')
        setSuggestions(prev => [{ text: `Generated ${data.n_conformers} conformers (best: ${confList[0].energy} kcal/mol)`, type: 'good' }, ...prev.slice(0, 4)])
      }
    } catch (e) {
      console.error('Conformer error:', e)
    }
    setLoading(false)
  }

  const addToReaction = () => {
    if (!smiles) return
    if (reactionSmiles) {
      setReactionSmiles(reactionSmiles + '.' + smiles)
    } else {
      setReactionSmiles(smiles)
    }
    setSuggestions(prev => [{ text: `Added ${molName || 'molecule'} to reaction`, type: 'info' }, ...prev.slice(0, 4)])
  }

  const clearReaction = () => {
    setReactionSmiles('')
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
    setIupacName('')
    setInchi('')
    setInchiKey('')
    setMolBlock('')
    setPdb3d('')
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

  const downloadFile = (content: string, filename: string, type: string) => {
    const blob = new Blob([content], { type })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = filename
    a.click()
    URL.revokeObjectURL(url)
  }

  const exportMOL = async () => {
    if (!molBlock) await generateMolBlock()
    if (molBlock) downloadFile(molBlock, `${molName || 'molecule'}.mol`, 'chemical/x-mdl-molfile')
  }

  const exportSDF = async () => {
    if (!molBlock) await generateMolBlock()
    if (molBlock) {
      const sdf = molBlock + '\n> <Name>\n' + (molName || 'Molecule') + '\n\n> <SMILES>\n' + smiles + '\n\n$$$$\n'
      downloadFile(sdf, `${molName || 'molecule'}.sdf`, 'chemical/x-mdl-sdfile')
    }
  }

  const exportInChI = async () => {
    if (!inchi) await fetchInChI()
    const content = `InChI=${inchi}\nInChIKey=${inchiKey}\nSMILES=${smiles}\n`
    downloadFile(content, `${molName || 'molecule'}.txt`, 'text/plain')
  }

  const exportPNG = async () => {
    const canvas = canvas2dRef.current
    if (!canvas) return
    const link = document.createElement('a')
    link.download = `${molName || 'molecule'}.png`
    link.href = canvas.toDataURL('image/png')
    link.click()
  }

  const RuleCheck = ({ pass, label, value }: { pass: boolean; label: string; value?: string }) => (
    <div className={`flex justify-between items-center text-xs py-1 px-2 rounded ${pass ? 'text-green-600' : 'text-red-500'}`}>
      <span>{pass ? '✓' : '✗'} {label}</span>
      {value && <span className="font-mono">{value}</span>}
    </div>
  )

  const filteredTemplates = HETEROCYCLES.filter(t =>
    t.name.toLowerCase().includes(templateFilter.toLowerCase()) ||
    t.smiles.toLowerCase().includes(templateFilter.toLowerCase())
  )

  return (
    <div className="h-full flex flex-col bg-gray-100">
      <div className="bg-slate-800 text-white px-4 py-2 flex items-center gap-3">
        <span className="text-sm font-medium text-cyan-400">ChemDraw</span>
        <span className="text-xs text-gray-400">Design, analyze & dock molecules</span>
        <div className="flex-1" />
        <span className="text-xs text-gray-500">Q: Cleanup | A: Analyze</span>
        <Button variant="outline" size="sm" onClick={cleanupStructure}>Cleanup</Button>
        <Button variant="outline" size="sm" onClick={clearEditor}>Clear</Button>
        <Button variant="outline" size="sm" onClick={copySmiles}>Copy SMILES</Button>
        <div className="relative">
          <Button variant="outline" size="sm" onClick={() => setShowExport(!showExport)}>Export ▾</Button>
          {showExport && (
            <div className="absolute right-0 top-full mt-1 w-40 bg-slate-700 border border-slate-600 rounded-lg shadow-lg z-50">
              {[
                { label: 'MOL File', action: exportMOL },
                { label: 'SDF File', action: exportSDF },
                { label: 'InChI/Key', action: exportInChI },
                { label: 'PNG Image', action: exportPNG },
              ].map(item => (
                <button key={item.label} onClick={() => { item.action(); setShowExport(false) }}
                  className="w-full text-left px-3 py-2 text-xs text-white hover:bg-slate-600 rounded-t-lg first:rounded-t-lg last:rounded-b-lg">
                  {item.label}
                </button>
              ))}
            </div>
          )}
        </div>
        <Button variant="outline" size="sm" onClick={analyzeMolecule} disabled={loading}>
          {loading ? '...' : 'Analyze'}
        </Button>
        <Button variant="primary" size="sm" onClick={dockMolecule}>Dock</Button>
        <Button variant="secondary" size="sm" onClick={optimizeMolecule}>AI Optimize</Button>
      </div>

      <div className="flex-1 flex overflow-hidden">
        <div className="w-72 bg-white border-r border-gray-200 flex flex-col overflow-y-auto">
          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700 mb-2">📤 Upload Molecule File</div>
            <div className="space-y-2">
              <label className="flex items-center justify-center gap-2 px-3 py-2 bg-blue-50 hover:bg-blue-100 border border-blue-200 rounded-lg cursor-pointer transition-colors text-sm text-blue-700">
                <span>📊</span>
                <span>Upload SDF File</span>
                <input type="file" accept=".sdf" className="hidden" onChange={e => handleFileUpload(e, 'sdf')} />
              </label>
              <label className="flex items-center justify-center gap-2 px-3 py-2 bg-green-50 hover:bg-green-100 border border-green-200 rounded-lg cursor-pointer transition-colors text-sm text-green-700">
                <span>🧬</span>
                <span>Upload MOL2 File</span>
                <input type="file" accept=".mol2" className="hidden" onChange={e => handleFileUpload(e, 'mol2')} />
              </label>
              <label className="flex items-center justify-center gap-2 px-3 py-2 bg-purple-50 hover:bg-purple-100 border border-purple-200 rounded-lg cursor-pointer transition-colors text-sm text-purple-700">
                <span>🧪</span>
                <span>Upload PDB File</span>
                <input type="file" accept=".pdb,.ent" className="hidden" onChange={e => handleFileUpload(e, 'pdb')} />
              </label>
            </div>
            {loadedFile && (
              <div className="mt-2 px-2 py-1 bg-gray-100 rounded text-xs text-gray-600 truncate">
                Loaded: {loadedFile}
              </div>
            )}
          </div>

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
            <button onClick={() => setShowTemplates(!showTemplates)}
              className="w-full text-left text-xs font-semibold text-gray-700 flex justify-between items-center">
              <span>📚 Template Library ({FDA_DRUGS.length + HETEROCYCLES.length})</span>
              <span>{showTemplates ? '▲' : '▼'}</span>
            </button>
            {showTemplates && (
              <div className="mt-2">
                <input type="text" placeholder="Search templates..." value={templateFilter}
                  onChange={e => setTemplateFilter(e.target.value)}
                  className="w-full px-2 py-1 text-xs border border-gray-300 rounded mb-2" />
                <div className="text-xs font-medium text-gray-500 mb-1">FDA Drugs</div>
                <div className="grid grid-cols-2 gap-1 mb-2">
                  {FDA_DRUGS.map(drug => (
                    <button key={drug.name}
                      onClick={() => { loadExample(drug.name, drug.smiles); setMolName(drug.name); setLoadedFile(null); }}
                      className="text-xs px-2 py-1.5 bg-gray-50 hover:bg-cyan-50 hover:text-cyan-700 border border-gray-200 rounded text-left transition-colors"
                      title={drug.use}>
                      {drug.name}
                    </button>
                  ))}
                </div>
                <div className="text-xs font-medium text-gray-500 mb-1">Heterocycles</div>
                <div className="grid grid-cols-2 gap-1 max-h-40 overflow-y-auto">
                  {filteredTemplates.map(t => (
                    <button key={t.name}
                      onClick={() => { loadExample(t.name, t.smiles); setMolName(t.name); }}
                      className="text-xs px-2 py-1.5 bg-gray-50 hover:bg-purple-50 hover:text-purple-700 border border-gray-200 rounded text-left transition-colors"
                      title={t.smiles}>
                      {t.name}
                    </button>
                  ))}
                </div>
              </div>
            )}
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
            <div className="text-xs font-semibold text-gray-700 mb-2">Identifiers</div>
            <div className="space-y-2">
              <div>
                <button onClick={fetchIUPAC} className="text-xs text-blue-600 hover:underline">
                  {iupacName ? `IUPAC: ${iupacName}` : 'Get IUPAC Name →'}
                </button>
              </div>
              <div>
                <button onClick={fetchInChI} className="text-xs text-blue-600 hover:underline">
                  {inchiKey ? `InChIKey: ${inchiKey}` : 'Get InChI/InChIKey →'}
                </button>
              </div>
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

          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700 mb-2">3D Conformer Generator</div>
            <button onClick={generateConformers} disabled={loading || !smiles}
              className="w-full text-xs px-2 py-1.5 bg-purple-50 hover:bg-purple-100 border border-purple-200 rounded text-purple-700 transition-colors disabled:opacity-50">
              {loading ? 'Generating...' : 'Generate 5 Conformers (ETKDG)'}
            </button>
            {conformers.length > 0 && (
              <div className="mt-2 space-y-1">
                {conformers.map((c, i) => (
                  <button key={c.idx}
                    onClick={() => setSelectedConformer(c.idx)}
                    className={`w-full text-xs px-2 py-1 rounded text-left transition-colors ${
                      selectedConformer === c.idx
                        ? 'bg-purple-100 border border-purple-300 text-purple-700'
                        : 'bg-gray-50 border border-gray-200 text-gray-600'
                    }`}>
                    #{i+1} {c.energy.toFixed(1)} kcal/mol {i === 0 && '⭐'}
                  </button>
                ))}
              </div>
            )}
          </div>

          <div className="p-3 border-b border-gray-200">
            <div className="text-xs font-semibold text-gray-700 mb-2">Reaction Mode</div>
            <button onClick={() => setReactionMode(!reactionMode)}
              className={`w-full text-xs px-2 py-1.5 rounded border transition-colors ${
                reactionMode ? 'bg-orange-100 border-orange-300 text-orange-700' : 'bg-gray-50 border-gray-200 text-gray-600'
              }`}>
              {reactionMode ? '🧪 Reaction Active' : '🧪 Enter Reaction Mode'}
            </button>
            {reactionMode && (
              <div className="mt-2 space-y-2">
                <button onClick={addToReaction}
                  className="w-full text-xs px-2 py-1.5 bg-blue-50 hover:bg-blue-100 border border-blue-200 rounded text-blue-700">
                  + Add Current to Reaction
                </button>
                <button onClick={clearReaction}
                  className="w-full text-xs px-2 py-1.5 bg-red-50 hover:bg-red-100 border border-red-200 rounded text-red-700">
                  Clear Reaction
                </button>
                {reactionSmiles && (
                  <div className="text-xs font-mono bg-gray-50 p-2 rounded break-all">
                    {reactionSmiles}
                  </div>
                )}
              </div>
            )}
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
            {reactionMode && (
              <span className="px-2 py-0.5 text-xs bg-orange-600 text-white rounded-full">🧪 Reaction</span>
            )}
            <div className="flex-1" />
            <span className="text-xs text-cyan-400 font-medium">{molName}</span>
          </div>

          {activeTab === '2d' ? (
            <div className="flex-1 flex items-center justify-center bg-gray-900">
              <canvas
                ref={canvas2dRef}
                width={600}
                height={400}
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
              <span className="truncate max-w-md">SMILES: {smiles}</span>
              <span>Q: Cleanup | A: Analyze | Drag to rotate | Scroll to zoom</span>
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
