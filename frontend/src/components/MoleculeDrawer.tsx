import { useState, useRef, useEffect } from 'react'
import { Modal, Button } from '@/components/ui'

interface MoleculeDrawerProps {
  isOpen: boolean
  onClose: () => void
  onSave: (smiles: string, molFormat: string) => void
  initialSmiles?: string
}

export function MoleculeDrawer({ isOpen, onClose, onSave, initialSmiles }: MoleculeDrawerProps) {
  const [smiles, setSmiles] = useState(initialSmiles || '')
  const [molFile, setMolFile] = useState('')
  const [format, setFormat] = useState<'smiles' | 'mol'>('smiles')
  const ketcherRef = useRef<any>(null)
  const [ketcherLoaded, setKetcherLoaded] = useState(false)
  const [KetcherComponent, setKetcherComponent] = useState<any>(null)

  useEffect(() => {
    if (isOpen && !ketcherLoaded) {
      loadKetcher()
    }
  }, [isOpen])

  async function loadKetcher() {
    try {
      const editorMod = await import('ketcher-react')
      const standaloneMod = await import('ketcher-standalone')
      setKetcherComponent(() => editorMod.Editor)
      // Store service provider constructor for use in child
      ;(window as any).__ketcherServiceProvider = standaloneMod.StandaloneStructServiceProvider
      setKetcherLoaded(true)
    } catch (e) {
      console.warn('Ketcher could not be loaded:', e)
    }
  }

  async function handleSave() {
    if (ketcherRef.current) {
      try {
        const structure = await ketcherRef.current.getStructure()
        const smilesStr = structure?.molecular?.toSmiles()
        if (smilesStr) {
          onSave(smilesStr, 'smiles')
          onClose()
          return
        }
      } catch (e) {
        console.warn('Could not get SMILES from Ketcher:', e)
      }
    }
    if (format === 'smiles' && smiles.trim()) {
      onSave(smiles.trim(), 'smiles')
      onClose()
    } else if (format === 'mol' && molFile.trim()) {
      onSave(molFile.trim(), 'mol')
      onClose()
    }
  }

  async function handleImportMol() {
    if (!molFile.trim()) return
    if (ketcherRef.current) {
      try {
        await ketcherRef.current.setMolecule(molFile)
      } catch (e) {
        console.warn('Could not import MOL:', e)
      }
    }
  }

  return (
    <Modal isOpen={isOpen} onClose={onClose} title="Draw Molecule" size="xl">
      <div className="space-y-4">
        <div className="flex gap-2 border-b pb-3">
          <button
            onClick={() => setFormat('smiles')}
            className={`px-4 py-2 text-sm font-medium rounded ${
              format === 'smiles'
                ? 'bg-primary text-white'
                : 'bg-gray-700 text-gray-300 hover:bg-gray-600'
            }`}
          >
            SMILES Input
          </button>
          <button
            onClick={() => setFormat('mol')}
            className={`px-4 py-2 text-sm font-medium rounded ${
              format === 'mol'
                ? 'bg-primary text-white'
                : 'bg-gray-700 text-gray-300 hover:bg-gray-600'
            }`}
          >
            MOL File
          </button>
        </div>

        {format === 'smiles' ? (
          <div>
            <label className="block text-sm font-medium text-gray-300 mb-1">
              SMILES String
            </label>
            <input
              type="text"
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              placeholder="Enter SMILES, e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
              className="w-full bg-gray-800 border border-gray-600 rounded-lg px-3 py-2 text-white text-sm"
            />
            <p className="text-xs text-gray-500 mt-1">
              Or use the MOL File tab above to draw and import
            </p>
          </div>
        ) : (
          <div>
            <label className="block text-sm font-medium text-gray-300 mb-1">
              MOL File Content
            </label>
            <textarea
              value={molFile}
              onChange={(e) => setMolFile(e.target.value)}
              placeholder="Paste MOL file content here..."
              rows={8}
              className="w-full bg-gray-800 border border-gray-600 rounded-lg px-3 py-2 text-white text-sm font-mono"
            />
            <button
              onClick={handleImportMol}
              disabled={!molFile.trim()}
              className="mt-2 px-3 py-1.5 text-xs bg-gray-700 text-gray-300 rounded hover:bg-gray-600 disabled:opacity-50"
            >
              Load into Editor
            </button>
          </div>
        )}

        {ketcherLoaded && KetcherComponent ? (
          <div className="border border-gray-600 rounded-lg overflow-hidden">
            <KetcherComponent
              staticResourcesUrl="/ketcher"
              structServiceProvider={new (window as any).__ketcherServiceProvider()}
              onInit={(ketcher: any) => { ketcherRef.current = ketcher }}
            />
          </div>
        ) : (
          <div className="h-[400px] bg-gray-800 rounded-lg flex items-center justify-center text-gray-500">
            <div className="text-center">
              <p className="mb-2">Ketcher editor loading...</p>
              <p className="text-xs">You can still paste SMILES above</p>
            </div>
          </div>
        )}

        <div className="flex justify-end gap-3 pt-3 border-t">
          <Button variant="secondary" onClick={onClose}>
            Cancel
          </Button>
          <Button onClick={handleSave} disabled={format === 'smiles' && !smiles.trim()}>
            Use Structure
          </Button>
        </div>
      </div>
    </Modal>
  )
}
