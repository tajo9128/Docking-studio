import { apiClient } from '@/lib/apiClient'

export interface ModificationResult {
  parent_smiles: string
  modified_smiles: string
  applied_transform: string
  properties: {
    mw: number
    logp: number
    hbd: number
    hba: number
    rot_bonds: number
    tpsa: number
  }
  docking_score?: number
  delta_score?: number
  source: 'database' | 'transformation' | 'autonomous'
}

export interface LigandModifierJob {
  job_id: string
  status: 'queued' | 'parsing' | 'fetching' | 'transforming' | 'validating' | 'docking' | 'completed' | 'failed' | 'cancelled'
  progress: number
  results: ModificationResult[]
  error?: string
}

export interface LigandModifierRequest {
  parent_smiles: string
  receptor_pdb: string
  mode: 'similarity_search' | 'prompt_based' | 'autonomous'
  prompt?: string
  database?: 'pubchem' | 'chembl'
  similarity_threshold?: number
  max_variants?: number
  docking_exhaustiveness?: number
}

export const ligandModifierAPI = {
  optimize: async (req: LigandModifierRequest): Promise<{ job_id: string; status: string; message: string }> => {
    const { data } = await apiClient.post('/api/ligand-modifier/optimize', req)
    return data
  },

  getStatus: async (jobId: string): Promise<LigandModifierJob> => {
    const { data } = await apiClient.get(`/api/ligand-modifier/status/${jobId}`)
    return data
  },

  cancel: async (jobId: string): Promise<{ message: string }> => {
    const { data } = await apiClient.delete(`/api/ligand-modifier/cancel/${jobId}`)
    return data
  },

  pollUntilComplete: async (
    jobId: string,
    onProgress?: (status: LigandModifierJob) => void,
    intervalMs: number = 2000,
    timeoutMs: number = 300000
  ): Promise<LigandModifierJob> => {
    const startTime = Date.now()

    while (Date.now() - startTime < timeoutMs) {
      const status = await ligandModifierAPI.getStatus(jobId)
      onProgress?.(status)

      if (['completed', 'failed', 'cancelled'].includes(status.status)) {
        return status
      }

      await new Promise(resolve => setTimeout(resolve, intervalMs))
    }

    throw new Error('Polling timeout')
  }
}
