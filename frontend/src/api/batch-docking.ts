import { apiClient } from '@/lib/apiClient'

export interface BatchDockingRequest {
  receptor_path: string
  ligand_library: string[]
  exhaustiveness?: number
  num_modes?: number
  center_x?: number
  center_y?: number
  center_z?: number
  size_x?: number
  size_y?: number
  size_z?: number
}

export interface BatchDockingResult {
  job_id: string
  status: string
  total_compounds: number
  results?: Array<{
    smiles: string
    best_score: number
    best_pose?: any
  }>
}

export async function startBatchDocking(request: BatchDockingRequest): Promise<{ job_id: string }> {
  const { data } = await apiClient.post('/batch/docking', request)
  return data
}

export async function getBatchDockingProgress(jobId: string): Promise<{ progress: number; status: string }> {
  const { data } = await apiClient.get(`/batch/docking/${jobId}/progress`)
  return data
}
