import { apiClient } from '@/lib/apiClient'
import type { PoseAnalysis, PoseRequest, RMSDResponse, BindingSiteResidue } from '@/lib/types'

export async function analyzePose(receptor: string, ligand: string): Promise<PoseAnalysis> {
  const { data } = await apiClient.post<PoseAnalysis>('/analyze', { receptor, ligand } as PoseRequest)
  return data
}

export async function analyzeAdvanced(receptor: string, ligand: string): Promise<PoseAnalysis> {
  const { data } = await apiClient.post<PoseAnalysis>('/analyze/advanced', { receptor, ligand } as PoseRequest)
  return data
}

export async function calculateRMSD(pdb1: string, pdb2: string): Promise<RMSDResponse> {
  const { data } = await apiClient.post<RMSDResponse>('/rmsd', { pdb1, pdb2 })
  return data
}

export async function getBindingSite(
  receptor: string,
  ligand: string
): Promise<{ residues: BindingSiteResidue[] }> {
  const { data } = await apiClient.post('/binding-site', { receptor, ligand } as PoseRequest)
  return data
}

export interface DockingHit {
  ligand_id: string
  vina_score?: number
  gnina_score?: number
  rf_score?: number
  [key: string]: any
}

export interface ExportTopHitsResponse {
  format: string
  content?: string
  results?: DockingHit[]
  filename: string
  count: number
}

export async function exportTopHits(
  dockingResults: DockingHit[],
  topN: number = 10,
  sortBy: string = 'vina_score',
  format: 'csv' | 'json' = 'csv'
): Promise<ExportTopHitsResponse> {
  const { data } = await apiClient.post<ExportTopHitsResponse>('/analysis/export/top-hits', {
    docking_results: dockingResults,
    top_n: topN,
    sort_by: sortBy,
    format,
  })
  return data
}
