import { apiClient } from '@/lib/apiClient'

export interface RankRequest {
  compounds: Array<{ smiles: string; score?: number }>
  method?: 'vina' | 'gnina' | 'consensus' | 'ml'
}

export interface RankResult {
  ranked: Array<{ smiles: string; score: number; rank: number }>
  method: string
}

export interface ConsensusRequest {
  compounds: Array<{ smiles: string; vina_score?: number; gnina_score?: number; rf_score?: number }>
  weights?: { vina: number; gnina: number; rf: number }
}

export interface ConsensusResult {
  compounds: Array<{ smiles: string; consensus_score: number; rank: number }>
}

export interface AnalysisADMETFilterRequest {
  compounds: Array<{ smiles: string }>
  rules: {
    mw_max?: number
    logp_max?: number
    tpsa_max?: number
    hbd_max?: number
    hba_max?: number
  }
}

export interface AnalysisADMETFilterResult {
  passed: Array<{ smiles: string }>
  failed: Array<{ smiles: string; reason: string }>
  pass_rate: number
}

export interface ReportRequest {
  job_id: string
  format?: 'pdf' | 'html' | 'json'
  include_structures?: boolean
}

export async function rankLigands(request: RankRequest): Promise<RankResult> {
  const { data } = await apiClient.post('/analysis/rank', request)
  return data
}

export async function filterByADMETAnalysis(request: AnalysisADMETFilterRequest): Promise<AnalysisADMETFilterResult> {
  const { data } = await apiClient.post('/analysis/filter/admet', request)
  return data
}

export async function consensusScoring(request: ConsensusRequest): Promise<ConsensusResult> {
  const { data } = await apiClient.post('/analysis/consensus', request)
  return data
}

export async function generateReport(request: ReportRequest): Promise<{ report_url: string; format: string }> {
  const { data } = await apiClient.post('/analysis/report', request)
  return data
}

export async function getInteractionsSummary(jobId: string): Promise<any> {
  const { data } = await apiClient.post('/analysis/interactions/summary', { job_id: jobId })
  return data
}
