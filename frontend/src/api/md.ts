import { apiClient } from '@/lib/apiClient'

export interface MDDynamicsRequest {
  pdb_content: string
  steps: number
  temperature: number
  pressure: number
  frame_interval: number
  solvent_model: string
  ionic_strength: number
  name: string
  notify_on_start: boolean
  notify_on_complete: boolean
}

export interface MDJobResponse {
  job_id: string
  status: string
  message: string
}

export interface MDJobStatus {
  status: 'pending' | 'running' | 'completed' | 'failed'
  progress: number
  message: string
  updated_at: string
  result?: MDResult
  error?: string
}

export interface MDResult {
  trajectory_path: string
  final_frame_path: string
  energy_csv_path: string
  n_frames: number
  n_steps: number
  sim_time_ns: number
  temperature_K: number
  avg_energy_kj_mol: number
  solvent_model: string
}

export interface MDAnalysisRequest {
  job_id: string
  trajectory_path?: string
  energy_csv_path?: string
}

export interface MDAnalysisResult {
  success: boolean
  output_file: string
  plot_data?: PlotlyData
  [key: string]: any
}

export interface PlotlyData {
  data: any[]
  layout: any
}

export interface MDPublicationRequest {
  job_id: string
  project_name: string
  analysis_job_id?: string
  compress: boolean
  notify_on_complete: boolean
}

export interface MDNotifyStatus {
  telegram: boolean
  discord: boolean
  slack: boolean
  email: boolean
}

export async function runDynamics(request: MDDynamicsRequest): Promise<MDJobResponse> {
  const { data } = await apiClient.post('/md/dynamics', request)
  return data
}

export async function getMDJobStatus(jobId: string): Promise<MDJobStatus> {
  const { data } = await apiClient.get(`/md/job/${jobId}`)
  return data
}

export async function analyzeRMSD(request: MDAnalysisRequest): Promise<MDAnalysisResult> {
  const { data } = await apiClient.post('/md/analysis/rmsd', request)
  return data
}

export async function analyzeRMSF(request: MDAnalysisRequest): Promise<MDAnalysisResult> {
  const { data } = await apiClient.post('/md/analysis/rmsf', request)
  return data
}

export async function analyzeEnergy(request: MDAnalysisRequest): Promise<MDAnalysisResult> {
  const { data } = await apiClient.post('/md/analysis/energy', request)
  return data
}

export async function analyzeGyration(request: MDAnalysisRequest): Promise<MDAnalysisResult> {
  const { data } = await apiClient.post('/md/analysis/gyration', request)
  return data
}

export async function analyzeSASA(request: MDAnalysisRequest): Promise<MDAnalysisResult> {
  const { data } = await apiClient.post('/md/analysis/sasa', request)
  return data
}

export async function analyzeHBonds(request: MDAnalysisRequest): Promise<MDAnalysisResult> {
  const { data } = await apiClient.post('/md/analysis/hbonds', request)
  return data
}

export async function runFullAnalysis(request: MDAnalysisRequest): Promise<MDJobResponse> {
  const { data } = await apiClient.post('/md/analysis/all', request)
  return data
}

export async function createPublication(request: MDPublicationRequest): Promise<{ success: boolean; package_path: string }> {
  const { data } = await apiClient.post('/md/publication/package', request)
  return data
}

export async function getNotifyStatus(): Promise<MDNotifyStatus> {
  const { data } = await apiClient.get('/md/notify/status')
  return data
}

export async function testNotification(channel: string = 'discord'): Promise<{ sent_to: string[] }> {
  const { data } = await apiClient.post('/md/notify/test', null, { params: { channel } })
  return data
}

export async function sendNotification(event: string, title: string, message: string): Promise<{ sent_to: string[] }> {
  const { data } = await apiClient.post('/md/notify', { event, title, message })
  return data
}

export async function minimizeStructure(pdbContent: string): Promise<{ job_id: string; status: string }> {
  const { data } = await apiClient.post('/md/minimize', null, { params: { pdb_content: pdbContent } })
  return data
}

export async function getMDHealth(): Promise<{ status: string; engine: string }> {
  const { data } = await apiClient.get('/md/health')
  return data
}

export interface MDEquilibrationRequest {
  pdb_content: string
  temperature: number
  pressure: number
  solvent_model: string
  ionic_strength: number
  name: string
}

export interface MMEquilibriumResult {
  equilibrated_pdb: string
  checkpoint: string
  minimization_energy_kj_mol: number
  nvt_energy_kj_mol: number
  npt_energy_kj_mol: number
  n_atoms: number
  ready_for_production: boolean
}

export async function runEquilibration(_request: MDEquilibrationRequest): Promise<MDJobResponse> {
  // Note: Backend route not yet implemented
  console.warn('MD equilibration endpoint not available on backend')
  return { job_id: '', status: 'error', message: 'Equilibration not yet available' }
}

export async function resumeSimulation(_jobId: string, _steps: number = 50000, _frameInterval: number = 500): Promise<MDJobResponse> {
  // Note: Backend route not yet implemented
  console.warn('MD resume endpoint not available on backend')
  return { job_id: '', status: 'error', message: 'Resume not yet available' }
}

export async function calculateMMGBSA(_request: { trajectory_path: string; receptor_pdb: string; ligand_pdb: string }): Promise<any> {
  // Note: Backend route not yet implemented
  console.warn('MD MM-GBSA endpoint not available on backend')
  return { success: false, error: 'MM-GBSA not yet available' }
}

export async function getMDGPUStatus(): Promise<any> {
  // Note: Backend route not yet implemented, use /api/md/gpu-info instead
  try {
    const { data } = await apiClient.get('/api/md/gpu-info')
    return data
  } catch {
    return { available: false, gpus: [], error: 'GPU info not available' }
  }
}
