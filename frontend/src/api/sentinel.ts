import { apiClient } from '@/lib/apiClient'

export interface SentinelStatus {
  queue_length: number
  active_jobs: number
  failed_jobs: number
  retry_count: number
}

export async function getSentinelStatus(): Promise<SentinelStatus> {
  const { data } = await apiClient.get('/sentinel/queue/status')
  return data
}

export async function monitorJob(jobId: string): Promise<any> {
  const { data } = await apiClient.post('/sentinel/monitor', { job_id: jobId })
  return data
}

export async function retryJob(jobId: string): Promise<{ success: boolean }> {
  const { data } = await apiClient.post('/sentinel/retry', { job_id: jobId })
  return data
}

export async function fallbackJob(jobId: string, fallbackEngine: string): Promise<any> {
  const { data } = await apiClient.post('/sentinel/fallback', { job_id: jobId, fallback_engine: fallbackEngine })
  return data
}

export async function escalateJob(jobId: string, reason: string): Promise<any> {
  const { data } = await apiClient.post('/sentinel/escalate', { job_id: jobId, reason })
  return data
}

export async function validateResult(jobId: string, validationResult: any): Promise<any> {
  const { data } = await apiClient.post('/sentinel/validate/result', { job_id: jobId, validation: validationResult })
  return data
}
