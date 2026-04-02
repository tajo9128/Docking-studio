import { apiClient } from '@/lib/apiClient'

export interface BenchmarkRequest {
  benchmark_type: 'pdbbind' | 'rmsd' | 'enrichment'
  dataset?: string
  parameters?: Record<string, any>
}

export interface BenchmarkResult {
  job_id: string
  status: string
  metrics?: Record<string, number>
  report_url?: string
}

export async function runBenchmark(request: BenchmarkRequest): Promise<{ job_id: string }> {
  const { data } = await apiClient.post('/api/benchmark/run', request)
  return data
}

export async function getBenchmarkResults(jobId: string): Promise<BenchmarkResult> {
  const { data } = await apiClient.get(`/api/benchmark/results/${jobId}`)
  return data
}
