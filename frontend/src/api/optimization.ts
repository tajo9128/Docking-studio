import { apiClient } from '@/lib/apiClient'

export interface LeadOptimizationRequest {
  lead_smiles: string
  target_property?: string
  max_variants?: number
}

export interface LeadOptimizationResult {
  variants: Array<{
    smiles: string
    score: number
    modifications: string[]
  }>
  best_variant?: string
}

export async function optimizeLead(request: LeadOptimizationRequest): Promise<LeadOptimizationResult> {
  const { data } = await apiClient.post('/optimize/lead', request)
  return data
}

export async function mutateLead(leadSmiles: string, mutationType: string = 'random'): Promise<{ mutated_smiles: string }> {
  const { data } = await apiClient.post('/optimize/mutate', { smiles: leadSmiles, mutation_type: mutationType })
  return data
}

export async function scoreVariant(smiles: string): Promise<{ score: number; properties: any }> {
  const { data } = await apiClient.post('/optimize/variant/score', { smiles })
  return data
}
