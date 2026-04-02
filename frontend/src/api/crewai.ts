import { apiClient } from '@/lib/apiClient'

export interface CrewStatus {
  status: string
  active_crews: number
  available_agents: number
}

export interface CrewAgent {
  name: string
  role: string
  goal: string
  backstory: string
  tools: string[]
}

export interface CrewKickoffRequest {
  crew_name: string
  inputs: Record<string, any>
}

export async function getCrewStatus(): Promise<CrewStatus> {
  const { data } = await apiClient.get('/crew/status')
  return data
}

export async function listAgents(): Promise<{ agents: CrewAgent[] }> {
  const { data } = await apiClient.get('/crew/agents')
  return data
}

export async function listCrews(): Promise<{ crews: Array<{ name: string; agents: string[]; tasks: number }> }> {
  const { data } = await apiClient.get('/crew/crews')
  return data
}

export async function kickoffCrew(request: CrewKickoffRequest): Promise<{ job_id: string }> {
  const { data } = await apiClient.post('/crew/kickoff', request)
  return data
}

export async function getCrewJobStatus(jobId: string): Promise<any> {
  const { data } = await apiClient.get(`/crew/job/${jobId}`)
  return data
}

export async function crewChat(message: string, conversationHistory?: any[]): Promise<any> {
  const { data } = await apiClient.post('/crew/chat', { message, history: conversationHistory })
  return data
}

export async function getMemoryStats(): Promise<any> {
  const { data } = await apiClient.get('/crew/memory/stats')
  return data
}

export async function getMemoryByExperiment(expId: string): Promise<any> {
  const { data } = await apiClient.get(`/crew/memory/${expId}`)
  return data
}

export async function getMemoryFailures(): Promise<any> {
  const { data } = await apiClient.get('/crew/memory/failures')
  return data
}

export async function validateTool(toolName: string): Promise<{ valid: boolean; error?: string }> {
  const { data } = await apiClient.post('/crew/validate/tool', { tool_name: toolName })
  return data
}

export async function orchestrateCrew(inputs: Record<string, any>): Promise<{ job_id: string }> {
  const { data } = await apiClient.post('/crew/orchestrate', { inputs })
  return data
}
