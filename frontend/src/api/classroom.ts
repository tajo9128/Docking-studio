import { apiClient } from '@/lib/apiClient'

export interface Assignment {
  id: string
  title: string
  description: string
  deadline: string
  created_by: string
  status: 'open' | 'closed'
}

export interface Submission {
  id: string
  assignment_id: string
  student_name: string
  submitted_at: string
  score?: number
  feedback?: string
}

export async function createAssignment(data: {
  title: string
  description: string
  deadline: string
}): Promise<{ success: boolean; assignment_id: string }> {
  const { data: result } = await apiClient.post('/classroom/assignment/create', data)
  return result
}

export async function joinAssignment(code: string, studentName: string): Promise<{ success: boolean }> {
  const { data } = await apiClient.post('/classroom/assignment/join', { code, student_name: studentName })
  return data
}

export async function submitAssignment(assignmentId: string, studentName: string, submissionData: any): Promise<{ success: boolean }> {
  const { data } = await apiClient.post('/classroom/assignment/submit', {
    assignment_id: assignmentId,
    student_name: studentName,
    submission: submissionData,
  })
  return data
}

export async function getInstructorDashboard(instructorId: string): Promise<{ assignments: Assignment[]; submissions: Submission[] }> {
  const { data } = await apiClient.get(`/classroom/instructor/${instructorId}`)
  return data
}

export async function getRubrics(): Promise<{ rubrics: any[] }> {
  const { data } = await apiClient.get('/classroom/rubrics')
  return data
}
