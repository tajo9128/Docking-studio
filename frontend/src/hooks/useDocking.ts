import { useEffect, useState, useCallback } from 'react'
import { getDockingStatus } from '@/api/docking'
import type { DockingProgress } from '@/lib/types'

export function useDockingStream(jobId: string | null) {
  const [progress, setProgress] = useState<DockingProgress | null>(null)
  const [error, setError] = useState<string | null>(null)

  const refetch = useCallback(async () => {
    if (!jobId) return null
    try {
      setError(null)
      const status = await getDockingStatus(jobId)
      setProgress(status)
      return status
    } catch (e) {
      setError('Failed to fetch docking status')
      return null
    }
  }, [jobId])

  useEffect(() => {
    if (!jobId) return

    refetch()

    const interval = setInterval(() => {
      refetch()
    }, 3000)

    return () => clearInterval(interval)
  }, [jobId, refetch])

  return { progress, error, refetch }
}
