import { apiClient } from '@/lib/apiClient'

export interface NotificationConfig {
  channel: string
  webhook_url?: string
  email?: string
  enabled: boolean
}

export async function getNotificationStatus(): Promise<{ channels: NotificationConfig[] }> {
  const { data } = await apiClient.get('/notifications/status')
  return data
}

export async function configureNotification(config: NotificationConfig): Promise<{ success: boolean }> {
  const { data } = await apiClient.post('/notifications/configure', config)
  return data
}

export async function sendNotification(event: string, title: string, message: string): Promise<{ sent_to: string[] }> {
  const { data } = await apiClient.post('/notifications/send', { event, title, message })
  return data
}

export async function testNotification(channel: string = 'discord'): Promise<{ sent_to: string[] }> {
  const { data } = await apiClient.post('/notifications/test', { channel })
  return data
}
