import { apiClient } from '@/lib/apiClient'

export interface DescriptorGroupsResponse {
  groups: string[]
  descriptors: Record<string, string[]>
}

export interface DescriptorsResponse {
  success: boolean
  descriptors: Record<string, number>[]
  valid_smiles: string[]
  failed_smiles: string[]
  n_valid: number
  n_failed: number
}

export interface DatasetUploadResult {
  X: number[][]
  y: number[]
  feature_names: string[]
  n_compounds: number
  n_features: number
  activity_mean: number
  activity_std: number
  activity_min: number
  activity_max: number
  nan_count: number
  failed_smiles: string[]
  failed_count: number
}

export interface TrainJobResponse {
  job_id: string
  status: string
  message: string
}

export interface TrainJobStatus {
  status: 'pending' | 'running' | 'completed' | 'failed'
  updated_at: string
  result?: TrainResult
  error?: string
}

export interface TrainResult {
  model_id: string
  model_name: string
  model_type: string
  metrics: {
    cv_r2: number
    cv_rmse: number
    cv_mae: number
    train_r2: number
    cv_std: number
  }
  cv_scores: number[]
  cv_folds: number
  n_compounds: number
  n_features: number
  scatter_plot: any
  created_at: string
}

export interface PredictionSingle {
  success: boolean
  smiles: string
  predicted_activity: number
  ad_status: 'in_domain' | 'warning' | 'out_of_domain' | 'unknown'
  ad_leverage?: number
  ad_warning_threshold?: number
  ad_danger_threshold?: number
}

export interface PredictionBatch {
  success: boolean
  predictions: {
    smiles: string
    predicted_activity: number | null
    ad_status: string
    ad_leverage?: number
    error?: string
  }[]
  n_total: number
  n_failed: number
  n_in_domain: number
  n_warning: number
  n_out_of_domain: number
}

export interface SavedModel {
  model_id: string
  name: string
  model_type: string
  feature_names: string[]
  n_features: number
  metrics: Record<string, number>
  activity_column: string
  descriptor_groups: string[]
  created_at: string
}

export async function getDescriptorGroups(): Promise<DescriptorGroupsResponse> {
  const { data } = await apiClient.get('/qsar/descriptor-groups')
  return data
}

export async function calculateDescriptors(
  smiles: string[],
  groups?: string[]
): Promise<DescriptorsResponse> {
  const { data } = await apiClient.post('/qsar/descriptors', { smiles, groups })
  return data
}

export async function uploadDataset(
  file: File,
  smilesCol: string,
  activityCol: string,
  groups?: string
): Promise<DatasetUploadResult> {
  const formData = new FormData()
  formData.append('file', file)
  formData.append('smiles_col', smilesCol)
  formData.append('activity_col', activityCol)
  if (groups) formData.append('groups', groups)

  const { data } = await apiClient.post('/qsar/descriptors/upload', formData, {
    headers: { 'Content-Type': 'multipart/form-data' },
  })
  return data
}

export async function startTraining(
  X: number[][],
  y: number[],
  featureNames: string[],
  modelType: string,
  modelName: string,
  activityColumn: string,
  descriptorGroups: string[],
  cvFolds: number = 5,
  modelParams?: Record<string, any>
): Promise<TrainJobResponse> {
  const { data } = await apiClient.post('/qsar/train', {
    X,
    y,
    feature_names: featureNames,
    model_type: modelType,
    model_name: modelName,
    activity_column: activityColumn,
    descriptor_groups: descriptorGroups,
    cv_folds: cvFolds,
    model_params: modelParams,
  })
  return data
}

export async function getTrainingStatus(jobId: string): Promise<TrainJobStatus> {
  const { data } = await apiClient.get(`/qsar/train/${jobId}/status`)
  return data
}

export async function getTrainingResults(jobId: string): Promise<{ status: string; result?: TrainResult; error?: string }> {
  const { data } = await apiClient.get(`/qsar/train/${jobId}/results`)
  return data
}

export async function predictSingle(
  modelId: string,
  smiles: string
): Promise<PredictionSingle> {
  const { data } = await apiClient.post('/qsar/predict', { model_id: modelId, smiles })
  return data
}

export async function predictBatch(
  modelId: string,
  smilesList: string[]
): Promise<PredictionBatch> {
  const { data } = await apiClient.post('/qsar/predict/batch', {
    model_id: modelId,
    smiles_list: smilesList,
  })
  return data
}

export async function listModels(): Promise<{ models: SavedModel[] }> {
  const { data } = await apiClient.get('/qsar/models')
  return data
}

export async function getModel(modelId: string): Promise<SavedModel> {
  const { data } = await apiClient.get(`/qsar/models/${modelId}`)
  return data
}

export async function deleteModel(modelId: string): Promise<{ success: boolean }> {
  const { data } = await apiClient.delete(`/qsar/models/${modelId}`)
  return data
}
