import { apiClient } from '@/lib/apiClient'

export interface PrepareProteinResponse {
  success: boolean
  pdb_path: string
  original_atoms: number
  final_atoms: number
  waters_removed: boolean
  hydrogens_added: boolean
}

export interface PrepareLigandResponse {
  success: boolean
  pdbqt_path?: string
  pdb_path?: string
  num_atoms: number
  meeko_used: boolean
  message?: string
}

export interface Interaction {
  ligand_atom: string
  receptor_atom: string
  distance_A: number
  type: 'H-bond' | 'hydrophobic'
  ligand_pos?: number[]
  receptor_pos?: number[]
}

export interface DetectInteractionsResponse {
  success: boolean
  h_bonds: Interaction[]
  hydrophobic_contacts: Interaction[]
  total_h_bonds: number
  total_hydrophobic: number
}

export async function prepareProtein(
  pdbContent: string,
  name: string = 'protein',
  removeWaters: boolean = true,
  addHydrogens: boolean = true
): Promise<PrepareProteinResponse> {
  const { data } = await apiClient.post<PrepareProteinResponse>('/rdkit/prepare_protein', {
    pdb_content: pdbContent,
    name,
    remove_waters: removeWaters,
    add_hydrogens: addHydrogens,
  })
  return data
}

export async function prepareLigand(
  pdbContent: string,
  name: string = 'ligand'
): Promise<PrepareLigandResponse> {
  const { data } = await apiClient.post<PrepareLigandResponse>('/rdkit/prepare_ligand', {
    pdb_content: pdbContent,
    name,
  })
  return data
}

export async function detectInteractions(
  receptorPdbContent: string,
  ligandPdbContent: string
): Promise<DetectInteractionsResponse> {
  const { data } = await apiClient.post<DetectInteractionsResponse>('/rdkit/detect_interactions', {
    receptor_pdb_content: receptorPdbContent,
    ligand_pdb_content: ligandPdbContent,
  })
  return data
}
