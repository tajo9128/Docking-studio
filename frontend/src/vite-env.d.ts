/// <reference types="vite/client" />

// Global type declarations for 3Dmol.js and smiles-drawer
declare global {
  interface Window {
    $3Dmol: any
  }
}

// SmilesDrawer module declaration
declare module 'smiles-drawer' {
  export default class SmilesDrawer {
    static parse(smiles: string, onSuccess: (tree: any) => void, onError: (error: any) => void): void
    constructor(options?: any)
    draw(tree: any, canvas: HTMLCanvasElement, theme?: string, infoOnly?: boolean): void
  }
  export class Drawer {
    constructor(options?: any)
    draw(tree: any, canvas: HTMLCanvasElement, theme?: string, infoOnly?: boolean): void
  }
}

export {}
