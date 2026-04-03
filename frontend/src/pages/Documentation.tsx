import { useState } from 'react'
import { useTheme } from '@/contexts/ThemeContext'

const sections = [
  { id: 'intro', title: '1. Introduction', icon: '📘' },
  { id: 'getting-started', title: '2. Getting Started', icon: '🚀' },
  { id: 'ligand-designer', title: '3. Ligand Designer', icon: '⬡' },
  { id: 'protein-prep', title: '4. Protein Preparation', icon: '🧬' },
  { id: 'docking', title: '5. Docking Engine', icon: '🎯' },
  { id: 'results', title: '6. Results & Visualization', icon: '📊' },
  { id: 'intelligence', title: '7. Intelligence Panel', icon: '🧠' },
  { id: 'advanced', title: '8. Advanced Features', icon: '⚡' },
  { id: 'ai', title: '9. AI & Optimization', icon: '🤖' },
  { id: 'developer', title: '10. Developer Docs', icon: '💻' },
  { id: 'deployment', title: '11. Deployment', icon: '☁️' },
  { id: 'formats', title: '12. File Formats', icon: '📁' },
  { id: 'validation', title: '13. Validation & Limitations', icon: '⚠️' },
  { id: 'troubleshooting', title: '14. Troubleshooting', icon: '🔧' },
  { id: 'roadmap', title: '15. Roadmap', icon: '🗺️' },
]

const comparisons = [
  {
    feature: 'Cost',
    biodockify: 'Free (MIT License)',
    discoveryStudio: 'Commercial ($15K-30K/year)',
    schrodinger: 'Commercial ($20K-40K/year)',
  },
  {
    feature: 'Docking Engines',
    biodockify: 'Vina + GNINA (CNN)',
    discoveryStudio: 'CDOCKER, LibDock',
    schrodinger: 'Glide',
  },
  {
    feature: 'AI Integration',
    biodockify: 'Built-in LLM assistant',
    discoveryStudio: 'Limited',
    schrodinger: 'Advanced (separate)',
  },
  {
    feature: 'Source Code',
    biodockify: 'Fully Open',
    discoveryStudio: 'Proprietary',
    schrodinger: 'Proprietary',
  },
  {
    feature: 'Extensibility',
    biodockify: 'Unlimited (APIs)',
    discoveryStudio: 'Limited',
    schrodinger: 'Limited',
  },
]

export function Documentation() {
  const { theme } = useTheme()
  const isDark = theme === 'dark'
  const [activeSection, setActiveSection] = useState('intro')

  const scrollToSection = (id: string) => {
    setActiveSection(id)
    document.getElementById(id)?.scrollIntoView({ behavior: 'smooth' })
  }

  return (
    <div className={`min-h-full ${isDark ? 'bg-gray-900 text-white' : 'bg-white text-gray-800'}`}>
      {/* Header */}
      <div className={`border-b ${isDark ? 'border-gray-700 bg-gray-800' : 'border-gray-200 bg-gray-50'} px-6 py-4`}>
        <h1 className="text-2xl font-bold flex items-center gap-2">
          <span>📘</span>
          <span>BioDockify Complete Guide</span>
        </h1>
        <p className={`text-sm mt-1 ${isDark ? 'text-gray-400' : 'text-gray-600'}`}>
          Product Manual • Developer Reference • Scientific Guide • v4.2.6
        </p>
      </div>

      <div className="flex">
        {/* Sidebar Navigation */}
        <aside className={`w-64 h-[calc(100vh-80px)] overflow-y-auto border-r ${isDark ? 'border-gray-700 bg-gray-800' : 'border-gray-200 bg-gray-50'} sticky top-0`}>
          <nav className="p-4 space-y-1">
            {sections.map((section) => (
              <button
                key={section.id}
                onClick={() => scrollToSection(section.id)}
                className={`w-full text-left px-3 py-2 rounded-lg text-sm transition-colors flex items-center gap-2 ${
                  activeSection === section.id
                    ? isDark
                      ? 'bg-blue-600 text-white'
                      : 'bg-blue-100 text-blue-700'
                    : isDark
                      ? 'text-gray-300 hover:bg-gray-700'
                      : 'text-gray-700 hover:bg-gray-200'
                }`}
              >
                <span>{section.icon}</span>
                <span>{section.title}</span>
              </button>
            ))}
          </nav>
        </aside>

        {/* Main Content */}
        <main className="flex-1 p-8 max-w-4xl">
          {/* Section 1: Introduction */}
          <section id="intro" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              1. Introduction
            </h2>
            <p className={`mb-4 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              BioDockify Studio is a unified, open-source platform for molecular docking, ligand design, 
              and computational drug discovery. It integrates 2D molecule drawing (Ketcher), molecular 
              docking (AutoDock Vina + GNINA), molecular dynamics (OpenMM), and medicinal chemistry 
              analysis into a single web-based interface.
            </p>
            
            <h3 className={`text-lg font-semibold mt-6 mb-3 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Key Capabilities
            </h3>
            <ul className={`list-disc pl-5 space-y-1 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <li>Ligand Design — 2D/3D molecule drawing with 60+ templates</li>
              <li>Molecular Docking — AutoDock Vina + GNINA CNN scoring</li>
              <li>MD Simulation — OpenMM-based molecular dynamics</li>
              <li>Medicinal Chemistry — Real-time property calculation</li>
              <li>AI Assistant — Natural language job explanation</li>
              <li>Job History — Persistent SQLite database</li>
            </ul>

            <div className={`mt-6 p-4 rounded-lg ${isDark ? 'bg-blue-900/30 border border-blue-700' : 'bg-blue-50 border border-blue-200'}`}>
              <p className={`text-sm italic ${isDark ? 'text-blue-200' : 'text-blue-800'}`}>
                <strong>Positioning:</strong> BioDockify provides a unified, accessible platform integrating 
                molecular drawing, docking, and medicinal chemistry analysis. While commercial platforms 
                offer highly optimized proprietary algorithms, BioDockify emphasizes flexibility, 
                transparency, and accessibility through open-source tools.
              </p>
            </div>
          </section>

          {/* Comparison Table */}
          <section className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              Comparison with Standard Software
            </h2>
            <div className={`overflow-x-auto rounded-lg border ${isDark ? 'border-gray-700' : 'border-gray-300'}`}>
              <table className={`w-full text-sm ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
                <thead className={isDark ? 'bg-gray-800' : 'bg-gray-100'}>
                  <tr>
                    <th className="px-4 py-3 text-left font-semibold">Feature</th>
                    <th className="px-4 py-3 text-left font-semibold text-cyan-500">BioDockify</th>
                    <th className="px-4 py-3 text-left font-semibold">Discovery Studio</th>
                    <th className="px-4 py-3 text-left font-semibold">Schrödinger</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-700/30">
                  {comparisons.map((row, idx) => (
                    <tr key={idx} className={isDark ? 'hover:bg-gray-800/50' : 'hover:bg-gray-50'}>
                      <td className="px-4 py-3 font-medium">{row.feature}</td>
                      <td className="px-4 py-3 text-cyan-500">{row.biodockify}</td>
                      <td className="px-4 py-3">{row.discoveryStudio}</td>
                      <td className="px-4 py-3">{row.schrodinger}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <p className={`mt-3 text-sm ${isDark ? 'text-gray-400' : 'text-gray-600'}`}>
              <strong>Note:</strong> BioDockify offers comparable core workflows with advantages in 
              transparency, extensibility, and zero licensing costs.
            </p>
          </section>

          {/* Section 2: Getting Started */}
          <section id="getting-started" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              2. Getting Started
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              System Requirements
            </h3>
            <div className={`grid grid-cols-2 gap-4 mt-3 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <div className={`p-3 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium mb-2">Minimum (CPU)</p>
                <ul className="text-sm space-y-1">
                  <li>• CPU: 4 cores, 2.5 GHz</li>
                  <li>• RAM: 8 GB</li>
                  <li>• Storage: 10 GB</li>
                  <li>• Docker Desktop 4.0+</li>
                </ul>
              </div>
              <div className={`p-3 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium mb-2">Recommended (GPU)</p>
                <ul className="text-sm space-y-1">
                  <li>• GPU: NVIDIA CUDA 11.8+</li>
                  <li>• VRAM: 4 GB+</li>
                  <li>• CPU: 8+ cores</li>
                  <li>• RAM: 16 GB+</li>
                </ul>
              </div>
            </div>

            <h3 className={`text-lg font-semibold mt-6 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Quick Start (Docker)
            </h3>
            <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${isDark ? 'bg-gray-800 text-gray-300' : 'bg-gray-900 text-gray-100'}`}>
{`git clone https://github.com/tajo9128/BioDockify-Studio-AI.git
cd BioDockify-Studio-AI
docker-compose up --build
# Access at http://localhost:3000`}
            </pre>
          </section>

          {/* Section 3: Ligand Designer */}
          <section id="ligand-designer" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              3. Ligand Designer Module
            </h2>
            <p className={`mb-4 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              The Ligand Designer provides a professional 2D molecular drawing interface powered by 
              Ketcher React. It supports chemical structure drawing, template libraries, and real-time 
              property calculation.
            </p>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Template Library (60+)
            </h3>
            <div className={`grid grid-cols-3 gap-3 text-sm ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium">FDA Drugs (12)</p>
                <p className="text-xs opacity-75">Aspirin, Ibuprofen, Caffeine...</p>
              </div>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium">Heterocycles (16)</p>
                <p className="text-xs opacity-75">Pyridine, Imidazole, Indole...</p>
              </div>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium">Amino Acids (8)</p>
                <p className="text-xs opacity-75">Alanine, Glycine, Cysteine...</p>
              </div>
            </div>

            <div className={`mt-4 p-3 rounded-lg ${isDark ? 'bg-amber-900/20 border border-amber-700' : 'bg-amber-50 border border-amber-200'}`}>
              <p className={`text-sm ${isDark ? 'text-amber-200' : 'text-amber-800'}`}>
                <strong>Comparison:</strong> Advanced template libraries in commercial tools are broader 
                (500+ templates). BioDockify focuses on commonly used medicinal chemistry templates.
              </p>
            </div>
          </section>

          {/* Section 4: Protein Preparation */}
          <section id="protein-prep" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              4. Protein Preparation
            </h2>
            <p className={`mb-4 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              BioDockify supports PDB and PDBQT formats for protein structures. Automatic processing 
              includes water removal, hydrogen addition, and charge assignment.
            </p>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Automatic Processing Steps
            </h3>
            <ol className={`list-decimal pl-5 space-y-1 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <li>Remove waters — Delete HOH residues to reduce noise</li>
              <li>Remove ligands — Extract HETATM groups</li>
              <li>Add hydrogens — Protonation at pH 7.4</li>
              <li>Repair chains — Fix missing atoms</li>
              <li>Assign charges — Gasteiger method for docking</li>
            </ol>
          </section>

          {/* Section 5: Docking Engine */}
          <section id="docking" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              5. Docking Engine (Core)
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Supported Engines
            </h3>
            <div className={`space-y-3 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-semibold text-cyan-500">AutoDock Vina 1.2.5</p>
                <p className="text-sm mt-1">Iterated local search global optimizer. Empirical force field scoring. 
                Best for standard virtual screening and pose prediction. Speed: ~1-2 min per ligand.</p>
              </div>
              <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-semibold text-cyan-500">GNINA 1.0 (CNN Scoring)</p>
                <p className="text-sm mt-1">Monte Carlo + CNN rescoring. Convolutional Neural Network scoring. 
                Higher accuracy on diverse sets. Speed: ~2-5 min per ligand.</p>
              </div>
            </div>

            <div className={`mt-4 p-3 rounded-lg ${isDark ? 'bg-green-900/20 border border-green-700' : 'bg-green-50 border border-green-200'}`}>
              <p className={`text-sm ${isDark ? 'text-green-200' : 'text-green-800'}`}>
                <strong>Key Advantage:</strong> GNINA's CNN scoring provides ML advantage for challenging 
                cases, though it requires more computation time.
              </p>
            </div>
          </section>

          {/* Section 6: Results */}
          <section id="results" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              6. Results & Visualization
            </h2>
            <p className={`mb-4 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              The Results panel displays docking outputs in multiple formats: score tables, 
              3D interactive viewer, interaction diagrams, and download options.
            </p>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Energy Interpretation Guide
            </h3>
            <table className={`w-full text-sm ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <thead className={isDark ? 'bg-gray-800' : 'bg-gray-100'}>
                <tr>
                  <th className="px-3 py-2 text-left">Affinity (kcal/mol)</th>
                  <th className="px-3 py-2 text-left">Interpretation</th>
                  <th className="px-3 py-2 text-left">Quality</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-700/30">
                <tr><td className="px-3 py-2">&lt; -10</td><td className="px-3 py-2">Very strong binding</td><td className="px-3 py-2 text-green-500">Excellent</td></tr>
                <tr><td className="px-3 py-2">-8 to -10</td><td className="px-3 py-2">Strong binding</td><td className="px-3 py-2 text-green-400">Good</td></tr>
                <tr><td className="px-3 py-2">-6 to -8</td><td className="px-3 py-2">Moderate binding</td><td className="px-3 py-2 text-yellow-500">Marginal</td></tr>
                <tr><td className="px-3 py-2">&gt; -6</td><td className="px-3 py-2">Weak binding</td><td className="px-3 py-2 text-red-400">Poor</td></tr>
              </tbody>
            </table>
          </section>

          {/* Section 7: Intelligence Panel */}
          <section id="intelligence" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              7. Intelligence Panel
            </h2>
            <p className={`mb-4 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              The Intelligence Panel provides real-time medicinal chemistry analysis as you design molecules. 
              This is a unique integrated feature combining multiple cheminformatics tools.
            </p>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Drug-Likeness Filters
            </h3>
            <div className={`grid grid-cols-2 gap-4 text-sm ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium mb-2">Lipinski's Rule of Five</p>
                <ul className="space-y-1 text-xs">
                  <li>• MW &lt; 500 Da</li>
                  <li>• LogP &lt; 5</li>
                  <li>• HBD ≤ 5</li>
                  <li>• HBA ≤ 10</li>
                </ul>
              </div>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium mb-2">Veber Rules</p>
                <ul className="space-y-1 text-xs">
                  <li>• TPSA ≤ 140 Å²</li>
                  <li>• Rotatable bonds ≤ 10</li>
                </ul>
              </div>
            </div>

            <div className={`mt-4 p-3 rounded-lg ${isDark ? 'bg-cyan-900/20 border border-cyan-700' : 'bg-cyan-50 border border-cyan-200'}`}>
              <p className={`text-sm ${isDark ? 'text-cyan-200' : 'text-cyan-800'}`}>
                <strong>USP:</strong> BioDockify integrates multiple medicinal chemistry tools into a 
                single interactive panel with real-time updates—a workflow advantage over commercial tools.
              </p>
            </div>
          </section>

          {/* Section 8: Advanced Features */}
          <section id="advanced" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              8. Advanced Features
            </h2>
            
            <div className={`space-y-4 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-semibold">NMR Prediction (Rule-Based)</p>
                <p className="text-sm mt-1">Rule-based estimation for ¹H and ¹³C chemical shifts. 
                <strong className={isDark ? 'text-red-400' : 'text-red-600'}> For educational purposes only.</strong></p>
              </div>
              
              <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-semibold">3D Conformer Generation</p>
                <p className="text-sm mt-1">ETKDG algorithm with MMFF94 optimization. Generate multiple 
                conformers ranked by energy.</p>
              </div>
              
              <div className={`p-4 rounded-lg ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-semibold">PubChem Similarity Search</p>
                <p className="text-sm mt-1">Tanimoto fingerprint search against PubChem database. 
                Returns similar compounds with CID and name.</p>
              </div>
            </div>
          </section>

          {/* Section 9: AI */}
          <section id="ai" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              9. AI & Optimization
            </h2>
            <p className={`mb-4 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              BioDockify AI provides natural language explanations for docking results, error 
              troubleshooting, and scientific concepts.
            </p>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Supported AI Providers
            </h3>
            <ul className={`list-disc pl-5 space-y-1 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <li><strong>Ollama</strong> — Local models (llama3.2, mistral)</li>
              <li><strong>OpenAI</strong> — Cloud API (GPT-4, GPT-3.5)</li>
              <li><strong>DeepSeek</strong> — Alternative cloud provider</li>
            </ul>

            <div className={`mt-4 p-3 rounded-lg ${isDark ? 'bg-amber-900/20 border border-amber-700' : 'bg-amber-50 border border-amber-200'}`}>
              <p className={`text-sm ${isDark ? 'text-amber-200' : 'text-amber-800'}`}>
                <strong>Disclaimer:</strong> AI predictions are for research assistance only. 
                Always validate with experimental data and expert knowledge.
              </p>
            </div>
          </section>

          {/* Section 10: Developer */}
          <section id="developer" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              10. Developer Documentation
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Architecture
            </h3>
            <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${isDark ? 'bg-gray-800 text-gray-300' : 'bg-gray-900 text-gray-100'}`}>
{`Frontend (React 18 + TypeScript)
    ↓ HTTP/WebSocket
Backend (FastAPI + Python 3.9+)
    ↓
Core (RDKit, OpenMM, Vina)
    ↓
Storage (SQLite + Filesystem)`}
            </pre>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Key API Endpoints
            </h3>
            <table className={`w-full text-sm ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <thead className={isDark ? 'bg-gray-800' : 'bg-gray-100'}>
                <tr>
                  <th className="px-3 py-2 text-left">Endpoint</th>
                  <th className="px-3 py-2 text-left">Method</th>
                  <th className="px-3 py-2 text-left">Purpose</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-700/30">
                <tr><td className="px-3 py-2 font-mono text-xs">/api/docking/run</td><td className="px-3 py-2">POST</td><td className="px-3 py-2">Submit docking job</td></tr>
                <tr><td className="px-3 py-2 font-mono text-xs">{"/api/jobs/{id}/full"}</td><td className="px-3 py-2">GET</td><td className="px-3 py-2">Complete job data</td></tr>
                <tr><td className="px-3 py-2 font-mono text-xs">/api/ai/job-explain</td><td className="px-3 py-2">POST</td><td className="px-3 py-2">AI explanation</td></tr>
                <tr><td className="px-3 py-2 font-mono text-xs">/api/chem/properties</td><td className="px-3 py-2">POST</td><td className="px-3 py-2">Calculate properties</td></tr>
              </tbody>
            </table>
          </section>

          {/* Section 11: Deployment */}
          <section id="deployment" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              11. Deployment Guide
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Docker Commands
            </h3>
            <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${isDark ? 'bg-gray-800 text-gray-300' : 'bg-gray-900 text-gray-100'}`}>
{`# Build and start
docker-compose up --build

# GPU support (NVIDIA)
docker-compose -f docker-compose.yml up

# Reset environment
docker-compose down
docker system prune -a`}
            </pre>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Scaling Options
            </h3>
            <ul className={`list-disc pl-5 space-y-1 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <li><strong>Single user:</strong> Local Docker (personal research)</li>
              <li><strong>Lab (5-10 users):</strong> Single server deployment</li>
              <li><strong>Department:</strong> Kubernetes cluster</li>
              <li><strong>Enterprise:</strong> Cloud auto-scaling (AWS/GCP/Azure)</li>
            </ul>
          </section>

          {/* Section 12: File Formats */}
          <section id="formats" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              12. File Formats & Standards
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Supported Formats
            </h3>
            <table className={`w-full text-sm ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <thead className={isDark ? 'bg-gray-800' : 'bg-gray-100'}>
                <tr>
                  <th className="px-3 py-2 text-left">Format</th>
                  <th className="px-3 py-2 text-left">Type</th>
                  <th className="px-3 py-2 text-left">Input</th>
                  <th className="px-3 py-2 text-left">Output</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-700/30">
                <tr><td className="px-3 py-2">PDB</td><td className="px-3 py-2">Protein</td><td className="px-3 py-2 text-green-500">✅</td><td className="px-3 py-2 text-green-500">✅</td></tr>
                <tr><td className="px-3 py-2">PDBQT</td><td className="px-3 py-2">Docking</td><td className="px-3 py-2 text-green-500">✅</td><td className="px-3 py-2 text-green-500">✅</td></tr>
                <tr><td className="px-3 py-2">MOL/SDF</td><td className="px-3 py-2">Molecule</td><td className="px-3 py-2 text-green-500">✅</td><td className="px-3 py-2 text-green-500">✅</td></tr>
                <tr><td className="px-3 py-2">SMILES</td><td className="px-3 py-2">String</td><td className="px-3 py-2 text-green-500">✅</td><td className="px-3 py-2 text-green-500">✅</td></tr>
                <tr><td className="px-3 py-2">XYZ</td><td className="px-3 py-2">Coordinates</td><td className="px-3 py-2 text-green-500">✅</td><td className="px-3 py-2 text-red-400">❌</td></tr>
              </tbody>
            </table>
          </section>

          {/* Section 13: Validation */}
          <section id="validation" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              13. Validation & Limitations
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Docking Limitations
            </h3>
            <ul className={`list-disc pl-5 space-y-1 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <li>Pose prediction: ~70-80% success rate (validate visually)</li>
              <li>Scoring correlation: R² ~0.5-0.7 (use for ranking, not absolute)</li>
              <li>Water effects: Often ignored (consider explicit waters)</li>
              <li>Protein flexibility: Rigid receptor only</li>
            </ul>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Critical Disclaimers
            </h3>
            <div className={`space-y-2 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <div className={`p-3 rounded ${isDark ? 'bg-red-900/20 border border-red-700' : 'bg-red-50 border border-red-200'}`}>
                <p className="text-sm"><strong>NMR Predictions:</strong> Rule-based estimates only. 
                NOT for structural determination.</p>
              </div>
              <div className={`p-3 rounded ${isDark ? 'bg-red-900/20 border border-red-700' : 'bg-red-50 border border-red-200'}`}>
                <p className="text-sm"><strong>Retrosynthesis:</strong> Basic rules only. 
                NOT for production route design.</p>
              </div>
            </div>

            <div className={`mt-4 p-3 rounded-lg ${isDark ? 'bg-blue-900/20 border border-blue-700' : 'bg-blue-50 border border-blue-200'}`}>
              <p className={`text-sm ${isDark ? 'text-blue-200' : 'text-blue-800'}`}>
                <strong>Credibility Statement:</strong> BioDockify provides accessible entry-level 
                computational chemistry tools. For mission-critical applications, validate with 
                experiments and consider specialized commercial tools.
              </p>
            </div>
          </section>

          {/* Section 14: Troubleshooting */}
          <section id="troubleshooting" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              14. Troubleshooting
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Common Issues
            </h3>
            <div className={`space-y-3 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium">Invalid SMILES</p>
                <p className="text-sm">Check SMILES syntax using online validator. Ensure proper bracket notation.</p>
              </div>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium">Docking Failure</p>
                <p className="text-sm">Check grid box size (20+ Å), verify ligand charges, ensure Docker is running.</p>
              </div>
              <div className={`p-3 rounded ${isDark ? 'bg-gray-800' : 'bg-gray-100'}`}>
                <p className="font-medium">Slow Performance</p>
                <p className="text-sm">Reduce exhaustiveness (32→8), use fewer CPU threads, check available RAM.</p>
              </div>
            </div>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Reset Commands
            </h3>
            <pre className={`p-4 rounded-lg overflow-x-auto text-sm ${isDark ? 'bg-gray-800 text-gray-300' : 'bg-gray-900 text-gray-100'}`}>
{`# Reset Docker
docker-compose down
docker system prune -a
docker-compose up --build

# Reset database (WARNING: data loss)
rm backend/storage/jobs.db`}
            </pre>
          </section>

          {/* Section 15: Roadmap */}
          <section id="roadmap" className="mb-12">
            <h2 className={`text-2xl font-bold mb-4 ${isDark ? 'text-white' : 'text-gray-900'}`}>
              15. Roadmap
            </h2>
            
            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Version 4.3 (Q2 2026)
            </h3>
            <ul className={`list-disc pl-5 space-y-1 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <li>DiffDock integration (diffusion model docking)</li>
              <li>Enhanced 3D viewer with WebGL 2.0</li>
              <li>Batch job templates</li>
              <li>Improved AI assistant context</li>
            </ul>

            <h3 className={`text-lg font-semibold mt-4 mb-2 ${isDark ? 'text-gray-200' : 'text-gray-800'}`}>
              Version 5.0 (Q3 2026)
            </h3>
            <ul className={`list-disc pl-5 space-y-1 ${isDark ? 'text-gray-300' : 'text-gray-700'}`}>
              <li>Multi-user collaboration</li>
              <li>Cloud storage integration</li>
              <li>Advanced virtual screening workflows</li>
              <li>ML-based activity prediction</li>
            </ul>

            <div className={`mt-6 p-4 rounded-lg ${isDark ? 'bg-gray-800 border border-gray-700' : 'bg-gray-100 border border-gray-300'}`}>
              <p className={`text-sm ${isDark ? 'text-gray-400' : 'text-gray-600'}`}>
                <strong>Support:</strong> GitHub Issues • API Docs: http://localhost:8000/docs • License: MIT
              </p>
            </div>
          </section>
        </main>
      </div>
    </div>
  )
}
