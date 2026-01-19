# BioDockify Docking Studio
## Complete Software Documentation | v1.1.0

---

# 📋 Executive Summary

**BioDockify Docking Studio** is a unified desktop application for molecular docking, visualization, and drug discovery research. It combines docking engines, molecular analysis tools, and AI-powered insights into a single professional platform.

| Metric | Value |
|--------|-------|
| **Version** | 1.1.0 |
| **Platform** | Windows 10/11 (64-bit) |
| **Technology** | Python, PyQt6, Docker |
| **License** | MIT |

---

# 🏗️ System Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        BioDockify Docking Studio                        │
├─────────────────────────────────────────────────────────────────────────┤
│  ┌─────────────────────────────────────────────────────────────────┐   │
│  │                    PRESENTATION LAYER (PyQt6)                    │   │
│  │  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌────────┐ │   │
│  │  │  Upload  │ │  Config  │ │ Progress │ │ Results  │ │Agent 0 │ │   │
│  │  │  Widget  │ │  Widget  │ │  Widget  │ │  Widget  │ │ Widget │ │   │
│  │  └──────────┘ └──────────┘ └──────────┘ └──────────┘ └────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                    │                                    │
│  ┌─────────────────────────────────┴───────────────────────────────┐   │
│  │                      SERVICE LAYER (FastAPI)                     │   │
│  │  ┌────────────────┐  ┌────────────────┐  ┌────────────────────┐ │   │
│  │  │  Job Manager   │  │ Analysis Svc   │  │   Parsing Svc      │ │   │
│  │  └────────────────┘  └────────────────┘  └────────────────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                    │                                    │
│  ┌─────────────────────────────────┴───────────────────────────────┐   │
│  │                         CORE LAYER                               │   │
│  │  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌────────┐ │   │
│  │  │ Engines  │ │ Parsers  │ │Analyzers │ │   Math   │ │Validators│ │   │
│  │  └──────────┘ └──────────┘ └──────────┘ └──────────┘ └────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                    │                                    │
│  ┌─────────────────────────────────┴───────────────────────────────┐   │
│  │                    INFRASTRUCTURE LAYER                          │   │
│  │  ┌────────────┐  ┌────────────┐  ┌────────────┐  ┌────────────┐ │   │
│  │  │   Docker   │  │  Database  │  │ Checkpoint │  │  Recovery  │ │   │
│  │  │  Manager   │  │  (SQLite)  │  │  Manager   │  │  Manager   │ │   │
│  │  └────────────┘  └────────────┘  └────────────┘  └────────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
│                                    │                                    │
│  ┌─────────────────────────────────┴───────────────────────────────┐   │
│  │                    EXTERNAL ENGINES (Docker)                     │   │
│  │  ┌──────────────────┐  ┌──────────────────┐  ┌────────────────┐ │   │
│  │  │  AutoDock Vina   │  │      Gnina       │  │    RDKit       │ │   │
│  │  └──────────────────┘  └──────────────────┘  └────────────────┘ │   │
│  └─────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────┘
```

---

# ✨ Complete Feature List

## 1. Molecular Docking Engine
| Feature | Description | Status |
|---------|-------------|--------|
| AutoDock Vina 1.2 | High-precision empirical scoring | ✅ |
| Gnina 1.0 | CNN-based deep learning scoring | ✅ |
| Flexible Ligand Docking | Full torsional flexibility | ✅ |
| Rigid Receptor Docking | Fixed protein backbone | ✅ |
| Flexible Residue Docking | Flexible binding site sidechains | ✅ |
| Grid Box Auto-Detection | Automatic search space | ✅ |
| Batch Docking | Multiple ligands simultaneously | ✅ |
| Re-scoring | Score poses with alternative methods | ✅ |

## 2. File Format Support
| Format | Read | Write | Description |
|--------|------|-------|-------------|
| PDB | ✅ | ✅ | Protein Data Bank |
| PDBQT | ✅ | ✅ | AutoDock format |
| MOL2 | ✅ | ✅ | Tripos SYBYL |
| SDF | ✅ | ✅ | MDL Structure Data |
| SMILES | ✅ | ✅ | Canonical SMILES |
| XYZ | ✅ | ❌ | Cartesian coordinates |

## 3. Molecular Analysis (BioDockviz)
| Analyzer | Capability |
|----------|------------|
| **Bond Detector** | Covalent bonds, bond orders, aromaticity |
| **Interaction Analyzer** | H-bonds, hydrophobic, π-stacking, salt bridges |
| **Molecular Engine** | Atom types, charges, ring systems |
| **Spatial Hash** | Fast 3D neighbor searches |

## 4. Drug Discovery Tools
| Tool | Purpose |
|------|---------|
| Lipinski's Rule of 5 | Oral bioavailability prediction |
| Veber's Rules | Oral drug-likeness |
| PAINS Filter | Pan-assay interference detection |
| Molecular Weight | MW calculation |
| LogP/LogD | Lipophilicity estimation |
| PSA/TPSA | Polar surface area |
| HBA/HBD Count | Hydrogen bond acceptors/donors |
| Rotatable Bonds | Flexibility metric |
| ADMET Prediction | Absorption, distribution, metabolism, excretion, toxicity |

## 5. Visualization Features
| Feature | Description |
|---------|-------------|
| 3D Molecule Viewer | Interactive WebGL rendering |
| Protein Surface | Solvent-accessible/molecular surface |
| Binding Pocket | Cavity detection and display |
| Interaction Diagram | 2D ligand interaction map |
| Pose Animation | Animated docking trajectory |
| Overlay Mode | Compare multiple poses |
| Measurement Tool | Distance, angle, dihedral |
| Screenshot Export | High-resolution image export |

## 6. AI Assistant (Agent Zero)
| Capability | Description |
|------------|-------------|
| Natural Language Query | Ask questions about results |
| Result Interpretation | Automated binding analysis |
| Parameter Suggestion | Optimal docking settings |
| Literature Search | Related publications |
| SMILES Generation | Generate molecules from description |
| Report Generation | Automated summary reports |

## 7. Job Management
| Feature | Description |
|---------|-------------|
| Job Queue | FIFO priority queue |
| Checkpointing | Resume from failure |
| Progress Tracking | Real-time status |
| Error Recovery | Automatic retry logic |
| Logging | Detailed execution logs |
| Resource Management | CPU/memory limits |

## 8. Data Management
| Feature | Description |
|---------|-------------|
| SQLite Database | Local persistent storage |
| Project Organization | Folder-based projects |
| Result History | Track all docking runs |
| Export Options | CSV, JSON, PDF reports |
| Backup/Restore | Data portability |

---

# 📁 Project Structure

```
Docking-studio/
├── src/                          # Source Code
│   ├── biodockify_main.py        # Entry point
│   ├── config.py                 # Configuration
│   ├── database.py               # SQLite database
│   ├── docker_manager.py         # Container management
│   ├── vina_engine.py            # Vina integration
│   ├── oddt_analyzer.py          # ODDT analysis
│   ├── rdkit_calculator.py       # RDKit calculations
│   ├── checkpoint_manager.py     # Job checkpoints
│   ├── recovery_manager.py       # Error recovery
│   ├── agent_zero.py             # AI assistant
│   │
│   ├── ui/                       # Desktop UI (PyQt6)
│   │   ├── main_window.py        # Main window
│   │   ├── upload_widget.py      # File upload
│   │   ├── configuration_widget.py
│   │   ├── progress_widget.py
│   │   ├── results_widget.py
│   │   ├── agent_zero_widget.py
│   │   └── theme.py              # UI theming
│   │
│   ├── core/                     # BioDockviz Core
│   │   ├── engines/              # Molecular engines
│   │   ├── parsers/              # File parsers
│   │   ├── analyzers/            # Analysis algorithms
│   │   ├── math/                 # Numerical utilities
│   │   └── validators.py         # Input validation
│   │
│   ├── api/                      # REST API
│   ├── services/                 # Business logic
│   ├── models/                   # Data models
│   └── schemas/                  # API schemas
│
├── docs/                         # Documentation
├── tests/                        # Unit tests
├── scripts/                      # Build scripts
├── build_windows.spec            # PyInstaller
├── installer.nsi                 # NSIS installer
└── docker-compose.yml            # Docker config
```

---

# 🔧 Technical Stack

| Component | Technology | Version |
|-----------|------------|---------|
| UI Framework | PyQt6 | 6.4+ |
| Backend API | FastAPI | 0.100+ |
| Database | SQLite | 3.x |
| Docking | AutoDock Vina | 1.2.5 |
| ML Docking | Gnina | 1.0 |
| Chemistry | RDKit | 2023.03+ |
| Analysis | ODDT | 0.7+ |
| Containers | Docker | 24.x |
| Packaging | PyInstaller | 6.x |
| Installer | NSIS | 3.x |

---

# 🚀 Installation

## Windows Installer
```
1. Download BioDockify-Setup-1.1.0.exe
2. Run as Administrator
3. Install Docker Desktop when prompted
4. Launch from Start Menu
```

## From Source
```bash
git clone https://github.com/tajo9128/Docking-studio.git
cd Docking-studio
pip install -r requirements.txt
python src/biodockify_main.py
```

---

# 📊 Workflow

```
[Upload] → [Configure] → [Dock] → [Analyze] → [Export]
   │           │           │          │           │
   ▼           ▼           ▼          ▼           ▼
 PDB/MOL2   Grid Box    Vina/Gnina  Interactions  PDF/CSV
 SDF/SMILES  Params      Docker     Visualization Reports
```

---

*BioDockify Docking Studio v1.1.0 - Unified Molecular Docking Platform*
