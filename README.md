# BioDockify Docking Studio

![Version](https://img.shields.io/badge/version-1.0.50-blue)
![License](https://img.shields.io/badge/license-Apache%202.0-green)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS-lightgrey)

**BioDockify Docking Studio** is a production-grade, desktop-based solution for molecular docking. By combining a modern PyQt6 interface with the robustness of containerized AutoDock Vina engines, it provides a seamless experience for computational chemists and drug discovery researchers.

![Screenshot](https://via.placeholder.com/800x450?text=BioDockify+Docking+Studio+v1.0.0)

## 🚀 Features

- **Intuitive GUI**: Drag-and-drop file inputs, 3D box configuration, and real-time progress monitoring.
- **Containerized Backend**: reproducible science using Docker-isolated environments.
- **Agent Zero™**: Self-healing AI that detects failures (timeouts, missed poses) and automatically recovers.
- **Advanced Analysis**: Integrated RDKit chemical descriptors and ODDT interaction scoring.

## 📦 Installation

### Windows 11
Download the [Setup Installer](https://github.com/tajo9128/Docking-studio/releases/download/v1.0.50/BioDockify-Docking-Studio.exe).

### macOS
Download the Disk Image for [Intel](https://github.com/tajo9128/Docking-studio/releases/download/v1.0.50/BioDockify-1.0.0-macos-intel.dmg) or [Apple Silicon](https://github.com/tajo9128/Docking-studio/releases/download/v1.0.50/BioDockify-1.0.0-macos-arm64.dmg).

*Requires Docker Desktop.*

## 📚 Documentation

- [Installation Guide](docs/installation.md)
- [User Guide](docs/user_guide.md)
- [Troubleshooting](docs/troubleshooting.md)
- [FAQ](docs/faq.md)

## 🏗️ Development

### Prerequisites
- Python 3.9+
- Docker Desktop

### Setup
```bash
git clone https://github.com/tajo9128/Docking-studio.git
cd Docking-studio
pip install -r requirements.txt
python src/main.py
```

## 📄 License
Released under the [Apache 2.0 License](LICENSE).

## ✍️ Citation
If you use BioDockify in your research, please cite:
> BioDockify Team. (2025). BioDockify Docking Studio (Version 1.0.0) [Computer software]. https://github.com/tajo9128/Docking-studio
