# Changelog

All notable changes to this project will be documented in this file.

## [1.0.3] - 2025-12-31

### Fixed
- **Build**: Switched to `collect_all('PyQt6')` in PyInstaller spec to force inclusion of all Qt6 binaries and plugins. This resolves persistent startup crashes.

## [1.0.2] - 2025-12-31

### Added
- **Automation**: Enabled fully automated CI/CD pipeline for Windows builds using GitHub Actions.

## [1.0.2] - 2025-12-31

### Added
- **Automation**: Enabled fully automated CI/CD pipeline for Windows builds using GitHub Actions.

## [1.0.1] - 2025-12-31

### Fixed
- **Critical**: Resolved `ModuleNotFoundError: No module named 'PyQt6'` by explicitly including UI modules in Windows build specification.
- Updated documentation links to point to the correct executable name.

## [1.0.0] - 2025-12-30

### Added
- Initial production release of BioDockify Docking Studio.
- **UI**: Modern Dark Theme PyQt6 interface with real-time job monitoring.
- **Backend**: Containerized AutoDock Vina engine managed via Docker SDK.
- **AI**: "Agent Zero" self-healing system for automatic error recovery (timeouts, bad poses).
- **Security**: Robust path traversal protection and input sanitization.
- **Data**: SQLite database with transactional integrity and checkpointing.

### Fixed
- Addressed Docker container leakage ensuring clean shutdown of worker containers.
- Fixed database connection leaks on Windows platforms.
- Resolved potential path traversal vulnerabilities in file upload handlers.

### Changed
- Optimized Docker image size using multi-stage builds.
- Improved error handling with distinct exception types for file and permission errors.
