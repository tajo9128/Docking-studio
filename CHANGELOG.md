# Changelog

All notable changes to this project will be documented in this file.

## [1.0.11] - 2025-12-31

### Fixed
- **Build**: Implemented "Force Copy" strategy. The build script now inspects `PyQt6.QtCore.__file__` to find the exact installation location on the runner and explicitly copies the entire `PyQt6` folder into the executable's distribution. This bypasses PyInstaller's discovery mechanism completely.

## [1.0.10] - 2025-12-31

### Fixed
- **Build**: Implemented `collect_submodules('PyQt6')` strategy. This dynamically finds all installed Qt modules in the build environment and forces them into `hidden_imports`, solving the namespace discovery issue. Updated `pathex` to `['.']`.

## [1.0.9] - 2025-12-31

### Fixed
- **Build**: completely rewrote `build_windows.spec` to use a clean, explicit `hidden_imports` list. Removed all fragile `collect_all` logic. This is the "Clean & Professional" configuration recommended for stability.

## [1.0.8] - 2025-12-31

### Fixed
- **Build**: Implemented "Aggressive Collection" strategy in build spec. Explicitly running `collect_all` on `PyQt6`, `QtWidgets`, `QtCore`, and `QtGui` individually to force-bundle missing namespace components.

## [1.0.7] - 2025-12-31

### Fixed
- **Build**: Fixed `NameError` in build spec by correctly initializing empty `datas` and `binaries` lists before use.

## [1.0.6] - 2025-12-31

### Fixed
- **Build**: Reverted manual `collect_all` hooks which caused namespace package errors. Relying on upgraded `PyInstaller 6.17.0` built-in hooks for clean PyQt6 bundling.

## [1.0.5] - 2025-12-31

### Fixed
- **Critical Build Fix**: Upgraded `PyInstaller` to v6.3.0+ and `PyQt6` to v6.6.1+. Added `copy_metadata` to build spec. This addresses the stubborn `ModuleNotFoundError` by using modern packaging hooks.

## [1.0.4] - 2025-12-31

### Added
- **UI**: Added custom professional application icon (`icon.ico`) to the executable.

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
