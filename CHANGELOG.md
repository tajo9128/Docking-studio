# Changelog

All notable changes to this project will be documented in this file.

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
