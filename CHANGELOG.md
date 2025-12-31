# Changelog

All notable changes to this project will be documented in this file.

## [1.0.31] - 2025-12-31

### Fixed
- **Build**: Fixed `SyntaxError: f-string expression part cannot include a backslash` in `src/biodockify_main.py`. This was caused by using a backslash inside an f-string expression (`{"\n".join(...)}`), which is invalid illegal in Python 3.11. The logic was refactored to pre-calculate the string key variable.

## [1.0.30] - 2025-12-31

### Fixed
- **Runtime**: Added import fallback logic to `src/biodockify_main.py`. The application now attempts to import `src.ui.main_window` if the standard `ui.main_window` import fails. This addresses potential package structure flattening or nesting issues within the PyInstaller bundle.
- **Debug**: Enhanced the error `MessageBox` to list the actual contents of the frozen bundle directory (`sys._MEIPASS`). This will transparently show whether the `ui` folder exists at the root or under a `src` subdirectory, allowing for rapid diagnosis of packaging errors.

## [1.0.29] - 2025-12-31

### Fixed
- **CI/CD**: Fixed `ERROR: script not found` by explicitly adding `src/biodockify_main.py` to the git repository. The previous rename operation resulted in the new file being untracked, causing it to be missing from the CI checkout.

## [1.0.28] - 2025-12-31

### Fixed
- **Build**: Renamed entry point from `src/main.py` to `src/biodockify_main.py`. The persistence of the error on line 38 (despite code changes in v1.0.27 moving that line) indicated that PyInstaller or the CI pipeline was using a stale/cached version of the script. Renaming the file forces a fresh analysis and bundling of the application logic.

## [1.0.27] - 2025-12-31

### Fixed
- **Debug**: Wrapped `ui` imports in `src/main.py` with a `try...except` block to catch `ImportError` and display the native Windows debug `MessageBox` (showing `sys.path` and CWD). This will provide definitive diagnostics if the module load fails again.
- **Build**: Reverted `pathex` in `build_windows.spec` to `['.', 'src']` (relative path). This configuration was previously proven to work in v1.0.20 (where the app successfully loaded `ui` modules before hitting a code typo), suggesting the absolute path change in v1.0.26 might have been counterproductive.

## [1.0.26] - 2025-12-31

### Fixed
- **Build**: Resolved persistent `ModuleNotFoundError` for `ui` modules by implementing a comprehensive fix in `build_windows.spec`. Updated `pathex` to use `os.path.abspath('src')` to guarantee the source directory is in the analysis path. Added both standard (`ui.*`) and fully qualified (`src.ui.*`) module names to `hidden_imports` to cover all potential import resolution strategies by PyInstaller.

## [1.0.25] - 2025-12-31

### Fixed
- **Build**: Fixed `SyntaxError: '(' was never closed` in `build_windows.spec`. This error was caused by a malformed code block and code duplication introduced in v1.0.24. Also restored the missing `tmp_ret = collect_all('PyQt6')` definition, ensuring that `datas` and `binaries` are correctly populated.

## [1.0.24] - 2025-12-31

### Fixed
- **Build**: Reverted to explicit `hidden_imports` list in `build_windows.spec`. The dynamic glob strategy in v1.0.22 failed to prevent `ModuleNotFoundError: No module named 'ui.main_window'`. By explicitly listing all UI widgets (`ui.main_window`, `ui.upload_widget`, etc.), we ensure PyInstaller bundles them correctly.

## [1.0.23] - 2025-12-31

### Fixed
- **Runtime**: Fixed `ModuleNotFoundError: No module named 'ui.main_window'` by emptying `src/ui/__init__.py`. The `__init__` file was attempting to bundle exports using relative imports (`from .main_window import ...`), which was failing in the frozen environment. Since `main.py` imports `MainWindow` explicitly from `ui.main_window`, these package-level exports were redundant and causing the crash.

## [1.0.22] - 2025-12-31

### Fixed
- **Build**: Fixed persistent `ModuleNotFoundError` for `ui` submodules by implementing dynamic discovery in `build_windows.spec`. Instead of manually listing files, the script now globs `src/ui/*.py` and adds all of them to `hidden_imports`. This ensures `upload_widget`, `configuration_widget`, and others are not missing.
- **Build**: Updated `pathex` to use `os.path.abspath('src')` for unambiguous path resolution.

## [1.0.21] - 2025-12-31

### Fixed
- **Code**: Fixed `ImportError: cannot import name 'QDockWidget' from 'PyQt6.QtGui'`. Correctly moved `QDockWidget` import to `PyQt6.QtWidgets` in `src/ui/main_window.py`. This runtime error was revealed after the build pipeline was successfully fixed in v1.0.20.

## [1.0.20] - 2025-12-31

### Fixed
- **Build**: Fixed `NameError: name 'tmp_ret' is not defined`. Explicitly imported `collect_all` and called `tmp_ret = collect_all('PyQt6')` before unpacking variables. This completes the `collect_all` logic intended in v1.0.18.

## [1.0.19] - 2025-12-31

### Fixed
- **Build**: Fixed `ModuleNotFoundError: No module named 'ui.main_window'`. Added `src` to PyInstaller's `pathex` and explicitly listed `ui` submodules in `hidden_imports`. This ensures the application can resolve its own internal modules when bundled.

## [1.0.18] - 2025-12-31

### Fixed
- **Build**: Switched purely to `collect_all('PyQt6')` strategy. Now that the CI environment is guaranteed to have `PyQt6` installed (via v1.0.17 fix), this standard hook should correctly gather all split packages, DLLs, and plugins that manual copying missed.

## [1.0.17] - 2025-12-31

### Fixed
- **Build**: Fixed `NameError: pyqt_path` by properly initializing variable at the start of spec.
- **CI/CD**: Added explicit `pip install PyQt6 --force-reinstall` and `pip list` debug info to GitHub Actions workflow. This addresses the root cause where `PyQt6` appeared to be missing from the build environment.

## [1.0.16] - 2025-12-31

### Fixed
- **Build**: Fixed critical bug in v1.0.15 where `import PyQt6` was outside the `try/except` block in the spec file, causing the build to crash. Moved it inside the safety block to allow fallback detection.

## [1.0.15] - 2025-12-31

### Fixed
- **Build**: Switched `PyQt6` copy logic in `build_windows.spec` to use `import PyQt6; os.path.dirname(PyQt6.__file__)`. This matches exactly what the Python interpreter sees, guaranteeing the correct folder is copied.

## [1.0.13] - 2025-12-31

### Fixed
- **Build**: Fixed `NameError: name 'block_cipher' is not defined` by restoring `block_cipher = None` to the top of `build_windows.spec`.

## [1.0.12] - 2025-12-31

### Fixed
- **Build**: Reverted `import PyQt6.QtCore` from spec file to fix build-time crash. Restored clean `hidden_imports` list (Option A) but explicitly added `PyQt6.sip` and ensured `pathex` includes root.

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
