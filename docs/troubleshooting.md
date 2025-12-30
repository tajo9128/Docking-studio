# Troubleshooting Guide

## Common Issues

### 1. "Docker is not connected"
**Symptoms**: Status bar shows "Docker: Disconnected" or Agent Zero reports `DOCKER_NOT_RUNNING`.
**Solution**:
- Open Docker Desktop.
- Wait for the engine to start (green icon).
- Restart BioDockify.

### 2. "Execution failed: No poses found"
**Symptoms**: Job fails with `DOCKING_NO_POSES`.
**Solution**:
- The search box might be misplaced. Verify X/Y/Z coordinates.
- box size might be too small. Try increasing dimensions by 5A.
- Agent Zero usually attempts to fix this by expanding the box automatically.

### 3. "Permission Denied" on File Upload
**Symptoms**: Error uploading receptor/ligand.
**Solution**:
- Ensure the files are not open in another program (e.g., PyMOL).
- Check file permissions.
- Move files to a simple path (e.g., Desktop) to avoid path length issues.

### 4. "Database Locked"
**Symptoms**: Application freezes or logs DB errors.
**Solution**:
- Restart the application.
- (Advanced) Delete `BioDockify.db` to reset history (Data loss warning).

## Logs
Log files are stored in:
- **Windows**: `%LOCALAPPDATA%\BioDockify\logs\`
- **macOS**: `~/.config/BioDockify/logs/`
