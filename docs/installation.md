# Installation Guide - BioDockify Docking Studio v1.0.0

## System Requirements

### Windows
- **OS**: Windows 11 or Windows 10 (64-bit)
- **RAM**: Minimum 8GB (16GB recommended)
- **Storage**: 2GB free space
- **Dependencies**: Docker Desktop for Windows

### macOS
- **OS**: macOS Sonoma (14.0) or newer
- **Chip**: Intel or Apple Silicon (M1/M2/M3)
- **RAM**: Minimum 8GB (16GB recommended)
- **Dependencies**: Docker Desktop for Mac

---

## Installation Steps

### Windows 11 (x64)
1. Download `BioDockify-Docking-Studio.exe` from the [Latest Release](https://github.com/tajo9128/Docking-studio/releases/download/v1.0.23/BioDockify-Docking-Studio.exe).
2. Double-click the installer to start the setup wizard.
3. Follow the on-screen instructions to complete the installation.
4. Ensure Docker Desktop is running before launching the application.

> [!NOTE]
> **Microsoft SmartScreen Warning**: Since this is a new release, Windows may display a blue window saying *"Windows protected your PC"*. This is normal for new software.
> Click **"More info"** and then **"Run anyway"** to install.

### macOS (Intel & Apple Silicon)
1. Download the appropriate `.dmg` file for your architecture:
   - Intel: `BioDockify-1.0.0-macos-intel.dmg`
   - Apple Silicon: `BioDockify-1.0.0-macos-arm64.dmg`
2. Open the disk image.
3. Drag `BioDockify.app` to your `Applications` folder.
4. **First Launch**: You may need to right-click the app and select "Open" to bypass security checks if the app is not signed with an Apple Developer ID.

---

## Post-Installation Verification
1. Launch BioDockify Docking Studio.
2. The status bar at the bottom should show "Docker: Connected".
3. If valid, you are ready to start docking!
