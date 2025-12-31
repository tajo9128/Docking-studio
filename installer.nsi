; BioDockify Docking Studio - Complete NSIS Installer
; Full setup with Docker Desktop dependency checking

!define APP_NAME "BioDockify Docking Studio"
!define APP_SHORT_NAME "BioDockify"
!define APP_VERSION "1.0.43"
!define APP_PUBLISHER "BioDockify Team"
!define APP_WEBSITE "https://github.com/tajo9128/Docking-studio"
!define DOCKER_URL "https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe"

; Installer Configuration
Name "${APP_NAME} ${APP_VERSION}"
OutFile "BioDockify-Setup-${APP_VERSION}.exe"
InstallDir "$PROGRAMFILES64\${APP_SHORT_NAME}"
InstallDirRegKey HKLM "Software\${APP_SHORT_NAME}" "InstallPath"
RequestExecutionLevel admin
SetCompressor /SOLID lzma
CRCCheck on

; Variables
Var DockerInstalled
Var DockerRunning
Var DockerVersion
Var DockerInstallPath
Var InstallDockerChoice
Var DockerDownloadPath

; Interface Settings
!include "MUI2.nsh"
!include "LogicLib.nsh"

; MUI Settings / Custom Branding
!define MUI_ABORTWARNING
!define MUI_ICON "biodockify_icon.ico"
!define MUI_UNICON "biodockify_icon.ico"

; Header Image (Top right of inner pages)
!define MUI_HEADERIMAGE
!define MUI_HEADERIMAGE_BITMAP "header.bmp" ; Relative to where makensis runs (dist/installer)
!define MUI_HEADERIMAGE_RIGHT

; Welcome/Finish Image (Left side)
!define MUI_WELCOMEFINISHPAGE_BITMAP "welcome.bmp"

; Custom Branding Text
BrandingText "BioDockify Docking Studio - Enterprise Edition"

!define MUI_COMPONENTSPAGE_NODESC

; Pages
!define MUI_WELCOMEPAGE_TITLE "Welcome to the BioDockify Setup Wizard"
!define MUI_WELCOMEPAGE_TEXT "This assistant will guide you through the installation of $(^Name).$\r$\n$\r$\nBioDockify is the industry standard for cloud-native molecular docking.$\r$\n$\r$\nIt is recommended that you close all other applications before starting Setup. This will make it possible to update relevant system files without having to reboot your computer.$\r$\n$\r$\n$_CLICK"

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "LICENSE"
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY

; Docker Dependency Check Page
Page custom DockerPageCreate DockerPageLeave

!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

; Languages
!insertmacro MUI_LANGUAGE "English"

; License Data
LicenseData "LICENSE"

; Icons
UninstallIcon "biodockify_icon.ico"

; Registry entries for Add/Remove Programs
!define REG_UNINSTALL "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}"

; ============================================================================
; INITIALIZATION - Check Docker Status
; ============================================================================
Function .onInit
    ; Clear errors
    ClearErrors
    
    ; Check if Docker is installed via registry
    ReadRegStr $DockerInstallPath HKLM "Software\Docker\Docker Desktop" "InstallLocation"
    ${If} ${Errors}
        ; Try alternate registry key if primary fails
        ReadRegStr $DockerVersion HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Docker Desktop" "DisplayVersion"
        ${If} ${Errors}
            StrCpy $DockerInstalled "0"
            StrCpy $DockerVersion "Not Installed"
            DetailPrint "Docker Desktop NOT found"
        ${Else}
            ReadRegStr $DockerInstallPath HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Docker Desktop" "InstallLocation"
            StrCpy $DockerInstalled "1"
            DetailPrint "Docker Desktop $DockerVersion found (alternate key)"
        ${EndIf}
    ${Else}
        StrCpy $DockerInstalled "1"
        ; Try to get Docker version
        ReadRegStr $DockerVersion HKLM "Software\Docker\Docker Desktop" "Version"
        ${If} ${Errors}
             DetailPrint "Docker Desktop installed at: $DockerInstallPath"
             StrCpy $DockerVersion "Unknown"
        ${Else}
             DetailPrint "Docker Desktop $DockerVersion found at: $DockerInstallPath"
        ${EndIf}
    ${EndIf}
    
    ; Check if Docker is running
    nsExec::ExecToLog 'docker ps' $R0
    Pop $R0
    ${If} $R0 == "0"
        StrCpy $DockerRunning "1"
        DetailPrint "Docker daemon is RUNNING"
    ${Else}
        StrCpy $DockerRunning "0"
        DetailPrint "Docker daemon is NOT running"
    ${EndIf}
FunctionEnd

; ============================================================================
; DOCKER DEPENDENCY PAGE
; ============================================================================
Function DockerPageCreate
    !insertmacro MUI_HEADER_TEXT "System Requirements Check"
    !insertmacro MUI_TEXT_TOP "Checking required dependencies for BioDockify Docking Studio..."
    
    ; Create custom page
    nsDialogs::Create 1018 -1 -1 -1 -1
    Pop $R0
    
    ; Header - Docker Status
    ${If} $DockerInstalled == "1"
        ${If} $DockerRunning == "1"
            ${NSD_CreateLabel} 1018 100u 30u 90% 12u "✅ Docker Desktop is INSTALLED and RUNNING"
            ${NSD_CreateLabel} 1018 100u 45u 90% 12u "Version: $DockerVersion"
            ${NSD_CreateLabel} 1018 100u 60u 90% 12u "Location: $DockerInstallPath"
            ${NSD_CreateLabel} 1018 100u 80u 90% 10u "All dependencies satisfied. You can proceed with installation."
            
            ; Show success image
            ${NSD_CreateIcon} 1018 0u 30u 40u 105 "SHELL32.dll" 0
            
            ; Allow to continue
            GetDlgItem $0 1 $R1
            EnableWindow $R1 1
        ${Else}
            ${NSD_CreateLabel} 1018 100u 30u 90% 12u "⚠️ Docker Desktop is INSTALLED but NOT RUNNING"
            ${NSD_CreateLabel} 1018 100u 45u 90% 12u "Version: $DockerVersion"
            ${NSD_CreateLabel} 1018 100u 65u 90% 10u "Docker Desktop must be started before using BioDockify."
            ${NSD_CreateLabel} 1018 100u 80u 90% 10u "Please start Docker Desktop from your Start Menu after installation."
            
            ; Show warning image
            ${NSD_CreateIcon} 1018 0u 30u 40u 101 "SHELL32.dll" 0
            
            ; Allow to continue with warning
            GetDlgItem $0 1 $R1
            EnableWindow $R1 1
        ${EndIf}
    ${Else}
        ${NSD_CreateLabel} 1018 100u 30u 90% 12u "❌ Docker Desktop is NOT INSTALLED"
        ${NSD_CreateLabel} 1018 100u 50u 90% 12u "BioDockify requires Docker Desktop to run molecular docking simulations."
        ${NSD_CreateLabel} 1018 100u 65u 90% 10u "Please install Docker Desktop before using BioDockify."
        
        ; Show error image
        ${NSD_CreateIcon} 1018 0u 30u 40u 106 "SHELL32.dll" 0
        
        ; Action buttons
        ${NSD_CreateButton} 1018 100u 100u 160u 25u "Download Docker Desktop" DownloadDocker
        ${NSD_CreateButton} 1018 270u 100u 160u 25u "Install from Disk" InstallDockerDisk
        ${NSD_CreateButton} 1018 100u 135u 160u 25u "Skip & Continue" SkipDocker
        
        ; Info text
        ${NSD_CreateLabel} 1018 100u 175u 90% 10u "⬇️ Click below to download and install Docker Desktop"
        ${NSD_CreateLabel} 1018 100u 190u 90% 10u "Docker will be installed before BioDockify."
        
        ; Disable Next button until action
        GetDlgItem $0 1 $R1
        EnableWindow $R1 0
    ${EndIf}
    
    ; Show the page
    nsDialogs::Show
    
    ; Destroy dialog
    ${NSD_DestroyDialog} 1018
FunctionEnd

; ============================================================================
; DOCKER PAGE ACTIONS
; ============================================================================

; Download Docker Desktop
Function DownloadDocker
    DetailPrint "Starting Docker Desktop download..."
    
    ; Set download path to temp directory
    StrCpy $DockerDownloadPath "$TEMP\DockerDesktopInstaller.exe"
    
    ; Disable buttons during download
    GetDlgItem $0 1 $R1
    EnableWindow $R1 0
    
    ; Update button text
    SendMessage $R0 ${WM_SETTEXT} 0 "STR:Downloading..."
    
    ; Download with progress
    inetc::get /CANCEL \
        "${DOCKER_URL}" \
        "$DockerDownloadPath" \
        /END
    
    Pop $R0
    ${If} $R0 == "OK"
        DetailPrint "Docker Desktop downloaded successfully to: $DockerDownloadPath"
        
        ; Offer to install immediately
        MessageBox MB_YESNO|MB_ICONQUESTION \
            "Docker Desktop has been downloaded.$\n$\nDo you want to install it now?" \
            IDYES InstallDownloadedDocker \
            IDNONE \
            "Download Complete"
    ${Else}
        MessageBox MB_OK|MB_ICONEXCLAMATION \
            "Failed to download Docker Desktop.$\n$\nPlease download it manually from:$\n$\n${DOCKER_URL}" \
            "Download Failed"
        
        ; Re-enable Next button to allow skipping
        GetDlgItem $0 1 $R1
        EnableWindow $R1 1
    ${EndIf}
FunctionEnd

; Install from Disk
Function InstallDockerDisk
    ; Open file dialog to select Docker installer
    nsDialogs::OpenFileDialog \
        "Select Docker Desktop Installer" \
        "$DESKTOP" \
        "Docker Installer (*.exe)|*.exe|All Files (*.*)|*.*" \
        "" \
        "" \
        $DockerInstallPath
    
    Pop $R0
    ${If} $R0 != "error"
        ${If} $DockerInstallPath != ""
            StrCpy $InstallDockerChoice "1"
            
            ; Enable Next button
            GetDlgItem $0 1 $R1
            EnableWindow $R1 1
            
            MessageBox MB_OK|MB_ICONINFORMATION \
                "Docker Desktop installer selected:$\n$DockerInstallPath$\n$\nClick Next to install it before BioDockify." \
                "Docker Installer Selected"
        ${EndIf}
    ${EndIf}
FunctionEnd

; Skip Docker Installation
Function SkipDocker
    MessageBox MB_YESNO|MB_ICONWARNING \
        "You chose to skip Docker Desktop installation.$\n$\nBioDockify will not work without Docker Desktop.$\n$\nAre you sure you want to continue?" \
        IDYES AllowSkip \
        IDNONE \
        "Warning"
    
    AllowSkip:
    ; Enable Next button
    GetDlgItem $0 1 $R1
    EnableWindow $R1 1
FunctionEnd

; Install downloaded Docker
Function InstallDownloadedDocker
    DetailPrint "Installing Docker Desktop from: $DockerDownloadPath"
    
    MessageBox MB_OK|MB_ICONINFORMATION \
        "BioDockify installer will now launch Docker Desktop setup.$\n$\nPlease complete the Docker installation and return to this installer.$\n$\nClick OK to continue." \
        "Launch Docker Setup"
    
    ; Execute Docker installer
    ExecWait '"$DockerDownloadPath"' $R0
    
    ${If} $R0 == "0"
        DetailPrint "Docker Desktop installation completed successfully"
        StrCpy $DockerInstalled "1"
        StrCpy $DockerInstallPath "$PROGRAMFILES\Docker\Docker"
    ${Else}
        DetailPrint "Docker Desktop installation failed with code: $R0"
        MessageBox MB_OK|MB_ICONEXCLAMATION \
            "Docker Desktop installation may not have completed successfully.$\n$\nError code: $R0" \
            "Installation Warning"
        
        ; Allow user to skip
        MessageBox MB_YESNO|MB_ICONWARNING \
            "Docker installation may have failed. Do you still want to continue?" \
            IDYES AllowSkip \
            IDNONE \
            "Warning"
    ${EndIf}
    
    AllowSkip:
    ; Enable Next button
    GetDlgItem $0 1 $R1
    EnableWindow $R1 1
FunctionEnd

; Docker page leave handler
Function DockerPageLeave
    ${If} $DockerInstalled == "0"
        ${If} $InstallDockerChoice == "1"
            ; User selected custom installer
            DetailPrint "Will install Docker from: $DockerInstallPath"
            ExecWait '"$DockerInstallPath"' /SILENT $R0
            
            ${If} $R0 == "0"
                DetailPrint "Docker Desktop installed successfully"
                StrCpy $DockerInstalled "1"
            ${Else}
                DetailPrint "Docker Desktop installation returned: $R0"
                MessageBox MB_OK|MB_ICONEXCLAMATION \
                    "Docker Desktop installation did not complete successfully.$\n$\nYou can install it manually later." \
                    "Installation Warning"
            ${EndIf}
        ${EndIf}
    ${Else}
        DetailPrint "Docker dependency satisfied, proceeding with installation..."
    ${EndIf}
FunctionEnd

; ============================================================================
; INSTALLATION SECTIONS
; ============================================================================

Section "BioDockify Application" SEC01
    SectionIn RO
    SetOutPath "$INSTDIR"
    
    ; Create installation directory
    CreateDirectory "$INSTDIR"
    
    ; Copy main executable (built by PyInstaller)
    File "BioDockify-Docking-Studio.exe"
    
    ; Copy templates directory
    CreateDirectory "$INSTDIR\templates"
    File /r "templates\*.*"
    
    ; Copy UI styles
    CreateDirectory "$INSTDIR\ui\styles"
    File /r "ui\styles\*.*"
    
    ; Copy icon
    File "biodockify_icon.ico"
    
    ; Copy license
    File "LICENSE"
    
    ; Copy README
    File "README.md"
    
    ; Create Start Menu shortcuts
    CreateDirectory "$SMPROGRAMS\${APP_SHORT_NAME}"
    CreateShortCut "$SMPROGRAMS\${APP_SHORT_NAME}\${APP_NAME}.lnk" \
        "$INSTDIR\BioDockify-Docking-Studio.exe" \
        "" \
        "" \
        "" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "Launch BioDockify Docking Studio"
    
    CreateShortCut "$SMPROGRAMS\${APP_SHORT_NAME}\Docker Desktop.lnk" \
        "$PROGRAMFILES64\Docker\Docker\Docker Desktop.exe" \
        "" \
        "" \
        "" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "Launch Docker Desktop"
    
    ; Create Desktop shortcut
    CreateShortCut "$DESKTOP\${APP_NAME}.lnk" \
        "$INSTDIR\BioDockify-Docking-Studio.exe" \
        "" \
        "" \
        "" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "BioDockify"
    
    ; Write registry entries
    WriteRegStr HKLM "Software\${APP_SHORT_NAME}" "InstallPath" "$INSTDIR"
    WriteRegStr HKLM "Software\${APP_SHORT_NAME}" "Version" "${APP_VERSION}"
    WriteRegStr HKLM "Software\${APP_SHORT_NAME}" "DockerRequired" "1"
    
    ; Add/Remove Programs entry
    WriteRegStr HKLM "${REG_UNINSTALL}" "DisplayName" "${APP_NAME}"
    WriteRegStr HKLM "${REG_UNINSTALL}" "DisplayVersion" "${APP_VERSION}"
    WriteRegStr HKLM "${REG_UNINSTALL}" "Publisher" "${APP_PUBLISHER}"
    WriteRegStr HKLM "${REG_UNINSTALL}" "UninstallString" "$INSTDIR\uninstall.exe"
    WriteRegStr HKLM "${REG_UNINSTALL}" "InstallLocation" "$INSTDIR"
    WriteRegStr HKLM "${REG_UNINSTALL}" "URLInfoAbout" "${APP_WEBSITE}"
    WriteRegDWORD HKLM "${REG_UNINSTALL}" "NoModify" 1
    WriteRegDWORD HKLM "${REG_UNINSTALL}" "NoRepair" 1
    WriteRegDWORD HKLM "${REG_UNINSTALL}" "EstimatedSize" "150000" ; 150 MB
    WriteRegStr HKLM "${REG_UNINSTALL}" "DisplayIcon" "$INSTDIR\BioDockify-Docking-Studio.exe"
    
    ; Write uninstaller
    WriteUninstaller "$INSTDIR\uninstall.exe"
    
    DetailPrint "BioDockify application files installed successfully"
SectionEnd

Section "Start Menu Shortcuts" SEC02
    SectionIn 1  ; Checked by default
    SetOutPath "$INSTDIR"
    
    ; Already created in main section, just ensure it exists
    CreateDirectory "$SMPROGRAMS\${APP_SHORT_NAME}"
    
    ; Create uninstall shortcut
    CreateShortCut "$SMPROGRAMS\${APP_SHORT_NAME}\Uninstall ${APP_NAME}.lnk" \
        "$INSTDIR\uninstall.exe" \
        "" \
        "" \
        "" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "Uninstall BioDockify"
    
    DetailPrint "Start menu shortcuts created"
SectionEnd

Section "Desktop Shortcut" SEC03
    SectionIn 1
    SetOutPath "$INSTDIR"
    
    ; Already created in main section
    DetailPrint "Desktop shortcut created"
SectionEnd

; ============================================================================
; POST-INSTALLATION ACTIONS
; ============================================================================
Section -PostInstall
    ; Check Docker status after installation
    ${If} $DockerInstalled == "1"
        ${If} $DockerRunning == "0"
            MessageBox MB_OK|MB_ICONINFORMATION \
                "Docker Desktop is installed but not running.$\n$\nPlease start Docker Desktop from your Start Menu before launching BioDockify." \
                "Docker Desktop Not Running"
        ${EndIf}
    ${Else}
        MessageBox MB_OK|MB_ICONWARNING \
            "Docker Desktop is not installed.$\n$\nBioDockify requires Docker Desktop to function.$\n$\nPlease install Docker Desktop from: ${DOCKER_URL}" \
            "Docker Desktop Required"
    ${EndIf}
    
    ; Offer to launch application
    SetShellVarContext all
    SetOverwrite on
    CreateShortCut "$SMSTARTUP\BioDockify-Docking-Studio.lnk" \
        "$INSTDIR\BioDockify-Docking-Studio.exe" \
        "" \
        "" \
        "" \
        0 \
        SW_SHOWNORMAL \
        "" \
        "Launch at startup"
    
    ; Ask if user wants to start with Windows
    MessageBox MB_YESNO|MB_ICONQUESTION \
        "Do you want BioDockify to start automatically when Windows starts?" \
        IDYES AddToStartup \
        IDNONE \
        "Auto-Start"
    
    AddToStartup:
    ; Shortcut already created above
    
    ; Ask to launch now
    MessageBox MB_YESNO|MB_ICONQUESTION \
        "BioDockify Docking Studio has been installed.$\n$\nDo you want to launch it now?" \
        IDYES LaunchApp \
        IDNONE \
        "Installation Complete"
    
    LaunchApp:
    Exec '"$INSTDIR\BioDockify-Docking-Studio.exe"'
SectionEnd

; ============================================================================
; UNINSTALLATION
; ============================================================================
Section "Uninstall"
    ; Close BioDockify if running
    nsExec::ExecToLog 'taskkill /F /IM BioDockify-Docking-Studio.exe' $R0
    Pop $R0
    
    ; Remove files
    Delete "$INSTDIR\BioDockify-Docking-Studio.exe"
    Delete "$INSTDIR\uninstall.exe"
    Delete "$INSTDIR\biodockify_icon.ico"
    Delete "$INSTDIR\LICENSE"
    Delete "$INSTDIR\README.md"
    
    RMDir /r "$INSTDIR\templates"
    RMDir /r "$INSTDIR\ui\styles"
    RMDir "$INSTDIR\ui"
    
    ; Remove startup shortcut
    Delete "$SMSTARTUP\BioDockify-Docking-Studio.lnk"
    
    ; Remove shortcuts
    Delete "$SMPROGRAMS\${APP_SHORT_NAME}\*.lnk"
    RMDir "$SMPROGRAMS\${APP_SHORT_NAME}"
    
    ; Remove desktop shortcut
    Delete "$DESKTOP\${APP_NAME}.lnk"
    
    ; Remove installation directory if empty
    RMDir "$INSTDIR"
    
    ; Remove registry entries
    DeleteRegKey HKLM "Software\${APP_SHORT_NAME}"
    DeleteRegKey HKLM "${REG_UNINSTALL}"
    
    ; Clean up Docker directory if we installed it (optional - commented out)
    ; Delete "$PROGRAMFILES64\Docker\Docker"
    
    DetailPrint "BioDockify has been uninstalled"
SectionEnd
