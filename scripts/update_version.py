import sys
import re
from pathlib import Path

def update_file(file_path, patterns, new_version):
    path = Path(file_path)
    if not path.exists():
        print(f"Warning: {file_path} not found.")
        return

    content = path.read_text(encoding='utf-8')
    original_content = content
    
    for pattern, replacement in patterns:
        # Replace using regex, standardizing on the new version
        # replacement format string expects {version}
        final_replacement = replacement.format(version=new_version)
        content = re.sub(pattern, final_replacement, content)
        
    if content != original_content:
        path.write_text(content, encoding='utf-8')
        print(f"Updated {file_path}")
    else:
        print(f"No changes needed for {file_path}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python scripts/update_version.py <new_version>")
        print("Example: python scripts/update_version.py 1.0.3")
        sys.exit(1)

    new_version = sys.argv[1]
    
    # Define update rules
    project_root = Path(__file__).parent.parent
    
    # 1. VERSION file
    (project_root / "VERSION").write_text(new_version, encoding='utf-8')
    print(f"Updated VERSION to {new_version}")

    # 2. README.md
    # Updates badges and download links
    readme_patterns = [
        (r"badge/version-[\d\.]+-blue", "badge/version-{version}-blue"),
        (r"/releases/download/v[\d\.]+/BioDockify", "/releases/download/v{version}/BioDockify")
    ]
    update_file(project_root / "README.md", readme_patterns, new_version)

    # 3. docs/installation.md
    # Updates download links
    install_patterns = [
        (r"/releases/download/v[\d\.]+/BioDockify", "/releases/download/v{version}/BioDockify")
    ]
    update_file(project_root / "docs/installation.md", install_patterns, new_version)

    # 4. CITATION.cff
    # Updates version field
    citation_patterns = [
        (r"version: [\d\.]+", "version: {version}")
    ]
    update_file(project_root / "CITATION.cff", citation_patterns, new_version)
    
    # 5. installer.nsi
    # Updates installer version
    nsis_patterns = [
        (r'!define APP_VERSION "[\d\.]+"', '!define APP_VERSION "{version}"')
    ]
    update_file(project_root / "installer.nsi", nsis_patterns, new_version)
    
    # 6. src/ui/main_window.py
    # Updates UI version labels
    ui_patterns = [
        (r'QLabel\("Version [\d\.]+"\)', 'QLabel("Version {version}")'),
        (r'QLabel\("v[\d\.]+"\)', 'QLabel("v{version}")'),
        (r"Version [\d\.]+</p>", "Version {version}</p>")
    ]
    update_file(project_root / "src/ui/main_window.py", ui_patterns, new_version)
    
    print(f"\nSUCCESS: All files updated to version {new_version}.")
    print("Next steps:")
    print(f"1. Update CHANGELOG.md manually.")
    print(f"2. git commit -am 'Bump version to v{new_version}'")
    print(f"3. git tag v{new_version}")
    print(f"4. git push origin main --tags")

if __name__ == "__main__":
    main()
