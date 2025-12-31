from PIL import Image
import sys
import os

def convert_to_ico(source, target):
    try:
        img = Image.open(source)
        # Resize to standard icon sizes
        icon_sizes = [(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (16, 16)]
        img.save(target, format='ICO', sizes=icon_sizes)
        print(f"Successfully converted {source} to {target}")
    except Exception as e:
        print(f"Error converting image: {e}")
        sys.exit(1)

if __name__ == "__main__":
    convert_to_ico("biodockify_icon.png", "src/ui/styles/icon.ico")
