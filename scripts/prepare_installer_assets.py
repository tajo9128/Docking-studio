from PIL import Image, ImageOps
import sys
import os

def prepare_assets(source_path, dist_dir):
    """
    Converts and resizes the source branding image into NSIS-compatible BMPs.
    """
    if not os.path.exists(source_path):
        print(f"Error: Source image not found at {source_path}")
        return

    if not os.path.exists(dist_dir):
        os.makedirs(dist_dir)

    try:
        # Load original image
        original = Image.open(source_path).convert('RGB')
        
        # 1. Create Header Image (150x57) - Logo on Right
        # NSIS MUI_HEADERIMAGE usually is just the small logo
        header_size = (150, 57)
        header_img = Image.new('RGB', header_size, (255, 255, 255))
        
        # Resize original to fit height
        ratio = min(header_size[0] / original.width, header_size[1] / original.height)
        new_size = (int(original.width * ratio), int(original.height * ratio))
        resized = original.resize(new_size, Image.Resampling.LANCZOS)
        
        # Center in the header box
        pos = ((header_size[0] - new_size[0]) // 2, (header_size[1] - new_size[1]) // 2)
        header_img.paste(resized, pos)
        
        header_path = os.path.join(dist_dir, 'header.bmp')
        header_img.save(header_path, 'BMP')
        print(f"Created Header Image: {header_path}")

        # 2. Create Welcome/Finish Sidebar (164x314)
        # We'll place the logo at the top on a white background
        side_size = (164, 314)
        side_img = Image.new('RGB', side_size, (255, 255, 255))
        
        # Resize to fit width with some padding
        target_width = 140
        ratio = target_width / original.width
        new_size = (int(original.width * ratio), int(original.height * ratio))
        resized = original.resize(new_size, Image.Resampling.LANCZOS)
        
        # Place at top with padding
        side_img.paste(resized, (12, 20))
        
        side_path = os.path.join(dist_dir, 'welcome.bmp')
        side_img.save(side_path, 'BMP')
        print(f"Created Welcome Image: {side_path}")
        
    except Exception as e:
        print(f"Error processing images: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Hardcoded path from user metadata for this session
    # In a real CI, this would be an env var or arg, but here we bridge the gap
    # The user file is in: C:/Users/tajo9/.gemini/antigravity/brain/560e6093-5f6d-4ea2-bb56-663d9107c913/uploaded_image_1767194043413.jpg
    
    # We will assume release.yml passes arguments: [script] [source_image] [dest_dir]
    if len(sys.argv) < 3:
        print("Usage: prepare_installer_assets.py <source_image> <dest_dir>")
        sys.exit(1)
        
    prepare_assets(sys.argv[1], sys.argv[2])
