from PIL import Image
import os
from tqdm import tqdm

# bin_types = ['mergerz_A', 'mergerz_B']
bin_types = ['NFWconc']
# bin_types = ['mergerz_mass_30_A', 'mergerz_mass_30_B', 'mergerz_mass_35_A', 'mergerz_mass_35_B', 'mergerz_mass_40_A', 'mergerz_mass_40_B', 'mergerz_mass_45']
snaps = [214]  # MTNG snapshots

def collect_all_pngs(folder):
    all_pngs = []
    for dirpath, _, filenames in os.walk(folder):
        for filename in filenames:
            if filename.lower().endswith('.png'):
                all_pngs.append(os.path.join(dirpath, filename))
    return sorted(all_pngs)

def resize_image(path, max_width=400):
    img = Image.open(path)
    if img.width > max_width:
        ratio = max_width / img.width
        new_size = (max_width, int(img.height * ratio))
        return img.resize(new_size, resample=Image.LANCZOS)
    return img

def pngs_to_gif(image_paths, out_path='output.gif', fps=60):
    duration = int(1000 / fps)  # 每帧毫秒
    print(f"[🧩 逐帧处理] 共 {len(image_paths)} 帧，目标帧率 {fps} FPS")
    
    images = []
    for path in tqdm(image_paths):
        img = resize_image(path)
        img = img.convert("P", palette=Image.ADAPTIVE)
        images.append(img)
        
    total_duration = (duration * len(images)) / 1000
    print(f"[⏱️ 动画时长] 共 {len(images)} 帧，每帧 {duration}ms, 动画总时长 ≈ {total_duration:.2f} 秒")

    print(f"[💾 保存动画] 保存为 {out_path}")
    images[0].save(out_path, save_all=True, append_images=images[1:], duration=duration, loop=0, optimize=True)
    print(f"[✅ 完成] 文件大小: {os.path.getsize(out_path)/(1024**2):.2f} MB")

save_dir = f'result/bootstrap_videos/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
if __name__ == "__main__":
    for bin_type in bin_types:
        for snap in snaps:
            root_folder = f"result/bootstrap_plots/MTNG/Hydro-Arepo/MTNG-L500-4320-A/with_{bin_type}/snap_{snap}"
            image_paths = collect_all_pngs(root_folder)
            pngs_to_gif(image_paths, out_path=f"{save_dir}/{bin_type}_{snap}_profiles.gif", fps=60)