"""Section 4.2: draw the reconstructed images of Figure 8 (b)-(d).

Each reconstruction is rescaled to [0, 255] (same as rescale_255() in
run_reconstruction.R) and shown as a 32 x 32 grey-scale image.
Run from the repository root:  python 05_compressed_sensing/plot_reconstruction.py
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

results_dir = os.path.join("05_compressed_sensing", "results")
figures_dir = os.path.join("05_compressed_sensing", "figures")
os.makedirs(figures_dir, exist_ok=True)

# (input csv, output pdf)
files_and_titles = [
    ('reconst_def_N=700.csv',  'fig8b_glmnet_default.pdf'),
    ('reconst_set_N=700.csv',  'fig8c_glmnet_proposed.pdf'),
    ('reconst_lars_N=700.csv', 'fig8d_lars.pdf'),
]

for file_name, save_name in files_and_titles:
    file_path = os.path.join(results_dir, file_name)

    # 1. データの読み込み
    df = pd.read_csv(file_path, header=0)
    raw_values = df.apply(pd.to_numeric, errors='coerce').values
    raw_values = raw_values.flatten()

    # 2. Rの関数と同じロジックで0-255にスケーリング
    v_min = raw_values.min()
    v_max = raw_values.max()
    if v_max != v_min:
        scaled_values = (raw_values - v_min) / (v_max - v_min) * 255
    else:
        scaled_values = raw_values - v_min

    # 3. 32x32にリシェイプ
    image_data = scaled_values.reshape(32, 32)

    # 4. 画像としてプロット（スケーリング済みなのでvmin=0, vmax=255で固定）
    plt.figure(figsize=(4, 4))
    plt.imshow(image_data, cmap='gray', vmin=0, vmax=255)
    plt.axis('off')

    # 5. 保存
    save_path = os.path.join(figures_dir, save_name)
    plt.savefig(save_path, format='pdf', bbox_inches='tight', pad_inches=0)
    plt.close()

    print(f"Saved (Scaled 0-255): {save_path}")
