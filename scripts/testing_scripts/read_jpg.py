import scanpy as sc
import pandas as pd
import os 
import matplotlib.pyplot as plt
import json
from PIL import Image
import numpy as np
Image.MAX_IMAGE_PIXELS = None  # Remove image size limit

img_path = "/home/ptruong/git/st_phylo/test/"
img_file = "H1_4.jpg"
image_path = os.path.join(img_path, img_file)

def read_jpg_as_np_array(image_path, return_hires = True):
    file_size = os.path.getsize(image_path)
    max_size = 1000 * 1024 * 1024  # 1000 MB, adjust as neede

    if file_size > max_size:
        print(f"Image file is too large ({file_size / (1024*1024):.2f} MB). Skipping.")
    else:
        with Image.open(image_path) as img:
            # Calculate new size (e.g., reduce to 25% of original size)
            new_size = tuple(dim // 4 for dim in img.size)
            img_resized = img.resize(new_size, Image.LANCZOS)
            img_array_lowres = np.array(img_resized)
        if return_hires == True:
            img = Image.open(image_path)
            img_array_hires = np.array(img)
            return img_array_lowres, img_array_hires
        else:
            return img_array_lowres, False

img_array_lowres, _ =  read_jpg_as_np_array(image_path, return_hires=False)       