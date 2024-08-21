import os
from PIL import Image
Image.MAX_IMAGE_PIXELS = None  # Remove image size limit
os.chdir("/home/ptruong/git/st_phylo/test")


img_path = "/home/ptruong/git/st_phylo/data/Histological_images/Patient_2/"
img_file = "H2_1.jpg"
image_path = os.path.join(img_path, img_file)

img = Image.open(image_path)
img.save("tissue_hires_image.png")
lowres_img = img.resize((img.width // 2, img.height // 2))
lowres_img.save("tissue_lowres_image.png")


