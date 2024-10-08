import os
import squidpy as sq
from PIL import Image
img_path = "/home/ptruong/git/st_phylo/test/"
img_file = "H3_2.jpg"
image_path = os.path.join(img_path, img_file)

with Image.open(image_path) as img:
    max_size = (3000, 3000)    
    img.thumbnail(max_size, Image.LANCZOS)
    resized_image_path = os.path.join(img_path, "H2_1_resized.jpg")
    img.save(resized_image_path)



img = sq.im.ImageContainer(resized_image_path)
#img = sq.im.ImageContainer(image_path)
sq.im.process(img_container, layer="image", method="smooth")


adata.uns["spatial"]["{sample_id}"]["images"]["hires"] = img

