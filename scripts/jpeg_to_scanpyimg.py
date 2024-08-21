import os
import cv2
from PIL import Image
import numpy as np
Image.MAX_IMAGE_PIXELS = None  # Remove image size limit

img_path = "/home/ptruong/git/st_phylo/test/"
img_file = "H3_2.jpg"
image_path = os.path.join(img_path, img_file)

def read_jpg_to_cv_img(image_path, return_hires = True):
    file_size = os.path.getsize(image_path)
    max_size = 1000 * 1024 * 1024  # 1000 MB, adjust as neede

    if file_size > max_size:
        print(f"Image file is too large ({file_size / (1024*1024):.2f} MB). Skipping.")
    else:
        with Image.open(image_path) as img:import os

import cv2
from PIL import Image
import numpy as np
Image.MAX_IMAGE_PIXELS = None  # Remove image size limit

img_path = "/home/ptruong/git/st_phylo/test/"
img_file = "H3_2.jpg"
image_path = os.path.join(img_path, img_file)

def read_jpg_to_cv_img(image_path, return_hires = True):
    file_size = os.path.getsize(image_path)
    max_size = 1000 * 1024 * 1024  # 1000 MB, adjust as neede

    if file_size > max_size:
        print(f"Image file is too large ({file_size / (1024*1024):.2f} MB). Skipping.")
    else:
        with Image.open(image_path) as img:
            # Calculate new size (e.g., reduce to 25% of original size)
            new_size = tuple(dim // 4 for dim in img.size)
            img_resized = img.resize(new_size, Image.LANCZOS)
            #img_array_lowres = np.array(img_resized)
        if return_hires == True:
            img = Image.open(image_path)
            #img_array_hires = np.array(img)
            return img_resized, img
        else:
            return img_resized, False


#def jpg_to_lowres_hires_img(image_path, tissue_hires_scalef, tissue_lowres_scalef, output_lowres, output_hires):

def jpg_to_lowres_hires_img(image_path, output_lowres, output_hires):
    img_resized, _ =  read_jpg_to_cv_img(image_path, return_hires=False)       

    max_size = (3000, 3000)  # Adjust this to a size that's manageable
    img_resized.thumbnail(max_size, Image.LANCZOS)
    #img_resized.show()  # To display the image

    # Convert the PIL image to a NumPy array and then to OpenCV format
    image = np.array(img_resized)

    # Convert RGB to BGR if necessary (since OpenCV uses B3GR by default)
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)

    # Load the original high-resolution histology image
    original_image = cv2.imread(image_path)

    # Scale factors from your data from json file
    #tissue_hires_scalef = 0.04111842
    #tissue_lowres_scalef = 0.012335527

    # Calculate the new dimensions for hires and lowres images
    #hires_dims = (int(original_image.shape[1] * tissue_hires_scalef), int(original_image.shape[0] * tissue_hires_scalef))
    #lowres_dims = (int(original_image.shape[1] * tissue_lowres_scalef), int(original_image.shape[0] * tissue_lowres_scalef))

    # Create the hires and lowres images
    #hires_image = cv2.resize(original_image, hires_dims, interpolation=cv2.INTER_AREA)
    #lowres_image = cv2.resize(original_image, lowres_dims, interpolation=cv2.INTER_AREA)
    #hires_image = cv2.resize(original_image, hires_dims, interpolation=cv2.INTER_LINEAR)
    #lowres_image = cv2.resize(original_image, lowres_dims, interpolation=cv2.INTER_LINEAR)

    # Save the images
    #os.chdir("/home/ptruong/git/st_phylo/test")
    #cv2.imwrite('hires_image.png', hires_image)
    #cv2.imwrite('lowres_image.png', lowres_image)
    cv2.imwrite(output_hires, image)
    cv2.imwrite(output_lowres, image)


os.chdir("/home/ptruong/git/st_phylo/test")
#jpg_to_lowres_hires_img(image_path, tissue_hires_scalef =  0.04111842, tissue_lowres_scalef = 0.012335527, 
#                        output_lowres = "tissue_lowres_image.png", 
#                        output_hires = "tissue_hires_image.png")

jpg_to_lowres_hires_img(image_path,  
                        output_lowres = "tissue_lowres_image.png", 
                        output_hires = "tissue_hires_image.png")
            # Calculate new size (e.g., reduce to 25% of original size)
            new_size = tuple(dim // 4 for dim in img.size)
            img_resized = img.resize(new_size, Image.LANCZOS)
            #img_array_lowres = np.array(img_resized)
        if return_hires == True:
            img = Image.open(image_path)
            #img_array_hires = np.array(img)
            return img_resized, img
        else:
            return img_resized, False


#def jpg_to_lowres_hires_img(image_path, tissue_hires_scalef, tissue_lowres_scalef, output_lowres, output_hires):

def jpg_to_lowres_hires_img(image_path, output_lowres, output_hires):
    img_resized, _ =  read_jpg_to_cv_img(image_path, return_hires=False)       

    max_size = (3000, 3000)  # Adjust this to a size that's manageable
    img_resized.thumbnail(max_size, Image.LANCZOS)
    #img_resized.show()  # To display the image

    # Convert the PIL image to a NumPy array and then to OpenCV format
    image = np.array(img_resized)

    # Convert RGB to BGR if necessary (since OpenCV uses B3GR by default)
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)

    # Load the original high-resolution histology image
    original_image = cv2.imread(image_path)

    # Scale factors from your data from json file
    #tissue_hires_scalef = 0.04111842
    #tissue_lowres_scalef = 0.012335527

    # Calculate the new dimensions for hires and lowres images
    #hires_dims = (int(original_image.shape[1] * tissue_hires_scalef), int(original_image.shape[0] * tissue_hires_scalef))
    #lowres_dims = (int(original_image.shape[1] * tissue_lowres_scalef), int(original_image.shape[0] * tissue_lowres_scalef))

    # Create the hires and lowres images
    #hires_image = cv2.resize(original_image, hires_dims, interpolation=cv2.INTER_AREA)
    #lowres_image = cv2.resize(original_image, lowres_dims, interpolation=cv2.INTER_AREA)
    #hires_image = cv2.resize(original_image, hires_dims, interpolation=cv2.INTER_LINEAR)
    #lowres_image = cv2.resize(original_image, lowres_dims, interpolation=cv2.INTER_LINEAR)

    # Save the images
    #os.chdir("/home/ptruong/git/st_phylo/test")
    #cv2.imwrite('hires_image.png', hires_image)
    #cv2.imwrite('lowres_image.png', lowres_image)
    cv2.imwrite(output_hires, image)
    cv2.imwrite(output_lowres, image)


os.chdir("/home/ptruong/git/st_phylo/test")
#jpg_to_lowres_hires_img(image_path, tissue_hires_scalef =  0.04111842, tissue_lowres_scalef = 0.012335527, 
#                        output_lowres = "tissue_lowres_image.png", 
#                        output_hires = "tissue_hires_image.png")

jpg_to_lowres_hires_img(image_path,  
                        output_lowres = "tissue_lowres_image.png", 
                        output_hires = "tissue_hires_image.png")