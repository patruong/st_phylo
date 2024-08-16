import cv2

# Load the original high-resolution histology image
original_image = cv2.imread('path_to_your_histology_image.jpg')

# Scale factors from your data
tissue_hires_scalef = 0.04111842
tissue_lowres_scalef = 0.012335527

# Calculate the new dimensions for hires and lowres images
hires_dims = (int(original_image.shape[1] * tissue_hires_scalef), int(original_image.shape[0] * tissue_hires_scalef))
lowres_dims = (int(original_image.shape[1] * tissue_lowres_scalef), int(original_image.shape[0] * tissue_lowres_scalef))

# Create the hires and lowres images
hires_image = cv2.resize(original_image, hires_dims, interpolation=cv2.INTER_AREA)
lowres_image = cv2.resize(original_image, lowres_dims, interpolation=cv2.INTER_AREA)

# Save the images
cv2.imwrite('hires_image.jpg', hires_image)
cv2.imwrite('lowres_image.jpg', lowres_image)