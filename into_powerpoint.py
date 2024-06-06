from pptx import Presentation
from pptx.util import Inches
import os

# Initialize a presentation object
prs = Presentation()

# Path to the directory containing the images
image_dir = 'img'

# Get the list of images in the directory
images = os.listdir(image_dir)

# Group images by gene name
image_groups = {}
for image in images:
    if image.endswith(('.png', '.jpg', '.jpeg')):
        # Extract the gene name from the image filename
        if 'Average' in image:
            gene_name = image.split('for ')[-1].split(' ')[0]
            if gene_name not in image_groups:
                image_groups[gene_name] = {}
            image_groups[gene_name]['average'] = image
        elif 'for each' in image:
            gene_name = image.split('Log transformed Intensity for ')[-1].split(' ')[0]
            if gene_name not in image_groups:
                image_groups[gene_name] = {}
            image_groups[gene_name]['log'] = image


# Create slides and add images
for gene_name, images in image_groups.items():
    slide = prs.slides.add_slide(prs.slide_layouts[5])  # Using a blank slide layout

    # Add title
    title_placeholder = slide.shapes.title
    title_placeholder.text = gene_name

    # Add log transformed intensity image
    if images["average"] is not None:
        log_image_path = os.path.join(image_dir, images["average"])
        slide.shapes.add_picture(log_image_path, Inches(0), Inches(2.5), width=Inches(5))

    # Add average of log transformed intensity image
    if images["log"] is not None:
        avg_image_path = os.path.join(image_dir, images["log"])
        slide.shapes.add_picture(avg_image_path, Inches(5), Inches(2.5), width=Inches(5))

# Save the presentation
prs.save('gene_images_presentation.pptx')
