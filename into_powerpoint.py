from pptx import Presentation
from pptx.util import Inches
import os

prs = Presentation()
image_dir = 'img'

images = os.listdir(image_dir)

image_groups = {}
for image in images:
    if image.endswith(('.png')):
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


for gene_name, images in image_groups.items():
    slide = prs.slides.add_slide(prs.slide_layouts[5])  # Using a blank slide layout

    title_placeholder = slide.shapes.title
    title_placeholder.text = gene_name

    if images["average"] is not None:
        log_image_path = os.path.join(image_dir, images["average"])
        slide.shapes.add_picture(log_image_path, Inches(0), Inches(2.5), width=Inches(5))
    if images["log"] is not None:
        avg_image_path = os.path.join(image_dir, images["log"])
        slide.shapes.add_picture(avg_image_path, Inches(5), Inches(2.5), width=Inches(5))

prs.save('gene_images_presentation.pptx')
