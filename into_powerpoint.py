from pptx import Presentation
from pptx.util import Inches
import os


def powerpoint_generate_for_julia(image_folder):
    prs = Presentation()
    image_dir = image_folder

    images = os.listdir(image_dir)

    image_groups = {}
    for image in images:
        if image.endswith(('.png', '.jpg')):
            #  I have two types of images naming as : Average of Log Transformed Intensity for ABL1 .png
            #                                         Log transformed Intensity for RAD9A for each DMSO and whel run.jpg
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
    # Got a dictionary of the set of files each with the average comparison and runs individual comparison.


    for gene_name, images in image_groups.items():
        # Using a blank slide layout
        slide = prs.slides.add_slide(prs.slide_layouts[5])

        title_placeholder = slide.shapes.title
        title_placeholder.text = gene_name

        # Set the position of each image.
        if images["average"] is not None:
            log_image_path = os.path.join(image_dir, images["average"])
            slide.shapes.add_picture(log_image_path, Inches(0), Inches(2.5), width=Inches(5))
        if images["log"] is not None:
            avg_image_path = os.path.join(image_dir, images["log"])
            slide.shapes.add_picture(avg_image_path, Inches(5), Inches(2.5), width=Inches(5))

    prs.save('gene_images_presentation.pptx')

if __name__ == '__main__':
    powerpoint_generate_for_julia("img")