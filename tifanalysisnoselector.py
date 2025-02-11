import numpy as np
import matplotlib.pyplot as plt
from tifffile import imread
from skimage.measure import label, regionprops
from scipy import ndimage
import os

inputfolder = 'nrcinput22'
outputfolder = 'nrcoutputadjustedthreshold'
# Process each image in the directory - automatically selects all tif files from input folder 'nrcinput22'
for filename in os.listdir(inputfolder):
    if filename.endswith('.tif'):
        print(f"Processing {filename}")

        # Load and prepare image by loading all files from path 'nrcinput22' is the input folder
        image_stack = imread(os.path.join(inputfolder, filename))

        # Calculate the mean of the first three slices to create a visually appealing output image
        mean_slice = np.mean(image_stack[0:3], axis=0)

        # Normalize the values to the range 0-1 to visualize HR and fluorescence in output
        mean_slice = (mean_slice - mean_slice.min()) / (mean_slice.max() - mean_slice.min())

        fo = image_stack[0].astype(float)  # Slice 1 that contains the values of Fo
        fm = image_stack[1].astype(float)  # Slice 2 that contains the values of Fm
        fv = fm - fo # calculate fv
        fvfm = np.divide(fv, fm, out=np.zeros_like(fv), where=fm > 15) #added threshold fm value to avoid background in calculations, in fm generally below 15, may be adjusted
        normalized_slice = (fvfm - np.min(fvfm)) / (np.max(fvfm) - np.min(fvfm)) # normalize fvfm values

        # Simple background removal using FM threshold
        leaf_mask = fm > 10

        # Clean up the mask
        leaf_mask = ndimage.binary_closing(leaf_mask, iterations=3) # Fill small gaps (holes) within the mask
        leaf_mask = ndimage.binary_opening(leaf_mask, iterations=1) # Remove small, isolated noisy regions for more accurate results

        # Calculate the average intensity of healthy tissue
        leaf_values = normalized_slice[leaf_mask]
        mean = np.mean(leaf_values)
        std = np.std(leaf_values)
        filtered_values = leaf_values[(leaf_values - mean) > -2 * std]
        healthy_intensity = np.mean(filtered_values)
        print(f"Healthy tissue baseline intensity: {healthy_intensity:.3f}")

        # Define structural element
        struct_element = np.array([[1, 1, 1],
                                 [1, 1, 1],
                                 [1, 1, 1]])  # 8-connected structural element for HR regions sampling

        # Create mask for unhealthy regions
        deviation_threshold = 0.05 # deviation above this percentile are considered unhealthy regions
        hr_mask = (normalized_slice < (healthy_intensity - deviation_threshold)) & leaf_mask # select regions with deviation above threshold

        # Adjust Region, use to grow with lower threshold may be adjusted to better detect edges
        neighbor_threshold = 0.05
        growing_mask = hr_mask.copy()
        # loop selecting neighboring regions meeting deviation neighbor_threshold
        while True:
            dilated = ndimage.binary_dilation(growing_mask, structure=struct_element)
            new_pixels = (dilated & ~growing_mask &
                        (normalized_slice < (healthy_intensity - neighbor_threshold)) &
                        leaf_mask)
            if not np.any(new_pixels):
                break
            growing_mask = growing_mask | new_pixels


        #replace hr_mask with calculated mask containing neighboring threshold deviating pixels
        hr_mask = growing_mask

        # Clean up the regions
        hr_mask = ndimage.binary_closing(hr_mask, structure=struct_element, iterations=5)
        hr_mask = ndimage.binary_opening(hr_mask, structure=struct_element, iterations=1)

        # Label and analyze regions
        labeled_regions = label(hr_mask, connectivity=2)
        regions = regionprops(labeled_regions, intensity_image=normalized_slice)

        # Filter regions based on size
        min_size = 200
        filtered_regions = [region for region in regions if region.area >= min_size]

        # Visualization
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Original image
        ax1.imshow(mean_slice, cmap='gray')
        ax1.set_title('Original Mean Image')
        ax1.text(10, 20, f'Healthy Avg: {healthy_intensity:.3f}',
                color='white', bbox=dict(facecolor='green', alpha=0.5))
        plt.colorbar(ax1.imshow(mean_slice, cmap='gray', vmin=0, vmax=1), ax=ax1)

        # Detected regions
        ax2.imshow(mean_slice, cmap='gray')
        ax2.text(10, 20, f'Healthy Avg: {healthy_intensity:.3f}',
                color='white', bbox=dict(facecolor='green', alpha=0.5))

        for region in filtered_regions:
            # Plot region boundary
            boundary = region.coords
            plt.plot(boundary[:, 1], boundary[:, 0],
                    color='red', linewidth=2, linestyle='-')

            # Add intensity information
            centroid = region.centroid
            mean_intensity = region.mean_intensity
            std_intensity = np.std(normalized_slice[region.coords[:, 0], region.coords[:, 1]])
            ax2.text(centroid[1], centroid[0], f'{mean_intensity:.3f}±{std_intensity:.3f}',
                    color='white', ha='center', va='center',
                    bbox=dict(facecolor='blue', alpha=0.5))

        ax2.set_title('Detected Low Intensity Regions')
        plt.tight_layout()

        # Save outputs
        output_filename = os.path.join(outputfolder, f'{os.path.splitext(filename)[0]}_output.png')
        plt.savefig(output_filename)
        plt.close()

        # Save region information
        txt_output_filename = os.path.join(outputfolder, f'{os.path.splitext(filename)[0]}_regions.txt')
        with open(txt_output_filename, 'w') as f:
            f.write(f"Analysis for {filename}\n")
            f.write(f"Healthy tissue baseline intensity: {healthy_intensity:.3f}\n\n")
            f.write("Region Details:\n")
            for i, region in enumerate(filtered_regions, 1):
                f.write(f"\nRegion {i}:\n")
                f.write(f"Area (pixels): {region.area}\n")
                f.write(f"Mean intensity: {region.mean_intensity:.3f}\n")
                f.write(f"Min intensity: {region.min_intensity:.3f}\n")
                f.write(f"Max intensity: {region.max_intensity:.3f}\n")
                f.write(f"Centroid location: ({region.centroid[0]:.1f}, {region.centroid[1]:.1f})\n")