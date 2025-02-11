import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.widgets import Button
from scipy import ndimage
from skimage.measure import label, regionprops
from tifffile import imread


class RegionSelector:
    def __init__(self):
        self.regions = []
        self.fig = None
        self.done = False
        self.current_points = []
        self.current_line = None
        self.current_polygon = None

    def onclick(self, event):
        if event.inaxes != self.ax:
            return
        if event.button == 1:  # Left click
            self.current_points.append((event.xdata, event.ydata))

            # Draw points
            self.ax.plot(event.xdata, event.ydata, 'ro')

            # Update line
            if len(self.current_points) > 1:
                xs, ys = zip(*self.current_points)
                if self.current_line:
                    self.current_line.remove()
                self.current_line, = self.ax.plot(xs, ys, 'r-')

            self.fig.canvas.draw()

    def onkey(self, event):
        if event.key == 'enter':  # Complete current rosette
            if len(self.current_points) > 2:  # Need at least 3 points
                # Close the polygon
                points = np.array(self.current_points)
                polygon = Polygon(points, fill=False, edgecolor='red')
                self.ax.add_patch(polygon)

                # Store the region
                self.regions.append(self.current_points.copy())

                # Clear current points for next rosette
                self.current_points = []
                if self.current_line:
                    self.current_line.remove()
                    self.current_line = None

                self.fig.canvas.draw()

    def done_callback(self, event):
        self.done = True
        plt.close(self.fig)

    def select_regions(self, image):
        self.fig, self.ax = plt.subplots()
        self.ax.imshow(image, cmap='gray')
        self.ax.set_title('Click to draw points, Enter to complete rosette')

        # Add Done button
        ax_button = plt.axes([0.8, 0.01, 0.1, 0.05])
        btn_done = Button(ax_button, 'Done')
        btn_done.on_clicked(self.done_callback)

        # Connect events
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)

        plt.show()


# Modified main loop
for filename in os.listdir('nrcinput22'):
    if filename.endswith('.tif'):
        print(f"Processing {filename}")

        # Load and prepare image as before
        image_stack = imread(os.path.join('nrcinput22', filename))
        mean_slice = np.mean(image_stack, axis=0)

        fo = image_stack[0].astype(float)  # Slice 1 (Fo)
        fm = image_stack[1].astype(float)  # Slice 2 (Fm)
        fv = fm - fo
        fvfm = np.divide(fv, fm, out=np.zeros_like(fv), where=fm > 15)
        normalized_slice = (fvfm - np.min(fvfm)) / (np.max(fvfm) - np.min(fvfm))

        # Create region selector and get user input
        region_selector = RegionSelector()
        region_selector.select_regions(mean_slice)

        # Create mask for selected regions
        mask = np.zeros_like(normalized_slice, dtype=bool)
        height, width = mask.shape

        # Create coordinate matrices
        x, y = np.meshgrid(np.arange(width), np.arange(height))

        for region_points in region_selector.regions:
            # Convert points to numpy arrays for mask creation
            region_points = np.array(region_points)
            from matplotlib.path import Path

            # Create path from points
            path = Path(region_points)

            # Create mesh grid for all image points
            coords = np.column_stack((x.ravel(), y.ravel()))

            # Test which points are inside the polygon
            mask_region = path.contains_points(coords).reshape(height, width)

            # Update the main mask
            mask = mask | mask_region

        # Apply mask to leaf_mask
        leaf_mask = (normalized_slice > 0.1) & mask
        leaf_mask = ndimage.binary_closing(leaf_mask, iterations=3)

        # Calculate the average intensity of healthy tissue (including selected regions)
        leaf_values = normalized_slice
        median = np.median(leaf_values)
        std = np.std(leaf_values)
        filtered_values = leaf_values[(leaf_values - median) > -2 * std]  # Only consider values above -2 standard deviations
        healthy_intensity = np.mean(filtered_values)
        healthy_intensity = np.mean(normalized_slice[leaf_mask])
        print(f"Healthy tissue baseline intensity: {healthy_intensity:.2f}")
        # Define structural element first
        struct_element = np.array([[1, 1, 1],
                                   [1, 1, 1],
                                   [1, 1, 1]])  # 8-connected structural element

        # Initial mask with higher threshold
        deviation_threshold = 0.15
        hr_mask = (normalized_slice < (healthy_intensity - deviation_threshold)) & leaf_mask

        # Region growing with lower threshold for neighboring pixels
        neighbor_threshold = 0.08
        growing_mask = hr_mask.copy()
        while True:
            dilated = ndimage.binary_dilation(growing_mask, structure=struct_element)
            new_pixels = (dilated & ~growing_mask &
                          (normalized_slice < (healthy_intensity - neighbor_threshold)) &
                          leaf_mask)
            if not np.any(new_pixels):
                break
            growing_mask = growing_mask | new_pixels

        hr_mask = growing_mask  # Replace original mask with grown regions

        # More aggressive diagonal closing
        hr_mask = ndimage.binary_closing(hr_mask, structure=struct_element, iterations=5)
        hr_mask = ndimage.binary_opening(hr_mask, structure=struct_element, iterations=1)

        # Label connected regions with diagonal connectivity
        labeled_regions = label(hr_mask, connectivity=2)  # 8-connectivity
        regions = regionprops(labeled_regions, intensity_image=normalized_slice)

        # Filter regions based on size
        min_size = 200
        filtered_regions = [region for region in regions if region.area >= min_size]

        # Visualization
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Original image
        ax1.imshow(normalized_slice, cmap='gray')
        ax1.set_title('Original Normalized Image')
        ax1.text(10, 20, f'Healthy Avg: {healthy_intensity:.2f}',
                 color='white', bbox=dict(facecolor='green', alpha=0.5))
        plt.colorbar(ax1.imshow(normalized_slice, cmap='gray'), ax=ax1)

        # Detected regions
        ax2.imshow(normalized_slice, cmap='gray')
        ax2.text(10, 20, f'Healthy Avg: {healthy_intensity:.2f}',
                 color='white', bbox=dict(facecolor='green', alpha=0.5))

        for region in filtered_regions:
            # Get the coordinates of the region boundary
            boundary = region.coords

            # Plot the region boundary
            contour = plt.plot(boundary[:, 1], boundary[:, 0],
                               color='red', linewidth=2, linestyle='-')

            # Add the mean intensity text at centroid
            centroid = region.centroid
            mean_intensity = region.mean_intensity
            std_intensity = np.std(normalized_slice[region.coords[:, 0], region.coords[:, 1]])
            ax2.text(centroid[1], centroid[0], f'{mean_intensity:.2f}±{std_intensity:.2f}',
                     color='white', ha='center', va='center',
                     bbox=dict(facecolor='blue', alpha=0.5))

        ax2.set_title('Detected Low Intensity Regions')
        plt.tight_layout()

        # Save the plot
        output_filename = os.path.join('nrcoutputadjustedthreshold', f'{os.path.splitext(filename)[0]}_output.png')
        plt.savefig(output_filename)
        plt.close()  # Close the figure to avoid memory overload
        # Create text file with region information
        txt_output_filename = os.path.join('nrcoutputadjustedthreshold', f'{os.path.splitext(filename)[0]}_regions.txt')
        with open(txt_output_filename, 'w') as f:
            f.write(f"Analysis for {filename}\n")
            f.write(f"Healthy tissue baseline intensity: {healthy_intensity:.2f}\n\n")
            f.write("Region Details:\n")
            for i, region in enumerate(filtered_regions, 1):
                f.write(f"\nRegion {i}:\n")
                f.write(f"Area (pixels): {region.area}\n")
                f.write(f"Mean intensity: {region.mean_intensity:.2f}\n")
                f.write(f"Min intensity: {region.min_intensity:.2f}\n")
                f.write(f"Max intensity: {region.max_intensity:.2f}\n")
                f.write(f"Centroid location: ({region.centroid[0]:.1f}, {region.centroid[1]:.1f})\n")