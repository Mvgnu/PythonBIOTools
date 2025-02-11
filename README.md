# PythonBIOTools

```markdown
# PythonBIOTools

## PAM Analysis Scripts

### Overview
This repository contains Python scripts for performing PAM (Position-Weight Matrix) analysis. These scripts are designed to analyze biological sequences and identify motifs based on PAM matrices.

### Script: `tifanalysisnoselector.py`

This script processes TIFF images to analyze and visualize regions of interest based on intensity deviations. It identifies unhealthy regions in the images by comparing pixel intensity values to a baseline intensity of healthy tissue.

### Installation
To use the PAM analysis scripts, clone this repository and install the required dependencies:
```bash
git clone https://github.com/Mvgnu/PythonBIOTools.git
cd PythonBIOTools
pip install -r requirements.txt
```

### Usage
Run the analysis script with the following command:
```bash
python tifanalysisnoselector.py
```

### Parameters
- `inputfolder`: Path to the input folder containing TIFF images (default: 'nrcinput22').
- `outputfolder`: Path to the output folder where results will be saved (default: 'nrcoutputadjustedthreshold').

### Script Workflow
1. **Loading Images**: Loads all TIFF images from the specified input folder.
2. **Image Processing**: Computes the mean of the first three slices, normalizes the values, and calculates the Fo, Fm, Fv, and Fv/Fm values.
3. **Background Removal**: Removes background using an FM threshold.
4. **Mask Cleaning**: Cleans up the mask using binary closing and opening operations.
5. **Intensity Analysis**: Calculates the average intensity of healthy tissue and identifies unhealthy regions based on deviations from this baseline.
6. **Region Analysis**: Labels and filters regions based on size.
7. **Visualization**: Generates and saves visualizations of the original mean image and detected low-intensity regions.
8. **Output**: Saves the results and region details to the specified output folder.

### Examples
To run the script, use the following command:
```bash
python tifanalysisnoselector.py
```

### Contributing
If you would like to contribute to this project, please fork the repository and submit a pull request.

### License
This project is licensed under the MIT License.

### Contact
For any questions or issues, please open an issue on this repository or contact the maintainer.
```
