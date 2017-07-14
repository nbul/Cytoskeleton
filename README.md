# Cytoskeleton

This is a script to analyse cytoskeleton organisation in 2D images

## Input data

The following three sub-folders within a single directory need to be present:
1. **borders** with the 8-bit tif images of cell outlines (e.g. E-cadherin staining) segmented in Packing Analyser. Each image will have a subfoder with the same name, that need to contain "tracked_bd.png" mask.

1. **cytoskeleton** with the corresponding tif projections (we use maximum projection in this case), which will be analysed for alignminet of cytoskeleton fibers within cells.

1. **cytoskeleton_average** with the corresponding tif projections (we use average projection in this case), which will be analysed for density and uniformity of cytoskeleton organisation.

Images should be numbered, but the numbers do not have to be sequencial, and the names might contain additional information. The numbers in image names will be used in the output files.

## Running the script
1. Start **actin_automatic_v1.m** script.

1. Select a folder with input data.

1. Select the parameters you would like to analyse: alignment of cytoskeleton (*MTSD*) or cytoskeleton distribution (*Density*).

1. In case you select density, choose the threshold method: Otsu's method with a threshold calculated on an adjusted image (*Otsu*), or calculating a threshold from intensity of edges detected in an image (*Edges*).

## Output files and format

The following folders will be automatically created depending on the selection at the start of a run.
1. **distribution** with cell-by-cell signal direction distributions for each image (**image number**_distribution.csv) in case *MTSD* option for analysis is selected.

1. **images_analysed** with adjusted images of cytoskeleton, projection of a cell mask on them and numbers of each analysed cell.<br>
![Example of analysed image](images/8_analysed_image.png)

1. **summary** with the summarysed data.

In the summary folder the subfolders will be created depending on the run choice: *edges*, *otsu* or *MTSD*. These subfolders will contain **image number**_summary**run choice**.csv files with cell-by-cell data for each analysed image. The cell numbers correspond to that in *images_analysed* files.

Additionally, in the summary folder itself will be files with data averaged per image: *Summary_edges.csv*, *Summary_otsu.csv* or *Summary_MTSD.csv* depending on the run choice.

## Output data explained

### MTSD
The following columns record cell-by-cell data:
1. **Cell** - cell number, which correspond to the cell number in *images_analysed*.
1. **Area** - area in pixels of each cell.
1. **Eccentricity** - eccentricity of each cell.
1. **Direction_cell** - direction of the main axis of each cell calculated from regionprops.
1. **SD** - a measure of cytoskeleton alignment, which is a standard deviation of signal direction distribution, fitted with Von-Misus distribution. For detailed explanation of how it is calculated see Gomez J.M., Chumakova L., Bulgakova N.A., Brown N.H. (2016) Microtubule organization is determined by the shape of epithelial cells. Nature Communications, 7. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/27779189).
1. **DEV** - a measure of direction of overall cytoskeleton orientation - mean of signal direction distribution, fitted with Von-Misus distribution.
1. **Elongation** - aspect ratio of each cell calculated from eccentricity.
1. **Alignment** - an alernative measure for cytoskeleton alignment, calculated as proportion of signal within 10° of mean signal direction.

### Density 
Independent on the chosen way to threshold, the resulting files will contain the following cell-by-cell data about the cytoskeleton signal:
1. **Cell** - cell number, which correspond to the cell number in *images_analysed*.
1. **Signal area** - proportion of pixels above threshold.
1. **Density** - mean intensity of pixels above threshold with subtracted mean intensity of pixels below threshold, and then normalized to mean intensity of pixels below threshold and signal area.
1. **Bundling** - mean intensity of pixels above threshold normalized to threshold.
1. **Uniformity** - standard deviation of signal in each cell without thresholding.
1. **Sparseness** - sparseness of the signal, calsulated as discribed in Lin A.C., Bygrave A.M., de Calignon A., Lee T., Miesenböck G. (2014). Sparse, decorrelated odor coding in the mushroom body enhances learned odor discrimination. Nat Neurosci. 17(4):559-68. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/24561998). 
1. **Skewness** - skewness of the signal above threshold.
1. **Kurtosis** - kurtosis of the signal above threshold.

## References
1. Gomez J.M., Chumakova L., Bulgakova N.A., Brown N.H. (2016) Microtubule organization is determined by the shape of epithelial cells. Nature Communications, 7. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/27779189).
1. Lin A.C., Bygrave A.M., de Calignon A., Lee T., Miesenböck G. (2014). Sparse, decorrelated odor coding in the mushroom body enhances learned odor discrimination. Nat Neurosci. 17(4):559-68. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/24561998).