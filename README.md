# A project to find check centromereic signal of Endogenous and Trans-gene proteins

This plugin takes as input a folder containing 4 colour .ims Z-stacks files with a channel order of DNA (DAPI), Endogenous protein (GFP), Transgene protein(Alexa 555), Centromere Marker (Alexa 647). 
The plugin also allows users to set thresholds for the centromeric/nuclear intensity for the 488 and 555 channels to define total or partial localisation to the centromeres.
The plugin:
- Finds the cell outlines using the 647 channel and cellpose with the cyto2 model with a 60 pixel cell diameter.
- For each Cell outline finds the nucleii in the DAPI channel using IJ 'Default' thresholding, a 250-Infinity pixel size and circularity=0.30-1.00
- For each Nucleus ROi finds the centrosomes in the 647 channel using the 'Triangle' thresholding and size=5-500 pixel size cut off.
- Outputs to a Results.csv in the original folder the areas and intensities of the nuclear and centromeric ROIs for each nucleus in the 488, 555 and 647 channels and whether the 488 and 555 show partial or full localisation based on the user defined thresholds.
- Outputs an Overlay image to the original folder of the Z-projected 488 (green),555 (red), 647 (blue) and nuclear and centromeric labels and outlines (white).

This repository is based off an example Maven project implementing an ImageJ2 command.
Visit [this link](https://github.com/imagej/example-imagej2-command/generate)
to view the original template repository.

### Pre-requisites:
To use this script you need to have Cellpose installed on your computer and the BIOP Cellpose wrapper installed in your IJ plugins.
To install Cellpose to use with ImageJ on a windows computer;
- Install Anaconda from https://www.anaconda.com/products/distribution
- Add Anaconda to path https://www.datacamp.com/tutorial/installing-anaconda-windows, you need to have admin rights on your computer.
- Install Cellpose https://github.com/MouseLand/cellpose#local-installation
- Add the ijl-utilities-wrapper jar to your ImageJ plugins folder (download the jar here: https://mvnrepository.com/artifact/ch.epfl.biop/ijl-utilities-wrappers/0.3.23)
- BIOP should appear as a plugin, navigate to BIOP -> Cellpose -> Define Env and prefs
- Set CellposeEnvDirectory to your Anaconda\envs\cellpose . On my computer this is C:\Users\username\Anaconda3\envs\cellpose (if this doesn’t exist you have not managed to install Cellpose correctly).
- Set EnvType to conda
- Only UseGpu should be ticked
- You can now use Cellpose through ImageJ