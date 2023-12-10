# ClinSig-V3
This repository contains version 3 of the "ClinSig" program, a Shiny application developed in RStudio that visualizes clinical variants from the NCBI ClinVar database in a lollipop plot format according to their clinical significance.

How to Use?
1. Query the NCBI ClinVar database by inputting the gene name to obtain clinical significance information.
2. Download the search output in text format using the "Downloads" button.
3. In the ClinSig_V3.R script, load the required packages/libraries and launch the application. Then on the file selection window, choose the text file downloaded in the previous step.

Version 3 Enhancements:
1- The two separate files for "manipulation and plotting" were removed. All the code was integrated into a single R script file.
2- A screen was added that allows the automatic loading of a text file when the program is run.
3- The color for each clinical significance information can now be customized.
4- Some errors encountered during plotting were resolved.
5- The ability to download plot images in "SVG, PNG, TIFF, and JPG" formats was added.

The key changes in Version 3 include:
1- Simplified the workflow by removing separate files for data handling and plotting.
2- Added a user-friendly front-end to load input files automatically.
3- Increased flexibility by making clinical significance colors user-defined.
4- Debugged issues causing errors during plot generation.
5- Enhanced usability with options to save plots in different image formats.

These improvements aim to streamline the process and provide more customization options to the end user.
