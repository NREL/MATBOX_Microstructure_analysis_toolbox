# MATBOX: Microstructure Analysis Toolbox documentation
**MATLAB open source microstructure analysis code embedded with a graphic user interface**\
*v0.91 March 18, 2021*\
 See release_notes.txt for a change log.

MATBOX is a MATLAB application for performing various microstructure-related tasks including **microstructure numerical generation**, **image filtering** and **microstructure segmentation**, **microstructure characterization**, result **three-dimensional visualization** and **result correlation**, and **microstructure meshing**. \
MATBOX was originally developed to analyse electrode microstructures for lithium ion batteries; however, the algorithms provided by the toolbox are widely applicable to other heterogeneous materials.

The toolbox provides a user-friendly experience thanks to a **Graphic-User Interface** and requires no coding to be used.

Installation and instructions are detailed in the pdf documentation. \
Run src/Main_menu/Main_menu.mlapp to start the toolbox (mlapp extension corresponds to MATLAB app created with app designer).

* MATBOX main menu and illustration of each module:
![Main menu and module illustrations](https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/blob/master/MATBOX.png)

* MATBOX modules inputs/outputs:
The inputs for most modules are stack tiff files (which once imported in MATLAB a 3D array).
![Main menu and module illustrations](https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/blob/master/IO.png)

## Authors
Contributions (excluding third-party software):
* Main developer and documentation writer: Francois Usseglio-Viretta (NREL)
* GUI development of the particle generation module: Prehit Patel
* Meshing module code adapted for monolithic mesh: Jeffery Allen (NREL)
* Contrast correction documentation/examples for adapthisteq: Elizabeth Bernhardt

## How to cite
If you produce results using the toolbox, or use some or parts of the algorithms contained within the toolbox, please quote them accordingly:
* For any results produced with the toolbox, please quote: F. L. E. Usseglio-Viretta et al., MATBOX: An Open-source Microstructure Analysis Toolbox for microstructure generation, segmentation, characterization, visualization, correlation, and meshing, SoftwareX, in preparation
* If you are generating additive phase with the energy criterion method, then please **also** quote: A. N. Mistry, K. Smith, and P. P. Mukherjee, Secondary Phase Stochastics in Lithium-Ion Battery Electrodes, ACS Appl. Mater. Interfaces 10(7) pp. 6317-6326 (2018), https://doi.org/10.1021/acsami.7b17771
* If you are calculating tortuosity factor, then please **also** quote: S.J. Cooper, A. Bertei, P.R. Shearing, J.A. Kilner, and N.P. Brandon, TauFactor: An open-source application for calculating tortuosity factors from tomographic data, SoftwareX, Volume 5, 2016, Pages 203-210
* If you are generating mesh, then please **also** quote: Q. Fang and D. A. Boas, Tetrahedral Mesh Generation From Volumetric Binary and Gray-scale Images, Proceedings of IEEE International Symposium on Biomedical Imaging 2009, 2009, Pages 1142-1145

## What's next?
- Publish in a journal article is next priority.
- Next release will focus on properly finish the generation module, dust off the the characterization module, and add a template to the correlation module to make it faster for the user.
If you have suggestions, please go to the discussion section of this repo.

## License
This toolbox uses BSD license. NREL Software Record number SWR-20-76. License file is in this folder. Third-party license are available in the third-party licences folder.

## Acknowledgments
This software was authored by the **National Renewable Energy Laboratory**, operated by **Alliance for Sustainable Energy, LLC**, for the **U.S. Department of Energy (DOE)** under Contract No. DE-AC36-08GO28308. Funding for algorithm development was pro-vided by the U.S. DOE Vehicle Technologies Office’s Computer-Aided Engineering of Batteries (CAEBAT) program (program manager Brian Cunningham). Application of the algorithm for fast-charge analysis was provided by the eXtreme Fast Charge Cell Evaluation of Lithium-Ion Batteries (XCEL) program (program manager Samuel Gillard).
