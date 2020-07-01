# MATBOX: Microstructure Analysis Toolbox documentation
**MATLAB open source microstructure analysis code embedded with a graphic user interface**\
*v0.9 06/29/2020*

MATBOX is a MATLAB application for performing various microstructure-related tasks including **microstructure numerical generation**, **image filtering** and **microstructure segmentation**, **microstructure characterization**, result **three-dimensional visualization** and **result correlation**, and **microstructure meshing**. \
MATBOX was originally developed to analyse electrode microstructures for lithium ion batteries; however, the algorithms provided by the toolbox are widely applicable to other heterogeneous materials.

The toolbox provides a user-friendly experience thanks to a **Graphic-User Interface** and requires no coding to be used.

Installation and instructions are detailed in the pdf documentation. \
Run src/Main_menu/Microstructure_analysis_toolbox_mainmenu.m to start the toolbox.

![Main menu and module illustrations](https://github.nrel.gov/fussegli/MATBOX_Microstructure_analysis_toolbox/blob/master/Image.png)

## Authors
Contributions (excluding third-party software):
* Main developer: Francois Usseglio-Viretta (NREL)
* GUI development of the particle generation module: Prehit Patel
* Meshing module code adapted for monolithic mesh: Jeffery Allen (NREL)

## License
This toolbox uses BSD license. NREL Software Record number SWR-20-76. License file is in this folder. Third-party license are available in the third-party licences folder.

## How to cite
If you produce results using the toolbox, or use some or parts of the algorithms contained within the toolbox, please quote them accordingly:
* For any results produced with the toolbox, please quote: F. L. E. Usseglio-Viretta et al., MATBOX: An Open-source Microstructure Analysis Toolbox for microstructure generation, segmentation, characterization, visualization, correlation, and meshing, SoftwareX, in preparation
* If you are calculating tortuosity factor, then please **also** quote: S.J. Cooper, A. Bertei, P.R. Shearing, J.A. Kilner, and N.P. Brandon, TauFactor: An open-source application for calculating tortuosity factors from tomographic data, SoftwareX, Volume 5, 2016, Pages 203-210
* If you are generating mesh, then please **also** quote: Q. Fang and D. A. Boas, Tetrahedral Mesh Generation From Volumetric Binary and Gray-scale Images, Proceedings of IEEE International Symposium on Biomedical Imaging 2009, 2009, Pages 1142-1145

## Acknowledgments
This software was authored by the **National Renewable Energy Laboratory**, operated by **Alliance for Sustainable Energy, LLC**, for the **U.S. Department of Energy (DOE)** under Contract No. DE-AC36-08GO28308. Funding for algorithm development was pro-vided by the U.S. DOE Vehicle Technologies Office’s Computer-Aided Engineering of Batteries (CAEBAT) program (program manager Brian Cunningham). Application of the algorithm for fast-charge analysis was provided by the eXtreme Fast Charge Cell Evaluation of Lithium-Ion Batteries (XCEL) program (program manager Samuel Gillard).
