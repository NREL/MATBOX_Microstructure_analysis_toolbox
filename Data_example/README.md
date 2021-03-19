# Data example
Several open-source segmented volumes are available directly from this repository.

## From the NREL Battery Microstructures open-source Library:
### In the folder "From computed tomography"
- nmc-1-cal-blackisporepluscbd.tif, a Nickel Manganese Cobalt oxide positive electrode (2-phase domain).
- graphite-5-cal-blackisporepluscbd, a graphite negative electrode (2-phase domain).

![3D view of computed tomography segmented microstructures](https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/blob/master/Data_example/From%20computed%20tomography/Image.png)

Source: https://www.nrel.gov/transportation/microstructure.html \
**More volumes are available from the above link, as well as gray level data set.**

## From the microstructure generation module of this toolbox:
### In the folder "From computed tomography with additives numerically generated"
- nmc-1-cal_generatedCBD_bridge.tif, bridge approach generation (distance based)
- nmc-1-cal_generatedCBD_Energy_w0.001.tif, energy approach generation with morphlogy parameter set to minimum 0.001
- nmc-1-cal_generatedCBD_Energy_w0.999.tif, energy approach generation with morphlogy parameter set to maximum 0.999
NMC electrode with additive (Carbon-binder domain) numerically generated.

![3D view of numerically generated additive phase](https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/blob/master/Data_example/From%20computed%20tomography%20with%20additives%20numerically%20generated/Bridge3D.png)

![3D view of numerically generated additive phase](https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/blob/master/Data_example/From%20computed%20tomography%20with%20additives%20numerically%20generated/Image2D_comparison.png)


### In the folder "Numerically generated"
- Ellipsoids.tif, unisize ellipsoid-based particles domain with a preferential alignment.
- Bi_layer.tif, a bi-layer material with variation of particle size and porosity.
- Bi_layer_poreformer.tif, a bi-layer material with pore former in one layer.
- Two_particle_size.tif, tri-domain material: pore (background), small and large spherical particles.

![3D view of numerically generated microstructures](https://github.com/NREL/MATBOX_Microstructure_analysis_toolbox/blob/master/Data_example/Numerically%20generated/Image.png)

*You can visualize these tif files using the Microstructure visualizatrion module.*
*Alternatively, you can use open source software ImageJ (or Fiji): https://imagej.net/Welcome*
