clearvars
%close all
clc

%% IMPORT
% E.g.:
% folderpath = 'C:\Users\fussegli\Desktop\Example\';
% Phaselabel = function_load_tif([folderpath 'phaselabel_run_1.tif']); % Only for visualization
% particle_data = load([folderpath 'Additionalinfo_particle_run_1.mat']);
% inputs = load([folderpath 'Inputs.mat']);
% folderpath = ...
% Phaselabel = ...
% particle_data = ...
% inputs = ...

folderpath = 'C:\Users\fussegli\Desktop\Bug Anudeep\mesh_NMC_LIB_ver6\mesh_NMC_LIB_ver6\';
%Phaselabel = function_load_tif([folderpath 'phaselabel_run_1.tif']); % Only for visualization
particle_data = load([folderpath 'Additionalinfo_particle_run_1.mat']);
inputs = load([folderpath 'Inputs.mat']);

%% PARAMTERS
scaling_factor = 4; % >1: upscaling, <1: downscaling
scale_diameterratio = 1.0; % Upscaling change slighly particle volume. You can scale particle diameter to correct it, or use porosity_target
Forceseparation = false; % True to separate particle (isolated particles)
Distanceseparation = 1; % If Forceseparation true, separation width in voxel
porosity_target = 0.0; % Set to 0 to ignore, otherwise pick a value between 0 and 1. Re-generation will stop once achieved porosity is less or equal to this target.
shuffleparticles = false; % If porosity target > 0, you can choose to shuffle particle order. Avoid discriminating the last pass/phase, and allow for some variations.

%% RE-GENERATE
% Unpack
particle_data = particle_data.particle_data;
domain_size = inputs.domain_size;
%forceodddiameter = inputs.forceodddiameter;
forceodddiameter = false;
cropparameters = zeros(3,10);

domain_size_rescaled = domain_size*scaling_factor;
particle_data_scaled = particle_data;
particle_data_scaled(:,3) = min(round( scaling_factor * particle_data(:,3) + 0.5 ), domain_size_rescaled(1));
particle_data_scaled(:,4) = min(round( scaling_factor * particle_data(:,4) + 0.5 ), domain_size_rescaled(2));
particle_data_scaled(:,5) = min(round( scaling_factor * particle_data(:,5) + 0.5 ), domain_size_rescaled(3));
particle_data_scaled(:,6) = round(particle_data(:,6) * scaling_factor);
particle_data_scaled(:,7) = round(particle_data(:,7) * scaling_factor);
particle_data_scaled(:,8) = round(particle_data(:,8) * scaling_factor);
% Force odd diameter
if forceodddiameter
    id_iseven = find(rem(particle_data_scaled(:,6), 2) == 0);
    particle_data_scaled(id_iseven,6) = particle_data_scaled(id_iseven,6)+1;
    id_iseven = find(rem(particle_data_scaled(:,7), 2) == 0);
    particle_data_scaled(id_iseven,7) = particle_data_scaled(id_iseven,7)+1;
    id_iseven = find(rem(particle_data_scaled(:,8), 2) == 0);
    particle_data_scaled(id_iseven,8) = particle_data_scaled(id_iseven,8)+1;
end
cropparameters(:,1) = round(cropparameters(:,1)*scaling_factor);
[Phaselabel_scaled, Particlelabel_scaled, overlapping_stat_scaled] = function_createvolume_fromparticledata(particle_data_scaled,domain_size_rescaled, cropparameters, scale_diameterratio, Forceseparation, Distanceseparation, porosity_target, shuffleparticles);

sum(sum(sum(Phaselabel_scaled==0)))/numel(Phaselabel_scaled)

%% VISUALIZE
Microstructure_comparison_visualization_interface(Phaselabel,Phaselabel_scaled);
Microstructure_basic_visualization_interface(Phaselabel_scaled)

%function_save_tif(uint8(Phaselabel_scaled),'C:\Users\fussegli\Documents\GitHub\MATBOX_Microstructure_analysis_toolbox\Nottoshare\Tutorial\RVE\Generated_volume\Macropore\Macropore.tif');
%function_save_tif(uint16(Particlelabel_scaled),'C:\Users\fussegli\Documents\GitHub\MATBOX_Microstructure_analysis_toolbox\Nottoshare\Tutorial\RVE\Generated_volume\Macropore\Macropore_particleid.tif');

Function_Volume_fractions(Phaselabel_scaled, 0.333, 'um')

% function_save_tif(uint8(Phaselabel_scaled),'C:\Users\fussegli\Documents\GitHub\MATBOX_Microstructure_analysis_toolbox\Nottoshare\Tutorial\RVE\Linear_variation_noRVE_inFOV\Slope.tif');
% function_save_tif(uint16(Particlelabel_scaled),'C:\Users\fussegli\Documents\GitHub\MATBOX_Microstructure_analysis_toolbox\Nottoshare\Tutorial\RVE\Linear_variation_noRVE_inFOV\Slope_particleid.tif');
% 
