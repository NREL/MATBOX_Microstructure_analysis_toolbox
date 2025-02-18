clearvars
close all
clc

%% IMPORT

folder = 'C:\Users\fussegli\Desktop\SCP_June24_generationforJulianne\Tube first\15nm\Transveral_anisotropy\';



foldersav = [folder 'Regen\'];
run = 1;

tube_id = function_load_tif([folder 'TubeId_run_' num2str(run) '.tif']);
tubes = load([folder 'TubeInfo_run_' num2str(run) '.mat']);
p = load([folder 'Inputs.mat']);

%% PARAMETERS

%volumefraction = 0.007324;
%volumefraction = 0.008545;
volumefraction = 0.009765;


% Scaling
scaling_factors.rad = 1;
scaling_factors.loc = 1;

% Unpack
tubes = tubes.tubes;
Bezier_pars = p.Bezier_pars;
cropp = p.cropp;
shape = p.p.shape;
n_tube = 0;

%% UPSCALE
domain_size = size(tube_id);
[tube_phase, tube_id, tube_skeleton] = function_upscale_tube(tube_id,tubes,domain_size,scaling_factors,Bezier_pars,cropp,shape,n_tube,volumefraction);

if ~exist(foldersav,'dir')
    mkdir(foldersav);
end 
function_save_tif(uint8(tube_phase),[foldersav 'Tubelabel_run_' num2str(run) '_scaled' num2str(scaling_factors.rad) '_' num2str(volumefraction*100,'%1.3f') 'p.tif']);
function_save_tif(uint8(tube_id),[foldersav 'TubeId_run_' num2str(run) '_scaled' num2str(scaling_factors.rad) '_' num2str(volumefraction*100,'%1.3f') 'p.tif']);

%Microstructure_basic_visualization_interface(tube_phase);