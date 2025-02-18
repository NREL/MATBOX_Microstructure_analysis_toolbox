clearvars
close all
clc

%% IMPORT
% E.g.:
% folderpath = 'C:\Users\fussegli\Desktop\Example\';
% particle_data = load([folderpath 'TubeInfo_run_1.mat']);
% inputs = load([folderpath 'Inputs.mat']);
% folderpath = ...
% tube_data = ...
% inputs = ...

% folderpath = 'C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\C45_CNR_CNT_comparison\Multiscale\CNR_only\Vf_0p25_voxelsize20p7nm\';
% tube_data = load([folderpath 'TubeInfo_run_1.mat']);
% inputs = load([folderpath 'Inputs.mat']);
% savefolder = 'C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\C45_CNR_CNT_comparison\Multiscale\CNR_only\Vf_0p0551_voxelsize20p7nm\';
% check_contiguity = false;

folderpath = 'C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\C45_CNR_CNT_comparison\Multiscale\SWCNT_only\High_curvature\';
tube_data = load([folderpath 'TubeInfo_run_1.mat']);
inputs = load([folderpath 'Inputs.mat']);
savefolder = 'C:\Users\fussegli\Desktop\SiHPC\Carbon nano tube\C45_CNR_CNT_comparison\Multiscale\SWCNT_only\High_curvature_vf0p0617\';
check_contiguity = true;

%% PARAMTERS
scaling_factors.loc = 1; % >1: upscaling, <1: downscaling
scaling_factors.rad = 1; % >1: upscaling, <1: downscaling
ntube = 0; % Set to 0 to ignore, otherwise pick a positive integer. Re-generation will stop once number of tube is equal to this target.
%volumefraction = 0.0290; % Set to 1 to ignore, otherwise pick a value between 0 and 1. Re-generation will stop once achieved volume fraction is higher or equal to this target.
volumefraction = 0.0617;

%% RE-GENERATE
% Unpack
tubes = tube_data.tubes;
domain_size = inputs.domain_size;
Bezier_pars = inputs.Bezier_pars;
cropp = inputs.cropp;
shape = inputs.p.shape;
[tube_phase, tube_id, tube_skeleton] = function_upscale_tube(tubes,domain_size,scaling_factors,Bezier_pars,cropp,shape,ntube,volumefraction);
nvoxel = numel(tube_phase);
volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel

%% CHECK CONTIGUITY
if check_contiguity
    ids = unique(tube_id);
    ids(ids==0)=[]; % remove background
    find_nocontiguous = false;
    for k=1:1:length(ids)
        BW = zeros(size(tube_phase));
        idx = find(tube_id==ids(k));
        BW(idx)=1;
        L = bwlabeln(BW,26);
        if length(unique(L))>2 % Non-contiguous
            tube_phase(idx)=0;
            tube_id(idx)=0;
            tube_skeleton(tube_skeleton==ids(k)) = 0;
            find_nocontiguous = true;
            fprintf('- Tube #%i is not continguous and has been removed\n',ids(k));            
        end
    end

    if find_nocontiguous
        disp 'Renumerotating...';
        % Recalculate volume fraction
        volume_fraction = sum(sum(sum(tube_phase~=0)))/nvoxel
        % Renumeroate
        old_ids = unique(tube_id);
        old_ids(old_ids==0)=[];
        k_tube = length(old_ids);
        new_ids = 1:1:k_tube;
        tmp = tube_id;
        for kid=1:1:k_tube
            tube_id(tmp == old_ids(kid)) = new_ids(kid);
        end
        clear tmp;
        tmp = tube_skeleton;
        for kid=1:1:k_tube
            tube_skeleton(tmp == old_ids(kid)) = new_ids(kid);
        end
        clear tmp;

        tmp = tubes;
        clear tubes;
        for kid=1:1:k_tube
            radius = tmp(old_ids(kid)).radius;
            points = tmp(old_ids(kid)).points;
            tubes(kid).radius = radius;
            tubes(kid).points = points;
        end

    else
        disp 'All tubes are contiguous';
    end


end

%% FORMAT
tube_phase = uint8(tube_phase);
if max(unique(tube_id))<=255
    tube_id = uint8(tube_id);
    tube_skeleton = uint8(tube_skeleton);
else
    tube_id = uint16(tube_id);
    tube_skeleton = uint16(tube_skeleton);
end

%% VISUALIZE
sum(sum(sum(tube_phase)))/numel(tube_phase)
Microstructure_basic_visualization_interface(tube_phase)
Microstructure_basic_visualization_interface(tube_id)

%% SAVE
function_save_tif( tube_phase, [savefolder 'Tubelabel_run_1.tif']);
function_save_tif( tube_id, [savefolder 'TubeId_run_1.tif']);
function_save_tif( tube_skeleton, [savefolder 'Tubeskeleton_run_1.tif']);
