function [] = Microstructure_generation_stochastic_frominputs(inputs, run_number, savefolder, save_progression, save_verification, scaling_factor, scaling_parameters)
% Stochastic generation algorithm from pre-existing input file
% In all cases, you must run the Microstructure_generation_stochastic app at least once to generate an input file.
% Then, please load this file: inputs=load('<your path>\Inputs.mat');
% - Default use, no scaling: Microstructure_generation_stochastic_frominputs(inputs, run_number, savefolder,save_progression, save_verification)
%   run_number = 2 % How many generation (for reproducibility)
%   savefolder = <where you want to save>. Path must end with '\'.
%   save_progression = true; % true of false
%   save_verification = true; % true of false
% - With scaling: Microstructure_generation_stochastic_frominputs(inputs, run_number, savefolder, save_progression, save_verification, scaling_factor, scaling_parameters)
%   Set run_number, savefolder,  save_progression and save_verification as done above, then:
%   scaling_factor = 2 % float, >0 (>1 upscaling, <1 downscaling)
%   - scaling_parameters.method = 1 % Re-create from partice information (recommended)
%       - scaling_parameters.scale_diameterratio = 1 % float, >0
%       - scaling_parameters.Forceseparation = false % true or false
%       - scaling_parameters.distance_separation = 1 % float, >0 in voxel length (used only if Forceseparation=true)
%   OR
%   - scaling_parameters.method = 2 % Classic up/down scaling
%       - scaling_parameters.scale_particlelabel = true; true (slow) or false (only phase label will be upscaled)
%

if nargin == 5 % Set default parameters
    scaling_factor = 1;
    scaling_parameters = [];
end

% Set saving options
inputs.saveoptions.folder = savefolder;
inputs.saveoptions.save_progression = save_progression;
inputs.saveoptions.save_verification  = save_verification ;
inputs.saveoptions.makevideo = save_progression;
for k_run = 1:1:run_number
    fprintf(['Run #'  num2str(k_run,'%i') '/' num2str(num2str(k_run,'%i'),'%i') ' ongoing... please wait. '])
    inputs.saveoptions.run_number = k_run;

    [microstructure3D, ~, particle_data, outcome, ~] = function_generate_ellipsoid_microstructure(inputs.domain_size,inputs.phase,inputs.tolerance,inputs.forceodddiameter,inputs.overlapping,inputs.acc,inputs.cropparameters,inputs.check_contiguity, inputs.pMask, inputs.stopingconditions, inputs.doverification, inputs.saveoptions);
    fprintf([outcome '\n']);

    if strcmp(outcome,'Success !')

        try
            save([savefolder 'Additionalinfo_run_' num2str(k_run) '.mat'],'microstructure3D','-mat');
        catch ME
            disp 'Additionalinfo_run_ not saved'
            ME.identifier
        end
        try
            save([savefolder 'Additionalinfo_particle_run_' num2str(k_run) '.mat'],'particle_data','-mat');
        catch
            disp 'Additionalinfo_particle_run_'
            ME.identifier
        end
        function_save_tif( uint8(microstructure3D.phase), [savefolder 'Phaselabel_run_' num2str(k_run) '.tif']);

        [microstructure3D.particle_id] = fct_intconvert(microstructure3D.particle_id);
        function_save_tif( microstructure3D.particle_id, [savefolder 'Particlelabel_run_' num2str(k_run) '.tif']);

        % Upscaling
        if scaling_factor~=1
            if scaling_parameters.method==1 % 'Re-create from partice information (recommended)'
                scaling_factor = 1/scaling_factor; % need to harmonize...
                domain_size_rescaled = inputs.domain_size*1/scaling_factor;
                particle_data_scaled = particle_data;
                particle_data_scaled(:,3) = min(round( 1/scaling_factor * particle_data(:,3) + 0.5 ), domain_size_rescaled(1));
                particle_data_scaled(:,4) = min(round( 1/scaling_factor * particle_data(:,4) + 0.5 ), domain_size_rescaled(2));
                particle_data_scaled(:,5) = min(round( 1/scaling_factor * particle_data(:,5) + 0.5 ), domain_size_rescaled(3));
                particle_data_scaled(:,6) = round(particle_data(:,6) * 1/scaling_factor);
                particle_data_scaled(:,7) = round(particle_data(:,7) * 1/scaling_factor);
                particle_data_scaled(:,8) = round(particle_data(:,8) * 1/scaling_factor);
                [microstructure3D_phaselabel_scaled, microstructure3D_particlelabel_scaled] = function_createvolume_fromparticledata(particle_data_scaled,domain_size_rescaled, scaling_parameters.scale_diameterratio, scaling_parameters.Forceseparation, scaling_parameters.distance_separation);
                save([savefolder 'Additionalinfo_particle_scaled_run_' num2str(k_run) '.mat'],'particle_data_scaled','-mat');
            else
                parameters_scaling.scaling_factor = 1/scaling_factor;
                parameters_scaling.label_or_greylevel = 'Label';
                parameters_scaling.background = 0;
                % Scale
                microstructure3D_phaselabel_scaled = function_scaling(microstructure3D.phase,parameters_scaling);
                if scaling_parameters.scale_particlelabel
                    microstructure3D_particlelabel_scaled = function_scaling(microstructure3D.particle_id,parameters_scaling);
                else
                    microstructure3D_particlelabel_scaled = [];
                end
            end
            function_save_tif( uint8(microstructure3D_phaselabel_scaled), [savefolder 'Phaselabel_scaled_run_' num2str(k_run) '.tif']);
            if ~isempty(microstructure3D_particlelabel_scaled)
                [microstructure3D_particlelabel_scaled] = fct_intconvert(microstructure3D_particlelabel_scaled);                
                function_save_tif( microstructure3D_particlelabel_scaled, [savefolder 'Particlelabel_scaled_run_' num2str(k_run) '.tif']);
            end
        end

    end

end

end