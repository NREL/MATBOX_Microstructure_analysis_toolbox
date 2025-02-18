function [Mseg,newtype,res] = fct_instancesegmentation(M,p)

if p.downscaling_factor>1 % To keep dimension coherent
    par.autoreshape = true;
    par.autoset_tomatch_downscaling = true;
    par.downscaling_factor = p.downscaling_factor;
    par.DeltaROI = false;
    [M,~,~] = fct_Crop(M,par);
end

sz = size(M);
dimension = length(sz);
newtype = 'Segmented (instance)';
res = [];

labels = unique(M);
BW = zeros(sz);
BW(M == p.labeltoinstance)=1;

if p.perslice && dimension == 3
    % Mseg = zeros(sz);
    % for z=1:1:sz(3)
    %     fprintf(['slice: ' num2str(z) '\n']);     
    %     slice_ = BW(:,:,z);
    %     if strcmp(p.instance_method,'Pseudo-Coulomb Repulsive Field (PCRF)')
    %         [~, slice_instance, ~, ~] = Function_Discrete_particle_size_PCRF_algorithm_v5(slice_,p);
    %     elseif strcmp(p.instance_method,'Watershed (custom)')
    %         [~,slice_instance] = Function_Discrete_particle_size_watershed_immersion_algoritv2(slice_,p);
    %     end
    % 
    %     pp.choice = 'Separate all non-background labels from each other';
    %     pp.sep_allwidth = 2;
    %     pp.background = 0;
    %     [~,~,idx_separation] = fct_separate_labels(slice_instance,pp);
    %     slice_(idx_separation) = 0;
    %     Mseg(:,:,z) = slice_;
    % 
    % end
else
    if strcmp(p.instance_method,'Pseudo-Coulomb Repulsive Field (PCRF)')
        [~, Mseg, ~, ~] = Function_Discrete_particle_size_PCRF_algorithm_v5(BW,p);
    elseif strcmp(p.instance_method,'Watershed (custom)')
        [~,Mseg] = Function_Discrete_particle_size_watershed_immersion_algoritv2(BW,p);
    elseif strcmp(p.instance_method,'Watershed (built-in)')
        [Mseg] = Segm_instance_watershed_builtin(BW,p);        
    end
end

if p.reapply_perinstance
    sav_downscaling_factor = p.downscaling_factor;
    % fix error induced by downscaling-upscaling
    sz = size(Mseg);

    Particle_Label_full = zeros(sz);
    unis = unique(Mseg);
    unis(unis==0)=[];

    p.display_downscalingerosion = false;
    p.downscaling_factor = 1; % No downscaling
    p.erosion_distance =0; % No erosion
    p.percluster = 0;
    p.details_convergence = false;
    p.visualize_2D = false;
    p.max_dist = round(p.max_dist * sav_downscaling_factor);
    p.max_dist = 40;

    for k=1:1:length(unis)
        k/length(unis)
        idx = find(Mseg==unis(k));
        [IX,IY,IZ]=ind2sub(size(Mseg),idx);
        x_minlabel = min(IX); x_maxlabel = max(IX);
        y_minlabel = min(IY); y_maxlabel = max(IY);
        if dimension==2
            z_minlabel = 1; z_maxlabel = 1;
            BW = Mseg(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel);
        else
            z_minlabel = min(IZ); z_maxlabel = max(IZ);
            BW = Mseg(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel,z_minlabel:z_maxlabel);
        end
        BW(BW~=unis(k))=0;
        BW(BW~=0)=1;

        if strcmp(p.instance_method,'Pseudo-Coulomb Repulsive Field (PCRF)')
            [~, Particle_Label_tmp, ~, ~] = Function_Discrete_particle_size_PCRF_algorithm_v5(BW,p);
        elseif strcmp(p.instance_method,'Watershed (custom)')
            [~,Particle_Label_tmp] = Function_Discrete_particle_size_watershed_immersion_algoritv2(BW,p);
        elseif strcmp(p.instance_method,'Watershed (built-in)')
            [~,Particle_Label_tmp] = Segm_instance_watershed_builtin(BW,p);            
        end

        m=max(max(max(Particle_Label_full)));
        Particle_Label_tmp = Particle_Label_tmp + m;
        Particle_Label_tmp(Particle_Label_tmp==m)=0;
        Particle_Label_full(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel,z_minlabel:z_maxlabel) = Particle_Label_full(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel,z_minlabel:z_maxlabel) + Particle_Label_tmp;
    end

    % Check C-PSD<D-PSD
    chess_pattern_removal = false;
    % Find index for the phase
    index_phase = find(M==1);
    % Get all coordinates
    [P1,P2,P3] = ind2sub(sz,index_phase);
    [Particle_Label_full] = Function_correct_DPSD_identification(M, Particle_Label_full, p.cpsd_refining, p.details_convergence, chess_pattern_removal, P1, P2, P3, p.visualize_2D);
    Mseg=Particle_Label_full;
end

if p.randomize
    [Mseg,~] = randomize_labels(Mseg);
end

% Convert to uint8 or uint16
[Mseg] = fct_intconvert(Mseg);

end