function [Mseg,newtype,res] = fct_instancesegmentation2(M,p)

sz = size(M);
dimension = length(sz);
newtype = 'Segmented (instance)';
res = [];

BW = zeros(sz);
BW(M == p.labeltoinstance)=1;
idlargestaxe = find(sz==max(sz));
idlargestaxe = idlargestaxe(1);

if p.perparts && strcmp(p.perparts_choice,'Auto')
    nvoxel = prod(sz);
    a = 1/(p.perpart_millions*1e6);
    b = 0;
    n = round(a*nvoxel+b);
    if n<2
        p.perparts = false;
    end
end


if p.perparts
    if strcmp(p.perparts_choice,'Manual')
        p.nparts = p.perpart_nparts;
    elseif strcmp(p.perparts_choice,'Auto')
        p.nparts = n;
    end

    fprintf('By part segmentation...\n')
    bounds = round(linspace(1,max(sz),p.nparts + p.nparts-1 + 2));
    bounds_independant = zeros(p.nparts,2);
    bounds_independant(1,1) = 1; bounds_independant(1,2) = bounds(3);
    for k=2:1:p.nparts
        bounds_independant(k,1) = bounds_independant(k-1,2) + 1;
        bounds_independant(k,2) = bounds(2*k+1);
    end

    bounds_overlapping = zeros(p.nparts-1,2);
    if strcmp(p.perpart_overlappinglength_choice,'Same as independant parts')
        for k=1:1:p.nparts-1
            bounds_overlapping(k,1) = bounds(2*k);
            bounds_overlapping(k,2) = bounds(2*k+2);
        end
    elseif strcmp(p.perpart_overlappinglength_choice,'Manual')
        l = bounds(3)-bounds(1)+1;
        if p.perpart_overlappinglength >= l
            for k=1:1:p.nparts-1
                bounds_overlapping(k,1) = bounds(2*k);
                bounds_overlapping(k,2) = bounds(2*k+2);
            end
        else
            for k=1:1:p.nparts-1
                center = bounds(2*k+1);
                bounds_overlapping(k,1) = max([1, center - round(p.perpart_overlappinglength/2)]);
                bounds_overlapping(k,2) = min([max(sz), center + round(p.perpart_overlappinglength/2)]);
            end
        end 
    end

    for k = 1:1:p.nparts
        fprintf('- Independent part,%i/%i\n',k,p.nparts)
        if dimension == 2
            if idlargestaxe==1
                sub = BW(bounds_independant(k,1):bounds_independant(k,2),:);
            elseif idlargestaxe==2
                sub = BW(:,bounds_independant(k,1):bounds_independant(k,2));
            end
        else
            if idlargestaxe==1
                sub = BW(bounds_independant(k,1):bounds_independant(k,2),:,:);
            elseif idlargestaxe==2
                sub = BW(:,bounds_independant(k,1):bounds_independant(k,2),:);
            elseif idlargestaxe==3
                sub = BW(:,:,bounds_independant(k,1):bounds_independant(k,2));
            end
        end
        BW_initial = sub;

        if p.erosion || p.downscaling
            % Check downscaling factor
            if k==1 && p.downscaling && (strcmp(p.downscaling_choice,'Auto') || strcmp(p.downscaling_choice,'Auto (and round)'))
                n0 = numel(sub);
                n1 = p.downscaling_millions*1e6;
                p.downscaling_factor = max([1, (n0/n1)^(1/dimension)]);
                % if strcmp(p.downscaling_choice,'Auto (and round)')
                %     p.downscaling_factor = round(p.downscaling_factor);
                % else
                %     p.downscaling_factor = round(p.downscaling_factor,2);
                % end
                p.downscaling_factor = round(p.downscaling_factor);
            end

            if k==1 && p.downscaling && p.downscaling_factor~=1
                fprintf('  Downscaling factor: %1.3f\n',p.downscaling_factor)
            end
            
            [sub, BW_beforeerosiondownscaling] = Segm_erosion_and_downscaling(sub,p);
        end

        if strcmp(p.instance_method,'Pseudo-Coulomb Repulsive Field (PCRF)')
            [Ind(k).Mseg] = Segm_PCRF(sub,p);
        elseif strcmp(p.instance_method,'Watershed (custom)')
            [Ind(k).Mseg] = Segm_customwatershed(sub,p);
        elseif strcmp(p.instance_method,'Watershed (built-in)')
            [Ind(k).Mseg] = Segm_builtinwatershed(sub,p);
        end
        
        if p.downscaling && p.downscaling_factor~=1 % Upscale
            [Ind(k).Mseg] = Segm_reupscaling(Ind(k).Mseg,BW_initial,BW_beforeerosiondownscaling,p);
        end
    end

    for k = 1:1:p.nparts-1
        fprintf('- Overlapping part,%i/%i\n',k,p.nparts-1)
        if dimension == 2
            if idlargestaxe==1
                sub = BW(bounds_overlapping(k,1):bounds_overlapping(k,2),:);
            elseif idlargestaxe==2
                sub = BW(:,bounds_overlapping(k,1):bounds_overlapping(k,2));
            end
        else
            if idlargestaxe==1
                sub = BW(bounds_overlapping(k,1):bounds_overlapping(k,2),:,:);
            elseif idlargestaxe==2
                sub = BW(:,bounds_overlapping(k,1):bounds_overlapping(k,2),:);
            elseif idlargestaxe==3
                sub = BW(:,:,bounds_overlapping(k,1):bounds_overlapping(k,2));
            end
        end
        BW_initial = sub;

        if p.erosion || p.downscaling
            [sub, BW_beforeerosiondownscaling] = Segm_erosion_and_downscaling(sub,p);
        end        

        if strcmp(p.instance_method,'Pseudo-Coulomb Repulsive Field (PCRF)')
            [Over(k).Mseg] = Segm_PCRF(sub,p);
        elseif strcmp(p.instance_method,'Watershed (custom)')
            [Over(k).Mseg] = Segm_customwatershed(sub,p);
        elseif strcmp(p.instance_method,'Watershed (built-in)')
            [Over(k).Mseg] = Segm_builtinwatershed(sub,p);
        end

        if p.downscaling && p.downscaling_factor~=1 % Upscale
            [Over(k).Mseg] = Segm_reupscaling(Over(k).Mseg,BW_initial,BW_beforeerosiondownscaling,p);
        end
    end

    % Combine
    Mseg = Segm_InstanceRecombineParts(idlargestaxe,bounds_independant,bounds_overlapping,Ind,Over);

    % Avoid duplicates (can happen if number of parts is too high)
    Mseg = uint32(Mseg).*uint32(BW);

else
    BW_initial = BW;

    if p.erosion || p.downscaling

        % Check downscaling factor
        if p.downscaling && (strcmp(p.downscaling_choice,'Auto') || strcmp(p.downscaling_choice,'Auto (and round)'))
            n0 = numel(BW);
            n1 = p.downscaling_millions*1e6;
            p.downscaling_factor = max([1, (n0/n1)^(1/dimension)]);
            % if strcmp(p.downscaling_choice,'Auto (and round)')
            %     p.downscaling_factor = round(p.downscaling_factor);
            % else
            %     p.downscaling_factor = round(p.downscaling_factor,2);
            % end
            p.downscaling_factor = round(p.downscaling_factor);
            if p.downscaling_factor~=1
                fprintf('  Downscaling factor: %1.3f\n',p.downscaling_factor)
            end
        end

        [BW, BW_beforeerosiondownscaling] = Segm_erosion_and_downscaling(BW,p);
    end

    if strcmp(p.instance_method,'Pseudo-Coulomb Repulsive Field (PCRF)')
        [Mseg] = Segm_PCRF(BW,p);
    elseif strcmp(p.instance_method,'Watershed (custom)')
        [Mseg] = Segm_customwatershed(BW,p);
    elseif strcmp(p.instance_method,'Watershed (built-in)')
        [Mseg] = Segm_builtinwatershed(BW,p);
    end

    if p.downscaling && p.downscaling_factor~=1 % Upscale
        [Mseg] = Segm_reupscaling(Mseg,BW_initial,BW_beforeerosiondownscaling,p);
    end


end

if p.randomize
    pp.basedondistance = false;
    pp.dist = [];
    pp.keepunchangedsomelabels = true;
    pp.cst = '0';
    [Mseg,~,~] = fct_randomize(Mseg,pp);
else
    % Convert to uint8 or uint16
    [Mseg] = fct_intconvert(Mseg); % Already done in fct_randomize
end







% 
% if p.reapply_perinstance
%     sav_downscaling_factor = p.downscaling_factor;
%     % fix error induced by downscaling-upscaling
%     sz = size(Mseg);
% 
%     Particle_Label_full = zeros(sz);
%     unis = unique(Mseg);
%     unis(unis==0)=[];
% 
%     p.display_downscalingerosion = false;
%     p.downscaling_factor = 1; % No downscaling
%     p.erosion_distance =0; % No erosion
%     p.percluster = 0;
%     p.details_convergence = false;
%     p.visualize_2D = false;
%     p.max_dist = round(p.max_dist * sav_downscaling_factor);
% 
%     for k=1:1:length(unis)
%         k/length(unis)
%         idx = find(Mseg==unis(k));
%         [IX,IY,IZ]=ind2sub(size(Mseg),idx);
%         x_minlabel = min(IX); x_maxlabel = max(IX);
%         y_minlabel = min(IY); y_maxlabel = max(IY);
%         if dimension==2
%             z_minlabel = 1; z_maxlabel = 1;
%             BW = Mseg(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel);
%         else
%             z_minlabel = min(IZ); z_maxlabel = max(IZ);
%             BW = Mseg(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel,z_minlabel:z_maxlabel);
%         end
%         BW(BW~=unis(k))=0;
%         BW(BW~=0)=1;
% 
%         if strcmp(p.instance_method,'Pseudo-Coulomb Repulsive Field (PCRF)')
%             [Particle_Label_tmp] = Segm_PCRF(BW,p);
%         elseif strcmp(p.instance_method,'Watershed (custom)')
%             [Particle_Label_tmp] = Segm_customwatershed(BW,p);
%         elseif strcmp(p.instance_method,'Watershed (built-in)')
%             [Particle_Label_tmp] = Segm_builtinwatershed(BW,p);            
%         end
% 
%         m=max(max(max(Particle_Label_full)));
%         Particle_Label_tmp = Particle_Label_tmp + m;
%         Particle_Label_tmp(Particle_Label_tmp==m)=0;
%         Particle_Label_full(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel,z_minlabel:z_maxlabel) = Particle_Label_full(x_minlabel:x_maxlabel,y_minlabel:y_maxlabel,z_minlabel:z_maxlabel) + Particle_Label_tmp;
%     end
% 
%     % Check C-PSD<D-PSD
%     chess_pattern_removal = false;
%     % Find index for the phase
%     index_phase = find(M==1);
%     % Get all coordinates
%     [P1,P2,P3] = ind2sub(sz,index_phase);
%     %[Particle_Label_full] = Function_correct_DPSD_identification(M, Particle_Label_full, p.cpsd_refining, p.details_convergence, chess_pattern_removal, P1, P2, P3, p.visualize_2D);
%     Mseg=Particle_Label_full;
% end


end