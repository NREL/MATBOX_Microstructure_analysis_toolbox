function [tube_phase, tube_id, tube_skeleton,k_tube] = function_upscale_tube_v2(tube_id,tubes,domain_size,scaling_factors,Bezier_pars,cropp,shape,n_tube,volumefraction)

nid = length(unique(tube_id))-1;

% (Up)downscale image
parameters_scaling.scaling_factor = 1/scaling_factors.loc;
parameters_scaling.label_or_greylevel = 'Label';
parameters_scaling.background = 0;
% Scale
tube_id_imgscaled = function_scaling(tube_id,parameters_scaling);
[D,IDX] = bwdist(tube_id_imgscaled);
tube_id_imgscaled_dilated = tube_id_imgscaled;
idx=find(D<=2);
tube_id_imgscaled_dilated(idx) = tube_id_imgscaled(IDX(idx));

new_domain_size = round(domain_size*scaling_factors.loc);
tube_phase = zeros(new_domain_size);
tube_id = zeros(new_domain_size);
tube_skeleton = zeros(new_domain_size);

Nq = max(new_domain_size)*10; % Set high to avoid non contiguity

if n_tube==0
    n_tube = length(tubes);
end
nvoxel = numel(tube_phase);
for k_tube = 1:1:n_tube
    current_volumefraction = sum(sum(sum(tube_phase==1)))/nvoxel;
    if current_volumefraction>=volumefraction
        break
    end
    points = tubes(k_tube).points;
    n_points = length(points);
    new_points = cell(1,n_points);
    for k_point = 1:1:n_points
        tmp = round(cell2mat(points(k_point))*scaling_factors.loc);
        for dir=1:1:3
            if tmp(dir)<1
                tmp(dir)=1;
            end
            if tmp(dir)>new_domain_size(dir)
                tmp(dir)=new_domain_size(dir);
            end
        end
        new_points(1,k_point) = {tmp};
    end

    dinit = tubes(k_tube).radius + 1 + tubes(k_tube).radius;
    dscal = dinit * scaling_factors.rad;
    new_radius = max([round(dscal-1)/2, 1]) ;
    %new_radius = tubes(k_tube).radius * scaling_factors.rad; 

    Q =[0 0 0;2 2 2];
    iterprecision = 0;
    while max(abs(diff( Q(:,1) ))) >1 || max(abs(diff( Q(:,2) ))) >1 || max(abs(diff( Q(:,3) ))) >1
        iterprecision = iterprecision+1;
        % % CALL THIRD PARTY CODE
        Q = custom_hobbysplines(new_points,'debug',false,'cycle',Bezier_pars.cycle,'tension',Bezier_pars.tension,'Nq',round(Nq*iterprecision));
        % % END OF THIRD PARTY CODE
        Q = round(Q);
    end

    for k=1:1:3
        idx = find(Q(:,k)>new_domain_size(k));
        if ~isempty(idx)
            Q(idx,k)=new_domain_size(k);
        end
        idx = find(Q(:,k)<1);
        if ~isempty(idx)
            Q(idx,k)=1;
        end
    end

    idQ = sub2ind(new_domain_size,Q(:,1),Q(:,2),Q(:,3));
    idQ(isnan(idQ))=[];

    tube_skeleton(idQ) = k_tube;
    tmpmask = zeros(new_domain_size);
    tmpmask(idQ)=1;
    if new_radius>0
        dmap = bwdist(tmpmask,shape);
        idx = dmap<=new_radius;
    else
        idx = idQ;
    end
    tube_phase(idx)=1;
    tube_id(idx)=k_tube;
end


% Very rarely, the tiny change in control point location will change the Bezier curve significantly resuling in overlapping
% with other tube or with mask, if any. So here we check if the upscaled tube is within the image-based upscaled-dilated tube.
% tube_id_sav = tube_id;
for kid = 1:1:nid
    imgscaled = zeros(new_domain_size);
    imgregen = zeros(new_domain_size);
    imgscaled(tube_id_imgscaled_dilated==kid)=1;
    idx_regen = tube_id==kid;
    imgregen(idx_regen)=1;
    res = sum(sum(sum(imgscaled .* imgregen)))  / sum(sum(sum(imgregen))); % find 0,1 and 1,0
    if res<0.75
        warning('Tube #%i upscaling modified its Bezier curve, so instead its image upscaled has been used',kid)
        % Replace with image rescaled
        tube_phase(idx_regen) = 0;
        tube_id(idx_regen) = 0;
        idx_scaled = tube_id_imgscaled==kid;
        tube_phase(idx_scaled) = 1;
        tube_id(idx_scaled) = kid;        
    end
end
% Microstructure_comparison_visualization_interface(tube_id_imgscaled,tube_id_sav,tube_id)

tube_phase(tube_phase~=0)=1;

% Crop
if sum(cropp.cut)>0
    l1 = cropp.cut(1)*scaling_factors.loc;
    l2 = cropp.cut(2)*scaling_factors.loc;
    l3 = cropp.cut(3)*scaling_factors.loc;
    tube_phase = tube_phase(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
    tube_id = tube_id(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
    tube_skeleton = tube_skeleton(l1:end-l1+1, l2:end-l2+1, l3:end-l3+1);
end

% Sometimes there is tiny overlap
tube_phase(tube_phase>1)=1;

end





