function [unstructured_microstructure, microstructure] = generate_rectangular_1D_channels(unstructured_microstructure, p)

p.thickness_axis = 3;
beta = 75;
Domain_size = size(unstructured_microstructure); % Dimension
Ttot = Domain_size(p.thickness_axis); % Total thickness
t = round(Ttot*p.t); % Thickness of the structured layer

if beta>=0
    half_w2b = p.half_w2 - t/tan(deg2rad(beta));
    if half_w2b>0
        alpha  = 90-beta;
        slope_rad = deg2rad(alpha);
    else
        half_w2b = p.half_w2 * p.w; % Half bottom channel width
        slope_rad = atan( (p.half_w2-half_w2b)/t ); % Slope angle
    end
else
    half_w2b = p.half_w2 * p.w; % Half bottom channel width
    slope_rad = atan( (p.half_w2-half_w2b)/t ); % Slope angle
end

%% INITIALIZE CROSS SECTION
nx = Domain_size(p.thickness_axis);
if p.thickness_axis==1
    if p.inplane_axis==2
        ny = Domain_size(3);
    elseif p.inplane_axis==3
        ny = Domain_size(2);
    end
elseif p.thickness_axis==2
    if p.inplane_axis==1
        ny = Domain_size(3);
    elseif p.inplane_axis==3
        ny = Domain_size(1);
    end  
elseif p.thickness_axis==3
    if p.inplane_axis==1
        ny = Domain_size(2);
    elseif p.inplane_axis==2
        ny = Domain_size(1);
    end   
end

%% CROSS SECTION PERIODIC PATTERN
periodic_pattern = zeros(nx,p.half_w2+p.half_w1); % Initialize
if t<Ttot % Pad layer if any
    periodic_pattern(1:Ttot-t,:)=1;
end
if slope_rad==0
    if p.crop_porous1_channel0
        periodic_pattern(Ttot-t+1:end,1:p.half_w1)=1;
    else
        periodic_pattern(Ttot-t+1:end,p.half_w2+1:end)=1;
    end
else
    for x = Ttot-t+1:1:Ttot
        local_x = x - (Ttot-t);
        local_w = round(tan(slope_rad)*local_x);
        if p.crop_porous1_channel0
            periodic_pattern(x,1:p.half_w1+p.half_w2-half_w2b-local_w+1)=1;
        else
            periodic_pattern(x,half_w2b+local_w+1:end)=1;
        end
    end
end
flipped_periodic_pattern = flip(periodic_pattern,2);

%% SECTION MASK FOR WHOLE IN-PLANE AXIS
[~,ny_pattern] = size(periodic_pattern);
n_complete_pattern = floor(ny/ny_pattern);
if p.crop_periodicity
    section = zeros(nx,n_complete_pattern*ny_pattern);
else
    section = zeros(nx,ny);
end
is_flipped=false;
y0 = 1;
for k=1:1:n_complete_pattern
    y1 = y0 + ny_pattern-1;
    if is_flipped
        section(:,y0:y1) = flipped_periodic_pattern;
    else
        section(:,y0:y1) = periodic_pattern;
    end
    is_flipped = ~is_flipped;
    y0 = y1+1;
end
if ~p.crop_periodicity && ny/ny_pattern>n_complete_pattern
    r = ny - n_complete_pattern*ny_pattern;
    if is_flipped
        section(:,y0:end) = flipped_periodic_pattern(:,1:r);
    else
        section(:,y0:end) = periodic_pattern(:,1:r);
    end
end
section_size = size(section);

%% 3D MASK
channel_mask = zeros(Domain_size);

if p.thickness_axis==1
    if p.inplane_axis==2
        foo=1;
    elseif p.inplane_axis==3
        foo=1;
    end
elseif p.thickness_axis==2
    if p.inplane_axis==1
        foo=1;
    elseif p.inplane_axis==3
        foo=1;
    end  
elseif p.thickness_axis==3
    if p.inplane_axis==1
        tmp = unstructured_microstructure(:,1:section_size(2),:); unstructured_microstructure=tmp; clear tmp;
        Domain_size = size(unstructured_microstructure);
        channel_mask = zeros(Domain_size);
        for k=1:1:Domain_size(1)
            channel_mask(k,:,:) = section';
        end
    elseif p.inplane_axis==2
        foo=1;
    end   
end

%% REMOVE PARTICLES
if strcmp(p.apply,'Homogenous medium')
    microstructure = channel_mask;
elseif strcmp(p.apply,'Heterogenous microstructure (cut through particles)')
    microstructure = unstructured_microstructure;
    microstructure(channel_mask==0)=0;
elseif strcmp(p.apply,'Heterogenous microstructure (remove particles)')
    microstructure = unstructured_microstructure;
    %microstructure = bwlabeln(microstructure,6);
    particles_id = unique(microstructure); particles_id(particles_id==0)=[]; n_particles = length(particles_id);
    idchannel = find(channel_mask==0);
    particles_id = microstructure(idchannel); particles_id=unique(particles_id);
    particles_id(particles_id==0)=[]; n_particles = length(particles_id);
    for k=1:1:n_particles
        current_id = particles_id(k);
        idx = find(microstructure==current_id);
        Id_on_mask = channel_mask(idx);
        check = unique(Id_on_mask);
        if length(check)==1
            microstructure(idx) = 0; % Remove this particle completely
        else % Particle is both within and outside the channel
            vol_particle = length(idx);
            vol_particle_withinchannel = sum(sum(sum(Id_on_mask==0)));
            percent_withinchannel = 100*vol_particle_withinchannel/vol_particle;
            if percent_withinchannel>=p.remove_threshold
                microstructure(idx) = 0; % Remove this particle completely
            else
                if rand <= percent_withinchannel/100
                    microstructure(idx) = 0; % Remove this particle completely
                end
            end
        end
    end
    microstructure(microstructure~=0)=1;
end

end

