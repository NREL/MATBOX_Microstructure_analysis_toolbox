function [unstructured_microstructure, microstructure] = generate_cylindrical_2D_channels(unstructured_microstructure, p)

Domain_size = size(unstructured_microstructure); % Dimension
Ttot = Domain_size(3); % Total thickness
t = round(Ttot*p.t); % Thickness of the structured layer
half_w2b = p.half_w2 * p.w; % Half bottom channel width
slope_rad = atan( (p.half_w2-half_w2b)/t ); % Slope angle

% Channel spacing (center to center)
Lcc = 2*p.half_w1 + 2*p.half_w2;
dx = Lcc;
dy = dx/2*tan(deg2rad(60));
dx=max(1,round(dx));
dy=max(1,round(dy));

% First channel
if p.channel_is_corner
    x=1;
    y=1;
else
    x=round(dx/2);
    y=round(dy/2);
end

% Other channels
xs1 = x:dx:Domain_size(2)+1;
xs2 = x-dx/2:dx:Domain_size(2)+1;
xs2(xs2<0)=[];
if xs2(1)==0
    xs2(1)=1;
end
if xs1(end)==Domain_size(2)+1
    xs1(end)=Domain_size(2);
end
if xs2(end)==Domain_size(2)+1
    xs2(end)=Domain_size(2);
end
ys = y:dy:Domain_size(1)+1;
if ys(end)==Domain_size(1)+1
    ys(end)=Domain_size(1);
end
channel_mask = zeros(Domain_size);
for ky=1:1:length(ys)
    y_=ys(ky);
    if mod(ky,2) == 1 
        xs = xs1;
    else
        xs = xs2;
    end
    for kx=1:1:length(xs)
        x_=xs(kx);
        channel_mask(y_,x_,1:t)=1;
    end
end

% Crop
if p.channel_crop
    xmin = min([min(xs1) min(xs2)]);
    xmax = max([max(xs1) max(xs2)]);
    ymin = min(ys);
    ymax = max(ys);
    channel_mask = channel_mask(ymin:ymax,xmin:xmax,:);
    unstructured_microstructure = unstructured_microstructure(ymin:ymax,xmin:xmax,:);
end

% Channel volume
for z=1:1:t
    r = round((t-z)*tan(slope_rad) + half_w2b);
    slice = channel_mask(:,:,z);
    dmap=bwdist(slice);
    slice(dmap<=r)=1;
    channel_mask(:,:,z)=slice;
end

channel_volume_ratio = sum(sum(sum(channel_mask)))/numel(channel_mask)

%% REMOVE PARTICLES
if strcmp(p.apply,'Homogenous medium')
    microstructure = channel_mask;
elseif strcmp(p.apply,'Heterogenous microstructure (cut through particles)')
    microstructure = unstructured_microstructure;
    microstructure(channel_mask==1)=0;
elseif strcmp(p.apply,'Heterogenous microstructure (remove particles)')
    microstructure = unstructured_microstructure;
    particles_id = unique(microstructure); particles_id(particles_id==0)=[]; n_particles = length(particles_id);
    idchannel = find(channel_mask==1);
    particles_id = microstructure(idchannel); particles_id=unique(particles_id);
    particles_id(particles_id==0)=[]; n_particles = length(particles_id);
    for k=1:1:n_particles
        current_id = particles_id(k);
        idx = find(microstructure==current_id);
        check = unique(channel_mask(idx));
        if length(check)==1
            microstructure(idx) = 0; % Remove this particle completely
        else % Particle is both within and outside the channel
            if rand < p.percentremove/100
                microstructure(idx) = 0; % Remove this particle completely
            end
        end
    end
    microstructure(microstructure~=0)=1;
end



