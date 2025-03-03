function [BW,BW_spheres,nodes_centroids_radius_id,chord_centroids_radius_id, sphere_connectivity] = function_generate_web_microstructure(centroids,p)

[n_particle,dimension] = size(centroids);

%% RE-CENTER AND FOV
x0 = min(centroids(:,1)); x1 = max(centroids(:,1));
y0 = min(centroids(:,2)); y1 = max(centroids(:,2));
if dimension == 3
    z0 = min(centroids(:,3)); z1 = max(centroids(:,3));
end

if p.recenter
    centroids(:,1) = centroids(:,1) - x0 + 1 ;
    centroids(:,2) = centroids(:,2) - y0 + 1 ;
    if dimension == 2
        sz = [x1-x0+1 y1-y0+1];
    else
        sz = [x1-x0+1 y1-y0+1 z1-z0+1];
        centroids(:,3) = centroids(:,3) - z0 + 1 ;
    end
else
    if dimension == 2
        sz = [x1 y1];
    else
        sz = [x1 y1 z1];
    end
end

if p.extend_FOV~=0
    centroids = centroids+p.extend_FOV;
    sz = sz + 2*p.extend_FOV;
end

%% WEB SKELETON

% Pick coordination target
[x_coordinations_target, ~, cdf_coordinations_target] = distributions(p.equiprobability_coordination, p.mean_coordination, p.sigma_coordination, p.maxdeviation_coordination);
coordinations_target = pickfromdistriubtion(p.maxdeviation_coordination,p.equiprobability_coordination,p.mean_coordination,cdf_coordinations_target,x_coordinations_target,n_particle);
coordinations_target = round(coordinations_target);
% Set distribution for connection radius
[x_connectionradius, ~, cdf_connectionradius] = distributions(p.equiprobability_connectionradius, p.mean_connectionradius, p.sigma_connectionradius, p.maxdeviation_connectionradius);

id_particles = [1:1:n_particle]';
connections = zeros(n_particle,max(x_coordinations_target));
coordination_current = zeros(n_particle,1);

Web_skeleton = zeros(sz);
all_connections = [];
for k=1:1:n_particle-1
    % All other particles
    remaining_centroids = centroids(k+1:end,:);
    remaining_ids = id_particles(k+1:end,:);
    [n_remaining,~] = size(remaining_centroids);

    % Current particle
    current_centroid = centroids(k,:);
    current_centroids = ones(n_remaining,dimension).*current_centroid;

    % Distance with all other particles sorted from close to far
    dists = (sum((current_centroids-remaining_centroids).^2,2)).^0.5;
    dists = [dists remaining_ids];
    dists = sortrows(dists,1,'ascend');

    % Number of new connection to establish
    r = coordinations_target(k) - coordination_current(k);

    if r>0
        for kk=1:1:n_remaining
            if dists(kk,1)>p.max_distance % Too far
                break
            else
                other_particle_id = dists(kk,2);
                if ismember(k,connections(other_particle_id,:)) % Already connected
                    continue
                else
                    if coordination_current(other_particle_id)<coordinations_target(other_particle_id)
                        % Establish connection
                        xa = current_centroid;
                        xb = centroids(other_particle_id,:);
                        Linf = max(abs(xb-xa))*2;

                        x = round(linspace(xa(1),xb(1),Linf));
                        y = round(linspace(xa(2),xb(2),Linf));
                        if dimension == 2
                            idx = sub2ind(sz,x,y);
                        else
                            z = round(linspace(xa(3),xb(3),Linf));
                            idx = sub2ind(sz,x,y,z);
                        end
                        
                        rad = pickfromdistriubtion(p.maxdeviation_connectionradius,p.equiprobability_connectionradius,p.mean_connectionradius,cdf_connectionradius,x_connectionradius,1);
                        Web_skeleton(idx)= round(rad);

                        coordination_current(other_particle_id) = coordination_current(other_particle_id)+1;
                        connections(other_particle_id,coordination_current(other_particle_id)) = k;

                        coordination_current(k) = coordination_current(k) + 1;
                        connections(k,coordination_current(k)) = other_particle_id;

                        tmp = [k other_particle_id rad];
                        all_connections = [all_connections; tmp];

                        r = r -1;
                        if r==0 % Enough connection established
                            break
                        end
                    end
                end
            end
        end
    end

end

%% CHORDS
chords = zeros(sz);
unis_chordradius = unique(Web_skeleton);
unis_chordradius(unis_chordradius==0)=[];
for k = 1:1:length(unis_chordradius)
    tmp = zeros(sz);
    tmp(Web_skeleton == unis_chordradius(k))=1;
    tmp_dmap=bwdist(tmp);
    chords(tmp_dmap<=unis_chordradius(k)) = 1;
end

%% NODES
if strcmp(p.node_choice,'Controlled with a distribution')
    % Pick coordination target
    [x_noderadius, ~, cdf_noderadius] = distributions(p.equiprobability_noderadius, p.mean_noderadius, p.sigma_noderadius, p.maxdeviation_noderadius);
    noderadius = pickfromdistriubtion(p.maxdeviation_noderadius,p.equiprobability_noderadius,p.mean_noderadius,cdf_noderadius,x_noderadius,n_particle);
elseif strcmp(p.node_choice,'Scaled with chord radius')
    noderadius = zeros(n_particle,1);
    x0s = max([ones(n_particle,1), centroids(:,1)-5],[],2); x1s = min([ones(n_particle,1)*sz(1) centroids(:,1)+5],[],2);
    y0s = max([ones(n_particle,1) centroids(:,2)-5],[],2); y1s = min([ones(n_particle,1)*sz(2) centroids(:,2)+5],[],2);
    if dimension == 3
        z0s = max([ones(n_particle,1) centroids(:,3)-5],[],2); z1s = min([ones(n_particle,1)*sz(3) centroids(:,3)+5],[],2);
    end
    for k=1:1:n_particle
        if dimension == 2
            sub = Web_skeleton(x0s(k):x1s(k),y0s(k):y1s(k));
        else
            sub = Web_skeleton(x0s(k):x1s(k),y0s(k):y1s(k),z0s(k):z1s(k));
        end
        noderadius(k) = p.node_radius_scalingfactor * max(max(max(sub)));
    end
elseif strcmp(p.node_choice,'Scaled with number of connection')
    noderadius = zeros(n_particle,1);
    x0s = max([ones(n_particle,1), centroids(:,1)-5],[],2); x1s = min([ones(n_particle,1)*sz(1) centroids(:,1)+5],[],2);
    y0s = max([ones(n_particle,1) centroids(:,2)-5],[],2); y1s = min([ones(n_particle,1)*sz(2) centroids(:,2)+5],[],2);
    if dimension == 3
        z0s = max([ones(n_particle,1) centroids(:,3)-5],[],2); z1s = min([ones(n_particle,1)*sz(3) centroids(:,3)+5],[],2);
    end
    for k=1:1:n_particle
        if dimension == 2
            sub = Web_skeleton(x0s(k):x1s(k),y0s(k):y1s(k));
        else
            sub = Web_skeleton(x0s(k):x1s(k),y0s(k):y1s(k),z0s(k):z1s(k));
        end
        b = max(max(max(sub))) - p.node_radius_a;
        noderadius(k) = (p.node_radius_a * coordination_current(k) + b)*p.node_radius_scalingfactor;
    end
end

nodes = zeros(sz);
if dimension == 2
    idx = sub2ind(sz,centroids(:,1),centroids(:,2));
else
    idx = sub2ind(sz,centroids(:,1),centroids(:,2),centroids(:,3));
end
unis_noderadius = unique(noderadius);
for k = 1:1:length(unis_noderadius)
    idx2 = find(noderadius == unis_noderadius(k));
    tmp = zeros(sz);
    tmp(idx(idx2))=1;
    dmap=bwdist(tmp);
    nodes(dmap<=unis_noderadius(k))=1;
end

%% JUNCTIONS
junctions = zeros(sz);
dmap2chord = bwdist(chords);
dmap2node = bwdist(nodes);
dmapcomb_1 = dmap2chord .* dmap2node;

BW = chords + nodes + junctions;
for k = 1:1:p.junction_dilationiteration
        
    if k==1
        dmapcomb = dmapcomb_1;
    else
        dmap2junc = bwdist(junctions);
        dmapcomb = dmapcomb_1 .* dmap2junc;
    end

    dmapcomb(BW~=0)=1e9;
    idx = dmapcomb == min(min(min(dmapcomb)));
    junctions(idx)=1;

    BW = chords + nodes + junctions;
    BW(BW~=0)=1;    
end

BW(BW~=0)=1; 
BW = uint8(BW);

%% SPHERES-ONLY REPRESENTATION
if p.convert_spheres
    sphere_connectivity = zeros(1e8,2);
    connn = 0;
    id_max = n_particle;
    BW_spheres = nodes;
    nodes_centroids_radius_id = [centroids noderadius [1:1:id_max]']; 
    chord_centroids_radius_id = [];
    [n_connection,~] = size(all_connections);
    for k=1:1:n_connection
        kp = all_connections(k,1);
        otherkp = all_connections(k,2);
        radius_chord = all_connections(k,3);
        xn0 = centroids(kp,:);
        rn0 = noderadius(kp);
        id0 = kp;
        
        xn1 = centroids(otherkp,:);
        rn1 = noderadius(otherkp);
        id1 = otherkp;

        dist = (sum((xn1-xn0).^2))^0.5;
        r_dist = dist-rn0-rn1;
        nsphere = round(r_dist/(2*radius_chord));
        rsphere = ceil(r_dist/(2*nsphere))+1;

        for ks=1:1:nsphere
            ds = rn0 + (2*(ks-1)+1)*rsphere;
            ds_norm = min([ds/dist 1]);
            xs = round(ds_norm*(xn1(1)-xn0(1)) + xn0(1));
            ys = round(ds_norm*(xn1(2)-xn0(2)) + xn0(2));
            xs = min([sz(1) xs]); ys = min([sz(2) ys]);
            xs = max([1 xs]); ys = max([1 ys]);            
            if dimension == 3
                zs = round(ds_norm*(xn1(3)-xn0(3)) + xn0(3));
                zs = min([sz(3) zs]);
                zs = max([1 zs]);
            end
            id_max = id_max+1;
            if dimension == 2
                chord_centroids_radius_id = [chord_centroids_radius_id; [xs ys rsphere id_max]];
            else
                chord_centroids_radius_id = [chord_centroids_radius_id; [xs ys zs rsphere id_max]];
            end

            connn = connn+1;
            if ks==1
                sphere_connectivity(connn,:) = [id0 id_max];
            elseif ks==nsphere
                sphere_connectivity(connn,:) = [id_max id1];
            else
                sphere_connectivity(connn,:) = [id_max-1 id_max];
            end

        end
    end

    unis_r = unique(chord_centroids_radius_id(:,dimension+1));
    for k=1:1:length(unis_r)
        idx = chord_centroids_radius_id(:,dimension+1)==unis_r(k);
        coors = chord_centroids_radius_id(idx,1:dimension);
        if dimension == 2
            idx = sub2ind(sz,coors(:,1),coors(:,2));
        else
            idx = sub2ind(sz,coors(:,1),coors(:,2),coors(:,3));
        end
        BWtmp = zeros(sz);
        BWtmp(idx)=1;
        dmap = bwdist(BWtmp);
        BW_spheres(dmap<=unis_r(k))=1;
    end
    BW_spheres = uint8(BW_spheres);

    % Connectivity
    sphere_connectivity(connn+1:end,:)=[];   
      
else
    BW_spheres = [];
    nodes_centroids_radius_id = [];
    chord_centroids_radius_id = [];
    sphere_connectivity = [];


end


%% FUNCTIONS
    function [x, pdf_, cdf_] = distributions(equiprobability, mean_, sigma_, maxdeviation)
        x = []; pdf_ = []; cdf_ =[];
        % Univariate Gaussian Distribution
        % Probability density function of a gaussian distribution
        pdf_GD = @(x,mu,sigma) 1./(2*pi*sigma.^2).^(0.5).*exp(-(x-mu).^2 ./ (2*sigma.^2));

        if ~equiprobability && maxdeviation~=0
            x=linspace(mean_-maxdeviation, mean_+maxdeviation, 1000);
            pdf_=pdf_GD(x,mean_,sigma_);
            cdf_ = pdf2cdf(x,pdf_);
        end

        function c = pdf2cdf(x,pdf)
            n=length(x);
            c=zeros(1,n);
            for i=2:1:n
                c(i) = trapz(x(1:i),pdf(1:i));
            end
            if c(end)<1
                c(end)=1;
            end
        end
    end

    function [val] = pickfromdistriubtion(maxdeviation,equiprobability,mean_val,cdf_,x_val,n)
        if maxdeviation==0
            val = ones(n,1)*mean_val;
        else
            if equiprobability
                val = rand(n,1)*(2*maxdeviation) + (mean_val-maxdeviation);
            else
                val = interp1(cdf_,x_val,rand(n,1));
            end
        end
    end



end