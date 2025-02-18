function [BW_grains] = function_add_grains(BW,p)


%% GENERATE SEEDS
sz = size(BW);
dimension = length(sz);

n_cubes = round(sz./p.average_size);

bounds_x = round(linspace(1,sz(1),n_cubes(1)+1));
bounds_y = round(linspace(1,sz(2),n_cubes(2)+1));
if dimension == 3
    bounds_z = round(linspace(1,sz(3),n_cubes(3)+1));
end

Seeds = zeros(sz);
for kx=1:1:n_cubes(1)
    x0 = bounds_x(kx)+1;
    x1 = bounds_x(kx+1)-1;
    mean_x = round((x0+x1)/2);
    maxdeviation_x = min([abs(mean_x-x0) abs(mean_x-x1)]);
    [x, ~, cdfx_] = distributions(p.equiprobability(1), mean_x, p.sigma_(1), maxdeviation_x);
    for ky=1:1:n_cubes(2)
        y0 = bounds_y(ky)+1;
        y1 = bounds_y(ky+1)-1;
        mean_y = round((y0+y1)/2);
        maxdeviation_y = min([abs(mean_y-y0) abs(mean_y-y1)]);
        [y, ~, cdfy_] = distributions(p.equiprobability(2), mean_y, p.sigma_(2), maxdeviation_y);
        if dimension == 3
            for kz=1:1:n_cubes(3)
                z0 = bounds_z(kz)+1;
                z1 = bounds_z(kz+1)-1;
                mean_z = round((z0+z1)/2);
                maxdeviation_z = min([abs(mean_z-z0) abs(mean_z-z1)]);
                [z, ~, cdfz_] = distributions(p.equiprobability(3), mean_z, p.sigma_(3), maxdeviation_z);
                [pos_x] = pickfromdistriubtion(maxdeviation_x,p.equiprobability(1),mean_x,cdfx_,x,1);
                [pos_y] = pickfromdistriubtion(maxdeviation_y,p.equiprobability(2),mean_y,cdfy_,y,1);
                [pos_z] = pickfromdistriubtion(maxdeviation_z,p.equiprobability(3),mean_z,cdfz_,z,1);
                Seeds(round(pos_x),round(pos_y),round(pos_z)) = 1;
            end
        else
            [pos_x] = pickfromdistriubtion(maxdeviation_x,p.equiprobability(1),mean_x,cdfx_,x,1);
            [pos_y] = pickfromdistriubtion(maxdeviation_y,p.equiprobability(2),mean_y,cdfy_,y,1);
            Seeds(round(pos_x),round(pos_y)) = 1;
        end

    end
end

%% ZONE OF INFLUENCE
[~,IDX] = bwdist(Seeds);
Zones = zeros(sz);

% Renumerotation
unis = unique(IDX);
n_unis = length(unis);
r = randperm(n_unis);
for k=1:1:n_unis
    Zones(IDX==unis(k))=r(k);
end

[counts_Zones, groupnames_Zones] = groupcounts(reshape(Zones,[numel(Zones) 1]));
Inter = double(BW).*Zones;
[counts_Inter, groupnames_Inter] = groupcounts(reshape(Inter,[numel(Inter) 1]));
id0 = find(groupnames_Inter==0);
if ~isempty(id0)
    groupnames_Inter(id0)=[];
    counts_Inter(id0)=[];
end

BW_grains = zeros(sz);
id = 0;
for k=1:1:length(groupnames_Inter)
    kk = groupnames_Inter(k);
    if counts_Inter(k) >= 0.5 * counts_Zones(kk)
        id = id+1;
        BW_grains(Zones==groupnames_Zones(kk)) = id;
    end
end

if max(max(max(BW_grains)))<=255
    BW_grains=uint8(BW_grains);
elseif max(max(max(BW_grains)))<=65535
    BW_grains=uint16(BW_grains);
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