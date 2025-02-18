%% ORIENTATION PARAMETERS
mean_azimuth = 0;
sigma_azimuth = [];
maxdeviation_azimuth = 180;
equiprobability_azimuth = true;

mean_elevation = 0;
maxdeviation_elevation = 90;
equiprobability_elevation = false;

%% NUMERIC PARAMETERS
n_elevation = 10;
sigmas_elevation = linspace(5,100,n_elevation);
n_point = 500;
D = 1;

Lzs = zeros(1,n_elevation);
Lxys = zeros(1,n_elevation);

%% ALGORITHM

azimuth_lim=rad2deg(atan(2));

n_elevation = length(sigmas_elevation);
[x_azimuth, ~, cdf_azimuth] = distributions(equiprobability_azimuth, mean_azimuth, sigma_azimuth, maxdeviation_azimuth);

ytmp = 1e-6*[1:1:1e4];

for k=1:1:n_elevation
    k
    sigma_elevation = sigmas_elevation(k);
    [x_elevation, ~, cdf_elevation] = distributions(equiprobability_elevation, mean_elevation, sigma_elevation, maxdeviation_elevation);

    Lzs_tmp = zeros(1,n_point);
    Lxys_tmp = zeros(1,n_point);
    for k_sample = 1:1:n_point
        cdf_elevation = cdf_elevation+ytmp;
        elevation = pickfromdistriubtion(maxdeviation_elevation,equiprobability_elevation,mean_elevation,cdf_elevation,x_elevation);
        azimuth= pickfromdistriubtion(maxdeviation_azimuth,equiprobability_azimuth,mean_azimuth,cdf_azimuth,x_azimuth);
        elevation = abs(elevation);
        azimuth = abs(azimuth);

        %Lz
        y = D/2*tan(deg2rad(elevation));
        if y>=D
            x = D/tan(deg2rad(elevation));
            Lzs_tmp(k_sample) = (x^2+y^2)^0.5;
        else
            Lzs_tmp(k_sample) = NaN;
        end

        % Lxy
        % Azimuth
        if azimuth>=azimuth_lim && azimuth <=180-azimuth_lim
            azimuth = abs(90-azimuth);
            t=D/cos(deg2rad(azimuth));
            % elevation
            y = t*tan(deg2rad(elevation));
            if y<=D/2
                Lxys_tmp(k_sample) = (t^2+y^2)^0.5;
            else
                Lxys_tmp(k_sample) = NaN;
            end
        else
            Lxys_tmp(k_sample) = NaN;
        end

    end

    Lzs(k) = nanmean(Lzs_tmp);
    Lxys(k) = nanmean(Lxys_tmp);




end

Fig = figure;
ax_=axes;
hold on
plot(sigmas_elevation,Lzs)
plot(sigmas_elevation,Lxys)



% [x_elevationchange, ~, p.cdf_elevationchange] = distributions(equiprobability_elevationchange, p.mean_elevationchange, p.sigma_elevationchange, p.maxdeviation_elevationchange);



function [x, pdf_, cdf_] = distributions(equiprobability, mean_, sigma_, maxdeviation)
x = []; pdf_ = []; cdf_ =[];
% Univariate Gaussian Distribution
% Probability density function of a gaussian distribution
pdf_GD = @(x,mu,sigma) 1./(2*pi*sigma.^2).^(0.5).*exp(-(x-mu).^2 ./ (2*sigma.^2));

if ~equiprobability && maxdeviation~=0
    x=linspace(mean_-maxdeviation, mean_+maxdeviation, 1e4);
    pdf_=pdf_GD(x,mean_,sigma_);
    cdf_ = pdf2cdf(x,pdf_);
end

    function c = pdf2cdf(x,pdf)
        n=length(x);
        c=zeros(1,n);
        for k=2:1:n
            c(k) = trapz(x(1:k),pdf(1:k));
        end
        if c(end)<1
            c(end)=1;
        end
    end
end

function [val] = pickfromdistriubtion(maxdeviation,equiprobability,mean_val,cdf_,x_val)
if maxdeviation==0
    val = mean_val;
else
    if equiprobability
        val = rand*(2*maxdeviation) + (mean_val-maxdeviation);
    else
        val = interp1(cdf_,x_val,rand);
    end
end
end
