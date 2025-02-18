function [distance_transform, fitted_diameter, dist_fct, analytical, error_radius] = Function_particle_size_distancemap_Algorithm(M, label, removeinclusion, distancelabel, voxel_size, density_fct_parameters)
% Fit distance map cumulative function with distance map cumulative
% function of a sphere of diameter D
% See article: DOI 10.1149/1945-7111/ab913b

%% NOTATION
% Binary phase
% =1 phase investigated
% =0 complementary phase
Domain_size=size(M); % Domain size of the microstructure
dim = length(Domain_size);
binary_phase=zeros(Domain_size); % Initialization

%% CREATE BINARY PHASE
if distancelabel==-1
    binary_phase(M==label)=1;
else
    binary_phase(M~=distancelabel)=1;
end

if removeinclusion
    L = bwlabeln(~binary_phase,6);
    [C,~,ic] = unique(L);
    if length(C)>2
        counts = accumarray(ic,1);
        value_counts = [C, counts];
        idx_binary = find(value_counts(:,1)==0);
        value_counts(idx_binary,:)=[];
        value_counts = sortrows(value_counts,-2);
        for k=2:1:length(C)-1
            binary_phase(L==value_counts(k,1)) = 1;
        end
    end
end

%% EUCLIDEAN DISTANCE MAP
distance_transform=double(bwdist(~binary_phase,'euclidean')); % Euclidean distance map
%distance_transform=distance_transform-0.5; % From voxel center to surface
distance_transform(distance_transform<=0)=0; % Complementary volume
distance_transform(M~=label)=0; % Complementary volume
% Unit conversion
distance_transform = distance_transform*voxel_size;

distance_transform_vector=reshape(distance_transform,[prod(Domain_size),1]); % Convert in a vector
distance_transform_vector(distance_transform_vector<=0)=[]; % % Remove all 0

%% CUMULATIVE AND DISTRIBUTION FUNCTIONS
[dist_fct, ~] = Function_probability_density(distance_transform_vector,[],density_fct_parameters);
distance = dist_fct.cumulative_fct(:,1);
cumulative= dist_fct.cumulative_fct(:,2);

%% FIT DIAMETER FOR THE ANALYTICAL CUMULATIVE FUNCTIONS FOR A SPHERE OF DIAMETER 2*R
% Parameters
gridres = 1200;
Rmin=0.5*voxel_size;
Rmax=2*max(distance);
nR = 400;
Rs = linspace(Rmin,Rmax,nR);

error_radius=zeros(nR,2);
if length(distance)>=2
    for kR=1:1:nR
        R=Rs(kR);
        d=linspace(0,R,gridres);
        if dim==3
            C_edm_sphere = (R-d).^3./R.^3; % cf. eq. 12
        else
            C_edm_sphere = (R-d).^2./R.^2;
        end
        if R<max(distance)
            d_tmp=[d max(distance)];
            C_edm_sphere_tmp = [C_edm_sphere 0];
            distance_tmp = distance;
            cumulative_tmp = cumulative;
        elseif R>max(distance)
            distance_tmp = [distance; R];
            cumulative_tmp = [cumulative; 0];
            d_tmp = d;
            C_edm_sphere_tmp = C_edm_sphere;
        else
            d_tmp = d;
            C_edm_sphere_tmp = C_edm_sphere;
            distance_tmp = distance;
            cumulative_tmp = cumulative;
        end
        common_xaxis=linspace(0,max(distance_tmp),gridres);
        analytical_ = interp1(d_tmp,C_edm_sphere_tmp,common_xaxis,'linear',1);
        numerical_ = interp1(distance_tmp,cumulative_tmp,common_xaxis,'linear',1);

        difference = analytical_-numerical_;
        integral_=trapz(common_xaxis,difference); % integral

        error_radius(kR,2)=abs(integral_); % Absolute error...
        error_radius(kR,1)=R; % ...obtained with this radius
    end
    % Find the mimimum
    index_=find(error_radius(:,2)==min(error_radius(:,2)));
    fitted_radius=error_radius(index_(1),1);
    fitted_diameter=2*fitted_radius;

    % Recalculate the fitted cumulative fct.
    R=Rs(index_(1));
    d=linspace(0,R,gridres);
    if dim==3
        C_edm_sphere = (R-d).^3./R.^3; % cf. eq. 12
    else
        C_edm_sphere = (R-d).^2./R.^2;
    end
    if R<max(distance)
        d_tmp=[d max(distance)];
        C_edm_sphere_tmp = [C_edm_sphere 0];
        distance_tmp = distance;
        cumulative_tmp = cumulative;
    elseif R>max(distance)
        distance_tmp = [distance; R];
        cumulative_tmp = [cumulative; 0];
        d_tmp = d;
        C_edm_sphere_tmp = C_edm_sphere;
    else
        d_tmp = d;
        C_edm_sphere_tmp = C_edm_sphere;
        distance_tmp = distance;
        cumulative_tmp = cumulative;
    end
    common_xaxis=linspace(0,max(distance_tmp),gridres);
    analytical=zeros(gridres,2);
    analytical(:,1) = common_xaxis;
    analytical(:,2) = interp1(d_tmp,C_edm_sphere_tmp,common_xaxis,'linear',1);
else
    fitted_diameter = NaN;
    analytical=NaN;
    error_radius=NaN;
end

end