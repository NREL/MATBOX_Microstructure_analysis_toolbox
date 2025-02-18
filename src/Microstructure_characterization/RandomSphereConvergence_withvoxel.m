function [] = RandomSphereConvergence_withvoxel(maxdiameter)
    arguments
        maxdiameter {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(maxdiameter,9)}
    end

disp('For randomly distributed isolated solid spherical particles. FOV = 15 diameters.');

filename = 'particle_data_sphereexample.mat'; % FOV: 357x357x357 voxels, sphere with 21 voxels, porosity 0.50. Then cropped to 15 spheres.
diameter_baseline = 17;
domain_size = [357 357 357];
if ispc
    separation_folder = '\';
else
    separation_folder = '/';
end
path_app = mfilename('fullpath');
higherlevelfolder = extractBetween(path_app,path_app(1),'src','Boundaries','inclusive');
filepath_path = [char(higherlevelfolder) '\Microstructure_generation' separation_folder filename];
if exist(filepath_path,'file')
    particle_data = load(filepath_path);
    particle_data = particle_data.particle_data;
else
    warning(['MATLAB did not find the file ' filename ' for the tortuosity dependence with mesh resolution']);
    warning('Default location is \MATBOX_Microstructure_analysis_toolbox\src\Microstructure_generation\');
end


porosity_target = 0.65; %  Before particle start to overlap, to test algorithm on 'isolated spheres' case for which analytical solution is known.
shuffleparticles = false;
scale_diameterratio = 1.0;
Forceseparation = false;
Distanceseparation = 1;

% Diameters
if (rem(maxdiameter, 2) == 0) % Even
    list_odd_voxelperdiameter = 7:2:maxdiameter+1; % n+1+n
    list_even_voxelperdiameter = 6:2:maxdiameter; % n+n
else
    list_odd_voxelperdiameter = 7:2:maxdiameter; % n+1+n
    list_even_voxelperdiameter = 6:2:maxdiameter-1; % n+n
end
list_odd_scalingfactor = list_odd_voxelperdiameter/diameter_baseline;
list_even_scalingfactor = list_even_voxelperdiameter/diameter_baseline;

% Ground truth
truth_bruggeman = 1.5; % Isolated spheres

% Initialization
n_odd = length(list_odd_voxelperdiameter);
n_even = length(list_even_voxelperdiameter);
numerical_odd_sphere_p = zeros(n_odd,1);
numerical_even_sphere_p = zeros(n_even,1); 

Fig=figure;
ax_=axes('Parent',Fig);
hold(ax_,'on');
col = gray(n_odd);

% % Calculation
% Odd diameter number
for k=n_odd:1:n_odd
    k/n_odd
    voxelperdiameter = list_odd_voxelperdiameter(k);
    fprintf('Voxel per diameter = %i\n',voxelperdiameter);
    % Create geometry
    scaling_factor = list_odd_scalingfactor(k);
    if scaling_factor~=1
        domain_size_rescaled = round(domain_size*scaling_factor);
        particle_data_scaled = particle_data;
        particle_data_scaled(:,3) = min(round( scaling_factor * particle_data(:,3) + 0.5 ), domain_size_rescaled(1));
        particle_data_scaled(:,4) = min(round( scaling_factor * particle_data(:,4) + 0.5 ), domain_size_rescaled(2));
        particle_data_scaled(:,5) = min(round( scaling_factor * particle_data(:,5) + 0.5 ), domain_size_rescaled(3));
        particle_data_scaled(:,6) = round(particle_data(:,6) * scaling_factor);
        particle_data_scaled(:,7) = round(particle_data(:,7) * scaling_factor);
        particle_data_scaled(:,8) = round(particle_data(:,8) * scaling_factor);
    else
        domain_size_rescaled = domain_size;
        particle_data_scaled = particle_data;
    end
    currentdiameter = particle_data_scaled(1,6);
    [Phaselabel_scaled, ~, ~] = function_createvolume_fromparticledata(particle_data_scaled,domain_size_rescaled, scale_diameterratio, Forceseparation, Distanceseparation, porosity_target, shuffleparticles);
    % Crop edges to remove edge effect
    Phaselabel_scaled = Phaselabel_scaled(currentdiameter:end-currentdiameter,currentdiameter:end-currentdiameter,currentdiameter:end-currentdiameter);

%     % Calculate Bruggeman exponent
%     porosity = sum(sum(sum(Phaselabel_scaled==0)))/numel(Phaselabel_scaled);
%     Tau_factor_result = TauFactor('InLine',1,0,0,~Phaselabel_scaled,[0 0 0;0 0 0;1 0 0],[1 1 1]);
%     tau =  Tau_factor_result.Tau_W1.Tau;
%     [results] = Function_TortuosityPorosity_fitting(porosity,tau, 'Bruggeman');
%     numerical_odd_sphere_p(k) = results.Bruggeman;

    % Calculate covariance function
    sz = size(Phaselabel_scaled);
    cov_orth = zeros(sz(1),2);

    cov_orth(:,1) = [1:1:sz(1)]/sz(1);
    for kk = 1:1:sz(1)
        covariogram = sum(sum(sum( Phaselabel_scaled(1:end-kk,:,:).*Phaselabel_scaled(kk+1:end,:,:) )));
        cov_orth(kk,2) = covariogram / ((sz(1)-kk)*sz(2)*sz(3));
    end
%     % centered derivate
%     cov_der = zeros(sz(1)-2,2);
%     kkk=0;
%     for kk=2:1:sz(1)-1
%         kkk=kkk+1;
%         cov_der(kkk,1) = cov_orth(kk,1);
%         cov_der(kkk,2) = (cov_orth(kk+1,2) - cov_orth(kk-1,2)) / (cov_orth(kk+1,1) - cov_orth(kk-1,1));
%     end
%     cov_der(:,2) = cov_der(:,2)/max(cov_der(:,2)); % Normalize

    

    %yyaxis left
    %plot(cov_orth(:,1),cov_orth(:,2),'Color',col(k,:));
    %yyaxis right
    %plot(cov_der(:,1),cov_der(:,2),'Color',col(k,:));
    
    

    TFmax = islocalmax(cov_orth(:,2));
    TFmin = islocalmin(cov_orth(:,2));
    plot(cov_orth(:,1),cov_orth(:,2),cov_orth(TFmax,1),cov_orth(TFmax,2),'r*',cov_orth(TFmin,1),cov_orth(TFmin,2),'ro')
    %pause(1)
    keyboard


end


lll

% Even diameter number
for k=1:1:n_even
    voxelperdiameter = list_even_voxelperdiameter(k);
    fprintf('Voxel per diameter = %i\n',voxelperdiameter);
    % Create geometry
    scaling_factor = list_even_scalingfactor(k);
    if scaling_factor~=1
        domain_size_rescaled = round(domain_size*scaling_factor);
        particle_data_scaled = particle_data;
        particle_data_scaled(:,3) = min(round( scaling_factor * particle_data(:,3) + 0.5 ), domain_size_rescaled(1));
        particle_data_scaled(:,4) = min(round( scaling_factor * particle_data(:,4) + 0.5 ), domain_size_rescaled(2));
        particle_data_scaled(:,5) = min(round( scaling_factor * particle_data(:,5) + 0.5 ), domain_size_rescaled(3));
        particle_data_scaled(:,6) = round(particle_data(:,6) * scaling_factor);
        particle_data_scaled(:,7) = round(particle_data(:,7) * scaling_factor);
        particle_data_scaled(:,8) = round(particle_data(:,8) * scaling_factor);
    else
        domain_size_rescaled = domain_size;
        particle_data_scaled = particle_data;
    end
    currentdiameter = particle_data_scaled(1,6);
    [Phaselabel_scaled, ~, ~] = function_createvolume_fromparticledata(particle_data_scaled,domain_size_rescaled, scale_diameterratio, Forceseparation, Distanceseparation, porosity_target, shuffleparticles);
    % Crop edges to remove edge effect
    Phaselabel_scaled = Phaselabel_scaled(currentdiameter:end-currentdiameter,currentdiameter:end-currentdiameter,currentdiameter:end-currentdiameter);

    % Calculate Bruggeman exponent
    porosity = sum(sum(sum(Phaselabel_scaled==0)))/numel(Phaselabel_scaled);
    Tau_factor_result = TauFactor('InLine',1,0,0,~Phaselabel_scaled,[0 0 0;0 0 0;1 0 0],[1 1 1]);
    tau =  Tau_factor_result.Tau_W1.Tau;
    [results] = Function_TortuosityPorosity_fitting(porosity,tau, 'Bruggeman');
    numerical_even_sphere_p(k) = results.Bruggeman;
end

% Plot
col = colororder;
Fig = figure; % Create figure
Fig.Name= 'Pore transport convergence';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.8*scrsz(3) 0.6*scrsz(4)]);
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
title(t,'Pore transport convergence','FontWeight','normal','Fontname','Times new roman','Fontsize',16)
t.Subtitle.String = 'For randomly distributed isolated solid spherical particles';
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = 'Times new roman';
t.Subtitle.FontSize = 14;
for k=1:1:2
    nexttile
    if k==1
        hold on
        title('TauFactor (finite difference, Dirichlet boundary condition)')
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_p(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameterl")
        plot(list_even_voxelperdiameter,numerical_even_sphere_p(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_bruggeman truth_bruggeman],"Color",'k',"LineStyle","--","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Bruggeman exponent')
    elseif k==2
        hold on
        title('TauFactor (finite difference, Dirichlet boundary condition)')
        plot(list_odd_voxelperdiameter,abs(100*(truth_bruggeman-numerical_odd_sphere_p(:,1))/truth_bruggeman),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_bruggeman-numerical_even_sphere_p(:,1))/truth_bruggeman),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        ylabel('Bruggeman exponent relative error (%)')
    end
    xlim([0 list_odd_voxelperdiameter(end)+9])
    xlabel('Number of voxel per diameter')
    set(gca,'Fontname','Times new roman','Fontsize',12)
    grid(gca,"on"); % Display grid
    set(gca,'XMinorGrid',"on",'YMinorGrid',"on"); % Display grid for minor thicks
    h_legend = legend(gca,'Location','best');
    hold off
end


end