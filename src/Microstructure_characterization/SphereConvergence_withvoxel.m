function [] = SphereConvergence_withvoxel(maxdiameter)
    arguments
        maxdiameter {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(maxdiameter,9)}
    end

disp( 'For unique sphere');

%% GROUND TRUTH
truth_diameter = 1; % 1 um (unit diameter)

% Sphere
truth_spheresurface = 4*pi*(truth_diameter/2)^2;
truth_spherevolume = 4/3*pi*(truth_diameter/2)^3;
turth_specificsurfacearea = truth_spheresurface/truth_spherevolume;

% Disc
truth_discperimeter = 2*pi*(truth_diameter/2);
truth_discarea = pi*(truth_diameter/2)^2;
turth_specificperimeterlength = truth_discperimeter/truth_discarea;

%% DIAMETER LIST
if (rem(maxdiameter, 2) == 0) % Even
    list_odd_voxelperdiameter = 7:2:maxdiameter+1; % n+1+n
    list_even_voxelperdiameter = 6:2:maxdiameter; % n+n
else
    list_odd_voxelperdiameter = 7:2:maxdiameter; % n+1+n
    list_even_voxelperdiameter = 6:2:maxdiameter-1; % n+n
end
%list_odd_voxelperdiameter = 7:2:201; % n+1+n
%list_even_voxelperdiameter = 6:2:200; % n+n

%% CORRECTIVE FACTOR
% Voxel representation induces a surface and perimeter overestimation
corrective_factor_spherearea = 2/3;
corrective_factor_discperimeter = pi/4;

%% INITIALIZATION
n_odd = length(list_odd_voxelperdiameter);
n_even = length(list_even_voxelperdiameter);

% Volume (area)
numerical_odd_spherevolume = zeros(n_odd,2);
numerical_even_spherevolume = zeros(n_even,2); 
numerical_odd_discarea = zeros(n_odd,2);
numerical_even_discarea = zeros(n_even,2); 

% Surface (perimeter)
numerical_odd_spheresurface = zeros(n_odd,2);
numerical_even_spheresurface = zeros(n_even,2); 
numerical_odd_discperimeter = zeros(n_odd,2);
numerical_even_discperimeter = zeros(n_even,2); 

% Specific surface area (specific perimeter length)
numerical_odd_spheresp = zeros(n_odd,2);
numerical_even_spheresp = zeros(n_even,2); 
numerical_odd_discsp = zeros(n_odd,2);
numerical_even_discsp = zeros(n_even,2); 

% Diameter
numerical_odd_sphere_diameterVOL = zeros(n_odd,2);
numerical_even_sphere_diameterVOL = zeros(n_even,2); 
numerical_odd_disc_diameterVOL = zeros(n_odd,2);
numerical_even_disc_diameterVOL = zeros(n_even,2); 

numerical_odd_sphere_diameterCPSD = zeros(n_odd,2);
numerical_even_sphere_diameterCPSD = zeros(n_even,2); 
numerical_odd_disc_diameterCPSD = zeros(n_odd,2);
numerical_even_disc_diameterCPSD = zeros(n_even,2); 

numerical_odd_sphere_diameterEDMF = zeros(n_odd,2);
numerical_even_sphere_diameterEDMF = zeros(n_even,2); 
numerical_odd_disc_diameterEDMF = zeros(n_odd,2);
numerical_even_disc_diameterEDMF = zeros(n_even,2); 

numerical_odd_sphere_diameterCHORD = zeros(n_odd,2);
numerical_even_sphere_diameterCHORD = zeros(n_even,2); 
numerical_odd_disc_diameterCHORD = zeros(n_odd,2);
numerical_even_disc_diameterCHORD = zeros(n_even,2); 

%% CALCULATION

% Odd diameter number
for k=1:1:n_odd
    voxelperdiameter = list_odd_voxelperdiameter(k);
    fprintf('Voxel per diameter = %i\n',voxelperdiameter);
        
    % Create geometry
    dim = voxelperdiameter+2;
    BW_sphere = zeros(dim,dim,dim);
    BW_disc = zeros(dim,dim);
    center = round(dim/2);
    BW_sphere(center,center,center) = 1;
    BW_disc(center,center) = 1;
    dmap_sphere = bwdist(BW_sphere);
    dmap_disc = bwdist(BW_disc);
    radius = floor(voxelperdiameter/2);
    BW_sphere(dmap_sphere<=radius) = 1;
    BW_disc(dmap_disc<=radius) = 1;
    if voxelperdiameter==15 % For illustration
        BW_disc_15 = BW_disc;
        BW_disc_15(center,center) = 2;
    end

    % % Calculate metrics (sphere)
    % Volume
    numbervoxel_sphere = sum(sum(sum(BW_sphere)));
    % Surface
    [~,n_voxel_phasesurface, ~] = Function_Specificsurface_direct_Algorithm(BW_sphere);
    % Diameter CPSD
    [Particle_size_CPSD] = Function_particle_size_CPSD_Algorithm(BW_sphere);
    all_diameters_CPSD = Particle_size_CPSD(BW_sphere==1);
    mean_diameter_CPSD = mean(all_diameters_CPSD);
    % Diameter EDMF
    density_fct_parameters.round_value = 3;
    density_fct_parameters.smooth_cumulative_fct = true;
    [~, mean_diameter_EDMF, ~, ~, ~] = Function_particle_size_distancemap_Algorithm(BW_sphere, 1, false, -1, 1, density_fct_parameters);
    % Diameter CHORD
    [mean_diameter_CHORD, ~, ~, ~, ~,~,~] = Function_particle_size_ChordLength_volume_Algorithm(BW_sphere,0);

    % Sphere approximation: thought voxel
    voxelsize = truth_diameter/(voxelperdiameter-1);
    numerical_odd_spherevolume(k,1) = numbervoxel_sphere * voxelsize^3;
    numerical_odd_spheresurface(k,1) = n_voxel_phasesurface * voxelsize^2 * corrective_factor_spherearea;
    numerical_odd_spheresp(k,1) = numerical_odd_spheresurface(k,1)/numerical_odd_spherevolume(k,1);
    numerical_odd_sphere_diameterVOL(k,1) = 2*((3*numerical_odd_spherevolume(k,1) / (4*pi))^(1/3));
    numerical_odd_sphere_diameterCPSD(k,1) = mean_diameter_CPSD * voxelsize;
    numerical_odd_sphere_diameterEDMF(k,1) = mean_diameter_EDMF * voxelsize;
    numerical_odd_sphere_diameterCHORD(k,1) = mean_diameter_CHORD(1) * voxelsize;
    
    % Sphere approximation: voxel edge
    voxelsize = truth_diameter/voxelperdiameter;
    numerical_odd_spherevolume(k,2) = numbervoxel_sphere * voxelsize^3;
    numerical_odd_spheresurface(k,2) = n_voxel_phasesurface * voxelsize^2 * corrective_factor_spherearea;
    numerical_odd_spheresp(k,2) = numerical_odd_spheresurface(k,2)/numerical_odd_spherevolume(k,2);
    numerical_odd_sphere_diameterVOL(k,2) = 2*((3*numerical_odd_spherevolume(k,2) / (4*pi))^(1/3));
    numerical_odd_sphere_diameterCPSD(k,2) = mean_diameter_CPSD * voxelsize;
    numerical_odd_sphere_diameterEDMF(k,2) = mean_diameter_EDMF * voxelsize;
    numerical_odd_sphere_diameterCHORD(k,2) = mean_diameter_CHORD(1) * voxelsize;

    % % Calculate metrics (disc)
    % Area
    numberpixel_disc = sum(sum(BW_disc));
    % Perimeter
    [~,n_pixel_phasesurface, ~] = Function_Specificsurface_direct_Algorithm(BW_disc);
    % Diameter CPSD
    [Particle_size_CPSD] = Function_particle_size_CPSD_Algorithm(BW_disc);
    all_diameters_CPSD = Particle_size_CPSD(BW_disc==1);
    mean_diameter_CPSD = mean(all_diameters_CPSD);   
    % Diameter EDMF
    [~, mean_diameter_EDMF, ~, ~, ~] = Function_particle_size_distancemap_Algorithm(BW_disc, 1, false, -1, 1, density_fct_parameters);
    % Diameter CHORD
    %[mean_diameter_CHORD, ~, ~, ~, ~,~,~] = Function_particle_size_ChordLength_volume_Algorithm(BW_disc,0);
    [mean_diameter_CHORD] = Function_particle_size_ChordLength_section_Algorithm(BW_disc,[0 0 1],0);

    % Disc approximation: though pixel
    pixelsize = truth_diameter/(voxelperdiameter-1);
    numerical_odd_discarea(k,1) = numberpixel_disc * pixelsize^2;
    numerical_odd_discperimeter(k,1) = n_pixel_phasesurface * pixelsize * corrective_factor_discperimeter;
    numerical_odd_discsp(k,1) = numerical_odd_discperimeter(k,1)/numerical_odd_discarea(k,1);
    numerical_odd_disc_diameterVOL(k,1) = 2*((numerical_odd_discarea(k,1) / pi)^(1/2));
    numerical_odd_disc_diameterCPSD(k,1) = mean_diameter_CPSD * pixelsize;
    numerical_odd_disc_diameterEDMF(k,1) = mean_diameter_EDMF * pixelsize;
    numerical_odd_disc_diameterCHORD(k,1) = mean_diameter_CHORD(3).diameter_eqdisc * pixelsize;

    % Disc approximation: pixel edge
    pixelsize = truth_diameter/voxelperdiameter;
    numerical_odd_discarea(k,2) = numberpixel_disc * pixelsize^2;    
    numerical_odd_discperimeter(k,2) = n_pixel_phasesurface * pixelsize * corrective_factor_discperimeter;
    numerical_odd_discsp(k,2) = numerical_odd_discperimeter(k,2)/numerical_odd_discarea(k,2);
    numerical_odd_disc_diameterVOL(k,2) = 2*((numerical_odd_discarea(k,2) / pi)^(1/2));
    numerical_odd_disc_diameterCPSD(k,2) = mean_diameter_CPSD * pixelsize;
    numerical_odd_disc_diameterEDMF(k,2) = mean_diameter_EDMF * pixelsize;
    numerical_odd_disc_diameterCHORD(k,2) = mean_diameter_CHORD(3).diameter_eqdisc * pixelsize;
    
end

% Even diameter number
for k=1:1:n_even
    voxelperdiameter = list_even_voxelperdiameter(k);
    fprintf('Voxel per diameter = %i\n',voxelperdiameter);
        
    % Create geometry
    dim = voxelperdiameter+2;
    BW_sphere = zeros(dim,dim,dim);
    BW_disc = zeros(dim,dim);
    center1 = dim/2;
    center2 = center1+1;
    BW_sphere(center1,center1,center1) = 1;
    BW_sphere(center2,center1,center1) = 1;
    BW_sphere(center1,center2,center1) = 1;
    BW_sphere(center2,center2,center1) = 1;
    BW_sphere(center1,center1,center2) = 1;
    BW_sphere(center2,center1,center2) = 1;
    BW_sphere(center1,center2,center2) = 1;
    BW_sphere(center2,center2,center2) = 1;
    BW_disc(center1,center1) = 1;
    BW_disc(center2,center1) = 1;
    BW_disc(center1,center2) = 1;
    BW_disc(center2,center2) = 1;    
    dmap_sphere = bwdist(BW_sphere);
    dmap_disc = bwdist(BW_disc);
    radius = (voxelperdiameter-2)/2;
    BW_sphere(dmap_sphere<=radius) = 1;
    BW_disc(dmap_disc<=radius) = 1;
    if voxelperdiameter==16 % For illustration
        BW_disc_16 = BW_disc;
        BW_disc_16(center1,center1) = 2;
        BW_disc_16(center2,center1) = 2;
        BW_disc_16(center1,center2) = 2;
        BW_disc_16(center2,center2) = 2;
    end

    % % Calculate metrics (sphere)
    % Volume
    numbervoxel_sphere = sum(sum(sum(BW_sphere)));
    % Surface
    [~,n_voxel_phasesurface, ~] = Function_Specificsurface_direct_Algorithm(BW_sphere);
    % Diameter CPSD
    [Particle_size_CPSD] = Function_particle_size_CPSD_Algorithm(BW_sphere);
    all_diameters_CPSD = Particle_size_CPSD(BW_sphere==1);
    mean_diameter_CPSD = mean(all_diameters_CPSD);
    % Diameter EDMF
    [~, mean_diameter_EDMF, ~, ~, ~] = Function_particle_size_distancemap_Algorithm(BW_sphere, 1, false, -1, 1, density_fct_parameters);
    % Diameter CHORD
    [mean_diameter_CHORD, ~, ~, ~, ~,~,~] = Function_particle_size_ChordLength_volume_Algorithm(BW_sphere,0);

    % Sphere approximation: thought voxel
    voxelsize = truth_diameter/(voxelperdiameter-1);
    numerical_even_spherevolume(k,1) = numbervoxel_sphere * voxelsize^3;
    numerical_even_spheresurface(k,1) = n_voxel_phasesurface * voxelsize^2 * corrective_factor_spherearea;
    numerical_even_spheresp(k,1) = numerical_even_spheresurface(k,1)/numerical_even_spherevolume(k,1);
    numerical_even_sphere_diameterVOL(k,1) = 2*(3*numerical_even_spherevolume(k,1) / (4*pi))^(1/3);
    numerical_even_sphere_diameterCPSD(k,1) = mean_diameter_CPSD * voxelsize;
    numerical_even_sphere_diameterEDMF(k,1) = mean_diameter_EDMF * voxelsize;
    numerical_even_sphere_diameterCHORD(k,1) = mean_diameter_CHORD(1) * voxelsize;

    % Sphere approximation: voxel edge
    voxelsize = truth_diameter/voxelperdiameter;
    numerical_even_spherevolume(k,2) = numbervoxel_sphere * voxelsize^3;
    numerical_even_spheresurface(k,2) = n_voxel_phasesurface * voxelsize^2 * corrective_factor_spherearea;
    numerical_even_spheresp(k,2) = numerical_even_spheresurface(k,2)/numerical_even_spherevolume(k,2);
    numerical_even_sphere_diameterVOL(k,2) = 2*(3*numerical_even_spherevolume(k,2) / (4*pi))^(1/3);
    numerical_even_sphere_diameterCPSD(k,2) = mean_diameter_CPSD * voxelsize;
    numerical_even_sphere_diameterEDMF(k,2) = mean_diameter_EDMF * voxelsize;
    numerical_even_sphere_diameterCHORD(k,2) = mean_diameter_CHORD(1) * voxelsize;

    % % Calculate metrics (disc)
    % Area
    numberpixel_disc = sum(sum(BW_disc));
    % Perimeter
    [~,n_pixel_phasesurface, ~] = Function_Specificsurface_direct_Algorithm(BW_disc);
    % Diameter CPSD
    [Particle_size_CPSD] = Function_particle_size_CPSD_Algorithm(BW_disc);
    all_diameters_CPSD = Particle_size_CPSD(BW_disc==1);
    mean_diameter_CPSD = mean(all_diameters_CPSD);     
    % Diameter EDMF
    [~, mean_diameter_EDMF, ~, ~, ~] = Function_particle_size_distancemap_Algorithm(BW_disc, 1, false, -1, 1, density_fct_parameters);
    % Diameter CHORD
    %[mean_diameter_CHORD, ~, ~, ~, ~,~,~] = Function_particle_size_ChordLength_volume_Algorithm(BW_disc,0);
    [mean_diameter_CHORD] = Function_particle_size_ChordLength_section_Algorithm(BW_disc,[0 0 1],0);

    % Disc approximation: though pixel
    pixelsize = truth_diameter/(voxelperdiameter-1);
    numerical_even_discarea(k,1) = numberpixel_disc * pixelsize^2;
    numerical_even_discperimeter(k,1) = n_pixel_phasesurface * pixelsize * corrective_factor_discperimeter;
    numerical_even_discsp(k,1) = numerical_even_discperimeter(k,1)/numerical_even_discarea(k,1);
    numerical_even_disc_diameterVOL(k,1) = 2*(numerical_even_discarea(k,1) / pi)^(1/2);
    numerical_even_disc_diameterCPSD(k,1) = mean_diameter_CPSD * pixelsize;
    numerical_even_disc_diameterEDMF(k,1) = mean_diameter_EDMF * pixelsize;
    numerical_even_disc_diameterCHORD(k,1) = mean_diameter_CHORD(3).diameter_eqdisc * pixelsize;

    % Disc approximation: pixel edge
    pixelsize = truth_diameter/voxelperdiameter;
    numerical_even_discarea(k,2) = numberpixel_disc * pixelsize^2;  
    numerical_even_discperimeter(k,2) = n_pixel_phasesurface * pixelsize * corrective_factor_discperimeter;
    numerical_even_discsp(k,2) = numerical_even_discperimeter(k,2)/numerical_even_discarea(k,2);
    numerical_even_disc_diameterVOL(k,2) = 2*(numerical_even_discarea(k,2) / pi)^(1/2);
    numerical_even_disc_diameterCPSD(k,2) = mean_diameter_CPSD * pixelsize;
    numerical_even_disc_diameterEDMF(k,2) = mean_diameter_EDMF * pixelsize;
    numerical_even_disc_diameterCHORD(k,2) = mean_diameter_CHORD(3).diameter_eqdisc * pixelsize;

end

%% PLOT SPECIFIC SURFACE
col = colororder;

Fig = figure; % Create figure
Fig.Name= 'Specific surface convergence';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.9*scrsz(3) 0.8*scrsz(4)]);
t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
title(t,'Specific surface convergence','FontWeight','normal','Fontname','Times new roman','Fontsize',16)
t.Subtitle.String = 'Assuming spherical particles wiht unit diameter (1 \mum) and surface corrective factor of 2/3';
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = 'Times new roman';
t.Subtitle.FontSize = 14;
for k=1:1:6
    nexttile
    if k==1
        hold on
        plot(list_odd_voxelperdiameter,numerical_odd_spherevolume(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_spherevolume(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_spherevolume(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_spherevolume(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_spherevolume truth_spherevolume],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere volume (\mum3)')
    elseif k==2
        hold on
        plot(list_odd_voxelperdiameter,numerical_odd_spheresurface(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_spheresurface(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_spheresurface(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_spheresurface(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_spheresurface truth_spheresurface],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere surface (\mum2)')
    elseif k==3
        hold on
        plot(list_odd_voxelperdiameter,numerical_odd_spheresp(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_spheresp(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_spheresp(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_spheresp(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[turth_specificsurfacearea turth_specificsurfacearea],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere specific surface area (\mum^{-1})')
    elseif k==4
        hold on
        plot(list_odd_voxelperdiameter,abs(100*(truth_spherevolume-numerical_odd_spherevolume(:,1))/truth_spherevolume),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_spherevolume-numerical_odd_spherevolume(:,2))/truth_spherevolume),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_spherevolume-numerical_even_spherevolume(:,1))/truth_spherevolume),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_spherevolume-numerical_even_spherevolume(:,2))/truth_spherevolume),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Sphere volume relative error (%)')
    elseif k==5
        hold on
        plot(list_odd_voxelperdiameter,abs(100*(truth_spheresurface-numerical_odd_spheresurface(:,1))/truth_spheresurface),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_spheresurface-numerical_odd_spheresurface(:,2))/truth_spheresurface),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_spheresurface-numerical_even_spheresurface(:,1))/truth_spheresurface),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_spheresurface-numerical_even_spheresurface(:,2))/truth_spheresurface),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        ylabel('Sphere surface relative error (%)')
    elseif k==6
        hold on
        plot(list_odd_voxelperdiameter,abs(100*(turth_specificsurfacearea-numerical_odd_spheresp(:,1))/turth_specificsurfacearea),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(turth_specificsurfacearea-numerical_odd_spheresp(:,2))/turth_specificsurfacearea),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(turth_specificsurfacearea-numerical_even_spheresp(:,1))/turth_specificsurfacearea),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(turth_specificsurfacearea-numerical_even_spheresp(:,2))/turth_specificsurfacearea),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        ylabel('Sphere specific surface area (%)')
    end
    xlim([0 list_odd_voxelperdiameter(end)+9])
    xlabel('Number of voxel per diameter')
    set(gca,'Fontname','Times new roman','Fontsize',12)
    grid(gca,"on"); % Display grid
    set(gca,'XMinorGrid',"on",'YMinorGrid',"on"); % Display grid for minor thicks
    h_legend = legend(gca,'Location','best');
    hold off
end

Fig = figure; % Create figure
Fig.Name= 'Specific perimeter convergence';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.9*scrsz(3) 0.8*scrsz(4)]);
t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
title(t,'Specific perimeter convergence','FontWeight','normal','Fontname','Times new roman','Fontsize',16)
t.Subtitle.String = 'Assuming disc wiht unit diameter (1 \mum) and perimeter corrective factor of \pi/4';
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = 'Times new roman';
t.Subtitle.FontSize = 14;
for k=1:1:6
    nexttile
    if k==1
        hold on
        plot(list_odd_voxelperdiameter,numerical_odd_discarea(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_discarea(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_discarea(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_discarea(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_discarea truth_discarea],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Disc area (\mum2)')
    elseif k==2
        hold on
        plot(list_odd_voxelperdiameter,numerical_odd_discperimeter(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_discperimeter(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_discperimeter(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_discperimeter(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_discperimeter truth_discperimeter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Disc perimeter (\mum)')
    elseif k==3
        hold on
        plot(list_odd_voxelperdiameter,numerical_odd_discsp(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_discsp(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_discsp(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_discsp(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[turth_specificperimeterlength turth_specificperimeterlength],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Disc specific perimeter (\mum^{-1})')
    elseif k==4
        hold on
        plot(list_odd_voxelperdiameter,abs(100*(truth_discarea-numerical_odd_discarea(:,1))/truth_discarea),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_discarea-numerical_odd_discarea(:,2))/truth_discarea),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_discarea-numerical_even_discarea(:,1))/truth_discarea),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_discarea-numerical_even_discarea(:,2))/truth_discarea),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Disc area relative error (%)')
    elseif k==5
        hold on
        plot(list_odd_voxelperdiameter,abs(100*(truth_discperimeter-numerical_odd_discperimeter(:,1))/truth_discperimeter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_discperimeter-numerical_odd_discperimeter(:,2))/truth_discperimeter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_discperimeter-numerical_even_discperimeter(:,1))/truth_discperimeter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_discperimeter-numerical_even_discperimeter(:,2))/truth_discperimeter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Disc perimeter relative error (%)')   
    elseif k==6
        hold on
        plot(list_odd_voxelperdiameter,abs(100*(turth_specificperimeterlength-numerical_odd_discsp(:,1))/turth_specificperimeterlength),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(turth_specificperimeterlength-numerical_odd_discsp(:,2))/turth_specificperimeterlength),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(turth_specificperimeterlength-numerical_even_discsp(:,1))/turth_specificperimeterlength),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(turth_specificperimeterlength-numerical_even_discsp(:,2))/turth_specificperimeterlength),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Disc specific perimeter relative error (%)')  
    end
    xlim([0 list_odd_voxelperdiameter(end)+9])
    xlabel('Number of pixel per diameter')
    set(gca,'Fontname','Times new roman','Fontsize',12)
    grid(gca,"on"); % Display grid
    set(gca,'XMinorGrid',"on",'YMinorGrid',"on"); % Display grid for minor thicks
    h_legend = legend(gca,'Location','best');
    hold off
end

%% PLOT DIAMETER
Fig = figure; % Create figure
Fig.Name= 'Diameter convergence';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.9*scrsz(3) 0.8*scrsz(4)]);
t = tiledlayout(2,4,'TileSpacing','Compact','Padding','Compact');
title(t,'Diameter convergence','FontWeight','normal','Fontname','Times new roman','Fontsize',16)
t.Subtitle.String = 'Assuming spherical particles wiht unit diameter (1 \mum)';
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = 'Times new roman';
t.Subtitle.FontSize = 14;
for k=1:1:8
    nexttile
    if k==1
        hold on
        title('Volume-based')
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterVOL(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterVOL(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterVOL(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterVOL(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere diameter (mum)')
    elseif k==2
        hold on
        title('C-PSD')
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterCPSD(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterCPSD(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterCPSD(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterCPSD(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere diameter (mum)')
    elseif k==3
        hold on
        title('EDMF')
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterEDMF(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterEDMF(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterEDMF(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterEDMF(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere diameter (mum)')
    elseif k==4
        hold on
        title('CLMF')
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterCHORD(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_sphere_diameterCHORD(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterCHORD(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_sphere_diameterCHORD(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere diameter (mum)')
    elseif k==5
        hold on
        title('Volume-based')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterVOL(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterVOL(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterVOL(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterVOL(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Sphere diameter relative error (%)')        
    elseif k==6
        hold on
        title('C-PSD')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterCPSD(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterCPSD(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterCPSD(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterCPSD(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Sphere diameter relative error (%)')   
    elseif k==7
        hold on
        title('EDMF')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterEDMF(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterEDMF(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterEDMF(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterEDMF(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Sphere diameter relative error (%)')   
    elseif k==8
        hold on
        title('CLMF')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterCHORD(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_sphere_diameterCHORD(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterCHORD(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_sphere_diameterCHORD(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Sphere diameter relative error (%)')       
    end
    xlim([0 list_odd_voxelperdiameter(end)+9])
    xlabel('Number of voxel per diameter')
    set(gca,'Fontname','Times new roman','Fontsize',12)
    grid(gca,"on"); % Display grid
    set(gca,'XMinorGrid',"on",'YMinorGrid',"on"); % Display grid for minor thicks
    h_legend = legend(gca,'Location','best');
    hold off
end

Fig = figure; % Create figure
Fig.Name= 'Diameter convergence';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.9*scrsz(3) 0.8*scrsz(4)]);
t = tiledlayout(2,4,'TileSpacing','Compact','Padding','Compact');
title(t,'Diameter convergence','FontWeight','normal','Fontname','Times new roman','Fontsize',16)
t.Subtitle.String = 'Assuming disc wiht unit diameter (1 \mum)';
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = 'Times new roman';
t.Subtitle.FontSize = 14;
for k=1:1:8
    nexttile
    if k==1
        hold on
        title('Volume-based')
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterVOL(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterVOL(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterVOL(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterVOL(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Disc diameter (mum)')
    elseif k==2
        hold on
        title('C-PSD')
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterCPSD(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterCPSD(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterCPSD(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterCPSD(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Disc diameter (mum)')
    elseif k==3
        hold on
        title('EDMF')
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterEDMF(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterEDMF(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterEDMF(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterEDMF(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere diameter (mum)')
    elseif k==4
        hold on
        title('CLMF')
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterCHORD(:,1),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,numerical_odd_disc_diameterCHORD(:,2),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterCHORD(:,1),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,numerical_even_disc_diameterCHORD(:,2),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")
        plot([list_even_voxelperdiameter(1),list_odd_voxelperdiameter(end)],[truth_diameter truth_diameter],"Color",'k',"LineStyle","-","LineWidth",2,"DisplayName","Ground truth")
        ylabel('Sphere diameter (mum)')
    elseif k==5
        hold on
        title('Volume-based')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterVOL(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterVOL(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterVOL(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterVOL(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Disc diameter relative error (%)')
    elseif k==6
        hold on
        title('C-PSD')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterCPSD(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterCPSD(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterCPSD(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterCPSD(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Disc diameter relative error (%)')
    elseif k==7
        hold on
        title('EDMF')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterEDMF(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterEDMF(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterEDMF(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterEDMF(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Disc diameter relative error (%)')
    elseif k==8        
        hold on
        title('CLMF')
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterCHORD(:,1))/truth_diameter),"Color",col(1,:),"LineStyle","-","LineWidth",2,"DisplayName","Odd diameter, through voxel")
        plot(list_odd_voxelperdiameter,abs(100*(truth_diameter-numerical_odd_disc_diameterCHORD(:,2))/truth_diameter),"Color",col(1,:),"LineStyle","--","LineWidth",2,"DisplayName","Odd diameter, voxel edge")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterCHORD(:,1))/truth_diameter),"Color",col(2,:),"LineStyle","-","LineWidth",2,"DisplayName","Even diameter, through voxel")
        plot(list_even_voxelperdiameter,abs(100*(truth_diameter-numerical_even_disc_diameterCHORD(:,2))/truth_diameter),"Color",col(2,:),"LineStyle","--","LineWidth",2,"DisplayName","Even diameter, voxel edge")        
        ylabel('Disc diameter relative error (%)') 
    end
    xlim([0 list_odd_voxelperdiameter(end)+9])
    xlabel('Number of pixel per diameter')
    set(gca,'Fontname','Times new roman','Fontsize',12)
    grid(gca,"on"); % Display grid
    set(gca,'XMinorGrid',"on",'YMinorGrid',"on"); % Display grid for minor thicks
    h_legend = legend(gca,'Location','best');
    hold off
end

end