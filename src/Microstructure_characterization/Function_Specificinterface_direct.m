function [] = Function_Specificinterface_direct(Phase_microstructure,infovol,opt,p, foo1, foo2)
% Calculate specific interface area
% Function_Specificinterface_direct(Phase_microstructure, infovol, opt, p, foo, foo) - when used with the MATBOX toolbox
% or
% Function_Specificinterface_direct(Phase_microstructure, doublet, voxelsize, unit, corrective_factor) - when used as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       doublet: a 1x2 array: the two labels for which the interface area is calculated
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       corrective factor, a scalar: surface area overestimation correction
%       si(phase i, phase j) = Interface phase i - phase j / domain volume
%       e.g.: Function_Specificinterface_direct(<your_3d_array>, [1,2],  0.4, 'um', 2/3); % for a volume with 1,2 being labels

%% DEFAULT VALUES
expected_number_argument = 6;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 5 % Case for function called as: Function_Specificinterface_direct(Phase_microstructure, doublet, voxelsize, unit, corrective_factor). Standalone use.
        doublet = infovol; clear infovol;
        voxelsize = opt; clear opt;
        unit = p; clear p;
        corrective_factor = foo1; clear foo1; clear foo2;
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['Intdirect_' char(t)];
        if ispc
            infovol.mainfolder=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
            separator = '\';
        else
            infovol.mainfolder = pwd;
            separator = '/';
        end
        infovol.volpath = [infovol.mainfolder separator infovol.volumesubfolder separator];    

        % Set default phase information
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        infovol.phaselabel = unique(Phase_microstructure);
        for k=1:1:length(infovol.phaselabel) % Loop over all unique values
            infovol.phasename(k,1) = {['Phase ' num2str(infovol.phaselabel(k))]}; % Assign phase name based on code value
            infovol.phasecolor(k,:) = allcolor(k,:);
        end
        infovol.voxelsize = voxelsize;
        infovol.unit = unit;
        
        % Set default direction information
        infovol.directionname = {'in-plane direction 1';'in-plane direction 2';'through-plane direction'};
        
        % Set default options
        opt.save.mat = true;
        opt.save.xls = true;
        opt.save.savefig = true;
        opt.save.fig_infig = true;
        opt.save.fig_format = {'png'};
        
        % Set display options
        opt.format.fontname = 'Times New Roman';
        opt.format.axefontsize =  12;
        opt.format.legendfontsize =  10;
        opt.format.titlefontsize =  14;
        opt.format.sgtitlefontsize =  16;
        opt.format.linewidth = 2;
        opt.format.grid = 'on'; opt.format.minorgrid = 'on';
        opt.format.autoclosefig = false;
        
        % No Voxel size dependence analysis
        p.scaling = 1;
        % No Representative Volume Element analysis
        p.RVE.number_RVE = 0;
        % Fractal analysis
        p.fractal_bertei.todo = false;     
        p.fractal_boxcounting = false;             

    else % Incorrect number of argument
        disp 'Error calling Function_Specificinterface_direct. Wrong number of argument.'
        help Function_Specificinterface_direct
        return
    end
    
else
    % Read parameters;
    corrective_factor = p.correctivefactor;
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    separator = '\';
else
    separator = '/';
end
Current_folder = [infovol.volpath 'Int_direct' separator];
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% VOLUME INFORMATION
Domain_size = size(Phase_microstructure);
if length(Domain_size)==2 % Add third column with unit value
    Domain_size = [Domain_size 1];
    number_dimension = 2;
    Area_str = 'Length';
    Volume_str = 'Area';
    Sp_str = 'Specific length';
else
    number_dimension =3; % 3D case
    Area_str = 'Surface';    
    Volume_str = 'Volume';
    Sp_str = 'Specific interface area';
end
number_phase = length(infovol.phaselabel); % Number of phase
voxel_number = prod(Domain_size); % Number of voxel
voxel_size = infovol.voxelsize;
voxel_unit = infovol.unit;

% Unit
Area_unit = [voxel_unit num2str(number_dimension-1)];
Volume_unit = [voxel_unit num2str(number_dimension)];
Int_unit = [voxel_unit '-1'];

%%
%% ALGORITHM ON WHOLE VOLUME
%%

colorinterface = [colororder; rand(1000,3)];

disp '    SPECIFIC INTERFACE AREA - DIRECT METHOD';
disp '    ---------------------------------------';
disp ' ';

%% CALCULATION
time_cpu_start_volume = cputime; % CPU start
time_clock_start_volume = tic; % Stopwatch start

% Name interfaces
[number_interface,~] = size(p.doublets);
interface_name = cell(number_interface,1);
p.todo = ones(1,number_interface);

for current_interface=1:1:number_interface
    index_phase_A = find([infovol.phaselabel] == p.doublets(current_interface,1));
    index_phase_B = find([infovol.phaselabel] == p.doublets(current_interface,2));
    %label_phase_A = p.doublets(current_interface,1);
    %label_phase_B = p.doublets(current_interface,2);
    name_phase_A = char(infovol.phasename(index_phase_A));
    name_phase_B = char(infovol.phasename(index_phase_B));
    %interface_label(current_interface) = [label_phase_A label_phase_B];
    interface_name(current_interface,1) = {[name_phase_A '-' name_phase_B]};
    results_correlation(current_interface).name = interface_name(current_interface,1);
    results_visualization(current_interface).name = interface_name(current_interface,1);
end

infovol_bis = infovol;
for current_interface=1:1:number_interface
    infovol_bis.phasename(current_interface) = interface_name(current_interface);
    infovol_bis.phasecolor(current_interface,:) = colorinterface(current_interface,:);
end

% Initialization (generic)
timedata = zeros(number_interface+1,3); timedata(1,1) = voxel_number;
timedata_domain = cell(number_interface+1,1);

% Initialization (algorithm-specific)
% Colunm 1: phase A
% Column 2: phase B
% Colunm 3: Interface area
% Column 4: Volume of the domain
% Column 5: Specific interface area (Interface area/Volume domain)
% Column 6: Specific interface area * corrective_factor
Specificinterface=zeros(number_interface,6); % Initialization
Specificinterface(:,4) = voxel_number*(voxel_size^number_dimension); % Domain volume
for current_interface=1:1:number_interface % Loop over all interfaces
    time_cpu_start_phase = cputime; % CPU start
    time_clock_start_phase = tic; % Stopwatch start

    % % Algorithm
    % Label
    label_phase_A = p.doublets(current_interface,1);
    label_phase_B = p.doublets(current_interface,2);
    Specificinterface(current_interface,1)=label_phase_A;
    Specificinterface(current_interface,2)=label_phase_B;
    % Create a binary microstructure : 1 = current phase A, 10 = current phase B, 0 = all the other phases
    binary_phase=zeros(Domain_size); % Initialization
    binary_phase(Phase_microstructure == label_phase_A) = 1;
    binary_phase(Phase_microstructure == label_phase_B) = 10;
    n_voxel_phases=sum(sum(sum(binary_phase~=0))); % Number of voxel of the two phases
    [Position_surface.interface(current_interface),surface_,LOC_surface.interface(current_interface).array] = Function_Specificinterface_direct_Algorithm(binary_phase); % Interface area
    Specificinterface(current_interface,3)=surface_*(voxel_size^(number_dimension-1));
    Specificinterface(current_interface,5)=Specificinterface(current_interface,3)/Specificinterface(current_interface,4);
    Specificinterface(current_interface,6)=Specificinterface(current_interface,5)*corrective_factor;
 
    % % Correlation
    results_correlation(current_interface).Specific_interface_directcounting = Specificinterface(current_interface,6);
    % % Visualization
    results_visualization(current_interface).Specific_interface_directcounting = LOC_surface.interface(current_interface).array;    

    % % Time
    timedata_domain(current_interface+1,1) = interface_name(current_interface);
    timedata(current_interface+1,1) = n_voxel_phases;
    timedata(current_interface+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
    timedata(current_interface+1,3) = toc(time_clock_start_phase); % Stopwatch elapsed time
end

% CPU and stopwatch time - end
timedata_domain(1,1) = {'Full volume'};
timedata(1,2) = cputime-time_cpu_start_volume; % CPU elapsed time
timedata(1,3) = toc(time_clock_start_volume); % Stopwatch elapsed time
timedata_pervolume = timedata(1,:);
timedata_perphase = timedata(2:end,:);

%% TABLES
% Time
Table_time = table(timedata_domain(:,1), timedata(:,1),timedata(:,2),timedata(:,3),...
    'VariableNames',{'Domain', 'Number of voxel','CPU time s' 'Stopwatch s'});
Results_specificinterfacearea.Table_time = Table_time; % Save in main table result
% Result calculated on whole volume
Table_Int_direct = table(interface_name,Specificinterface(:,1),Specificinterface(:,2),Specificinterface(:,3),Specificinterface(:,4),Specificinterface(:,5),ones(number_interface,1)*corrective_factor,Specificinterface(:,6),...
    'VariableNames',{'Interface' 'Label A' 'Label B' [Area_str ' ' Area_unit] [Volume_str ' ' Volume_unit] [Sp_str ' ' Int_unit] 'Corrective factor' ['Corrected ' Sp_str ' ' Int_unit]});%
Results_specificinterfacearea.Table_Int_direct = Table_Int_direct; % Save in main table result
 

%% SAVE TABLES
if opt.save.xls
    filename = 'Interface_direct'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Phaseinterface_over_Domainvolume';
    DATA_writetable.sheet(1).table=Table_Int_direct; 
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
disp([Sp_str ' = interface ' Area_str ' / domain ' Volume_str])
disp(Table_Int_direct)
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME
%%

%% ALONG DIRECTIONS
% README

% Specific surface area Sp =  sum(g(x)), with g(x)=n(x)*face area/phase volume
% and n(x) the number of face at the position x along a given direction
% This discrete sum defintion will be converted in a continous (integral)
% sum definition: Sp = 1/D*integral(f(x),dx,0,D) = sum(g(x))
% To get f(x) [um-1], sum(g(x)) is derived.

% This is analogous with calculating a probability density function.
% (Although, the integral calculated here is not a probability but [um-1])

% Pros: obtain correct shape of the curve Sp=f(x).
% For a randomly distributed phase with no gradation, f(x)=cst=<Sp>

%% CALCULATION OF g(x), with g(x)=n(x)*face area/volume
for current_interface=1:1:number_interface
    for direction=1:1:number_dimension
        % Set position in real unit
        Position_surface.interface(current_interface).direction(direction).value(:,1)=Position_surface.interface(current_interface).direction(direction).value(:,1)*voxel_size;
        % Set face surface in real unit
        Position_surface.interface(current_interface).direction(direction).value(:,3)=Position_surface.interface(current_interface).direction(direction).value(:,2)*(voxel_size^(number_dimension-1));
        %  Set specific surface area in real unit
        Position_surface.interface(current_interface).direction(direction).value(:,3)=Position_surface.interface(current_interface).direction(direction).value(:,3)/Specificinterface(current_interface,4);
        % Corrective factor
        Position_surface.interface(current_interface).direction(direction).value(:,3)=Position_surface.interface(current_interface).direction(direction).value(:,3)*corrective_factor;
    end
end

%% DISTRIBUTION
% Cumulative function sum(g(x))
for current_interface=1:1:number_interface
    for direction=1:1:number_dimension
        Cumulative_surface.interface(current_interface).direction(direction).value=zeros(1,2*Domain_size(direction)-1);
        for i=1:1:2*Domain_size(direction)-1
            for j=1:1:i
                Cumulative_surface.interface(current_interface).direction(direction).value(1,j)=Cumulative_surface.interface(current_interface).direction(direction).value(1,j)+Position_surface.interface(current_interface).direction(direction).value(i,3);
            end
        end
    end
end
% Derivate function (i.e. 1/D*f(x))
for current_interface=1:1:number_interface
    for direction=1:1:number_dimension
        c_length = 2*Domain_size(direction)-1;
        for i=2:1:(c_length-1)
            Derivate_Cumulative_surface.interface(current_interface).direction(direction).value(1,i)=-( Cumulative_surface.interface(current_interface).direction(direction).value(1,i+1)-Cumulative_surface.interface(current_interface).direction(direction).value(1,i-1) ) / (2*0.5*voxel_size);
        end
        Derivate_Cumulative_surface.interface(current_interface).direction(direction).value(1,1)=-(Cumulative_surface.interface(current_interface).direction(direction).value(1,2)-Cumulative_surface.interface(current_interface).direction(direction).value(1,1))/(0.5*voxel_size);
        Derivate_Cumulative_surface.interface(current_interface).direction(direction).value(1,c_length)=-(Cumulative_surface.interface(current_interface).direction(direction).value(1,c_length)-Cumulative_surface.interface(current_interface).direction(direction).value(1,c_length-1))/(0.5*voxel_size);
    end
end
% The specific surface area f(x)
for current_interface=1:1:number_interface
    for direction=1:1:number_dimension
        Specific_surface.interface(current_interface).direction(direction).value = zeros(1,2*Domain_size(direction)-1);
        Specific_surface.interface(current_interface).direction(direction).value = Derivate_Cumulative_surface.interface(current_interface).direction(direction).value*Domain_size(direction)*(voxel_size);
    end
end

%% VERIFICATION
% the relation Sp = 1/D*integral(f(x),dx,0,D) is checked
Specificsurface_integralvalue=zeros(number_interface,7);
for current_interface=1:1:number_interface
    for direction=1:1:number_dimension
        X_= Position_surface.interface(current_interface).direction(direction).value(:,1);
        Y_= Specific_surface.interface(current_interface).direction(direction).value(1,:);
        Specificsurface_integralvalue(current_interface,1)=Specificinterface(current_interface,6); % Reference value
        Specificsurface_integralvalue(current_interface,direction+1)=(1/(Domain_size(direction)*(voxel_size)))*trapz(X_,Y_); % Integral value
        % Relative error in percent
        Specificsurface_integralvalue(current_interface,direction+1+3)=100*(Specificsurface_integralvalue(current_interface,direction+1)-Specificsurface_integralvalue(current_interface,1))/Specificsurface_integralvalue(current_interface,1);
    end
end

%% MANAGING RESULTS
% Results are saved in a table

clear Variable_name_table;
Variable_name_table(1)={['Position ' voxel_unit]};
for current_interface=1:1:number_interface
    Variable_name_table(current_interface+1)=interface_name(current_interface);
end

% g(x)
for direction=1:1:number_dimension
    array_g = zeros(2*Domain_size(direction)-1,number_interface+1);
    array_g(:,1)= Position_surface.interface(current_interface).direction(direction).value(:,1);
    for current_interface=1:1:number_interface
        array_g(:,current_interface+1)= Position_surface.interface(current_interface).direction(direction).value(:,3);
    end
    Gx.direction(direction).table = array2table(array_g,'VariableNames',Variable_name_table);
end
% f(x)
for direction=1:1:number_dimension
    array_f = zeros(2*Domain_size(direction)-1,number_interface+1);
    array_f(:,1)= Position_surface.interface(current_interface).direction(direction).value(:,1);
    for current_interface=1:1:number_interface
        array_f(:,current_interface+1)= Specific_surface.interface(current_interface).direction(direction).value(1,:);
    end
    Fx.direction(direction).table = array2table(array_f,'VariableNames',Variable_name_table);
end
 
% Verification
Table_verification_integral = table(interface_name,Specificsurface_integralvalue(:,1),Specificsurface_integralvalue(:,2),Specificsurface_integralvalue(:,3),Specificsurface_integralvalue(:,4),Specificsurface_integralvalue(:,5),Specificsurface_integralvalue(:,6),Specificsurface_integralvalue(:,7),...
    'VariableNames',{'Interface' 'Reference' 'Integration direction 1' 'Integration direction 2' 'Integration direction 3' 'Relative error percent direction 1' 'Relative error percent direction 2' 'Relative error percent direction 3'});
fprintf('> Calculated along direction:\n\n');
fprintf('  Integral verification: relative error should be close to 0 percent:\n\n');
disp(Table_verification_integral)

Results_specificinterfacearea.Table_evolution_Gx = Gx; % Save in main table result
Results_specificinterfacearea.Table_evolution_Fx = Fx;
Results_specificinterfacearea.Table_verification_integral=Table_verification_integral;

%% SAVE TABLES
if opt.save.xls
    filename = 'Int_direct_along_directions'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    sheet=0;
    for current_direction=1:1:number_dimension % Loop over all directions
        sheet=sheet+1;
        DATA_writetable.sheet(sheet).name=['Along direction ' num2str(current_direction) ' g'];
        DATA_writetable.sheet(sheet).table=Gx.direction(current_direction).table;
    end
    for current_direction=1:1:number_dimension % Loop over all directions
        sheet=sheet+1;
        DATA_writetable.sheet(sheet).name=['Along direction ' num2str(current_direction) ' f'];
        DATA_writetable.sheet(sheet).table=Fx.direction(current_direction).table;
    end    
    sheet=sheet+1;
    % Data : Integral summation
    DATA_writetable.sheet(sheet).name='Integration_f';
    DATA_writetable.sheet(sheet).table=Table_verification_integral;  
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% FIGURES
strunit = voxel_unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)';
    Intunit = '\mum^{-1}';
else
    axisunit = ['(' strunit ')'];
    Intunit = ['' voxel_unit ' ^{-1}'];
end

scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Specific interface area, direct method (integral definition)'; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
for current_direction=1:1:number_dimension % Iterate over axe
    sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    h_title=title ({'Specific interface area, direct method',sprintf('Int_{p} = phase-phase interface %s / domain %s',Area_str,Volume_str)}); % Set title font
    % Plot graphs
    for current_interface=1:1:number_interface % Loop over phases
        x_ = Position_surface.interface(current_interface).direction(current_direction).value(:,1);
        y_ = Specific_surface.interface(current_interface).direction(current_direction).value(1,:);
        plot(x_,y_,'Color', colorinterface(current_interface,:),'LineWidth',opt.format.linewidth,'DisplayName',char(interface_name(current_interface)));
    end
    % Axis label
    t_ = xlabel(' ');
    t_1 = sprintf('Position along %s ',char(infovol.directionname(current_direction)));
    t_2 = axisunit;
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = ylabel(' '); 
    str1 = 'Specific interface distribution function f(x), such as';
    str2 = ['<Int_{p}> = 1/L \times \int_{0}^{L} f(x) dx (' Intunit ')'];
    t_.String= {str1,str2};
    % Legend
    for current_interface=1:1:number_interface
        str_legend(current_interface).name = [char(interface_name(current_interface)) ', <Int_{p}>=' num2str(Specificinterface(current_interface,6),'%1.3f') Intunit];
    end
    h_legend = legend(sub_axes,str_legend.name,'Location','best');
    % - Grid
    grid(sub_axes,opt.format.grid); % Display grid
    set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
    set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
    h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
    h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
    hold(sub_axes,'off'); % Relase figure    
end
sgtitle(Fig,'Specific interface area along directions','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
if opt.save.savefig % Save figure
    filename= 'Specific_interface_direct_along_direction';
    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
end
if opt.format.autoclosefig
    close(Fig); % Do not keep open figures
end

%% FRACTAL DIMENSION (COUNTING BOX)
if p.fractal_boxcounting.todo
    fprintf(['> ' Sp_str ' fractal dimension (box-counting method)\n']);
    p.fractal_boxcounting.topology_dimension = number_dimension-1;
    p.fractal_boxcounting.plot = false;
    for current_interface=1:1:number_interface % Loop over line and call box counting algorithm
        [N(:,current_interface),box_lengths,fractal_dimension(:,current_interface),fractal_dimension_convergence(:,:,current_interface)] = Function_fractaldimension_boxcounting(LOC_surface.interface(current_interface).array,p.fractal_boxcounting);
    end

    % Table
    Table_fractaldimension = table(interface_name(:,1),fractal_dimension(1,:)',fractal_dimension(2,:)',fractal_dimension(3,:)',fractal_dimension(4,:)',...
        'VariableNames',{'Interface' 'Fit from 1 to' 'Fractal dimension' 'Topology dimension' 'Fractal propensity'});
    for current_interface=1:1:number_interface
        Table_boxlength(current_interface).t = table(fractal_dimension_convergence(:,1,current_interface),fractal_dimension_convergence(:,2,current_interface),fractal_dimension_convergence(:,3,current_interface),...
        'VariableNames',{'Fit from 1 to' 'Fractal dimension' 'Fit norm error'});
    end
    Results_specificinterfacearea.Table_fractaldimension = Table_fractaldimension; % Save in main table result
    disp(Table_fractaldimension)

    % Save table
    if opt.save.xls
        filename = 'PhaseInterface_Fractaldimension_boxcounting'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name = 'Fractal dimension';
        DATA_writetable.sheet(1).table = Table_fractaldimension;
        for current_interface=1:1:number_interface
            DATA_writetable.sheet(1+current_interface).name = char(interface_name(current_interface,1));
            DATA_writetable.sheet(1+current_interface).table = Table_boxlength(current_interface).t;
        end           
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end

    % Figure
    if p.fractal_boxcounting.topology_dimension==2
        propertyname= 'Phase interface (area)'; % Figure name
    else
        propertyname = 'Phase interface (line)'; % Figure name
    end
    filename = 'PhaseInterface_fractal_dimension_boxcounting';
    function_fractalfig(propertyname,filename,Current_folder, box_lengths,N,fractal_dimension_convergence,p.fractal_boxcounting.topology_dimension,number_interface,p,opt,infovol_bis);

    % Correlation
    results_correlation(current_interface).PhaseInterface_fractaldimension_boxcounting = fractal_dimension(2,current_interface);
    results_correlation(current_interface).PhaseInterface_fractalpropensity_boxcounting = abs(p.fractal_boxcounting.topology_dimension - fractal_dimension(2,current_interface));
end

%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%

if length(p.scaling)>=2 % Check if voxel size analysis is asked
    size_choice = p.scaling;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed
    
    %% CALCULATION
    % Initialization
    property_voxelsizedependence = zeros(number_resize+1,number_interface+1);
    property_voxelsizedependence(1,1)=voxel_size;
    property_voxelsizedependence(1,2:end)=Specificinterface(:,6)';

    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fractaldimension_voxelsizedependence = zeros(number_resize+1,number_interface+1);
        fractaldimension_voxelsizedependence(1,1) = voxel_size;
        fractaldimension_voxelsizedependence(1,2:end) = fractal_dimension(2,:);
    end    

    % Loop on each voxel size
    for current_iteration=1:1:number_resize

        % % Microstructure scaling
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        property_voxelsizedependence(current_iteration+1,1)=current_voxel_size;
        % Set parameters
        parameters_scaling.scaling_factor = size_choice(current_iteration);
        parameters_scaling.label_or_greylevel = 'Label';
        parameters_scaling.background = min(infovol.phaselabel);
        % Scale
        Phase_microstructure_resized = function_scaling(Phase_microstructure,parameters_scaling);

        % CPU and stopwatch time - start
        time_cpu_start_volume = cputime; % CPU start
        time_clock_start_volume = tic; % Stopwatch start
        % Number of voxel of the current resized microstructure
        voxel_number_tmp=numel(Phase_microstructure_resized);
        Current_domain_size = size(Phase_microstructure_resized);
        n_=prod(Current_domain_size); % Number of voxel of the domain
        Volume_=n_*(current_voxel_size^number_dimension);

        for current_interface=1:1:number_interface % Loop over all phases
            time_cpu_start_phase = cputime; % CPU start
            time_clock_start_phase = tic; % Stopwatch start

            % % Algorithm: SPECIFIC FOR EACH FILE
            % Label
            label_phase_A = p.doublets(current_interface,1);
            label_phase_B = p.doublets(current_interface,2);
            % Create a binary microstructure : 1 = current phase A, 10 = current phase B, 0 = all the other phases
            binary_phase=zeros(Current_domain_size); % Initialization
            binary_phase(Phase_microstructure_resized == label_phase_A) = 1;
            binary_phase(Phase_microstructure_resized == label_phase_B) = 10;
            n_voxel_phases=sum(sum(sum(binary_phase~=0))); % Number of voxel of the two phases
            [~,surface_,LOC_surface_resized.interface(current_interface).array] = Function_Specificinterface_direct_Algorithm(binary_phase); % Interface area
            surface_=surface_*(current_voxel_size^(number_dimension-1));
            Int=surface_/Volume_;
            property_voxelsizedependence(current_iteration+1,current_interface+1)=Int*corrective_factor;

            % % Time
            timedata_perphase = [timedata_perphase; [n_voxel_phases (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

        end
        % CPU and stopwatch time - end
        timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];

        % % Fractal
        if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
            for current_interface=1:1:number_interface % Loop over line and call box counting algorithm
                [~,~,fractal_dimension_tmp, ~] = Function_fractaldimension_boxcounting(LOC_surface_resized.interface(current_interface).array,p.fractal_boxcounting);
                fractaldimension_voxelsizedependence(current_iteration+1,1)=size_choice(current_iteration)*voxel_size;
                fractaldimension_voxelsizedependence(current_iteration+1,current_interface+1) = fractal_dimension_tmp(2,1);
            end      
        end

    end
    clear Phase_microstructure_resized;

    % Sort per voxel size
    property_voxelsizedependence = sortrows(property_voxelsizedependence,1);   
    
    %% EXTRAPOLATION TO 0 nm
    fprintf('> Specific interface area dependence with the voxel size (with corrective factor)\n');
    tmp = zeros(number_resize+2,number_interface+1); % + 0 nm and + initial voxel size
    x=property_voxelsizedependence(:,1);
    if strcmp(p.scaling_extrapolation,'Linear')
        interpolation_voxelsize_order=1;
    elseif strcmp(p.scaling_extrapolation,'Quadratic')
        interpolation_voxelsize_order=2;
    elseif strcmp(p.scaling_extrapolation,'Cubic')
        interpolation_voxelsize_order=3;
    end
    max_order = length(p.scaling)-1;
    interpolation_voxelsize_order = min(interpolation_voxelsize_order,max_order);
    fprintf('  Extrapolation to zero voxel size: polynomial of order %i\n\n',interpolation_voxelsize_order);
    for current_interface=1:1:number_interface
        y=property_voxelsizedependence(:,current_interface+1);
        pi = polyfit(x,y,interpolation_voxelsize_order);
        vq = polyval(pi,0);
        tmp(1,current_interface+1)=vq;
        interpolation_voxelsize(current_interface).pi=pi;
    end
    tmp(2:end,:) = property_voxelsizedependence;
    property_voxelsizedependence = tmp; clear tmp;

    if p.fractal_bertei.todo
        fractal_dimension_Bertei = zeros(number_interface,3);
        % Fractal dimension according to Bertei et al., https://doi.org/10.1016/j.nanoen.2017.06.028 (Richardson/Mandelbrot formula)
        % Log(property) = m + (2-fractal dimension)*log(voxel size)
        logV = log(property_voxelsizedependence(2:end,1));
        logPs = zeros(length(logV),number_interface);
        for current_interface=1:1:number_interface
            logP = log(property_voxelsizedependence(2:end,current_interface+1));
            logPs(:,current_interface)=logP;
            pf = polyfit(logV,logP,1);
            fractal_dimension_Bertei(current_interface,1) = 2-pf(1);
            fractal_dimension_Bertei(current_interface,2) = 2; % Topology dimension, surface
            fractal_dimension_Bertei(current_interface,3) = abs( fractal_dimension_Bertei(current_interface,2) - fractal_dimension_Bertei(current_interface,1) ); % Fractal propensity
            % Correlation
            results_correlation(current_interface).PhaseInterface_fractaldimension_Mandelbrot = fractal_dimension_Bertei(current_interface,1);
            results_correlation(current_interface).PhaseInterface_fractalpropensity_Mandelbrot = fractal_dimension_Bertei(current_interface,3);
        end
    end    
    
    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
    for current_interface=1:1:number_interface
        Variable_name_table(1+current_interface)=interface_name(current_interface);
    end
    % Table
    Table_Intdirect_voxelsizedependence = array2table(property_voxelsizedependence,...
        'VariableNames',Variable_name_table);
    if p.fractal_bertei.todo
        Variable_name_table={'Voxel size log'}; % Columns name
        for current_interface=1:1:number_interface
            Variable_name_table(1+current_interface)={[char(interface_name(current_interface,1)) ' log']};
        end
        Table_Intdirect_voxelsizedependence_loglog = array2table([logV logPs],'VariableNames',Variable_name_table);
        Table_Fractaldimension_Bertei = table(interface_name(:,1),fractal_dimension_Bertei(:,1),fractal_dimension_Bertei(:,2),fractal_dimension_Bertei(:,3),...
            'VariableNames',{'Interface','Fractal dimension','Topology dimension','Fractal propensity'});   
    end

    
    %% DISPLAY TEXT RESULTS
    disp(Table_Intdirect_voxelsizedependence)
    if p.fractal_bertei.todo
        disp(Table_Intdirect_voxelsizedependence_loglog)
        fprintf('Richardson/Mandelbrot formula: Log(property) = m + (2-fractal dimension)*log(voxel size))\n');
        disp(Table_Fractaldimension_Bertei)
    end    
    
    %% SAVE RESULTS
    Results_specificinterfacearea.voxelsizedependence = Table_Intdirect_voxelsizedependence; % Save in main table result
    if opt.save.xls
        filename = 'Int_direct_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Specific interface area';
        DATA_writetable.sheet(1).table=Table_Intdirect_voxelsizedependence;
        if p.fractal_bertei.todo
            DATA_writetable.sheet(2).name='Log log';
            DATA_writetable.sheet(2).table=Table_Intdirect_voxelsizedependence_loglog;
            DATA_writetable.sheet(3).name='Fractal dimension Mandelbrot';
            DATA_writetable.sheet(3).table=Table_Fractaldimension_Bertei;
        end        
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    %% FIGURES
    parameters_figure.plotlog = true;
    parameters_figure.propertyname = Sp_str;
    parameters_figure.method = 'direct counting';
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence;
    parameters_figure.number_phase = number_interface;
    parameters_figure.str_ylabel = [sprintf('Int_{p} = phase-phase interface %s/domain %s',Area_str,Volume_str) '(' Intunit ')'];
    parameters_figure.propertynameunit = Intunit;
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'Int_direct_voxel_size_dependence';
    parameters_figure.infovol = infovol_bis;
    parameters_figure.opt = opt;
    parameters_figure.todo = ones(number_interface,1);
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    

    %% FRACTAL DIMENSION
    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fprintf(['> ' Sp_str ' fractal dimension (box counting) dependence with the voxel size\n']);
        fprintf('  Extrapolation to zero voxel size: polynomial of order 1\n\n');

        % Sort by voxel size
        fractaldimension_voxelsizedependence = sortrows(fractaldimension_voxelsizedependence,1);

        % Extrapolation to 0 voxel size
        tmp = zeros(number_resize+2,number_interface+1); % + 0 nm and + initial voxel size
        for current_interface=1:1:number_interface
            y=fractaldimension_voxelsizedependence(:,current_interface+1);
            pi = polyfit(x,y,1);
            vq = polyval(pi,0);
            tmp(1,current_interface+1)=vq;
            interpolation_voxelsize(current_interface).pi=pi;
        end
        tmp(2:end,:) = fractaldimension_voxelsizedependence;
        fractaldimension_voxelsizedependence = tmp; clear tmp;

        % Managing result
        Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
        for current_interface=1:1:number_interface
            Variable_name_table(1+current_interface)=interface_name(current_interface);
        end
        % Table
        Table_Fractaldimension_voxelsizedependence = array2table(fractaldimension_voxelsizedependence,...
            'VariableNames',Variable_name_table);
        fprintf('    fitted from s=1 to s=2\n');
        disp(Table_Fractaldimension_voxelsizedependence)

        % Save result
        Results_specificinterfacearea.Fractaldimension_voxelsizedependence = Table_Fractaldimension_voxelsizedependence; % Save in main table result
        if opt.save.xls
            filename = 'PhaseInterface_Fractaldimension_boxcounting_voxel_size_dependence'; % Filename without extension
            % Prepare the data
            clear DATA_writetable
            DATA_writetable.sheet(1).name='Fit from 1 to 2';
            DATA_writetable.sheet(1).table=Table_Fractaldimension_voxelsizedependence;
            % Save function
            Function_Writetable(Current_folder,filename,DATA_writetable)
        end

        % Correlation
        for current_interface=1:1:number_interface
            results_correlation(current_interface).PhaseSurface_fractaldimension_extrapolated = fractaldimension_voxelsizedependence(1,current_interface+1) ;
            results_correlation(current_interface).PhaseSurface_fractalpropensity_extrapolated = abs(p.fractal_boxcounting.topology_dimension - fractaldimension_voxelsizedependence(1,current_interface+1));
        end

        % Figure
        parameters_figure.plotlog = false;
        parameters_figure.figname = 'Phase Interface fractal dimension';
        parameters_figure.propertyname = 'Phase Interface fractal dimension';
        parameters_figure.method = 'Box counting';
        parameters_figure.property_voxelsizedependence = fractaldimension_voxelsizedependence;
        parameters_figure.number_phase = number_interface;
        parameters_figure.str_ylabel = 'Fractal dimension';
        parameters_figure.propertynameunit = [];
        parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
        parameters_figure.Current_folder = Current_folder;
        parameters_figure.filename = 'PhaseInterface_Fractaldimension_voxel_size_dependence';
        parameters_figure.infovol = infovol_bis;
        parameters_figure.opt = opt;
        parameters_figure.todo = p.todo;
        Function_create_figure_voxelsizedependence(parameters_figure) % Figures
    end

end

%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%

if p.RVE.number_RVE>0
    for k_RVE = 1:1:p.RVE.number_RVE % Loop over all RVE asked
        RVEparameters = p.RVE.RVE(k_RVE);
        RVE(k_RVE).RVEparameters = RVEparameters; % For result structure
        if opt.save.xls || opt.save.savefig
            Sub_folder_RVE = [Current_folder RVEparameters.savename separator];
            while exist(Sub_folder_RVE,'dir')
                RVEparameters.savename = [RVEparameters.savename '_bis'];
                Sub_folder_RVE = [Current_folder RVEparameters.savename separator];
            end
            mkdir(Sub_folder_RVE);
        end
        if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
            thresholds = RVEparameters.threshold_std;
        else
            thresholds = RVEparameters.threshold_reldiff;
        end
        n_threshold = length(thresholds);

        % Nested analysis ?
        if RVEparameters.donested % Yes
            [n_nestedRVE, xcrop] = Function_NestedRVE(RVEparameters,Domain_size);
            Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_interface,3,3); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
        else % No
            n_nestedRVE = 0;
        end        

        for k_nestedRVE = 0:1:n_nestedRVE

            %% NESTED ANALYSIS
            if k_nestedRVE == 0
                Domain_size_nested = Domain_size;
            else
                x0 = xcrop(k_nestedRVE,1); x1 = xcrop(k_nestedRVE,2);
                y0 = xcrop(k_nestedRVE,3); y1 = xcrop(k_nestedRVE,4);
                z0 = xcrop(k_nestedRVE,5); z1 = xcrop(k_nestedRVE,6);
                Phase_microstructure_nested = Phase_microstructure(x0:x1,y0:y1,z0:z1);
                Domain_size_nested = size(Phase_microstructure_nested);
            end
            
            %% SUBDOMAINS
            [All_subdomain,GROUP_SUBDOMAIN, Wholevolume_size] = Function_get_Subdomains(RVEparameters, Domain_size_nested, Domain_size); % Location of all subdomains
            Wholevolume_size(1:2) = Wholevolume_size(1:2)*voxel_size;
            if RVEparameters.donested
                Result_nestedRVE(k_nestedRVE+1,1,:,:,1) = Wholevolume_size(1);
                if k_nestedRVE ==1
                    fprintf('       Nested analysis\n');
                end
                if k_nestedRVE > 0
                    fprintf('          Cropping iteration: #%i/%i (cubic root volume=%1.1f %s)\n',k_nestedRVE,n_nestedRVE,Wholevolume_size(1),infovol.unit);
                end                
                if strcmp(RVEparameters.type,'C')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'D')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3) = Wholevolume_size(2);                    
                elseif strcmp(RVEparameters.type,'G')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'H')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2) = Wholevolume_size(2);                    
                end
            end

            % Information about subdomains
            RVE(k_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),All_subdomain(:,12),...
                'VariableNames',{'Subdomain Id' 'Group Id' 'Number subdomain' 'Equivalent cubic length' 'Equivalent section length' 'Length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
            [number_subdomain,~] = size(All_subdomain); % The number of subdomain
            number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size

            %% ALGORITHM
            % Colunm 1 is the subdomain id
            % Colunm 2 and 3 are the sizes of the subdomain.
            Property_eachsubdomain = zeros(number_subdomain,number_interface+4);
            % Property calculated for each subdomain
            for subdomain_id = 1:1:number_subdomain
                % Boundary of the subdomain
                x0 = All_subdomain(subdomain_id,7); x1 = All_subdomain(subdomain_id,8);
                y0 = All_subdomain(subdomain_id,9); y1 = All_subdomain(subdomain_id,10);
                z0 = All_subdomain(subdomain_id,11); z1 = All_subdomain(subdomain_id,12);
                clear current_subdomain;
                % Crop volume
                if k_nestedRVE == 0
                    current_subdomain = Phase_microstructure(x0:x1,y0:y1,z0:z1);
                else
                    current_subdomain = Phase_microstructure_nested(x0:x1,y0:y1,z0:z1);
                end                    
                Current_domain_size = size(current_subdomain);                

                % CPU and stopwatch time - start
                time_cpu_start_volume = cputime; % CPU start
                time_clock_start_volume = tic; % Stopwatch start
                % Number of voxel of the current resized microstructure
                voxel_number_tmp=numel(current_subdomain);

                Property_eachsubdomain(subdomain_id,1)=subdomain_id;
                % Equivalent size of the subdomain
                Property_eachsubdomain(subdomain_id,2)=All_subdomain(subdomain_id,4)*voxel_size; % Cubic root length
                Property_eachsubdomain(subdomain_id,3)=All_subdomain(subdomain_id,5)*voxel_size; % Square root length
                Property_eachsubdomain(subdomain_id,4)=All_subdomain(subdomain_id,6)*voxel_size; % Length

                n_=prod(Current_domain_size); % Number of voxel of the domain
                Volume_=n_*(voxel_size^number_dimension);

                for current_interface=1:1:number_interface % Loop over all phases
                    time_cpu_start_phase = cputime; % CPU start
                    time_clock_start_phase = tic; % Stopwatch start

                    % % Algorithm: SPECIFIC FOR EACH FILE
                    % Label
                    label_phase_A = p.doublets(current_interface,1);
                    label_phase_B = p.doublets(current_interface,2);
                    % Create a binary microstructure : 1 = current phase A, 10 = current phase B, 0 = all the other phases
                    binary_phase=zeros(Current_domain_size); % Initialization
                    binary_phase(current_subdomain == label_phase_A) = 1;
                    binary_phase(current_subdomain == label_phase_B) = 10;
                    n_voxel_phases=sum(sum(sum(binary_phase~=0))); % Number of voxel of the two phases
                    [~,surface_,~] = Function_Specificinterface_direct_Algorithm(binary_phase); % Interface area
                    surface_=surface_*(voxel_size^(number_dimension-1));
                    Int=surface_/Volume_;
                    Property_eachsubdomain(subdomain_id,current_interface+4)=Int*corrective_factor;

                    % % Time
                    timedata_perphase = [timedata_perphase; [n_voxel_phases (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

                end
                % CPU and stopwatch time - end
                timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
            end

            %% STATISTICAL ANALYSIS and RVE SIZE
            [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis(number_group_size,number_interface,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters);
            
            %% SAVE FOR CORRELATION
            if k_nestedRVE == 0
                for k_threshold=1:1:n_threshold
                    for current_interface=1:1:number_interface
                        if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                            if Size_RVE(1,current_interface,2,1)~=0
                                str_ = ['Intdirect_RVE_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                str_(str_=='.')='p';
                                results_correlation(current_interface).(str_) = Size_RVE(1,current_interface,2,1);
                            end
                            if strcmp(RVEparameters.type,'C')
                                if Size_RVE(1,current_interface,2,2)~=0
                                    str_ = ['Intdirect_RVE_squarerootFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_interface).(str_) = Size_RVE(1,current_interface,2,2);
                                end
                            elseif strcmp(RVEparameters.type,'D')
                                if Size_RVE(1,current_interface,2,2)~=0
                                    str_ = ['Intdirect_RVE_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_interface).(str_) = Size_RVE(1,current_interface,2,2);
                                end
                            end
                        else
                            if Size_convergence(1,current_interface,2,1)~=0
                                str_ = ['Intdirect_conv_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                str_(str_=='.')='p';
                                results_correlation(current_interface).(str_) = Size_convergence(1,current_interface,2,1);
                            end
                            if strcmp(RVEparameters.type,'G')
                                if Size_convergence(1,current_interface,2,2)~=0
                                    str_ = ['Intdirect_conv_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_interface).(str_) = Size_convergence(1,current_interface,2,2);
                                end
                            elseif strcmp(RVEparameters.type,'H')
                                if Size_convergence(1,current_interface,2,2)~=0
                                    str_ = ['Intdirect_conv_areaFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_interface).(str_) = Size_convergence(1,current_interface,2,2);
                                end
                            end
                        end
                    end
                end
            end            

            %% MANAGING RESULTS
            p.todo = ones(number_interface,1);
            [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,RVE,k_RVE,number_interface,number_interface,infovol_bis,p);

            if RVEparameters.donested
                for k_threshold=1:1:n_threshold
                    for current_interface=1:1:number_interface
                        if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,1,1) = Size_RVE(k_threshold,current_interface,1,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,2,1) = Size_RVE(k_threshold,current_interface,2,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,3,1) = Size_RVE(k_threshold,current_interface,3,1);
                        else
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,1,1) = Size_convergence(k_threshold,current_interface,1,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,2,1) = Size_convergence(k_threshold,current_interface,2,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,3,1) = Size_convergence(k_threshold,current_interface,3,1);
                        end
                        if strcmp(RVEparameters.type,'C')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,1,2) = Size_RVE(k_threshold,current_interface,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,2,2) = Size_RVE(k_threshold,current_interface,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,3,2) = Size_RVE(k_threshold,current_interface,3,2);
                        elseif strcmp(RVEparameters.type,'D')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,1,3) = Size_RVE(k_threshold,current_interface,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,2,3) = Size_RVE(k_threshold,current_interface,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,3,3) = Size_RVE(k_threshold,current_interface,3,2);
                        elseif strcmp(RVEparameters.type,'G')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,1,3) = Size_convergence(k_threshold,current_interface,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,2,3) = Size_convergence(k_threshold,current_interface,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,3,3) = Size_convergence(k_threshold,current_interface,3,2);
                        elseif strcmp(RVEparameters.type,'H')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,1,2) = Size_convergence(k_threshold,current_interface,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,2,2) = Size_convergence(k_threshold,current_interface,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_interface,3,2) = Size_convergence(k_threshold,current_interface,3,2);
                        end
                    end
                end
            end

            %% TEXT DISPLAY AND SAVE RESULTS
            if k_nestedRVE == 0
                propertyname=Sp_str;
                RVEparameters.disp_parRVE = true;
                Function_subdomains_display_and_save(RVE,k_RVE,RVEparameters,number_interface,number_interface,propertyname,Sub_folder_RVE,opt,infovol_bis,p);
            end

            %% FIGURES
            if k_nestedRVE == 0
                parameters_figure.propertyname = propertyname;
                parameters_figure.propertynameunit = Intunit;
                parameters_figure.RVE = RVEparameters;
                parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                parameters_figure.savefolder = Sub_folder_RVE;
                parameters_figure.number_phase = number_interface;
                parameters_figure.number_phase_todo = number_interface;
                parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
                parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
                parameters_figure.derivative_convergence = derivative_convergence;
                parameters_figure.relativedifference_convergence = relativedifference_convergence;
                parameters_figure.Size_RVE = Size_RVE;
                parameters_figure.convergence_criterion = RVEparameters.threshold_reldiff;
                parameters_figure.Size_convergence = Size_convergence;
                parameters_figure.Wholevolume_size = Wholevolume_size;
                parameters_figure.Wholevolume_results = Specificinterface(:,6);
                parameters_figure.infovol = infovol_bis;
                parameters_figure.todo = ones(number_interface,1);
                parameters_figure.opt = opt;
                Function_create_figures_RVE(parameters_figure) % Figures
            end
        end

        %% NESTED ANALYSIS RESULT
        if RVEparameters.donested
            % Table
            [RVE] = Function_nestedtable(RVE,k_RVE,RVEparameters,number_interface,propertyname,Result_nestedRVE,Sub_folder_RVE,opt,infovol_bis,p);
            % Figure
            parameters_figure.Result_nestedRVE = Result_nestedRVE;
            Function_create_figures_nestedRVE(parameters_figure) % Figures
            % Save
            RVE(k_RVE).nestedanalysis = Result_nestedRVE;
        end

        Results_specificinterfacearea.RVE.Sp = RVE; % Save in main table result
    end
end

%%
%% ENDING FUNCTION
%%

%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_specificinterfacearea.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_specificinterfacearea.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_specificinterfacearea.Table_date = Table_date;

if opt.save.xls
    % Prepare the data
    clear DATA_writetable
    % Time per volume
    DATA_writetable.sheet(1).name='Time per volume';
    DATA_writetable.sheet(1).table=Table_time_pervolume;
    % Time per phase
    DATA_writetable.sheet(2).name='Time per phase';
    DATA_writetable.sheet(2).table=Table_time_perphase;    
    % Data : Date
    DATA_writetable.sheet(3).name='Date';
    DATA_writetable.sheet(3).table=Table_date;
    % Save function
    Function_Writetable(Current_folder,'Int_direct_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'Int_direct_calculation_time', 'Specific interface area (direct)', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_Int_surface_area'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_Int_direct'],'Results_specificinterfacearea')
end

end