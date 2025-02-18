function [] = Function_Specificsurface_direct(Phase_microstructure,infovol,opt,p, foo, foo2)
% Calculate specific surface area
% Function_Specificsurface_direct(Phase_microstructure, infovol, opt, p, foo) - when used with the MATBOX toolbox
% or
% Function_Specificsurface_direct(Phase_microstructure, voxelsize, unit, corrective_factor, donotconsiderisolatedcomvol) - when used as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       corrective factor, a scalar: surface area overestimation correction (2/3 for sphere)
%       donotconsiderisolatedcomvol, boolean: if true connectivity analysis is performed on the complementary volume: isolated and unknows clusters will not be considered for the interface calculation.
%       by default: sp(phase) = surface area phase / domain volume
%       e.g.: Function_Specificsurface_direct(<your_3d_array>, 0.4, 'um', 2/3, false);

%% DEFAULT VALUES
expected_number_argument = 6;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 5 % Case for function called as: Function_Specificsurface_direct(Phase_microstructure, voxelsize, unit, corrective_factor). Standalone use.
        voxelsize = infovol; clear infovol;
        unit = opt; clear opt;
        corrective_factor = p; clear p;
        donotconsiderisolatedcomvol = foo; clear foo;
        Sp_definition = 'Phase surface divided by domain volume';
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['Spdirect_' char(t)];
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
            p.todo(k)=1;
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
        p.fractal_boxcounting.todo = false;       

    else % Incorrect number of argument
        disp 'Error calling Function_Specificsurface_direct. Wrong number of argument.'
        help Function_Specificsurface_direct
        return
    end
    
else
    % Read parameters;
    Sp_definition = p.definition;
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
Current_folder = [infovol.volpath 'Sp_direct' separator];
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
    Sp_str = 'Specific surface area';
end
number_phase = length(infovol.phaselabel); % Number of phase
voxel_number = prod(Domain_size); % Number of voxel
voxel_size = infovol.voxelsize;
voxel_unit = infovol.unit;

% Unit
Area_unit = [voxel_unit num2str(number_dimension-1)];
Volume_unit = [voxel_unit num2str(number_dimension)];
Sp_unit = [voxel_unit '-1'];

%% INITIALIZE RESULTS (USE FOR CORRELATION)
current_phase_todo = 0;
for current_phase=1:1:number_phase
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        results_correlation(current_phase_todo).name = infovol.phasename(current_phase,1);
        results_visualization(current_phase_todo).name = infovol.phasename(current_phase,1);
    end
end

%% PARAMETERS
if strcmp(Sp_definition,'Phase surface divided by phase volume')
    refvolume = 2;
    refvalue = 4;
    str_figure = 'Phase_{surface} / Phase_{volume}';
    %str_correlation = 'Sp_phase_surface_over_phase_volume_directcounting';
elseif strcmp(Sp_definition,'Phase surface divided by domain volume')
    refvolume = 5;
    refvalue = 7;    
    str_figure = 'Phase_{surface} / Domain_{volume}';
    %str_correlation = 'Sp_phase_surface_over_domain_volume_directcounting';
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

disp '    SPECIFIC SURFACE AREA - DIRECT METHOD';
disp '    -------------------------------------';
disp ' ';

%% CALCULATION
time_cpu_start_volume = cputime; % CPU start
time_clock_start_volume = tic; % Stopwatch start

% Initialization (generic)
number_phase_todo = sum(p.todo); % How many phase are we going to analyse ?
timedata = zeros(number_phase_todo+1,3); timedata(1,1) = voxel_number;
timedata_domain = cell(number_phase_todo+1,1);
phasename_todo = cell(number_phase_todo,1);
phaselabel = zeros(number_phase_todo,1);

% Initialization (algorithm-specific)
% Colunm 1: Surface area of the phase
% Column 2: Volume of the phase
% Column 3: Specific surface area (surface phase/volume phase)
% Column 4: Specific surface area (surface phase/volume phase) * corrective_factor
% Column 5: Volume of the domain
% Column 6: Specific surface area (surface phase/volume domain)
% Column 7: Specific surface area (surface phase/volume domain) * corrective_factor
Specificsurface_phase=zeros(number_phase_todo,7); % Initialization
Specificsurface_phase_removeisolatedcompvol=zeros(number_phase_todo,7); % Initialization
ratio_comp = zeros(number_phase_todo,3); % Ratio after/before removing complementary isolated and unknow volume: volume ratio, sp ratio (/volume phase), sp ratio (/volume domain)
current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        current_phase_todo=current_phase_todo+1;
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);

        % % Algorithm
        phaselabel(current_phase_todo,1) = infovol.phaselabel(current_phase);
        % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
        binary_phase=zeros(Domain_size); % Initialization
        binary_phase(Phase_microstructure == phaselabel(current_phase_todo,1)) = 1; % Binary phase
        n_voxel_phase=sum(sum(sum(binary_phase==1))); % Number of voxel of the phase        
        Specificsurface_phase(current_phase_todo,2)=n_voxel_phase*(voxel_size^number_dimension); % Phase volume

        [Position_surface.phase(current_phase_todo),n_voxel_phasesurface, LOC_surface.phase(current_phase_todo).array] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
        Specificsurface_phase(current_phase_todo,1)=n_voxel_phasesurface*(voxel_size^(number_dimension-1)); % Phase surface area
        Specificsurface_phase(current_phase_todo,3)=Specificsurface_phase(current_phase_todo,1)/Specificsurface_phase(current_phase_todo,2); % Specific surface area
        Specificsurface_phase(current_phase_todo,4)=Specificsurface_phase(current_phase_todo,3)*corrective_factor; % Specific surface area * corrective factor
        Specificsurface_phase(current_phase_todo,5)=voxel_number*(voxel_size^number_dimension); % Domain volume
        Specificsurface_phase(current_phase_todo,6)=Specificsurface_phase(current_phase_todo,1)/Specificsurface_phase(current_phase_todo,5); % Specific surface area
        Specificsurface_phase(current_phase_todo,7)=Specificsurface_phase(current_phase_todo,6)*corrective_factor; % Specific surface area * corrective factor
        % % Correlation
        results_correlation(current_phase_todo).Sp_phase_surface_over_phase_volume_directcount = Specificsurface_phase(current_phase_todo,4);
        results_correlation(current_phase_todo).Sp_phase_surface_over_domain_volume_directcount = Specificsurface_phase(current_phase_todo,7);
        % % Visualization
        results_visualization(current_phase_todo).Sp_phase_directcounting = LOC_surface.phase(current_phase_todo).array;

        if p.donotconsiderisolatedcompvol
            Specificsurface_phase_removeisolatedcompvol(current_phase_todo,2)=Specificsurface_phase(current_phase_todo,2);
            [binary_phase, ratio_comp(current_phase_todo,1)] = removeisolatedcompvol(binary_phase);
            [Position_surface.phase(current_phase_todo),n_voxel_phasesurface, LOC_surface.phase(current_phase_todo).array] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
            Specificsurface_phase_removeisolatedcompvol(current_phase_todo,1)=n_voxel_phasesurface*(voxel_size^(number_dimension-1)); % Phase surface area
            Specificsurface_phase_removeisolatedcompvol(current_phase_todo,3)=Specificsurface_phase_removeisolatedcompvol(current_phase_todo,1)/Specificsurface_phase_removeisolatedcompvol(current_phase_todo,2); % Specific surface area
            Specificsurface_phase_removeisolatedcompvol(current_phase_todo,4)=Specificsurface_phase_removeisolatedcompvol(current_phase_todo,3)*corrective_factor; % Specific surface area * corrective factor
            Specificsurface_phase_removeisolatedcompvol(current_phase_todo,5)=voxel_number*(voxel_size^number_dimension); % Domain volume
            Specificsurface_phase_removeisolatedcompvol(current_phase_todo,6)=Specificsurface_phase_removeisolatedcompvol(current_phase_todo,1)/Specificsurface_phase_removeisolatedcompvol(current_phase_todo,5); % Specific surface area
            Specificsurface_phase_removeisolatedcompvol(current_phase_todo,7)=Specificsurface_phase_removeisolatedcompvol(current_phase_todo,6)*corrective_factor; % Specific surface area * corrective factor
            ratio_comp(current_phase_todo,2) = Specificsurface_phase_removeisolatedcompvol(current_phase_todo,4) / Specificsurface_phase(current_phase_todo,4);
            ratio_comp(current_phase_todo,3) = Specificsurface_phase_removeisolatedcompvol(current_phase_todo,7) / Specificsurface_phase(current_phase_todo,7);

            % % Correlation
            results_correlation(current_phase_todo).Sp_phase_surface_over_phase_volume_directcount_removedcompvol = Specificsurface_phase_removeisolatedcompvol(current_phase_todo,4);
            results_correlation(current_phase_todo).Sp_phase_surface_over_domain_volume_directcount_removedcompvol = Specificsurface_phase_removeisolatedcompvol(current_phase_todo,7);
            results_correlation(current_phase_todo).Sp_phase_surface_over_phase_volume_directcount_compvolratio = ratio_comp(current_phase_todo,1);
            results_correlation(current_phase_todo).Sp_phase_surface_over_phase_volume_directcount_ratio = ratio_comp(current_phase_todo,2);
            results_correlation(current_phase_todo).Sp_phase_surface_over_domain_volume_directcount_ratio = ratio_comp(current_phase_todo,3);
            % % Visualization
            results_visualization(current_phase_todo).Sp_phase_directcounting_removedcompvol = LOC_surface.phase(current_phase_todo).array;
        end

        % % Time
        timedata_domain(current_phase_todo+1,1) = infovol.phasename(current_phase,1);
        timedata(current_phase_todo+1,1) = n_voxel_phase;
        timedata(current_phase_todo+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
        timedata(current_phase_todo+1,3) = toc(time_clock_start_phase); % Stopwatch elapsed time   

    end
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
Results_specificsurfacearea.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_Sp_phasevolume_direct = table(phasename_todo,phaselabel,Specificsurface_phase(:,1),Specificsurface_phase(:,2),Specificsurface_phase(:,3),ones(number_phase_todo,1)*corrective_factor,Specificsurface_phase(:,4),...
    'VariableNames',{'Phase' 'Label' [Area_str ' ' Area_unit] [Volume_str ' phase ' Volume_unit] [Sp_str ' ' Sp_unit] 'Corrective factor' ['Corrected ' Sp_str ' ' Sp_unit]});%
Results_specificsurfacearea.Table_Sp_phasevolume_direct = Table_Sp_phasevolume_direct; % Save in main table result
Table_Sp_domainvolume_direct = table(phasename_todo,phaselabel,Specificsurface_phase(:,1),Specificsurface_phase(:,5),Specificsurface_phase(:,6),ones(number_phase_todo,1)*corrective_factor,Specificsurface_phase(:,7),...
    'VariableNames',{'Phase' 'Label' [Area_str ' ' Area_unit] [Volume_str ' domain ' Volume_unit] [Sp_str ' ' Sp_unit] 'Corrective factor' ['Corrected ' Sp_str ' ' Sp_unit]});%
Results_specificsurfacearea.Table_Sp_domainvolume_direct = Table_Sp_domainvolume_direct; % Save in main table result
if p.donotconsiderisolatedcompvol
    Table_Sp_phasevolume_direct_removedcompvol = table(phasename_todo,phaselabel,Specificsurface_phase_removeisolatedcompvol(:,1),Specificsurface_phase_removeisolatedcompvol(:,2),Specificsurface_phase_removeisolatedcompvol(:,3),ones(number_phase_todo,1)*corrective_factor,Specificsurface_phase_removeisolatedcompvol(:,4),...
        'VariableNames',{'Phase' 'Label' [Area_str ' ' Area_unit] [Volume_str ' phase ' Volume_unit] [Sp_str ' ' Sp_unit] 'Corrective factor' ['Corrected ' Sp_str ' ' Sp_unit]});%
    Results_specificsurfacearea.Table_Sp_phasevolume_direct_removedcompvol = Table_Sp_phasevolume_direct_removedcompvol; % Save in main table result
    Table_Sp_domainvolume_direct_removedcompvol = table(phasename_todo,phaselabel,Specificsurface_phase_removeisolatedcompvol(:,1),Specificsurface_phase_removeisolatedcompvol(:,5),Specificsurface_phase_removeisolatedcompvol(:,6),ones(number_phase_todo,1)*corrective_factor,Specificsurface_phase_removeisolatedcompvol(:,7),...
        'VariableNames',{'Phase' 'Label' [Area_str ' ' Area_unit] [Volume_str ' domain ' Volume_unit] [Sp_str ' ' Sp_unit] 'Corrective factor' ['Corrected ' Sp_str ' ' Sp_unit]});%
    Results_specificsurfacearea.Table_Sp_domainvolume_direct_removedcompvol = Table_Sp_domainvolume_direct_removedcompvol; % Save in main table result
    Table_Sp_ratio = table(phasename_todo,phaselabel,ratio_comp(:,1),ratio_comp(:,2),ratio_comp(:,3),...
        'VariableNames',{'Phase' 'Label' 'Complementary volume ratio' 'Sp ratio (phase volume)' 'Sp ratio (domain volume)'});
    Results_specificsurfacearea.Table_Sp_ratio = Table_Sp_ratio; % Save in main table result
end

%% SAVE TABLES
if opt.save.xls
    filename = 'Sp_direct'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Phasesurface over Phasevolume';
    DATA_writetable.sheet(1).table=Table_Sp_phasevolume_direct;
    DATA_writetable.sheet(2).name='Phasesurface over Domainvolume';
    DATA_writetable.sheet(2).table=Table_Sp_domainvolume_direct;    
    if p.donotconsiderisolatedcompvol
        DATA_writetable.sheet(3).name='(Phase volume)no isolated comp';
        DATA_writetable.sheet(3).table=Table_Sp_phasevolume_direct_removedcompvol;
        DATA_writetable.sheet(4).name='(Domain volume)no isolated comp';
        DATA_writetable.sheet(4).table=Table_Sp_domainvolume_direct_removedcompvol;
        DATA_writetable.sheet(5).name='Ratio';
        DATA_writetable.sheet(5).table=Table_Sp_ratio;        
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
disp([Sp_str ' = phase ' Area_str ' / phase ' Volume_str])
disp(Table_Sp_phasevolume_direct)
if p.donotconsiderisolatedcompvol
    disp('   Isolated complementary volume removed')
    disp(Table_Sp_phasevolume_direct_removedcompvol)
end
disp([Sp_str ' = phase ' Area_str ' / domain ' Volume_str])
disp(Table_Sp_domainvolume_direct)
if p.donotconsiderisolatedcompvol
    disp('   Isolated complementary volume removed')
    disp(Table_Sp_domainvolume_direct_removedcompvol)
    disp('Ratio')
    disp(Table_Sp_ratio)    
end
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

if p.donotconsiderisolatedcompvol 
    Specificsurface_phase = Specificsurface_phase_removeisolatedcompvol;
end

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
for current_phase_todo=1:1:number_phase_todo
    for direction=1:1:number_dimension
        % Set position in real unit
        Position_surface.phase(current_phase_todo).direction(direction).value(:,1)=Position_surface.phase(current_phase_todo).direction(direction).value(:,1)*voxel_size;
        % Set face surface in real unit
        Position_surface.phase(current_phase_todo).direction(direction).value(:,3)=Position_surface.phase(current_phase_todo).direction(direction).value(:,2)*(voxel_size^(number_dimension-1));
        %  Set specific surface area in real unit
        Position_surface.phase(current_phase_todo).direction(direction).value(:,3)=Position_surface.phase(current_phase_todo).direction(direction).value(:,3)/Specificsurface_phase(current_phase_todo,refvolume);
        % Corrective factor
        Position_surface.phase(current_phase_todo).direction(direction).value(:,3)=Position_surface.phase(current_phase_todo).direction(direction).value(:,3)*corrective_factor;
    end
end

%% DISTRIBUTION
% Cumulative function sum(g(x))
for current_phase_todo=1:1:number_phase_todo
    for direction=1:1:number_dimension
        Cumulative_surface.phase(current_phase_todo).direction(direction).value=zeros(1,2*Domain_size(direction)-1);
        for i=1:1:2*Domain_size(direction)-1
            for j=1:1:i
                Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,j)=Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,j)+Position_surface.phase(current_phase_todo).direction(direction).value(i,3);
            end
        end
    end
end
% Derivate function (i.e. 1/D*f(x))
for current_phase_todo=1:1:number_phase_todo
    for direction=1:1:number_dimension
        c_length = 2*Domain_size(direction)-1;
        for i=2:1:(c_length-1)
            Derivate_Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,i)=-( Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,i+1)-Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,i-1) ) / (2*0.5*voxel_size);
        end
        Derivate_Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,1)=-(Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,2)-Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,1))/(0.5*voxel_size);
        Derivate_Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,c_length)=-(Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,c_length)-Cumulative_surface.phase(current_phase_todo).direction(direction).value(1,c_length-1))/(0.5*voxel_size);
    end
end
% The specific surface area f(x)
for current_phase_todo=1:1:number_phase_todo
    for direction=1:1:number_dimension
        Specific_surface.phase(current_phase_todo).direction(direction).value = zeros(1,2*Domain_size(direction)-1);
        Specific_surface.phase(current_phase_todo).direction(direction).value = Derivate_Cumulative_surface.phase(current_phase_todo).direction(direction).value*Domain_size(direction)*(voxel_size);
    end
end

%% VERIFICATION
% the relation Sp = 1/D*integral(f(x),dx,0,D) is checked
Specificsurface_integralvalue=zeros(number_phase_todo,7);
for current_phase_todo=1:1:number_phase_todo
    for direction=1:1:number_dimension
        X_= Position_surface.phase(current_phase_todo).direction(direction).value(:,1);
        Y_= Specific_surface.phase(current_phase_todo).direction(direction).value(1,:);
        Specificsurface_integralvalue(current_phase_todo,1)=Specificsurface_phase(current_phase_todo,refvalue); % Reference value
        Specificsurface_integralvalue(current_phase_todo,direction+1)=(1/(Domain_size(direction)*(voxel_size)))*trapz(X_,Y_); % Integral value
        % Relative error in percent
        Specificsurface_integralvalue(current_phase_todo,direction+1+3)=100*(Specificsurface_integralvalue(current_phase_todo,direction+1)-Specificsurface_integralvalue(current_phase_todo,1))/Specificsurface_integralvalue(current_phase_todo,1);
    end
end

%% MANAGING RESULTS
% Results are saved in a table

clear Variable_name_table;
Variable_name_table(1)={['Position ' voxel_unit]};
for current_phase_todo=1:1:number_phase_todo
    Variable_name_table(current_phase_todo+1)=phasename_todo(current_phase_todo);
end

% g(x)
for direction=1:1:number_dimension
    array_g = zeros(2*Domain_size(direction)-1,number_phase_todo+1);
    array_g(:,1)= Position_surface.phase(current_phase_todo).direction(direction).value(:,1);
    for current_phase_todo=1:1:number_phase_todo
        array_g(:,current_phase_todo+1)= Position_surface.phase(current_phase_todo).direction(direction).value(:,3);
    end
    Gx.direction(direction).table = array2table(array_g,'VariableNames',Variable_name_table);
end
% f(x)
for direction=1:1:number_dimension
    array_f = zeros(2*Domain_size(direction)-1,number_phase_todo+1);
    array_f(:,1)= Position_surface.phase(current_phase_todo).direction(direction).value(:,1);
    for current_phase_todo=1:1:number_phase_todo
        array_f(:,current_phase_todo+1)= Specific_surface.phase(current_phase_todo).direction(direction).value(1,:);
    end
    Fx.direction(direction).table = array2table(array_f,'VariableNames',Variable_name_table);
end
 
% Verification
Table_verification_integral = table(phasename_todo,Specificsurface_integralvalue(:,1),Specificsurface_integralvalue(:,2),Specificsurface_integralvalue(:,3),Specificsurface_integralvalue(:,4),Specificsurface_integralvalue(:,5),Specificsurface_integralvalue(:,6),Specificsurface_integralvalue(:,7),...
    'VariableNames',{'Phase' 'Reference' 'Integration direction 1' 'Integration direction 2' 'Integration direction 3' 'Relative error percent direction 1' 'Relative error percent direction 2' 'Relative error percent direction 3'});
fprintf('> Calculated along direction:\n\n');
fprintf('  Integral verification: relative error should be close to 0 percent:\n\n');
disp(Table_verification_integral)

Results_specificsurfacearea.Table_evolution_Gx = Gx; % Save in main table result
Results_specificsurfacearea.Table_evolution_Fx = Fx;
Results_specificsurfacearea.Table_verification_integral=Table_verification_integral;

%% SAVE TABLES
if opt.save.xls
    filename = 'Sp_direct_along_directions'; % Filename without extension
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
    Spunit = '\mum^{-1}';
else
    axisunit = ['(' strunit ')'];
    Spunit = ['' voxel_unit ' ^{-1}'];
end

scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Specific surface area, direct method (integral definition)'; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
for current_direction=1:1:number_dimension % Iterate over axe
    sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    h_title=title ({'Specific surface area, direct method',sprintf('S_{p} = %s ',str_figure)}); % Set title font
    % Plot graphs
    current_phase_todo = 0;
    for current_phase=1:1:number_phase % Loop over phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;
            x_ = Position_surface.phase(current_phase_todo).direction(current_direction).value(:,1);
            y_ = Specific_surface.phase(current_phase_todo).direction(current_direction).value(1,:);
            plot(x_,y_,'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
        end
    end
    % Axis label
    t_ = xlabel(' ');
    t_1 = sprintf('Position along %s ',char(infovol.directionname(current_direction)));
    t_2 = axisunit;
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = ylabel(' '); 
    str1 = 'Specific surface distribution function f(x), such as';
    str2 = ['<S_{p}> = 1/L \times \int_{0}^{L} f(x) dx (' Spunit ')'];
    t_.String= {str1,str2};
    % Legend
    for current_phase_todo=1:1:number_phase_todo
        str_legend(current_phase_todo).name = [char(phasename_todo(current_phase_todo)) ', <S_{p}>=' num2str(Specificsurface_phase(current_phase_todo,refvalue),'%1.3f') Spunit];
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
sgtitle(Fig,'Specific surface area along directions','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
if opt.save.savefig % Save figure
    filename= 'Specific_surface_direct_along_direction';
    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
end
if opt.format.autoclosefig
    close(Fig); % Do not keep open figures
end

%% FRACTAL DIMENSION (COUNTING BOX)
if p.fractal_boxcounting.todo
    fprintf('> Phase surface fractal dimension (box-counting method)\n');
    p.fractal_boxcounting.topology_dimension = number_dimension-1;
    p.fractal_boxcounting.plot = false;
    for current_phase_todo=1:1:number_phase_todo % Loop over phase and call box counting algorithm
        binary_array = zeros(size(Phase_microstructure));
        binary_array(Phase_microstructure==phaselabel(current_phase_todo,1))=1;
        [N(:,current_phase_todo),box_lengths,fractal_dimension(:,current_phase_todo),fractal_dimension_convergence(:,:,current_phase_todo)] = Function_fractaldimension_boxcounting(LOC_surface.phase(current_phase_todo).array,p.fractal_boxcounting);
    end

    % Table
    Table_fractaldimension = table(phasename_todo,fractal_dimension(1,:)',fractal_dimension(2,:)',fractal_dimension(3,:)',fractal_dimension(4,:)',...
        'VariableNames',{'Phase' 'Fit from 1 to' 'Fractal dimension' 'Topology dimension' 'Fractal propensity'});
    for current_phase_todo=1:1:number_phase_todo
        Table_boxlength(current_phase_todo).t = table(fractal_dimension_convergence(:,1,current_phase_todo),fractal_dimension_convergence(:,2,current_phase_todo),fractal_dimension_convergence(:,3,current_phase_todo),...
        'VariableNames',{'Fit from 1 to' 'Fractal dimension' 'Fit norm error'});
    end
    Results_specificsurfacearea.Table_fractaldimension = Table_fractaldimension; % Save in main table result
    disp(Table_fractaldimension)

    % Save table
    if opt.save.xls
        filename = 'PhaseSurface_Fractaldimension_boxcounting'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name = 'Fractal dimension';
        DATA_writetable.sheet(1).table = Table_fractaldimension;
        for current_phase_todo=1:1:number_phase_todo
            DATA_writetable.sheet(1+current_phase_todo).name = char(phasename_todo(current_phase_todo));
            DATA_writetable.sheet(1+current_phase_todo).table = Table_boxlength(current_phase_todo).t;
        end           
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end

    % Figure
    if p.fractal_boxcounting.topology_dimension==2
        propertyname= 'Phase surface (area)'; % Figure name
    else
        propertyname = 'Phase surface (line)'; % Figure name
    end
    filename = 'PhaseSurface_fractal_dimension_boxcounting';
    function_fractalfig(propertyname,filename,Current_folder, box_lengths,N,fractal_dimension_convergence,p.fractal_boxcounting.topology_dimension,number_phase,p,opt,infovol);

    % Correlation
    results_correlation(current_phase_todo).PhaseSurface_fractaldimension_boxcounting = fractal_dimension(2,current_phase_todo);
    results_correlation(current_phase_todo).PhaseSurface_fractalpropensity_boxcounting = abs(p.fractal_boxcounting.topology_dimension - fractal_dimension(2,current_phase_todo));
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
    property_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1);
    property_voxelsizedependence(1,1)=voxel_size;
    property_voxelsizedependence(1,2:end)=Specificsurface_phase(:,refvalue)';

    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fractaldimension_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1);
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

        current_phase_todo = 0;
        for current_phase=1:1:number_phase % Loop over all phases
            if p.todo(current_phase)
                time_cpu_start_phase = cputime; % CPU start
                time_clock_start_phase = tic; % Stopwatch start

                current_phase_todo=current_phase_todo+1;

                % % Algorithm: SPECIFIC FOR EACH FILE
                code_tmp = infovol.phaselabel(current_phase);
                Numbervoxel_phase_tmp= sum(sum(sum(Phase_microstructure_resized==code_tmp )));
                binary_phase=zeros(Current_domain_size); % Binary phase
                binary_phase(Phase_microstructure_resized == code_tmp) = 1;
                if strcmp(Sp_definition,'Phase surface divided by phase volume')
                    n_=sum(sum(sum(binary_phase==1))); % Number of voxel of the phase
                else
                    n_=prod(Current_domain_size); % Number of voxel of the domain
                end
                Volume_=n_*(current_voxel_size^number_dimension);

                if p.donotconsiderisolatedcompvol
                    [binary_phase, ~] = removeisolatedcompvol(binary_phase);
                end

                [~,surface_,LOC_surface_resized.phase(current_phase_todo).array] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
                surface_=surface_*(current_voxel_size^(number_dimension-1)); % Surface area in square micrometer
                Specificsurface_=surface_/Volume_; % Specific surface area
                property_voxelsizedependence(current_iteration+1,current_phase_todo+1)=Specificsurface_*corrective_factor;

                % % Time
                timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

            end
        end
        % CPU and stopwatch time - end
        timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];

        % % Fractal
        if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
            for current_phase_todo=1:1:number_phase_todo % Loop over phase and call box counting algorithm
                [~,~,fractal_dimension_tmp, ~] = Function_fractaldimension_boxcounting(LOC_surface_resized.phase(current_phase_todo).array,p.fractal_boxcounting);
                fractaldimension_voxelsizedependence(current_iteration+1,1)=size_choice(current_iteration)*voxel_size;
                fractaldimension_voxelsizedependence(current_iteration+1,current_phase_todo+1) = fractal_dimension_tmp(2,1);
            end      
        end


    end
    clear Phase_microstructure_resized;

    % Sort per voxel size
    property_voxelsizedependence = sortrows(property_voxelsizedependence,1);   
    
    %% EXTRAPOLATION TO 0 nm
    fprintf('> Specific surface area dependence with the voxel size (with corrective factor)\n');
    tmp = zeros(number_resize+2,number_phase_todo+1); % + 0 nm and + initial voxel size
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
    for current_phase_todo=1:1:number_phase_todo
        y=property_voxelsizedependence(:,current_phase_todo+1);
        pi = polyfit(x,y,interpolation_voxelsize_order);
        vq = polyval(pi,0);
        tmp(1,current_phase_todo+1)=vq;
        interpolation_voxelsize(current_phase_todo).pi=pi;
    end
    tmp(2:end,:) = property_voxelsizedependence;
    property_voxelsizedependence = tmp; clear tmp;

    if p.fractal_bertei.todo
        fractal_dimension_Bertei = zeros(number_phase_todo,3);
        % Fractal dimension according to Bertei et al., https://doi.org/10.1016/j.nanoen.2017.06.028 (Richardson/Mandelbrot formula)
        % Log(property) = m + (2-fractal dimension)*log(voxel size)
        logV = log(property_voxelsizedependence(2:end,1));
        logPs = zeros(length(logV),number_phase_todo);
        for current_phase_todo=1:1:number_phase_todo
            logP = log(property_voxelsizedependence(2:end,current_phase_todo+1));
            logPs(:,current_phase_todo)=logP;
            pf = polyfit(logV,logP,1);
            fractal_dimension_Bertei(current_phase_todo,1) = 2-pf(1);
            fractal_dimension_Bertei(current_phase_todo,2) = 2; % Topology dimension, surface
            fractal_dimension_Bertei(current_phase_todo,3) = abs( fractal_dimension_Bertei(current_phase_todo,2) - fractal_dimension_Bertei(current_phase_todo,1) ); % Fractal propensity
            % Correlation
            results_correlation(current_phase_todo).PhaseSurface_fractaldimension_Mandelbrot = fractal_dimension_Bertei(current_phase_todo,1);
            results_correlation(current_phase_todo).PhaseSurface_fractalpropensity_Mandelbrot = fractal_dimension_Bertei(current_phase_todo,3);
        end
    end    
    
    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
    for current_phase_todo=1:1:number_phase_todo
        Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
    end
    % Table
    Table_Spdirect_voxelsizedependence = array2table(property_voxelsizedependence,...
        'VariableNames',Variable_name_table);
    if p.fractal_bertei.todo
        Variable_name_table={'Voxel size log'}; % Columns name
        for current_phase_todo=1:1:number_phase_todo
            Variable_name_table(1+current_phase_todo)={[char(phasename_todo(current_phase_todo)) ' log']};
        end
        Table_Spdirect_voxelsizedependence_loglog = array2table([logV logPs],'VariableNames',Variable_name_table);
        Table_Fractaldimension_Bertei = table(phasename_todo(:,1),fractal_dimension_Bertei(:,1),fractal_dimension_Bertei(:,2),fractal_dimension_Bertei(:,3),...
            'VariableNames',{'Phase','Fractal dimension','Topology dimension','Fractal propensity'});   
    end        
    
    %% DISPLAY TEXT RESULTS
    disp(Table_Spdirect_voxelsizedependence)
    if p.fractal_bertei.todo
        disp(Table_Spdirect_voxelsizedependence_loglog)
        fprintf('Richardson/Mandelbrot formula: Log(property) = m + (2-fractal dimension)*log(voxel size))\n');
        disp(Table_Fractaldimension_Bertei)
    end        
    
    %% SAVE RESULTS
    Results_specificsurfacearea.voxelsizedependence = Table_Spdirect_voxelsizedependence; % Save in main table result
    if opt.save.xls
        filename = 'Sp_direct_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Specific surface area';
        DATA_writetable.sheet(1).table=Table_Spdirect_voxelsizedependence;
        if p.fractal_bertei.todo
            DATA_writetable.sheet(2).name='Log log';
            DATA_writetable.sheet(2).table=Table_Spdirect_voxelsizedependence_loglog;
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
    parameters_figure.number_phase = number_phase_todo;
    parameters_figure.str_ylabel = [sprintf('S_{p} = %s ',str_figure) '(' Spunit ')'];
    parameters_figure.propertynameunit = Spunit;
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'Sp_direct_voxel_size_dependence';
    parameters_figure.infovol = infovol;
    parameters_figure.opt = opt;
    parameters_figure.todo = p.todo;
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    

    %% FRACTAL DIMENSION
    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fprintf('> Phase surface fractal dimension (box counting) dependence with the voxel size\n');
        fprintf('  Extrapolation to zero voxel size: polynomial of order 1\n\n');

        % Sort by voxel size
        fractaldimension_voxelsizedependence = sortrows(fractaldimension_voxelsizedependence,1);

        % Extrapolation to 0 voxel size
        tmp = zeros(number_resize+2,number_phase_todo+1); % + 0 nm and + initial voxel size
        for current_phase_todo=1:1:number_phase_todo
            y=fractaldimension_voxelsizedependence(:,current_phase_todo+1);
            pi = polyfit(x,y,1);
            vq = polyval(pi,0);
            tmp(1,current_phase_todo+1)=vq;
            interpolation_voxelsize(current_phase_todo).pi=pi;           
        end
        tmp(2:end,:) = fractaldimension_voxelsizedependence;
        fractaldimension_voxelsizedependence = tmp; clear tmp;

        % Managing result
        Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
        for current_phase_todo=1:1:number_phase_todo
            Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
        end
        % Table
        Table_Fractaldimension_voxelsizedependence = array2table(fractaldimension_voxelsizedependence,...
            'VariableNames',Variable_name_table);
        fprintf('    fitted from s=1 to s=2\n');
        disp(Table_Fractaldimension_voxelsizedependence)

        % Save result
        Results_specificsurfacearea.Fractaldimension_voxelsizedependence = Table_Fractaldimension_voxelsizedependence; % Save in main table result
        if opt.save.xls
            filename = 'PhaseSurface_Fractaldimension_boxcounting_voxel_size_dependence'; % Filename without extension
            % Prepare the data
            clear DATA_writetable
            DATA_writetable.sheet(1).name='Fit from 1 to 2';
            DATA_writetable.sheet(1).table=Table_Fractaldimension_voxelsizedependence;
            % Save function
            Function_Writetable(Current_folder,filename,DATA_writetable)
        end

        % Correlation
        for current_phase_todo=1:1:number_phase_todo
            results_correlation(current_phase_todo).PhaseSurface_fractaldimension_extrapolated = fractaldimension_voxelsizedependence(1,current_phase_todo+1) ;
            results_correlation(current_phase_todo).PhaseSurface_fractalpropensity_extrapolated = abs(p.fractal_boxcounting.topology_dimension - fractaldimension_voxelsizedependence(1,current_phase_todo+1));
        end

        % Figure
        parameters_figure.plotlog = false; 
        parameters_figure.figname = 'Phase surface fractal dimension';
        parameters_figure.propertyname = 'Phase surface fractal dimension';
        parameters_figure.method = 'Box counting';
        parameters_figure.property_voxelsizedependence = fractaldimension_voxelsizedependence;
        parameters_figure.number_phase = number_phase_todo;
        parameters_figure.str_ylabel = 'Fractal dimension';
        parameters_figure.propertynameunit = [];
        parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
        parameters_figure.Current_folder = Current_folder;
        parameters_figure.filename = 'PhaseSurface_Fractaldimension_voxel_size_dependence';
        parameters_figure.infovol = infovol;
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
            Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_phase_todo,3,3); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
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
            Property_eachsubdomain = zeros(number_subdomain,number_phase_todo+4);
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

                current_phase_todo = 0;
                for current_phase=1:1:number_phase % Loop over all phases
                    if p.todo(current_phase)
                        time_cpu_start_phase = cputime; % CPU start
                        time_clock_start_phase = tic; % Stopwatch start
                        current_phase_todo=current_phase_todo+1;

                        % % Algorithm: SPECIFIC FOR EACH FILE
                        code_tmp = infovol.phaselabel(current_phase);
                        Numbervoxel_phase_tmp= sum(sum(sum(current_subdomain==code_tmp )));
                        binary_phase=zeros(Current_domain_size); % Initialization
                        binary_phase(current_subdomain == code_tmp) = 1;
                        if strcmp(Sp_definition,'Phase surface divided by phase volume')
                            n_=sum(sum(sum(binary_phase==1))); % Number of voxel of the phase
                        else
                            n_=prod(Current_domain_size); % Number of voxel of the domain
                        end
                        Volume_=n_*(voxel_size^number_dimension);

                        if p.donotconsiderisolatedcompvol
                            [binary_phase, ~] = removeisolatedcompvol(binary_phase);
                        end

                        [~,surface_,~] = Function_Specificsurface_direct_Algorithm(binary_phase); % Surface area
                        surface_=surface_*(voxel_size^(number_dimension-1)); % Surface area
                        Specificsurface_=surface_/Volume_; % Specific surface area
                        Property_eachsubdomain(subdomain_id,current_phase_todo+4)=Specificsurface_*corrective_factor;                        

                        % % Time
                        timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

                    end
                end
                % CPU and stopwatch time - end
                timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
            end

            %% STATISTICAL ANALYSIS and RVE SIZE
            [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis(number_group_size,number_phase_todo,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters);
            
            %% SAVE FOR CORRELATION
            if k_nestedRVE == 0
                for k_threshold=1:1:n_threshold
                    current_phase_todo = 0;
                    for current_phase=1:1:number_phase
                        if p.todo(current_phase)
                            current_phase_todo=current_phase_todo+1;
                            if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                if Size_RVE(1,current_phase_todo,2,1)~=0
                                    str_ = ['Spdirect_RVE_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,1);
                                end
                                if strcmp(RVEparameters.type,'C')
                                    if Size_RVE(1,current_phase_todo,2,2)~=0
                                        str_ = ['Spdirect_RVE_squarerootFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                    end
                               elseif strcmp(RVEparameters.type,'D')
                                    if Size_RVE(1,current_phase_todo,2,2)~=0
                                        str_ = ['Spdirect_RVE_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                    end
                                end   
                            else
                                if Size_convergence(1,current_phase_todo,2,1)~=0
                                    str_ = ['Spdirect_conv_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,1);
                                end
                                if strcmp(RVEparameters.type,'G')
                                    if Size_convergence(1,current_phase_todo,2,2)~=0
                                        str_ = ['Spdirect_conv_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,2);
                                    end
                                elseif strcmp(RVEparameters.type,'H')
                                    if Size_convergence(1,current_phase_todo,2,2)~=0
                                        str_ = ['Spdirect_conv_areaFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,2);
                                    end
                                end
                            end
                        end
                    end
                end
            end

            %% MANAGING RESULTS
            [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,RVE,k_RVE,number_phase,number_phase_todo,infovol,p);

            if RVEparameters.donested
                for k_threshold=1:1:n_threshold
                    current_phase_todo = 0;
                    for current_phase=1:1:number_phase
                        if p.todo(current_phase)
                            current_phase_todo=current_phase_todo+1;
                            if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1) = Size_RVE(k_threshold,current_phase_todo,1,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1) = Size_RVE(k_threshold,current_phase_todo,2,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1) = Size_RVE(k_threshold,current_phase_todo,3,1);
                            else
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1) = Size_convergence(k_threshold,current_phase_todo,1,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1) = Size_convergence(k_threshold,current_phase_todo,2,1);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1) = Size_convergence(k_threshold,current_phase_todo,3,1);
                            end
                            if strcmp(RVEparameters.type,'C')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2) = Size_RVE(k_threshold,current_phase_todo,3,2);
                            elseif strcmp(RVEparameters.type,'D')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3) = Size_RVE(k_threshold,current_phase_todo,3,2);                                    
                            elseif strcmp(RVEparameters.type,'G')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3) = Size_convergence(k_threshold,current_phase_todo,3,2);
                            elseif strcmp(RVEparameters.type,'H')
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2) = Size_convergence(k_threshold,current_phase_todo,3,2);                                
                            end
                        end
                    end
                end
            end

            %% TEXT DISPLAY AND SAVE RESULTS
            if k_nestedRVE == 0
                propertyname=Sp_str;
                RVEparameters.disp_parRVE = true;
                Function_subdomains_display_and_save(RVE,k_RVE,RVEparameters,number_phase,number_phase_todo,propertyname,Sub_folder_RVE,opt,infovol,p);
            end

            %% FIGURES
            if k_nestedRVE == 0
                parameters_figure.propertyname = propertyname;
                parameters_figure.propertynameunit = Spunit;
                parameters_figure.RVE = RVEparameters;
                parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                parameters_figure.savefolder = Sub_folder_RVE;
                parameters_figure.number_phase = number_phase;
                parameters_figure.number_phase_todo = number_phase_todo;
                parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
                parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
                parameters_figure.derivative_convergence = derivative_convergence;
                parameters_figure.relativedifference_convergence = relativedifference_convergence;
                parameters_figure.Size_RVE = Size_RVE;
                parameters_figure.convergence_criterion = RVEparameters.threshold_reldiff;
                parameters_figure.Size_convergence = Size_convergence;
                parameters_figure.Wholevolume_size = Wholevolume_size;
                parameters_figure.Wholevolume_results = Specificsurface_phase(:,refvalue);
                parameters_figure.infovol = infovol;
                parameters_figure.todo = p.todo;
                parameters_figure.opt = opt;
                Function_create_figures_RVE(parameters_figure) % Figures
            end
        end

        %% NESTED ANALYSIS RESULT
        if RVEparameters.donested
            % Table
            [RVE] = Function_nestedtable(RVE,k_RVE,RVEparameters,number_phase,propertyname,Result_nestedRVE,Sub_folder_RVE,opt,infovol,p);
            % Figure
            parameters_figure.Result_nestedRVE = Result_nestedRVE;
            Function_create_figures_nestedRVE(parameters_figure) % Figures
            % Save
            RVE(k_RVE).nestedanalysis = Result_nestedRVE;
        end

        Results_specificsurfacearea.RVE.Sp = RVE; % Save in main table result
    end
end

%%
%% ENDING FUNCTION
%%

%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_specificsurfacearea.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_specificsurfacearea.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_specificsurfacearea.Table_date = Table_date;

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
    Function_Writetable(Current_folder,'Sp_direct_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'Sp_direct_calculation_time', 'Specific surface area (direct)', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_Sp_surface_area'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_Sp_direct'],'Results_specificsurfacearea')
end

end

%% LOCAL FUNCTIONS
function [binary_phase, ratio] = removeisolatedcompvol(binary_phase)
n0 = sum(sum(sum(~binary_phase)));

[L,n] = bwlabeln(~binary_phase,6);
cluster_ids = unique(L);
cluster_ids(cluster_ids==0) = []; % Remove background
sz = size(binary_phase);
for kcluster = 1:1:n
    idx = find(L==cluster_ids(kcluster));
    [IX, IY, IZ] = ind2sub(sz,idx);
    if min(IX)>1 && min(IY)>1 && min(IZ)>1 && max(IX)<sz(1) && max(IY)<sz(2) && max(IZ)<sz(3)
        binary_phase(idx)=1;
    end
end

% if n~=1
%     % Number of voxel per cluster
%     [GC,GR] = groupcounts(reshape(L,[numel(L),1]));
%     % Remove background
%     tmp = [GC, GR];
%     id0 = find(GR==0);
%     tmp(id0,:)=[];
%     % Identify complemenatry of main cluster (ie either isolated or unknown clusters)
%     idmax = find(tmp(:,1)==max(tmp(:,1)));
%     tmp(idmax,:)=[];
%     idisolateds = tmp(:,2);
%     % Assign to binary phase
%     for kiso = 1:1:length(idisolateds)
%         binary_phase(L==idisolateds(kiso)) = 1;
%     end
% end

n1 = sum(sum(sum(~binary_phase)));
ratio = n1/n0;

end