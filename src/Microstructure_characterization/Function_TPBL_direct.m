function [] = Function_TPBL_direct(Phase_microstructure,infovol,opt,p, foo1, foo2)
% Calculate triple phase boundary length density
% Function_TPBL_direct(Phase_microstructure, infovol, opt, p, foo, foo) - when used with the MATBOX toolbox
% or
% Function_TPBL_direct(Phase_microstructure, triplet, voxelsize, unit, corrective_factor) - when used as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       triplet: a 1x3 array: the three labels for which the length is calculated
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       corrective factor, a scalar: surface area overestimation correction (1/1.455 is recommended)
%       TPBL(phase i, phase j, phase k) = length of contact line between phase i - phase j - phase k / domain volume
%       e.g.: Function_TPBL_direct(<your_3d_array>, [0,1,2],  0.4, 'um', 1/1.455); % for a volume with 0,1,2 being labels

%% DEFAULT VALUES
expected_number_argument = 6;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 5 % Case for function called as: Function_TPBL_direct(Phase_microstructure, triplet, voxelsize, unit, corrective_factor). Standalone use.
        triplet = infovol; clear infovol;
        voxelsize = opt; clear opt;
        unit = p; clear p;
        corrective_factor = foo1; clear foo1; clear foo2;
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['TPBLdirect_' char(t)];
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
        p.fractal_boxcounting.todo = false;  

    else % Incorrect number of argument
        disp 'Error calling Function_TPBL_direct. Wrong number of argument.'
        help Function_TPBL_direct
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
Current_folder = [infovol.volpath 'TPBL_direct' separator];
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% VOLUME INFORMATION
Domain_size = size(Phase_microstructure);
number_dimension = 3;
if length(Domain_size)==2 % Add third column with unit value
    disp 'Error calling Function_TPBL_direct. Array is 2D.'
    return
end
number_phase = length(infovol.phaselabel); % Number of phase
voxel_number = prod(Domain_size); % Number of voxel
voxel_size = infovol.voxelsize;
voxel_unit = infovol.unit;

% Unit
Length_unit = voxel_unit;
Volume_unit = [voxel_unit num2str(number_dimension)];
TPBL_unit = [voxel_unit '-2']; % density

%%
%% ALGORITHM ON WHOLE VOLUME
%%

colorline = [colororder; rand(1000,3)];

disp '    TRIPLE PHASE BOUNDARY LENGTH density - DIRECT METHOD';
disp '    ----------------------------------------------------';
disp ' ';

%% CALCULATION
time_cpu_start_volume = cputime; % CPU start
time_clock_start_volume = tic; % Stopwatch start

%% NAME INTERFACES
[number_line,~] = size(p.triplets);
line_name = cell(number_line,1);
p.todo = ones(1,number_line);

for current_line=1:1:number_line
    index_phase_A = find([infovol.phaselabel] == p.triplets(current_line,1));
    index_phase_B = find([infovol.phaselabel] == p.triplets(current_line,2));
    index_phase_C = find([infovol.phaselabel] == p.triplets(current_line,3));    
    name_phase_A = char(infovol.phasename(index_phase_A));
    name_phase_B = char(infovol.phasename(index_phase_B));
    name_phase_C = char(infovol.phasename(index_phase_C));
    line_name(current_line,1) = {[name_phase_A '-' name_phase_B '-' name_phase_C]};
    results_correlation(current_line).name = line_name(current_line,1);
    results_visualization(current_line).name = line_name(current_line,1);
end

infovol_bis = infovol;
for current_line=1:1:number_line
    infovol_bis.phasename(current_line) = line_name(current_line);
    infovol_bis.phasecolor(current_line,:) = colorline(current_line,:);
end

% Initialization (generic)
timedata = zeros(number_line+1,3); timedata(1,1) = voxel_number;
timedata_domain = cell(number_line+1,1);

% Initialization (algorithm-specific)
% Colunm 1: phase A
% Column 2: phase B
% Column 3: phase C
% Colunm 4: Contact line length
% Column 5: Volume of the domain
% Column 6: TPBL density (line length area/Volume domain)
% Column 7: TPBL density * corrective_factor
TPBLdensity=zeros(number_line,7); % Initialization
TPBLdensity(:,5) = voxel_number*(voxel_size^number_dimension); % Domain volume
for current_line=1:1:number_line % Loop over all interfaces
    time_cpu_start_phase = cputime; % CPU start
    time_clock_start_phase = tic; % Stopwatch start

    % % Algorithm
    % Label
    label_phase_A = p.triplets(current_line,1);
    label_phase_B = p.triplets(current_line,2);
    label_phase_C = p.triplets(current_line,3);
    TPBLdensity(current_line,1)=label_phase_A;
    TPBLdensity(current_line,2)=label_phase_B;
    TPBLdensity(current_line,3)=label_phase_C;
    % Create a binary microstructure : 1 = current phase A, 10 = current phase B, 0 = all the other phases
    binary_phase=zeros(Domain_size); % Initialization
    binary_phase(Phase_microstructure == label_phase_A) = 1;
    binary_phase(Phase_microstructure == label_phase_B) = 2;
    binary_phase(Phase_microstructure == label_phase_C) = 3;
    n_voxel_phases=sum(sum(sum(binary_phase~=0))); % Number of voxel of the three phases
    [Position_line.line(current_line),length_,LOC_tpbl.line(current_line).array] = Function_TPBL_direct_Algorithm(binary_phase, [1 2 3]); % TPBL
    TPBLdensity(current_line,4)=length_*voxel_size;
    TPBLdensity(current_line,6)=TPBLdensity(current_line,4)/TPBLdensity(current_line,5);
    TPBLdensity(current_line,7)=TPBLdensity(current_line,6)*corrective_factor;
 
    % % Correlation
    results_correlation(current_line).TPBL_directcounting = TPBLdensity(current_line,7);
    % % Visualization
    results_visualization(current_line).TPBL_directcounting = LOC_tpbl.line(current_line).array;

    % % Time
    timedata_domain(current_line+1,1) = line_name(current_line);
    timedata(current_line+1,1) = n_voxel_phases;
    timedata(current_line+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
    timedata(current_line+1,3) = toc(time_clock_start_phase); % Stopwatch elapsed time
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
Results_TPBL.Table_time = Table_time; % Save in main table result
% Result calculated on whole volume
Table_TPBL_direct = table(line_name,TPBLdensity(:,1),TPBLdensity(:,2),TPBLdensity(:,3),TPBLdensity(:,4),TPBLdensity(:,5),TPBLdensity(:,6),ones(number_line,1)*corrective_factor,TPBLdensity(:,7),...
    'VariableNames',{'Contact line' 'Label A' 'Label B' 'Label C' ['Length ' Length_unit] ['Volume ' Volume_unit] ['TPBL density ' TPBL_unit] 'Corrective factor' ['Corrected TPBL density ' TPBL_unit]});%
Results_TPBL.Table_TPBL_direct = Table_TPBL_direct; % Save in main table result
 

%% SAVE TABLES
if opt.save.xls
    filename = 'TPBL_direct'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='TPBL_over_Domainvolume';
    DATA_writetable.sheet(1).table=Table_TPBL_direct; 
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
disp('TPBL density = contact length / domain volume')
disp(Table_TPBL_direct)
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME
%%

%% ALONG DIRECTIONS
% README

% Triple Phase Boundary length density TPBLd =  sum(g(x)), with g(x)=n(x)*length/domain volume
% and n(x) the number of line at the position x along a given direction
% This discrete sum defintion will be converted in a continous (integral)
% sum definition: TPBLd = 1/D*integral(f(x),dx,0,D) = sum(g(x))
% To get f(x) [um-2], sum(g(x)) is derived.

% This is analogous with calculating a probability density function.
% (Although, the integral calculated here is not a probability but [um-2])

% Pros: obtain correct shape of the curve TPBLd=f(x).
% For a randomly distributed phase with no gradation, f(x)=cst=<TPBLd>

%% CALCULATION OF g(x), with g(x)=n(x)*line length/volume
for current_line=1:1:number_line
    for direction=1:1:number_dimension
        % Set position in real unit
        Position_line.line(current_line).direction(direction).value(:,1)=Position_line.line(current_line).direction(direction).value(:,1)*voxel_size;
        % Set line length in real unit
        Position_line.line(current_line).direction(direction).value(:,3)=Position_line.line(current_line).direction(direction).value(:,2)*voxel_size;
        %  Set TPBLd in real unit
        Position_line.line(current_line).direction(direction).value(:,3)=Position_line.line(current_line).direction(direction).value(:,3)/TPBLdensity(current_line,5);
        % Corrective factor
        Position_line.line(current_line).direction(direction).value(:,3)=Position_line.line(current_line).direction(direction).value(:,3)*corrective_factor;
    end
end

%% DISTRIBUTION
% Cumulative function sum(g(x))
for current_line=1:1:number_line
    for direction=1:1:number_dimension
        Cumulative_length.line(current_line).direction(direction).value=zeros(1,Domain_size(direction));
        for i=1:1:Domain_size(direction)
            for j=1:1:i
                Cumulative_length.line(current_line).direction(direction).value(1,j)=Cumulative_length.line(current_line).direction(direction).value(1,j)+Position_line.line(current_line).direction(direction).value(i,3);
            end
        end
    end
end
% Derivate function (i.e. 1/D*f(x))
for current_line=1:1:number_line
    for direction=1:1:number_dimension
        c_length = Domain_size(direction);
        for i=2:1:(c_length-1)
            Derivate_Cumulative_length.line(current_line).direction(direction).value(1,i)=-( Cumulative_length.line(current_line).direction(direction).value(1,i+1)-Cumulative_length.line(current_line).direction(direction).value(1,i-1) ) / (2*voxel_size);
        end
        Derivate_Cumulative_length.line(current_line).direction(direction).value(1,1)=-(Cumulative_length.line(current_line).direction(direction).value(1,2)-Cumulative_length.line(current_line).direction(direction).value(1,1))/(1*voxel_size);
        Derivate_Cumulative_length.line(current_line).direction(direction).value(1,c_length)=-(Cumulative_length.line(current_line).direction(direction).value(1,c_length)-Cumulative_length.line(current_line).direction(direction).value(1,c_length-1))/(1*voxel_size);
    end
end
% The specific surface area f(x)
for current_line=1:1:number_line
    for direction=1:1:number_dimension
        TPBL_density.line(current_line).direction(direction).value = zeros(1,Domain_size(direction));
        TPBL_density.line(current_line).direction(direction).value = Derivate_Cumulative_length.line(current_line).direction(direction).value*Domain_size(direction)*(voxel_size);
    end
end

%% VERIFICATION
% the relation Sp = 1/D*integral(f(x),dx,0,D) is checked
TPBL_integralvalue=zeros(number_line,7);
for current_line=1:1:number_line
    for direction=1:1:number_dimension
        X_= Position_line.line(current_line).direction(direction).value(:,1);
        Y_= TPBL_density.line(current_line).direction(direction).value(1,:);
        TPBL_integralvalue(current_line,1)=TPBLdensity(current_line,7); % Reference value
        TPBL_integralvalue(current_line,direction+1)=(1/(Domain_size(direction)*(voxel_size)))*trapz(X_,Y_); % Integral value
        % Relative error in percent
        TPBL_integralvalue(current_line,direction+1+3)=100*(TPBL_integralvalue(current_line,direction+1)-TPBL_integralvalue(current_line,1))/TPBL_integralvalue(current_line,1);
    end
end

%% MANAGING RESULTS
% Results are saved in a table

clear Variable_name_table;
Variable_name_table(1)={['Position ' voxel_unit]};
for current_line=1:1:number_line
    Variable_name_table(current_line+1)=line_name(current_line);
end

% g(x)
for direction=1:1:number_dimension
    array_g = zeros(Domain_size(direction),number_line+1);
    array_g(:,1)= Position_line.line(current_line).direction(direction).value(:,1);
    for current_line=1:1:number_line
        array_g(:,current_line+1)= Position_line.line(current_line).direction(direction).value(:,3);
    end
    Gx.direction(direction).table = array2table(array_g,'VariableNames',Variable_name_table);
end
% f(x)
for direction=1:1:number_dimension
    array_f = zeros(Domain_size(direction),number_line+1);
    array_f(:,1)= Position_line.line(current_line).direction(direction).value(:,1);
    for current_line=1:1:number_line
        array_f(:,current_line+1)= TPBL_density.line(current_line).direction(direction).value(1,:);
    end
    Fx.direction(direction).table = array2table(array_f,'VariableNames',Variable_name_table);
end
 
% Verification
Table_verification_integral = table(line_name,TPBL_integralvalue(:,1),TPBL_integralvalue(:,2),TPBL_integralvalue(:,3),TPBL_integralvalue(:,4),TPBL_integralvalue(:,5),TPBL_integralvalue(:,6),TPBL_integralvalue(:,7),...
    'VariableNames',{'Contact line' 'Reference' 'Integration direction 1' 'Integration direction 2' 'Integration direction 3' 'Relative error percent direction 1' 'Relative error percent direction 2' 'Relative error percent direction 3'});
fprintf('> Calculated along direction:\n\n');
fprintf('  Integral verification: relative error should be close to 0 percent:\n\n');
disp(Table_verification_integral)

Results_TPBL.Table_evolution_Gx = Gx; % Save in main table result
Results_TPBL.Table_evolution_Fx = Fx;
Results_TPBL.Table_verification_integral=Table_verification_integral;

%% SAVE TABLES
if opt.save.xls
    filename = 'TPBL_direct_along_directions'; % Filename without extension
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
    TPBLunit = '\mum^{-2}';
else
    axisunit = ['(' strunit ')'];
    TPBLunit = ['' voxel_unit ' ^{-2}'];
end

scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Triple phase boundary length density, direct method (integral definition)'; % Figure name
Fig.Color='white'; % Background colour
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*number_dimension/3 scrsz(4)*1/2]); % Full screen figure
for current_direction=1:1:number_dimension % Iterate over axe
    sub_axes=subplot(1,number_dimension,current_direction,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    h_title=title ({'Triple phase boundary length density, direct method','TPBL_{d} = 3-phase contact length / domain volume'}); % Set title font
    % Plot graphs
    for current_line=1:1:number_line % Loop over phases
        x_ = Position_line.line(current_line).direction(current_direction).value(:,1);
        y_ = TPBL_density.line(current_line).direction(current_direction).value(1,:);
        plot(x_,y_,'Color', colorline(current_line,:),'LineWidth',opt.format.linewidth,'DisplayName',char(line_name(current_line)));
    end
    % Axis label
    t_ = xlabel(' ');
    t_1 = sprintf('Position along %s ',char(infovol.directionname(current_direction)));
    t_2 = axisunit;
    t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
    t_ = ylabel(' '); 
    str1 = 'TPBL density distribution function f(x), such as';
    str2 = ['<TPBL_{d}> = 1/L \times \int_{0}^{L} f(x) dx (' TPBLunit ')'];
    t_.String= {str1,str2};
    % Legend
    for current_line=1:1:number_line
        str_legend(current_line).name = [char(line_name(current_line)) ', <TPBL_{p}>=' num2str(TPBLdensity(current_line,7),'%1.3f') TPBLunit];
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
sgtitle(Fig,'Triple phase boundary length density along directions','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
if opt.save.savefig % Save figure
    filename= 'TPBL_direct_along_direction';
    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
end
if opt.format.autoclosefig
    close(Fig); % Do not keep open figures
end

%% FRACTAL DIMENSION (COUNTING BOX)
if p.fractal_boxcounting.todo
    fprintf('> TPBL fractal dimension (box-counting method)\n');
    p.fractal_boxcounting.topology_dimension = 1;
    p.fractal_boxcounting.plot = false;
    for current_line=1:1:number_line % Loop over line and call box counting algorithm
        [N(:,current_line),box_lengths,fractal_dimension(:,current_line),fractal_dimension_convergence(:,:,current_line)] = Function_fractaldimension_boxcounting(LOC_tpbl.line(current_line).array,p.fractal_boxcounting);
    end

    % Table
    Table_fractaldimension = table(line_name(:,1),fractal_dimension(1,:)',fractal_dimension(2,:)',fractal_dimension(3,:)',fractal_dimension(4,:)',...
        'VariableNames',{'Contact line' 'Fit from 1 to' 'Fractal dimension' 'Topology dimension' 'Fractal propensity'});
    for current_line=1:1:number_line
        Table_boxlength(current_line).t = table(fractal_dimension_convergence(:,1,current_line),fractal_dimension_convergence(:,2,current_line),fractal_dimension_convergence(:,3,current_line),...
        'VariableNames',{'Fit from 1 to' 'Fractal dimension' 'Fit norm error'});
    end
    Results_Volumefraction.Table_fractaldimension = Table_fractaldimension; % Save in main table result
    disp(Table_fractaldimension)

    % Save table
    if opt.save.xls
        filename = 'TPBL_Fractaldimension_boxcounting'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name = 'Fractal dimension';
        DATA_writetable.sheet(1).table = Table_fractaldimension;
        for current_line=1:1:number_line
            DATA_writetable.sheet(1+current_line).name = char(line_name(current_line,1));
            DATA_writetable.sheet(1+current_line).table = Table_boxlength(current_line).t;
        end           
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end

    % Figure
    propertyname = 'TPBL surface'; % Figure name
    filename = 'TPBL_fractal_dimension_boxcounting';
    function_fractalfig(propertyname,filename,Current_folder, box_lengths,N,fractal_dimension_convergence,p.fractal_boxcounting.topology_dimension,number_line,p,opt,infovol_bis);

    % Correlation
    results_correlation(current_line).TPBL_fractaldimension_boxcounting = fractal_dimension(2,current_line);
    results_correlation(current_line).TPBL_fractalpropensity_boxcounting = abs(p.fractal_boxcounting.topology_dimension - fractal_dimension(2,current_line));
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
    property_voxelsizedependence = zeros(number_resize+1,number_line+1);
    property_voxelsizedependence(1,1)=voxel_size;
    property_voxelsizedependence(1,2:end)=TPBLdensity(:,7)';

    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fractaldimension_voxelsizedependence = zeros(number_resize+1,number_line+1);
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

        for current_line=1:1:number_line % Loop over all phases
            time_cpu_start_phase = cputime; % CPU start
            time_clock_start_phase = tic; % Stopwatch start

            % % Algorithm: SPECIFIC FOR EACH FILE
            % Label
            label_phase_A = p.triplets(current_line,1);
            label_phase_B = p.triplets(current_line,2);
            label_phase_C = p.triplets(current_line,3);
            % Create a binary microstructure : 1 = current phase A, 10 = current phase B, 0 = all the other phases
            binary_phase=zeros(Current_domain_size); % Initialization
            binary_phase(Phase_microstructure_resized == label_phase_A) = 1;
            binary_phase(Phase_microstructure_resized == label_phase_B) = 2;
            binary_phase(Phase_microstructure_resized == label_phase_C) = 3;
            n_voxel_phases=sum(sum(sum(binary_phase~=0))); % Number of voxel of the three phases
            [~,length_,LOC_tpbl_resized.line(current_line).array] = Function_TPBL_direct_Algorithm(binary_phase, [1 2 3]); % TPBL
            length_=length_*current_voxel_size;
            property_voxelsizedependence(current_iteration+1,current_line+1)=length_/Volume_*corrective_factor;

            % % Time
            timedata_perphase = [timedata_perphase; [n_voxel_phases (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

        end
        % CPU and stopwatch time - end
        timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];

        % % Fractal
        if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
            for current_line=1:1:number_line % Loop over line and call box counting algorithm
                [~,~,fractal_dimension_tmp, ~] = Function_fractaldimension_boxcounting(LOC_tpbl_resized.line(current_line).array,p.fractal_boxcounting);
                fractaldimension_voxelsizedependence(current_iteration+1,1)=size_choice(current_iteration)*voxel_size;
                fractaldimension_voxelsizedependence(current_iteration+1,current_line+1) = fractal_dimension_tmp(2,1);
            end      
        end

    end
    clear Phase_microstructure_resized;

    % Sort per voxel size
    property_voxelsizedependence = sortrows(property_voxelsizedependence,1);   
    
    %% EXTRAPOLATION TO 0 nm
    fprintf('> Triple phase boundary length density dependence with the voxel size (with corrective factor)\n');
    tmp = zeros(number_resize+2,number_line+1); % + 0 nm and + initial voxel size
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
    for current_line=1:1:number_line
        y=property_voxelsizedependence(:,current_line+1);
        pi = polyfit(x,y,interpolation_voxelsize_order);
        vq = polyval(pi,0);
        tmp(1,current_line+1)=vq;
        interpolation_voxelsize(current_line).pi=pi;
    end
    tmp(2:end,:) = property_voxelsizedependence;
    property_voxelsizedependence = tmp; clear tmp;

    if p.fractal_bertei.todo
        fractal_dimension_Bertei = zeros(number_line,3);
        % Fractal dimension according to Bertei et al., https://doi.org/10.1016/j.nanoen.2017.06.028 (Richardson/Mandelbrot formula)
        % Log(property) = m + (1-fractal dimension)*log(voxel size)
        logV = log(property_voxelsizedependence(2:end,1));
        logPs = zeros(length(logV),number_line);
        for current_line=1:1:number_line
            logP = log(property_voxelsizedependence(2:end,current_line+1));
            logPs(:,current_line)=logP;
            pf = polyfit(logV,logP,1);
            fractal_dimension_Bertei(current_line,1) = 1-pf(1);
            fractal_dimension_Bertei(current_line,2) = 1; % Topology dimension, line
            fractal_dimension_Bertei(current_line,3) = abs( fractal_dimension_Bertei(current_line,2) - fractal_dimension_Bertei(current_line,1) ); % Fractal propensity
            % Correlation
            results_correlation(current_line).TPBL_fractaldimension_Mandelbrot = fractal_dimension_Bertei(current_line,1);
            results_correlation(current_line).TPBL_fractalpropensity_Mandelbrot = fractal_dimension_Bertei(current_line,3);
        end
    end
    
    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
    for current_line=1:1:number_line
        Variable_name_table(1+current_line)=line_name(current_line);
    end
    % Table
    Table_TPBLdirect_voxelsizedependence = array2table(property_voxelsizedependence,...
        'VariableNames',Variable_name_table);

    if p.fractal_bertei.todo
        Variable_name_table={'Voxel size log'}; % Columns name
        for current_line=1:1:number_line
            Variable_name_table(1+current_line)={[char(line_name(current_line)) ' log']};
        end
        Table_TPBLdirect_voxelsizedependence_loglog = array2table([logV logPs],'VariableNames',Variable_name_table);
        Table_Fractaldimension_Bertei = table(line_name(:,1),fractal_dimension_Bertei(:,1),fractal_dimension_Bertei(:,2),fractal_dimension_Bertei(:,3),...
            'VariableNames',{'Contact line','Fractal dimension','Topology dimension','Fractal propensity'});   
    end

    %% DISPLAY TEXT RESULTS
    disp(Table_TPBLdirect_voxelsizedependence)
    if p.fractal_bertei.todo
        disp(Table_TPBLdirect_voxelsizedependence_loglog)
        fprintf('Richardson/Mandelbrot formula: Log(property) = m + (1-fractal dimension)*log(voxel size))\n');
        disp(Table_Fractaldimension_Bertei)
    end
    
    %% SAVE RESULTS
    Results_TPBL.voxelsizedependence = Table_TPBLdirect_voxelsizedependence; % Save in main table result
    if opt.save.xls
        filename = 'TPBL_direct_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Triple phase boundary length density';
        DATA_writetable.sheet(1).table=Table_TPBLdirect_voxelsizedependence;
        if p.fractal_bertei.todo
            DATA_writetable.sheet(2).name='Log log';
            DATA_writetable.sheet(2).table=Table_TPBLdirect_voxelsizedependence_loglog;
            DATA_writetable.sheet(3).name='Fractal dimension Mandelbrot';
            DATA_writetable.sheet(3).table=Table_Fractaldimension_Bertei;
        end
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end

    %% FIGURES
    parameters_figure.plotlog = true;
    parameters_figure.propertyname = 'Triple phase boundary length density';
    parameters_figure.method = 'direct counting';
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence;
    parameters_figure.number_phase = number_line;
    parameters_figure.str_ylabel = ['TPBL_{d} = 3-phase contact length / domain volume' '(' TPBLunit ')'];
    parameters_figure.propertynameunit = TPBLunit;
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'TPBL_direct_voxel_size_dependence';
    parameters_figure.infovol = infovol_bis;
    parameters_figure.opt = opt;
    parameters_figure.todo = ones(number_line,1);
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    

    %% FRACTAL DIMENSION
    if p.fractal_boxcounting.todo && p.fractal_boxcounting.voxelsize
        fprintf('> TPBL fractal dimension (box counting) dependence with the voxel size\n');
        fprintf('  Extrapolation to zero voxel size: polynomial of order 1\n\n');

        % Sort by voxel size
        fractaldimension_voxelsizedependence = sortrows(fractaldimension_voxelsizedependence,1);

        % Extrapolation to 0 voxel size
        tmp = zeros(number_resize+2,number_line+1); % + 0 nm and + initial voxel size
        for current_line=1:1:number_line
            y=fractaldimension_voxelsizedependence(:,current_line+1);
            pi = polyfit(x,y,1);
            vq = polyval(pi,0);
            tmp(1,current_line+1)=vq;
            interpolation_voxelsize(current_line).pi=pi;                
        end
        tmp(2:end,:) = fractaldimension_voxelsizedependence;
        fractaldimension_voxelsizedependence = tmp; clear tmp;
    
        % Managing result
        Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
        for current_line=1:1:number_line
            Variable_name_table(1+current_line)=line_name(current_line);
        end
        % Table
        Table_Fractaldimension_voxelsizedependence = array2table(fractaldimension_voxelsizedependence,...
            'VariableNames',Variable_name_table);
        fprintf('    fitted from s=1 to s=2\n');
        disp(Table_Fractaldimension_voxelsizedependence)

        % Save result
        Results_Volumefraction.Fractaldimension_voxelsizedependence = Table_Fractaldimension_voxelsizedependence; % Save in main table result
        if opt.save.xls
            filename = 'TPBL_Fractaldimension_boxcounting_voxel_size_dependence'; % Filename without extension
            % Prepare the data
            clear DATA_writetable
            DATA_writetable.sheet(1).name='Fit from 1 to 2';
            DATA_writetable.sheet(1).table=Table_Fractaldimension_voxelsizedependence;
            % Save function
            Function_Writetable(Current_folder,filename,DATA_writetable)
        end

        % Correlation
        for current_line=1:1:number_line
            results_correlation(current_line).TPBL_fractaldimension_boxcounting_extrapolated = fractaldimension_voxelsizedependence(1,current_line+1) ;
            results_correlation(current_line).TPBL_fractalpropensity_boxcounting_extrapolated = abs(p.fractal_boxcounting.topology_dimension - fractaldimension_voxelsizedependence(1,current_line+1));
        end

        % Figure
        parameters_figure.plotlog = false; 
        parameters_figure.figname = 'TPBL fractal dimension';
        parameters_figure.propertyname = 'TPBL fractal dimension';
        parameters_figure.method = 'Box counting';
        parameters_figure.property_voxelsizedependence = fractaldimension_voxelsizedependence;
        parameters_figure.number_phase = number_line;
        parameters_figure.str_ylabel = 'Fractal dimension';
        parameters_figure.propertynameunit = [];
        parameters_figure.interpolation_voxelsize = interpolation_voxelsize;
        parameters_figure.Current_folder = Current_folder;
        parameters_figure.filename = 'TPBL_Fractaldimension_voxel_size_dependence';
        parameters_figure.infovol = infovol_bis;
        parameters_figure.opt = opt;
        parameters_figure.todo = p.todo;
        Function_create_figure_voxelsizedependence(parameters_figure) % Figures
    end

end


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%

propertyname='Triple phase boundary length density';
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
            Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_line,3,3); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
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
            Property_eachsubdomain = zeros(number_subdomain,number_line+4);
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

                for current_line=1:1:number_line % Loop over all phases
                    time_cpu_start_phase = cputime; % CPU start
                    time_clock_start_phase = tic; % Stopwatch start

                    % % Algorithm: SPECIFIC FOR EACH FILE
                    % Label
                    label_phase_A = p.triplets(current_line,1);
                    label_phase_B = p.triplets(current_line,2);
                    label_phase_C = p.triplets(current_line,3);
                    % Create a binary microstructure : 1 = current phase A, 10 = current phase B, 0 = all the other phases
                    binary_phase=zeros(Current_domain_size); % Initialization
                    binary_phase(current_subdomain == label_phase_A) = 1;
                    binary_phase(current_subdomain == label_phase_B) = 2;
                    binary_phase(current_subdomain == label_phase_C) = 3;
                    n_voxel_phases=sum(sum(sum(binary_phase~=0))); % Number of voxel of the three phases
                    [~,length_,~] = Function_TPBL_direct_Algorithm(binary_phase, [1 2 3]); % TPBL
                    length_=length_*voxel_size;
                    Property_eachsubdomain(subdomain_id,current_line+4)=length_/Volume_*corrective_factor;

                    % % Time
                    timedata_perphase = [timedata_perphase; [n_voxel_phases (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

                end
                % CPU and stopwatch time - end
                timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
            end

            %% STATISTICAL ANALYSIS and RVE SIZE
            [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis(number_group_size,number_line,GROUP_SUBDOMAIN,Property_eachsubdomain,voxel_size,RVEparameters);
            
            %% SAVE FOR CORRELATION
            if k_nestedRVE == 0
                for k_threshold=1:1:n_threshold
                    for current_line=1:1:number_line
                        if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                            if Size_RVE(1,current_line,2,1)~=0
                                str_ = ['TPBLdirect_RVE_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                str_(str_=='.')='p';
                                results_correlation(current_line).(str_) = Size_RVE(1,current_line,2,1);
                            end
                            if strcmp(RVEparameters.type,'C')
                                if Size_RVE(1,current_line,2,2)~=0
                                    str_ = ['TPBLdirect_RVE_squarerootFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_line).(str_) = Size_RVE(1,current_line,2,2);
                                end
                            elseif strcmp(RVEparameters.type,'D')
                                if Size_RVE(1,current_line,2,2)~=0
                                    str_ = ['TPBLdirect_RVE_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_line).(str_) = Size_RVE(1,current_line,2,2);
                                end
                            end
                        else
                            if Size_convergence(1,current_line,2,1)~=0
                                str_ = ['TPBLdirect_conv_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                str_(str_=='.')='p';
                                results_correlation(current_line).(str_) = Size_convergence(1,current_line,2,1);
                            end
                            if strcmp(RVEparameters.type,'G')
                                if Size_convergence(1,current_line,2,2)~=0
                                    str_ = ['TPBLdirect_conv_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_line).(str_) = Size_convergence(1,current_line,2,2);
                                end
                            elseif strcmp(RVEparameters.type,'H')
                                if Size_convergence(1,current_line,2,2)~=0
                                    str_ = ['TPBLdirect_conv_areaFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                    str_(str_=='.')='p';
                                    results_correlation(current_line).(str_) = Size_convergence(1,current_line,2,2);
                                end
                            end
                        end
                    end
                end
            end            

            %% MANAGING RESULTS
            p.todo = ones(number_line,1);
            [RVE] = Function_subdomains_manage_results(Property_eachsubdomain, Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,RVE,k_RVE,number_line,number_line,infovol_bis,p);

            if RVEparameters.donested
                for k_threshold=1:1:n_threshold
                    for current_line=1:1:number_line
                        if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,1,1) = Size_RVE(k_threshold,current_line,1,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,2,1) = Size_RVE(k_threshold,current_line,2,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,3,1) = Size_RVE(k_threshold,current_line,3,1);
                        else
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,1,1) = Size_convergence(k_threshold,current_line,1,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,2,1) = Size_convergence(k_threshold,current_line,2,1);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,3,1) = Size_convergence(k_threshold,current_line,3,1);
                        end
                        if strcmp(RVEparameters.type,'C')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,1,2) = Size_RVE(k_threshold,current_line,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,2,2) = Size_RVE(k_threshold,current_line,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,3,2) = Size_RVE(k_threshold,current_line,3,2);
                        elseif strcmp(RVEparameters.type,'D')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,1,3) = Size_RVE(k_threshold,current_line,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,2,3) = Size_RVE(k_threshold,current_line,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,3,3) = Size_RVE(k_threshold,current_line,3,2);
                        elseif strcmp(RVEparameters.type,'G')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,1,3) = Size_convergence(k_threshold,current_line,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,2,3) = Size_convergence(k_threshold,current_line,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,3,3) = Size_convergence(k_threshold,current_line,3,2);
                        elseif strcmp(RVEparameters.type,'H')
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,1,2) = Size_convergence(k_threshold,current_line,1,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,2,2) = Size_convergence(k_threshold,current_line,2,2);
                            Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_line,3,2) = Size_convergence(k_threshold,current_line,3,2);
                        end
                    end
                end
            end

            %% TEXT DISPLAY AND SAVE RESULTS
            if k_nestedRVE == 0
                RVEparameters.disp_parRVE = true;
                Function_subdomains_display_and_save(RVE,k_RVE,RVEparameters,number_line,number_line,propertyname,Sub_folder_RVE,opt,infovol_bis,p);
            end

            %% FIGURES
            if k_nestedRVE == 0
                parameters_figure.propertyname = propertyname;
                parameters_figure.propertynameunit = TPBLunit;
                parameters_figure.RVE = RVEparameters;
                parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                parameters_figure.savefolder = Sub_folder_RVE;
                parameters_figure.number_phase = number_line;
                parameters_figure.number_phase_todo = number_line;
                parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
                parameters_figure.Property_eachsubdomain = Property_eachsubdomain;
                parameters_figure.derivative_convergence = derivative_convergence;
                parameters_figure.relativedifference_convergence = relativedifference_convergence;
                parameters_figure.Size_RVE = Size_RVE;
                parameters_figure.convergence_criterion = RVEparameters.threshold_reldiff;
                parameters_figure.Size_convergence = Size_convergence;
                parameters_figure.Wholevolume_size = Wholevolume_size;
                parameters_figure.Wholevolume_results = TPBLdensity(:,7);
                parameters_figure.infovol = infovol_bis;
                parameters_figure.todo = ones(number_line,1);
                parameters_figure.opt = opt;
                Function_create_figures_RVE(parameters_figure) % Figures
            end
        end

        %% NESTED ANALYSIS RESULT
        if RVEparameters.donested
            % Table
            [RVE] = Function_nestedtable(RVE,k_RVE,RVEparameters,number_line,propertyname,Result_nestedRVE,Sub_folder_RVE,opt,infovol_bis,p);
            % Figure
            parameters_figure.Result_nestedRVE = Result_nestedRVE;
            Function_create_figures_nestedRVE(parameters_figure) % Figures
            % Save
            RVE(k_RVE).nestedanalysis = Result_nestedRVE;
        end

        Results_TPBL.RVE.Sp = RVE; % Save in main table result
    end
end

%%
%% ENDING FUNCTION
%%

%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_TPBL.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_TPBL.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_TPBL.Table_date = Table_date;

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
    Function_Writetable(Current_folder,'TPBL_direct_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'Int_direct_calculation_time', 'Triple phase boundary length (direct)', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_TPBL_length_area'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_TPBL_direct'],'Results_TPBL')
end

%% SAVE VISUALIZATION
Current_folder = [infovol.volpath 'Visualization' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Visualization_TPBL_direct'],'results_visualization','-v7.3');    

end