function [] = Function_particle_size_CHORD(Phase_microstructure,infovol,opt,p,foo1, foo2)
% Calculate Particle size fitting the cumulative function of the chord lengths
% Function_particle_size_CHORD(Phase_microstructure, infovol, opt, p) - when used with the MATBOX toolbox
% or
% Function_particle_size_CHORD(Phase_microstructure, labels, directions, voxelsize, unit, size_filter) - when used as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       labels: either 'All' (all phased are characterized) or a 1D array listing the labels to characterize
%       directions, a 1D array: which directions to investigate. Set equal to 0 for volume analysis only.
%       voxelsize, a scalar: the voxel length
%       unit, a string: the unit name of the voxel length
%       size_filter, a scalar: filter applied on chord length, in voxel length. Set to 0 to ignore.
%       e.g.: Function_particle_size_CHORD(<your_3d_array>, 0, [1,2,3], 0.4, 'um', 0);

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 6 % Case for function called as: Function_particle_size_CHORD(Phase_microstructure, labels, directions, voxelsize, unit). Standalone use.
        labels = infovol; clear infovol;
        directions = opt; clear opt;
        voxelsize = p; clear p;
        unit = foo1;  clear foo1;
        size_filter = foo2; clear foo2;
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['DiameterChord_' char(t)];
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
            if isnumeric(labels)
                if ismember(infovol.phaselabel(k),labels)
                    p.todo(k)=1;
                else
                    p.todo(k)=0;
                end
            elseif strcmp(labels,'All')
                p.todo(k)=1;
            else
                disp 'Error calling Function_particle_size_CHORD. Incorrect labels argument.'
                help Function_particle_size_CHORD
                return
            end
        end
        if ~sum(p.todo)
            disp 'Error calling Function_particle_size_CHORD. Labels argument correspond to no existing label.'
            help Function_particle_size_CHORD
            return
        end

        % Directions investigated
        for k=1:1:3
            p.direction_todo(k)=0;
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

        % Set default colormap
        p.colormap = 'jet';
        p.colormap_background = 'white';        
        
        % No Voxel size dependence analysis
        p.scaling = 1;
        % No Representative Volume Element analysis
        p.RVE.number_RVE = 0;

    else % Incorrect number of argument
        disp 'Error calling Function_particle_size_CHORD. Wrong number of argument.'
        help Function_particle_size_CHORD
        return
    end
    
else
    % Read parameters
    size_filter = p.size_filter;
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    separator = '\';
else
    separator = '/';
end
Current_folder = [infovol.volpath 'Diameter_Chord' separator];
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% VOLUME INFORMATION
Domain_size = size(Phase_microstructure);
if length(Domain_size)==2 % Add third column with unit value
    Domain_size = [Domain_size 1];
    number_dimension = 2;
else
    number_dimension =3; % 3D case
end
number_phase = length(infovol.phaselabel); % Number of phase
voxel_number = prod(Domain_size); % Number of voxel
voxel_size = infovol.voxelsize;
voxel_unit = infovol.unit;

%% INITIALIZE RESULTS (USE FOR CORRELATION and VISUALIZATION)
current_phase_todo = 0;
for current_phase=1:1:number_phase
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        results_correlation(current_phase_todo).name = infovol.phasename(current_phase,1);
        results_visualization(current_phase_todo).name = infovol.phasename(current_phase,1);
    end
end

%% PARAMETERS
% Note: those are not algorithm paramter (c-PSD is parameter free), but post-processing smoothing parameter to get a better looking distribution function

% Cumulative and distribution functions parameters
density_fct_parameters.round_value = 3;
density_fct_parameters.smooth_cumulative_fct = true;

%%
%% ALGORITHM ON WHOLE VOLUME
%%

disp '    PARTICLE SIZE - CHORD LENGTH MAP FITTING (CLMF) METHOD';
disp '    ------------------------------------------------------';
disp ' ';

%% CALCULATION
time_cpu_start_volume = cputime; % CPU start
time_clock_start_volume = tic; % Stopwatch start

% Initialization (generic)
number_phase_todo = sum(p.todo); % How many phase are we going to analyse ?
number_direction_todo = sum(p.direction_todo);
Numbervoxel_phase=zeros(number_phase_todo,1); 
timedata = zeros(number_phase_todo+1,3); timedata(1,1) = voxel_number;
timedata_domain = cell(number_phase_todo+1,1);
phasename_todo = cell(number_phase_todo,1);
phaselabel = zeros(number_phase_todo,1);
% Initialization (algorithm-specific)
fitted_diameters = zeros(number_phase_todo,3); % mean

current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        current_phase_todo=current_phase_todo+1;
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);
        Numbervoxel_phase(current_phase_todo,1)= sum(sum(sum(Phase_microstructure==phaselabel(current_phase_todo,1) )));

        % % Algorithm
        phaselabel(current_phase_todo,1) = infovol.phaselabel(current_phase);
        % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
        binary_phase=zeros(Domain_size); % Initialization
        binary_phase(Phase_microstructure == phaselabel(current_phase_todo,1)) = 1; % Binary phase
        [fitted_diameters(current_phase_todo,:), Res(current_phase_todo).Rs, Res(current_phase_todo).difference, Res(current_phase_todo).CDF, ChordLength_alongaxe1,ChordLength_alongaxe2,ChordLength_alongaxe3] = Function_particle_size_ChordLength_volume_Algorithm(binary_phase,size_filter);
        % Equivalent mean diameter
        fitted_diameters(current_phase_todo,:) = fitted_diameters(current_phase_todo,:)*voxel_size;

        % % correlation
        results_correlation(current_phase_todo).Particle_diameter_alongaxe1_CHORD = fitted_diameters(current_phase_todo,1);
        results_correlation(current_phase_todo).Particle_diameter_alongaxe2_CHORD = fitted_diameters(current_phase_todo,2);
        results_correlation(current_phase_todo).Particle_diameter_alongaxe3_CHORD = fitted_diameters(current_phase_todo,3);

        % % Visualization
        results_visualization(current_phase_todo).ChordLength_alongaxe1 = ChordLength_alongaxe1;
        results_visualization(current_phase_todo).ChordLength_alongaxe2 = ChordLength_alongaxe2;
        results_visualization(current_phase_todo).ChordLength_alongaxe3 = ChordLength_alongaxe3;
        
        % % Time
        timedata_domain(current_phase_todo+1,1) = infovol.phasename(current_phase,1);
        timedata(current_phase_todo+1,1) = Numbervoxel_phase(current_phase_todo,1);
        timedata(current_phase_todo+1,2) = cputime-time_cpu_start_phase; % CPU elapsed time
        timedata(current_phase_todo+1,3) = toc(time_clock_start_phase); % CPU elapsed time  

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
Results_chord.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
Table_chord_d50 = table(phasename_todo,fitted_diameters(:,1),fitted_diameters(:,2),fitted_diameters(:,3),...
    'VariableNames',{'Phase' [char(infovol.directionname(1)) ' ' voxel_unit] [char(infovol.directionname(2)) ' ' voxel_unit] [char(infovol.directionname(3)) ' ' voxel_unit]});
for current_phase_todo=1:1:number_phase_todo
    for current_direction=1:1:3
        Chord_cumulative.phase(current_phase_todo).direction(current_direction).table = table(Res(current_phase_todo).CDF(current_direction).x*voxel_size, Res(current_phase_todo).CDF(current_direction).numerical, Res(current_phase_todo).CDF(current_direction).analytical,...
            'VariableNames',{['Chord length ' voxel_unit],'Numerical cumulative function','Analytical cumulative function',});
    end
end
Results_chord.Table_chord_d50 = Table_chord_d50;
Results_chord.Table_chord_cumulative = Chord_cumulative;

%% SAVE TABLES
if opt.save.xls
    filename = 'CHORD'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Fitted_diameter';
    DATA_writetable.sheet(1).table=Table_chord_d50;
    sheet_=1;    
    % Chord cumulative distribution
    for current_phase_todo = 1:1:number_phase_todo
        for current_direction=1:1:3
            sheet_=sheet_+1;
            DATA_writetable.sheet(sheet_).name=[char(phasename_todo(current_phase_todo,1)) ' along ' num2str(current_direction)];
            DATA_writetable.sheet(sheet_).table=Chord_cumulative.phase(current_phase_todo).direction(current_direction).table;
        end
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
disp 'Particle diameter';
disp(Table_chord_d50)
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%
scrsz = get(0,'ScreenSize'); % Screen resolution

strunit = voxel_unit;
if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
    axisunit = '(\mum)';
    Dunit = '\mum';
else
    axisunit = ['(' strunit ')'];
    Dunit = voxel_unit;
end

%% CUMULATIVE FUNCTION DIAMETER FITTING
col_ = colororder;
current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;

        Fig= figure;
        Fig.Name= sprintf('Chord length: fitted diameters of phase %s',char(phasename_todo(current_phase_todo,1)));
        Fig.Color='white'; % Background colour
        set(Fig, 'Position', [scrsz(1) scrsz(2) scrsz(3)*4/5 scrsz(4)*1/2]);

        sub_axes=subplot(1,2,1,'Parent',Fig); % Create axes
        hold(sub_axes,'on'); % Active subplot
        h_title=title (' ','FontName',opt.format.fontname);
        h_title.String= {'Analytical - numerical integral difference of','the chord length cumulative fct'};
        for k=1:1:3
            plot(2*Res(current_phase_todo).Rs ,Res(current_phase_todo).difference(:,k)','LineWidth',opt.format.linewidth, 'DisplayName',char(infovol.directionname(k)));
        end
        xlabel(['Ellipsoid diameters ' axisunit]); % Axis label
        ylabel('Difference');
        grid(sub_axes,opt.format.grid); % Display grid
        set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
        set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
        h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
        h_legend = legend(sub_axes,'Location','best');
        h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
        hold(sub_axes,'off'); % - Figure has been done

        sub_axes=subplot(1,2,2,'Parent',Fig); % Create axes
        hold(sub_axes,'on'); % Active subplot
        h_title=title (' ','FontName',opt.format.fontname);
        h_title.String= {'Chord length cumulative fct','calculated for the phase and a unique sphere with a fitted diameter'};
        for current_direction=1:1:3
            h_phase = plot(Res(current_phase_todo).CDF(current_direction).x*voxel_size, Res(current_phase_todo).CDF(current_direction).numerical);
            h_sphere=plot(Res(current_phase_todo).CDF(current_direction).x*voxel_size,Res(current_phase_todo).CDF(current_direction).analytical);
            set(h_phase,'LineStyle','-','Color',col_(current_direction,:),'LineWidth',opt.format.linewidth,'LineStyle','-','DisplayName',[char(infovol.directionname(current_direction)) ', ' char(phasename_todo(current_phase_todo,1))]);
            set(h_sphere,'LineStyle','--','Color',col_(current_direction,:),'LineWidth',opt.format.linewidth,'LineStyle','--','DisplayName',['Sphere with fitted diameter = ' num2str(fitted_diameters(current_phase_todo,current_direction),'%4.2f') ' ' Dunit] );
        end
        xlabel(['Chord length ' axisunit]); % - Axis label
        ylabel('Cumulative function');
        grid(sub_axes,opt.format.grid); % Display grid
        set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
        set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % - Fontname and fontsize
        h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
        h_legend = legend(sub_axes,'Location','best');
        h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
        hold(sub_axes,'off'); % - Figure has been done

        sgtitle(Fig,['Chord length cumulative function diameter fitting, ' char(phasename_todo(current_phase_todo,1))] ,'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
        if opt.save.savefig % Save figure
            filename= sprintf('Chordlength_fitted_%s',char(phasename_todo(current_phase_todo,1)));
            function_savefig(Fig, Current_folder, filename, opt.save); % Call function
        end
        if opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end
end


%% PARTICLE DISTANCE TO SURFACE MAP
% Color map
str = ['myColorMap = ' p.colormap '(256);'];
eval(str);
if strcmp(p.colormap_background,'white')
    myColorMap(1,:) = 1; 
elseif strcmp(p.colormap_background,'black')
    myColorMap(1,:) = 0; 
elseif strcmp(p.colormap_background,'grey')
    myColorMap(1,:) = [0.75 0.75 0.75];      
end

for current_phase_todo=1:1:number_phase_todo % Loop over phases
    Fig = figure; % Create figure
    Fig.Name= ['Chord length, phase ' char(phasename_todo(current_phase_todo))]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]); % Full screen figure
    idaxe = 0;
    for view_normal_to=1:1:number_dimension
        for chorddirection = 1:1:number_dimension
            if chorddirection==1
                L = results_visualization(current_phase_todo).ChordLength_alongaxe1;
            elseif chorddirection==2
                L = results_visualization(current_phase_todo).ChordLength_alongaxe2;
            else
                L = results_visualization(current_phase_todo).ChordLength_alongaxe3;
            end
            idaxe = idaxe+1;
            sub_axes=subplot(number_dimension,number_dimension,idaxe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            h_title=title ({['Chord length along ' char(infovol.directionname(chorddirection))],['View normal to ' char(infovol.directionname(view_normal_to))]}); % Set title font
            % - Plot graphs
            if view_normal_to==1
                tmp=squeeze(L(round(Domain_size(1)/2),:,:));
                h=image(tmp,'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',char(infovol.directionname(3)));
                t_1y = sprintf('Position along %s ',char(infovol.directionname(2)));
                set(h, 'XData', [0, Domain_size(3)*voxel_size]);
                set(h, 'YData', [0, Domain_size(2)*voxel_size]);
            elseif view_normal_to==2
                tmp=squeeze(L(:,round(Domain_size(2)/2),:));
                h=image(tmp,'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',char(infovol.directionname(3)));
                t_1y = sprintf('Position along %s ',char(infovol.directionname(1)));
                set(h, 'XData', [0, Domain_size(3)*voxel_size]);
                set(h, 'YData', [0, Domain_size(1)*voxel_size]);
            elseif view_normal_to==3
                h=image(L(:,:,round(Domain_size(3)/2)),'CDataMapping','scaled');
                t_1x = sprintf('Position along %s ',char(infovol.directionname(1)));
                t_1y = sprintf('Position along %s ',char(infovol.directionname(2)));
                set(h, 'XData', [0, Domain_size(1)*voxel_size]);
                set(h, 'YData', [0, Domain_size(2)*voxel_size]);
            end
            axis equal; axis tight;
            % - Axis label
            t_ = xlabel(' ');
            t_2 = axisunit;
            t_.String= [t_1x t_2]; % Sprintf does not accept greek characters
            t_ = ylabel(' ');
            t_2 = axisunit;
            t_.String= [t_1y t_2]; % Sprintf does not accept greek characters
            set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize-4); % Fontname and fontsize
            colormap(myColorMap);
            % Create colorbar
            %colorbar('peer',sub_axes);
            h=colorbar(sub_axes);
            ylabel(h, ['Chord length ' axisunit]);
            set(h,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize);
            h_title.FontSize = opt.format.titlefontsize-2; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
    end
    sgtitle(Fig,['Chord length, phase, ' char(phasename_todo(current_phase_todo))],'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
    if opt.save.savefig % Save figure
        filename= sprintf('Chord_map_%s', char(phasename_todo(current_phase_todo)));
        function_savefig(Fig, Current_folder, filename, opt.save); % Call function
    end
    if opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end
end

%% SECTION ALONG POSITION
if sum(p.direction_todo)
    fprintf('Diameters fitted slice per slice, calculation is on-going (slow)...\n\n');

    % Calculation
    for current_phase_todo=1:1:number_phase_todo % Loop over all phases
        binary_phase=zeros(Domain_size);
        binary_phase(Phase_microstructure == phaselabel(current_phase_todo,1)) = 1; % Binary phase
        [fitted2D(current_phase_todo).dir] = Function_particle_size_ChordLength_section_Algorithm(binary_phase,p.direction_todo,size_filter);
    end

    % Plot
    for current_phase_todo=1:1:number_phase_todo % Loop over all phases
        Fig = figure; % Create figure
        Fig.Name= ['Fitted section diameter along directions,' char(phasename_todo(current_phase_todo,1))]; % Figure name
        Fig.Color='white'; % Background colour
        set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]); % Full screen figure
        current_direction_todo = 0;
        for current_direction=1:1:3
            if p.direction_todo(current_direction)
                current_direction_todo = current_direction_todo+1;
                sub_axes=subplot(1,number_direction_todo,current_direction_todo,'Parent',Fig);
                hold(sub_axes,'on'); % Active subplot
                h_title=title(['Section diameters along ',char(infovol.directionname(current_direction))]);
                tmp = fitted2D(current_phase_todo).dir(current_direction);
                h1 = plot(tmp.x*voxel_size, tmp.diameter_ellipse(:,1)*voxel_size,'LineWidth',opt.format.linewidth);
                h2 = plot(tmp.x*voxel_size, tmp.diameter_ellipse(:,2)*voxel_size,'LineWidth',opt.format.linewidth);
                streq = sprintf('D_{eq,%s}',char(infovol.directionname(current_direction)));
                heq = plot(tmp.x*voxel_size, tmp.diameter_eqdisc(:,1)*voxel_size,'LineWidth',opt.format.linewidth,'DisplayName',streq);
                if current_direction==1
                    str1 = sprintf('D_{%s,%s}',char(infovol.directionname(2)),char(infovol.directionname(current_direction)));
                    str2 = sprintf('D_{%s,%s}',char(infovol.directionname(3)),char(infovol.directionname(current_direction)));
                elseif current_direction==2
                    str1 = sprintf('D_{%s,%s}',char(infovol.directionname(1)),char(infovol.directionname(current_direction)));
                    str2 = sprintf('D_{%s,%s}',char(infovol.directionname(3)),char(infovol.directionname(current_direction)));
                else
                    str1 = sprintf('D_{%s,%s}',char(infovol.directionname(1)),char(infovol.directionname(current_direction)));
                    str2 = sprintf('D_{%s,%s}',char(infovol.directionname(2)),char(infovol.directionname(current_direction)));
                end
                set(h1,'DisplayName',str1)
                set(h2,'DisplayName',str2)
                % Axis label
                t_ = xlabel(' ');
                t_1 = sprintf('Position along %s ',char(infovol.directionname(current_direction)));
                t_2 = axisunit;
                t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
                ylabel(['Fitted diameters ' axisunit]);

                h_legend = legend(sub_axes,'Location','best');
                h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize

                % - Grid
                grid(sub_axes,opt.format.grid); % Display grid
                set(sub_axes,'XMinorGrid',opt.format.minorgrid,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
                set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
                h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
                h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
                hold(sub_axes,'off'); % Relase figure
            end
        end
        sgtitle(Fig,{'Section diameters along directions', char(phasename_todo(current_phase_todo))},'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
        if opt.save.savefig % Save figure
            filename= ['Diameters_chord_along_direction_' char(phasename_todo(current_phase_todo))];
            function_savefig(Fig, Current_folder, filename, opt.save); % Call function
        end
        if opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end

    % Table
    for current_phase_todo=1:1:number_phase_todo
        for current_direction=1:1:3
            if p.direction_todo(current_direction)
                tmp = fitted2D(current_phase_todo).dir(current_direction);
                if current_direction==1
                    str1 = ['D_{' char(infovol.directionname(2)) ',' char(infovol.directionname(current_direction))    '}'];
                    str2 = ['D_{' char(infovol.directionname(3)) ',' char(infovol.directionname(current_direction))    '}'];
                elseif current_direction==2
                    str1 = ['D_{' char(infovol.directionname(1)) ',' char(infovol.directionname(current_direction))    '}'];
                    str2 = ['D_{' char(infovol.directionname(3)) ',' char(infovol.directionname(current_direction))    '}'];
                else
                    str1 = ['D_{' char(infovol.directionname(1)) ',' char(infovol.directionname(current_direction))    '}'];
                    str2 = ['D_{' char(infovol.directionname(2)) ',' char(infovol.directionname(current_direction))    '}'];
                end
                streq = ['D_{eq,' char(infovol.directionname(current_direction)) '}',];
                Chord_evolution.phase(current_phase_todo).direction(current_direction).table = table((tmp.x*voxel_size)', tmp.diameter_ellipse(:,1)*voxel_size, tmp.diameter_ellipse(:,2)*voxel_size, tmp.diameter_eqdisc(:,1)*voxel_size,...
                    'VariableNames',{['Position along ' char(infovol.directionname(current_direction)) ' ' voxel_unit],[str1 ' ' voxel_unit],[str2 ' ' voxel_unit],[streq ' ' voxel_unit]});
            end
        end
    end
    Results_chord.Alongdirections = Chord_evolution; % Save in main table result

    if opt.save.xls
        filename = 'Diameters_chord_along_directions'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        sheet_=0;
        for current_phase_todo=1:1:number_phase_todo
            for current_direction=1:1:3
                if p.direction_todo(current_direction)
                    sheet_=sheet_+1;
                    DATA_writetable.sheet(sheet_).name=[char(phasename_todo(current_phase_todo,1)) ' along ' char(infovol.directionname(current_direction))];
                    DATA_writetable.sheet(sheet_).table=Chord_evolution.phase(current_phase_todo).direction(current_direction).table;
                end
            end
        end
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end

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
    property_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1,number_dimension);
    property_voxelsizedependence(1,1,:)=voxel_size;
    for current_direction_todo=1:1:number_dimension
        property_voxelsizedependence(1,2:end,current_direction_todo)=fitted_diameters(:,current_direction_todo)';
    end

    % Loop on each voxel size
    for current_iteration=1:1:number_resize

        % % Microstructure scaling
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        property_voxelsizedependence(current_iteration+1,1,:)=current_voxel_size;
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
                % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
                binary_phase=zeros(Current_domain_size); % Initialization
                code_tmp = infovol.phaselabel(current_phase);
                Numbervoxel_phase_tmp= sum(sum(sum(Phase_microstructure_resized==code_tmp )));
                binary_phase(Phase_microstructure_resized == code_tmp) = 1; % Binary phase
                [fitted_diameters_resized(current_phase_todo,:), ~, ~, ~, ~, ~, ~] = Function_particle_size_ChordLength_volume_Algorithm(binary_phase,size_filter);
                % % Time
                timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase)/number_dimension toc(time_clock_start_phase)/number_dimension]];
            end
        end

        fitted_diameters_resized = fitted_diameters_resized*current_voxel_size;
        for current_direction_todo=1:1:number_dimension
            property_voxelsizedependence(current_iteration+1,2:end,current_direction_todo)=fitted_diameters_resized(:,current_direction_todo)';
        end

        % CPU and stopwatch time - end
        timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume)/number_dimension toc(time_clock_start_volume)/number_dimension]];

    end
    clear Phase_microstructure_resized;

    % Sort per voxel size
    for current_direction_todo=1:1:number_dimension
        property_voxelsizedependence(:,:,current_direction_todo) = sortrows(property_voxelsizedependence(:,:,current_direction_todo),1);  
    end

    %% EXTRAPOLATION TO 0 nm
    x=property_voxelsizedependence(:,1,1);

    tmp = zeros(number_resize+2,number_phase_todo+1,number_dimension); % + 0 nm and + initial voxel size
    if strcmp(p.scaling_extrapolation,'Linear')
        interpolation_voxelsize_order=1;
    elseif strcmp(p.scaling_extrapolation,'Quadratic')
        interpolation_voxelsize_order=2;
    elseif strcmp(p.scaling_extrapolation,'Cubic')
        interpolation_voxelsize_order=3;
    end
    max_order = length(p.scaling)-1;
    interpolation_voxelsize_order = min(interpolation_voxelsize_order,max_order);
    for current_phase_todo=1:1:number_phase_todo
        current_direction_todo = 0;
        for current_direction=1:1:number_dimension
                current_direction_todo=current_direction_todo+1;
                y=property_voxelsizedependence(:,current_phase_todo+1,current_direction_todo);
                pi = polyfit(x,y,interpolation_voxelsize_order);
                vq = polyval(pi,0);
                tmp(1,current_phase_todo+1,current_direction_todo)=vq;
                interpolation_voxelsize(current_phase_todo,current_direction_todo).pi=pi;
                % For correlation
                results_correlation(current_phase_todo).(['Particle_diameter_alongaxe' num2str(current_direction) 'CHORD_extrapolated']) = vq;
        end
    end
    tmp(2:end,:,:) = property_voxelsizedependence(:,:,:);
    property_voxelsizedependence = tmp; clear tmp;

    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
    for current_phase_todo=1:1:number_phase_todo
        Variable_name_table(1+current_phase_todo)=phasename_todo(current_phase_todo);
    end
    % Table
    for current_direction = 1:1:number_dimension
        Table_chord_d50_voxelsizedependence(current_direction).table = array2table(property_voxelsizedependence(:,:,current_direction),...
            'VariableNames',Variable_name_table);
    end 

    %% DISPLAY TEXT RESULTS
    fprintf('> Extrapolation to zero voxel size: polynomial of order %i\n\n',interpolation_voxelsize_order);
    for current_direction = 1:1:number_dimension
        fprintf('> Along %s\n\n',char(infovol.directionname(current_direction)));
        disp(Table_chord_d50_voxelsizedependence(current_direction).table)
    end
    
    %% SAVE RESULTS
    Results_chord.Table_chord_d50_voxelsizedependence = Table_chord_d50_voxelsizedependence; % Save in main table result
    if opt.save.xls
        filename = 'CHORD_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        sheet=0;
        for current_direction = 1:1:number_dimension
            sheet=sheet+1;
            DATA_writetable.sheet(sheet).name=char(infovol.directionname(current_direction));
            DATA_writetable.sheet(sheet).table=Table_chord_d50_voxelsizedependence(current_direction).table;
        end   
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    %% FIGURES
    parameters_figure.plotlog = false;
    parameters_figure.propertyname = 'Mean diameter';
    parameters_figure.method = 'Chord length';
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,1);
    parameters_figure.number_phase = number_phase_todo;
    parameters_figure.str_ylabel = ['D_{50} along ' char(infovol.directionname(1)) ' ' axisunit];
    parameters_figure.propertynameunit = Dunit;
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,1);
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = 'D1_Chord_voxel_size_dependence';
    parameters_figure.infovol = infovol;
    parameters_figure.opt = opt;
    parameters_figure.todo = p.todo;
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures    

    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,2);
    parameters_figure.str_ylabel = ['D_{50} along ' char(infovol.directionname(2)) ' ' axisunit];
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,2);
    parameters_figure.filename = 'D2_Chord_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures  

    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,3);
    parameters_figure.str_ylabel = ['D_{50} along ' char(infovol.directionname(3)) ' ' axisunit];
    parameters_figure.interpolation_voxelsize = interpolation_voxelsize(:,3);
    parameters_figure.filename = 'D3_Chord_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures  
end

%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%

if p.RVE.number_RVE>0
    n_property = 3;
    str_property(1).corrname = 'Particle_diameter_alongaxe1_CHORD';
    str_property(2).corrname = 'Particle_diameter_alongaxe2_CHORD';
    str_property(3).corrname = 'Particle_diameter_alongaxe3_CHORD';
    str_property(1).propertyname = ['Mean diameter along ' char(infovol.directionname(1))];
    str_property(2).propertyname = ['Mean diameter along ' char(infovol.directionname(2))];
    str_property(3).propertyname = ['Mean diameter along ' char(infovol.directionname(3))];
    res(1).Wholevolume_results = fitted_diameters(:,1);
    res(2).Wholevolume_results = fitted_diameters(:,2);
    res(3).Wholevolume_results = fitted_diameters(:,3);
    
    for k_RVE = 1:1:p.RVE.number_RVE % Loop over all RVE asked
        RVEparameters = p.RVE.RVE(k_RVE);
        for property_RVE=1:1:n_property
            Rk(property_RVE).RVE(k_RVE).RVEparameters = RVEparameters;
        end
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
            Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_phase_todo,3,3, n_property); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
            % 1 fitted diameter, 2 mean distance to surface
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
                Result_nestedRVE(k_nestedRVE+1,1,:,:,1,:) = Wholevolume_size(1);
                if k_nestedRVE ==1
                    fprintf('       Nested analysis\n');
                end
                if k_nestedRVE > 0
                    fprintf('          Cropping iteration: #%i/%i (cubic root volume=%1.1f %s)\n',k_nestedRVE,n_nestedRVE,Wholevolume_size(1),infovol.unit);
                end                
                if strcmp(RVEparameters.type,'C')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2,:) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'D')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3,:) = Wholevolume_size(2);                    
                elseif strcmp(RVEparameters.type,'G')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,3,:) = Wholevolume_size(2);
                elseif strcmp(RVEparameters.type,'H')
                    Result_nestedRVE(k_nestedRVE+1,1,:,:,2,:) = Wholevolume_size(2);                    
                end
            end

            % Information about subdomains
            for property_RVE=1:1:n_property
                Rk(property_RVE).RVE(k_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),All_subdomain(:,12),...
                    'VariableNames',{'Subdomain Id' 'Group Id' 'Number subdomain' 'Equivalent cubic length' 'Equivalent section length' 'Length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
            end
            [number_subdomain,~] = size(All_subdomain); % The number of subdomain
            number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size

            %% ALGORITHM
            % Colunm 1 is the subdomain id
            % Colunm 2 and 3 are the sizes of the subdomain.
            Property_eachsubdomain = zeros(number_subdomain,number_phase_todo+4,n_property);
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

                Property_eachsubdomain(subdomain_id,1,:)=subdomain_id;
                % Equivalent size of the subdomain
                Property_eachsubdomain(subdomain_id,2,:)=All_subdomain(subdomain_id,4)*voxel_size; % Cubic root length
                Property_eachsubdomain(subdomain_id,3,:)=All_subdomain(subdomain_id,5)*voxel_size; % Square root length
                Property_eachsubdomain(subdomain_id,4,:)=All_subdomain(subdomain_id,6)*voxel_size; % Length

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
                        [fitted_diameters_sub(current_phase_todo,:), ~, ~, ~, ~, ~, ~] = Function_particle_size_ChordLength_volume_Algorithm(binary_phase,size_filter);
                        Property_eachsubdomain(subdomain_id,current_phase_todo+4,1)=fitted_diameters_sub(current_phase_todo,1)* voxel_size;
                        Property_eachsubdomain(subdomain_id,current_phase_todo+4,2)=fitted_diameters_sub(current_phase_todo,2)* voxel_size;
                        Property_eachsubdomain(subdomain_id,current_phase_todo+4,3)=fitted_diameters_sub(current_phase_todo,3)* voxel_size;

                        % % Time
                        timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

                    end
                end

                % CPU and stopwatch time - end
                timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
            end

            for property_RVE=1:1:n_property
                %% STATISTICAL ANALYSIS and RVE SIZE
                [Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence] = Function_subdomains_statistical_analysis(number_group_size,number_phase_todo,GROUP_SUBDOMAIN,Property_eachsubdomain(:,:,property_RVE),voxel_size,RVEparameters);

                %% SAVE FOR CORRELATION
                if k_nestedRVE == 0
                    for k_threshold=1:1:n_threshold
                        current_phase_todo = 0;
                        for current_phase=1:1:number_phase
                            if p.todo(current_phase)
                                current_phase_todo=current_phase_todo+1;
                                if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                    if Size_RVE(1,current_phase_todo,2,1)~=0
                                        str_ = [str_property(property_RVE).corrname '_RVE_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,1);
                                    end
                                    if strcmp(RVEparameters.type,'C')
                                        if Size_RVE(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_RVE_squarerootFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                            str_(str_=='.')='p';
                                            results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                        end
                                    elseif strcmp(RVEparameters.type,'D')
                                        if Size_RVE(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_RVE_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                            str_(str_=='.')='p';
                                            results_correlation(current_phase_todo).(str_) = Size_RVE(1,current_phase_todo,2,2);
                                        end
                                    end
                                else
                                    if Size_convergence(1,current_phase_todo,2,1)~=0
                                        str_ = [str_property(property_RVE).corrname '_conv_cubicroot_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                        str_(str_=='.')='p';
                                        results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,1);
                                    end
                                    if strcmp(RVEparameters.type,'G')
                                        if Size_convergence(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_conv_lengthFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
                                            str_(str_=='.')='p';
                                            results_correlation(current_phase_todo).(str_) = Size_convergence(1,current_phase_todo,2,2);
                                        end
                                    elseif strcmp(RVEparameters.type,'H')
                                        if Size_convergence(1,current_phase_todo,2,2)~=0
                                            str_ = [str_property(property_RVE).corrname '_conv_areaFOV_thresh' num2str(thresholds(k_threshold),'%1.2f') '_' RVEparameters.savename];
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
                Rk(property_RVE).RVE = Function_subdomains_manage_results(Property_eachsubdomain(:,:,property_RVE), Property_subdomains_statistics, Size_RVE, derivative_convergence, relativedifference_convergence, Size_convergence, RVEparameters,Rk(property_RVE).RVE,k_RVE,number_phase,number_phase_todo,infovol,p);

                if RVEparameters.donested
                    for k_threshold=1:1:n_threshold
                        current_phase_todo = 0;
                        for current_phase=1:1:number_phase
                            if p.todo(current_phase)
                                current_phase_todo=current_phase_todo+1;
                                if strcmp(RVEparameters.type,'A') || strcmp(RVEparameters.type,'B') || strcmp(RVEparameters.type,'C') || strcmp(RVEparameters.type,'D')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1,property_RVE) = Size_RVE(k_threshold,current_phase_todo,1,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1,property_RVE) = Size_RVE(k_threshold,current_phase_todo,2,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1,property_RVE) = Size_RVE(k_threshold,current_phase_todo,3,1);
                                else
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,1,property_RVE) = Size_convergence(k_threshold,current_phase_todo,1,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,1,property_RVE) = Size_convergence(k_threshold,current_phase_todo,2,1);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,1,property_RVE) = Size_convergence(k_threshold,current_phase_todo,3,1);
                                end
                                if strcmp(RVEparameters.type,'C')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2,property_RVE) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2,property_RVE) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2,property_RVE) = Size_RVE(k_threshold,current_phase_todo,3,2);
                                elseif strcmp(RVEparameters.type,'D')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3,property_RVE) = Size_RVE(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3,property_RVE) = Size_RVE(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3,property_RVE) = Size_RVE(k_threshold,current_phase_todo,3,2);
                                elseif strcmp(RVEparameters.type,'G')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,3,property_RVE) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,3,property_RVE) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,3,property_RVE) = Size_convergence(k_threshold,current_phase_todo,3,2);
                                elseif strcmp(RVEparameters.type,'H')
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,1,2,property_RVE) = Size_convergence(k_threshold,current_phase_todo,1,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,2,2,property_RVE) = Size_convergence(k_threshold,current_phase_todo,2,2);
                                    Result_nestedRVE(k_nestedRVE+1,k_threshold+1,current_phase_todo,3,2,property_RVE) = Size_convergence(k_threshold,current_phase_todo,3,2);
                                end
                            end
                        end
                    end
                end

                %% TEXT DISPLAY AND SAVE RESULTS
                if k_nestedRVE == 0
                    if property_RVE==1
                        RVEparameters.disp_parRVE = true;
                    else
                        RVEparameters.disp_parRVE = false;
                    end
                    Function_subdomains_display_and_save(Rk(property_RVE).RVE,k_RVE,RVEparameters,number_phase,number_phase_todo,str_property(property_RVE).propertyname,Sub_folder_RVE,opt,infovol,p);
                end

                %% FIGURES
                if k_nestedRVE == 0
                    parameters_figure.propertyname = str_property(property_RVE).propertyname;
                    parameters_figure.propertynameunit = Dunit;
                    parameters_figure.RVE = RVEparameters;
                    parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                    parameters_figure.savefolder = Sub_folder_RVE;
                    parameters_figure.number_phase = number_phase;
                    parameters_figure.number_phase_todo = number_phase_todo;
                    parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
                    parameters_figure.Property_eachsubdomain = Property_eachsubdomain(:,:,property_RVE);
                    parameters_figure.derivative_convergence = derivative_convergence;
                    parameters_figure.relativedifference_convergence = relativedifference_convergence;
                    parameters_figure.Size_RVE = Size_RVE;
                    parameters_figure.convergence_criterion = RVEparameters.threshold_reldiff;
                    parameters_figure.Size_convergence = Size_convergence;
                    parameters_figure.Wholevolume_size = Wholevolume_size;
                    parameters_figure.Wholevolume_results = res(property_RVE).Wholevolume_results;
                    parameters_figure.infovol = infovol;
                    parameters_figure.todo = p.todo;
                    parameters_figure.opt = opt;
                    Function_create_figures_RVE(parameters_figure) % Figures
                end
            end
        end

        %% NESTED ANALYSIS RESULT
        if RVEparameters.donested
            for property_RVE=1:1:n_property
                % Table
                [Rk(property_RVE).RVE] = Function_nestedtable(Rk(property_RVE).RVE,k_RVE,RVEparameters,number_phase,str_property(property_RVE).propertyname,Result_nestedRVE(:,:,:,:,:,property_RVE),Sub_folder_RVE,opt,infovol,p);
                % Figure
                parameters_figure.propertyname = str_property(property_RVE).propertyname;
                parameters_figure.Result_nestedRVE = Result_nestedRVE(:,:,:,:,:,property_RVE);
                Function_create_figures_nestedRVE(parameters_figure) % Figures
                % Save
                Rk(property_RVE).RVE(k_RVE).nestedanalysis = Result_nestedRVE(:,:,:,:,:,property_RVE);
            end
        end

    end
    for property_RVE=1:1:n_property
        Results_chord.RVE.(str_property(property_RVE).corrname) = Rk(property_RVE).RVE; % Save in main table result
    end    

end

%%
%% ENDING FUNCTION
%%
        
%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_chord.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_chord.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_chord.Table_date = Table_date;

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
    Function_Writetable(Current_folder,'D50_Chord_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'D50_Chord_calculation_time', 'Diameter (Chord)', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_D50_chord'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_D50_chord'],'Results_chord')
end
%% SAVE VISUALIZATION
Current_folder = [infovol.volpath 'Visualization' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Visualization_particlediameter_chord'],'results_visualization','-v7.3');    

end