function [] = Function_Tortuosity_factor_taufactor(Phase_microstructure,infovol,opt,p)
% Calculate tortuosity factor (normalized effective diffusion coefficient)
% Function_Tortuosity_factor_taufactor(Phase_microstructure, infovol, opt, p) - when used with the MATBOX toolbox
% or
% Function_Tortuosity_factor_taufactor(Phase_microstructure, labels, directions) - when used as a standalone function
% with: Phase_microstructure, a 3D array: the 3D segmented volumes
%       labels, a 1D array: which phases to investigate
%       directions, a 1D array: which directions to investigate
%       e.g.: Function_Tortuosity_factor_taufactor(<your_3d_array>, 0, [1,2,3]);

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 3 % Case for function called as: Function_Tortuosity_factor_taufactor(Phase_microstructure, labels, directions). Standalone use.
        labels = infovol; clear infovol;
        directions = opt; clear opt;
                
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        infovol.volumesubfolder = ['Tortuosity_' char(t)];
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
            if sum(labels==infovol.phaselabel(k)) % Phase investigated
                p.todo(k)=1;
            else
                p.todo(k)=0;
            end
        end
        % Directions investigated
        for k=1:1:3
            if sum(directions==k)
                p.direction_todo(k)=1;
            else
                p.direction_todo(k)=0;
            end
        end        
        
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

    else % Incorrect number of argument
        disp 'Error calling Function_Tortuosity_factor_taufactor. Wrong number of argument.'
        Function_Tortuosity_factor_taufactor
        return
    end
    
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    separator = '\';
else
    separator = '/';
end
Current_folder = [infovol.volpath 'Tortuosity' separator];
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

%% INITIALIZE RESULTS (USE FOR CORRELATION)
current_phase_todo = 0;
for current_phase=1:1:number_phase
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        results_correlation(current_phase_todo).name = infovol.phasename(current_phase,1);
    end
end

%%
%% ALGORITHM ON WHOLE VOLUME
%%

disp '    TORTUOSITY FACTOR';
disp '    -----------------';
disp ' ';

%% CALCULATION
time_cpu_start_volume = cputime; % CPU start
time_clock_start_volume = tic; % Stopwatch start

% Initialization (generic)
number_phase_todo = sum(p.todo); % How many phase are we going to analyse ?
Numbervoxel_phase=zeros(number_phase_todo,1); 
timedata = zeros(number_phase_todo+1,3); timedata(1,1) = voxel_number;
timedata_domain = cell(number_phase_todo+1,1);
phasename_todo = cell(number_phase_todo,1);
phaselabel = zeros(number_phase_todo,1);
% Initialization (algorithm-specific)
number_direction_todo = sum(p.direction_todo);
Tortuosity_factor = zeros(number_phase_todo,number_direction_todo);
Bruggeman_exponent = zeros(number_phase_todo,number_direction_todo);
Normalized_effective_diffusion_coefficient = zeros(number_phase_todo,number_direction_todo);
Volumefraction = zeros(number_phase_todo,1);
if number_direction_todo==3
    Tortuosity_factor_anisotropy = zeros(number_phase_todo,3); % x/z, y/z, (x/z+y/z)/2
    Bruggeman_exponent_anisotropy = zeros(number_phase_todo,3);
    Normalized_effective_diffusion_coefficient_anisotropy = zeros(number_phase_todo,3);
end

current_phase_todo = 0;
for current_phase=1:1:number_phase % Loop over all phases
    if p.todo(current_phase)
        current_phase_todo=current_phase_todo+1;
        phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);
        phaselabel(current_phase_todo,1) = infovol.phaselabel(current_phase);
        % Create a binary microstructure : 1 = current analysed phase, 0 = complementay phase
        binary_phase=zeros(Domain_size); % Initialization
        binary_phase(Phase_microstructure == phaselabel(current_phase_todo,1)) = 1; % Binary phase
        % Volume fraction
        Numbervoxel_phase(current_phase_todo,1)= sum(sum(sum(Phase_microstructure==phaselabel(current_phase_todo,1) )));
        Volumefraction(current_phase_todo,1) = Numbervoxel_phase(current_phase_todo,1)/voxel_number; % Defintion of volume fraction

        time_cpu_start_phase = cputime; % CPU start
        time_clock_start_phase = tic; % Stopwatch start

        current_direction_todo = 0;
        for current_direction=1:1:number_dimension
            if p.direction_todo(current_direction)
                current_direction_todo=current_direction_todo+1;

                directionname_todo(current_direction_todo,1) = infovol.directionname(current_direction);
                % % Algorithm
                % Call tortuosity algorithm (cf. Tau factor manual)
                if current_direction==1
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                    if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
                        Tau_factor_result.Tau_W1.Tau = NaN;
                    end
                    Tortuosity_factor(current_phase_todo,current_direction_todo)=Tau_factor_result.Tau_W1.Tau;
                elseif current_direction==2
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                    if Tau_factor_result.Tau_W2.Tau == Inf
                        Tau_factor_result.Tau_W2.Tau = NaN;
                    end
                    Tortuosity_factor(current_phase_todo,current_direction_todo)=Tau_factor_result.Tau_W2.Tau;
                elseif current_direction==3
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                    if Tau_factor_result.Tau_W3.Tau == Inf
                        Tau_factor_result.Tau_W3.Tau = NaN;
                    end
                    Tortuosity_factor(current_phase_todo,current_direction_todo)=Tau_factor_result.Tau_W3.Tau;
                end
                % Effective diffusion coefficient and Bruggeman exponent
                Normalized_effective_diffusion_coefficient(current_phase_todo,current_direction_todo) = Volumefraction(current_phase_todo,1)/Tortuosity_factor(current_phase_todo,current_direction_todo);
                tau = Tortuosity_factor(current_phase_todo,current_direction_todo);
                eps = Volumefraction(current_phase_todo,1);
                Bruggeman_exponent(current_phase_todo,current_direction_todo) = 1 - log(tau)/log(eps);

                % % Correlation
                [str_direction] = function_remove_emptyandspecialcharacter_string(char(infovol.directionname(current_direction)));
                results_correlation(current_phase_todo).(['Tortuosity_' num2str(current_direction) '_' str_direction]) = Tortuosity_factor(current_phase_todo,current_direction_todo);
                results_correlation(current_phase_todo).(['Normalized_diffusivity_' num2str(current_direction) '_' str_direction]) = Normalized_effective_diffusion_coefficient(current_phase_todo,current_direction_todo);
                results_correlation(current_phase_todo).(['MacMullin_number_' num2str(current_direction) '_' str_direction]) = 1/Normalized_effective_diffusion_coefficient(current_phase_todo,current_direction_todo);
                results_correlation(current_phase_todo).(['Bruggeman_exponent_' num2str(current_direction) '_' str_direction]) = Bruggeman_exponent(current_phase_todo,current_direction_todo);
            end
        end

        % % Time
        timedata_domain(current_phase_todo+1,1) = infovol.phasename(current_phase,1);
        timedata(current_phase_todo+1,1) = Numbervoxel_phase(current_phase_todo,1);
        timedata(current_phase_todo+1,2) = (cputime-time_cpu_start_phase)/number_direction_todo; % CPU elapsed time
        timedata(current_phase_todo+1,3) = (toc(time_clock_start_phase))/number_direction_todo; % CPU elapsed time

        if number_direction_todo==3 % Anisotropy
            Tortuosity_factor_anisotropy(current_phase_todo,1) = Tortuosity_factor(current_phase_todo,1) / Tortuosity_factor(current_phase_todo,3);
            Tortuosity_factor_anisotropy(current_phase_todo,2) = Tortuosity_factor(current_phase_todo,2) / Tortuosity_factor(current_phase_todo,3);
            Tortuosity_factor_anisotropy(current_phase_todo,3) = (Tortuosity_factor_anisotropy(current_phase_todo,1) + Tortuosity_factor_anisotropy(current_phase_todo,2))/2;

            Bruggeman_exponent_anisotropy(current_phase_todo,1) = Bruggeman_exponent(current_phase_todo,1) / Bruggeman_exponent(current_phase_todo,3);
            Bruggeman_exponent_anisotropy(current_phase_todo,2) = Bruggeman_exponent(current_phase_todo,2) / Bruggeman_exponent(current_phase_todo,3);
            Bruggeman_exponent_anisotropy(current_phase_todo,3) = (Bruggeman_exponent_anisotropy(current_phase_todo,1) + Bruggeman_exponent_anisotropy(current_phase_todo,2))/2;

            Normalized_effective_diffusion_coefficient_anisotropy(current_phase_todo,1) = Normalized_effective_diffusion_coefficient(current_phase_todo,1) / Normalized_effective_diffusion_coefficient(current_phase_todo,3);
            Normalized_effective_diffusion_coefficient_anisotropy(current_phase_todo,2) = Normalized_effective_diffusion_coefficient(current_phase_todo,2) / Normalized_effective_diffusion_coefficient(current_phase_todo,3);
            Normalized_effective_diffusion_coefficient_anisotropy(current_phase_todo,3) = (Normalized_effective_diffusion_coefficient_anisotropy(current_phase_todo,1) + Normalized_effective_diffusion_coefficient_anisotropy(current_phase_todo,2))/2;            

            % Correlation
            results_correlation(current_phase_todo).Tortuosity_anisotropy = Tortuosity_factor_anisotropy(current_phase_todo,3);
            results_correlation(current_phase_todo).Normalized_effective_diffusion_coefficient_anisotropy = Normalized_effective_diffusion_coefficient_anisotropy(current_phase_todo,3);
            results_correlation(current_phase_todo).MacMullin_number = 1/Normalized_effective_diffusion_coefficient_anisotropy(current_phase_todo,3);
            results_correlation(current_phase_todo).Bruggeman_exponent_anisotropy = Bruggeman_exponent_anisotropy(current_phase_todo,3);     
        end

    end
end
% CPU and stopwatch time - end
timedata_domain(1,1) = {'Full volume'};
timedata(1,2) = (cputime-time_cpu_start_volume)/number_direction_todo; % CPU elapsed time
timedata(1,3) = (toc(time_clock_start_volume))/number_direction_todo; % Stopwatch elapsed time
timedata_pervolume = timedata(1,:);
timedata_perphase = timedata(2:end,:);

%% TABLES
% Time
Table_time = table(timedata_domain(:,1), timedata(:,1),timedata(:,2),timedata(:,3),...
    'VariableNames',{'Domain', 'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Tortuosity.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
for current_phase_todo=1:1:number_phase_todo
   Tortuositytaufactor.phase(current_phase_todo).table = table(infovol.directionname(p.direction_todo),ones(number_direction_todo,1)*Volumefraction(current_phase_todo,1),Tortuosity_factor(current_phase_todo,:)',Bruggeman_exponent(current_phase_todo,:)',Normalized_effective_diffusion_coefficient(current_phase_todo,:)',...
        'VariableNames',{'Direction','Volume fraction','Tortuosity factor','Bruggeman exponent','Normalized Deff',});
end
Results_Tortuosity.Tortuositytaufactor = Tortuositytaufactor; % Save in main table result

%% SAVE TABLES
if opt.save.xls
    filename = 'Tortuosity'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    for current_phase_todo=1:1:number_phase_todo
        DATA_writetable.sheet(current_phase_todo).name=char(phasename_todo(current_phase_todo,1));
        DATA_writetable.sheet(current_phase_todo).table=Tortuositytaufactor.phase(current_phase_todo).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
for current_phase_todo=1:1:number_phase_todo
    fprintf('  For the phase %s:\n\n',char(phasename_todo(current_phase_todo,1)));    
    disp(Tortuositytaufactor.phase(current_phase_todo).table)
end
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%% FIGURE
scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Transport properties';
Fig.Color='white'; % Background colour
set(Fig,'position',round([scrsz(1) scrsz(2) 2/3*scrsz(3)  2/3*scrsz(4)]) ); % Full screen figure
% - Create axes as a subplot
X = categorical(phasename_todo);
X = reordercats(X,phasename_todo);
for id_axe=1:1:4
    sub_axes=subplot(2,2,id_axe,'Parent',Fig);
    hold(sub_axes,'on'); % Active subplot
    if id_axe==1
        h_title=title ('Volume fraction \epsilon'); % Title
        h_=bar(X,Volumefraction,0.5); % Width=1
        ylabel('Volume fraction \epsilon'); % Axis label
    elseif id_axe==2
        h_title=title ('Tortuosity factor \tau'); % Title
        h_=bar(X,Tortuosity_factor,1); % Width=1
        ylabel('Tortuosity factor \tau'); % Axis label
        ylim([1 Inf]);
    elseif id_axe==3
        h_title=title ('Bruggeman exponent p, \tau=\epsilon^{1-p}'); % Title
        h_=bar(X,Bruggeman_exponent,1); % Width=1
        ylabel('Bruggeman exponent p'); % Axis label
        xls = xlim;
        plot([xls(1) xls(2)],[1 1],'Color','k','LineStyle','--','LineWidth',2);
        plot([xls(1) xls(2)],[1.5 1.5],'Color','k','LineStyle','--','LineWidth',2);
        plot([xls(1) xls(2)],[2.0 2.0],'Color','k','LineStyle','--','LineWidth',2);
        xt = [xls(2) xls(2) xls(2)];
        yt = [1 1.5 2.0];
        str = {'Rule of mixture','Spheres','Cylinders'};
        text(xt,yt,str,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize);
    elseif id_axe==4
        h_title=title ('Normalized effective diffusion coefficient'); % Title
        h_=bar(X,Normalized_effective_diffusion_coefficient,1); % Width=1
        ylabel('D_{eff}/D_{bulk}'); % Axis label
    end
    if id_axe>1
        h_legend = legend(sub_axes,infovol.directionname(p.direction_todo),'Location','best'); % Legend
    end
    set(sub_axes,'YGrid',opt.format.grid); % Display grid
    set(sub_axes,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
    set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
    h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
    if id_axe>1
        h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
    end
    hold(sub_axes,'off'); % Relase figure
end
sgtitle(Fig,'Transport properties','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
if opt.save.savefig % Save figure
    filename= 'Transport_properties';
    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
end
if opt.format.autoclosefig
    close(Fig); % Do not keep open figures
end

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%

%% ALONG DIRECTIONS

if sum(p.direction_section)>number_dimension
    fprintf('> Calculated section per section:\n\n');
    section_direction_todo=0;
    for current_sectiondirection=1:1:number_dimension
        if p.direction_section(current_sectiondirection)>1
            section_direction_todo = section_direction_todo+1;
            for current_transportdirection=1:1:number_dimension
                if p.direction_todo(current_transportdirection)
                    % Initialization
                    Section(current_sectiondirection).transport(1).porosity = zeros(p.direction_section(current_sectiondirection),number_phase_todo+2);
                    Section(current_sectiondirection).transport(current_transportdirection).tortuosity = zeros(p.direction_section(current_sectiondirection),number_phase_todo+2);
                    Section(current_sectiondirection).transport(current_transportdirection).Bruggeman = zeros(p.direction_section(current_sectiondirection),number_phase_todo+2);
                    Section(current_sectiondirection).transport(current_transportdirection).Deff = zeros(p.direction_section(current_sectiondirection),number_phase_todo+2);                        

                    % Section bounds
                    xs = round(linspace(1,Domain_size(current_sectiondirection),p.direction_section(current_sectiondirection)+1));
                    xstart = xs(1:end-1);
                    xend = [xstart(2:end)-1 xs(end)];
                    Section(current_sectiondirection).transport(1).porosity(:,1)=xstart; Section(current_sectiondirection).transport(1).porosity(:,2)=xend;
                    Section(current_sectiondirection).transport(current_transportdirection).tortuosity(:,1)=xstart; Section(current_sectiondirection).transport(current_transportdirection).tortuosity(:,2)=xend;
                    Section(current_sectiondirection).transport(current_transportdirection).Bruggeman(:,1)=xstart; Section(current_sectiondirection).transport(current_transportdirection).Bruggeman(:,2)=xend;
                    Section(current_sectiondirection).transport(current_transportdirection).Deff(:,1)=xstart; Section(current_sectiondirection).transport(current_transportdirection).Deff(:,2)=xend;

                    % Calculation
                    current_phase_todo = 0;
                    for current_phase=1:1:number_phase % Loop over all phases
                        if p.todo(current_phase)
                            current_phase_todo=current_phase_todo+1;
                            binary_phase=zeros(Domain_size); % Initialization
                            binary_phase(Phase_microstructure == phaselabel(current_phase_todo,1)) = 1; % Binary phase

                            for ksection = 1:1:p.direction_section(current_sectiondirection)
                                if current_sectiondirection==1
                                    BWsection = binary_phase(xstart(ksection):xend(ksection),:,:);
                                elseif current_sectiondirection==2
                                    BWsection = binary_phase(:,xstart(ksection):xend(ksection),:);
                                else
                                    BWsection = binary_phase(:,:,xstart(ksection):xend(ksection));
                                end

                                % Porosity
                                Section(current_sectiondirection).transport(1).porosity(ksection,current_phase_todo+2) = sum(sum(sum(BWsection)))/numel(BWsection);
                                
                                % Tortuosity
                                if current_transportdirection==1
                                    Tau_factor_result = TauFactor('InLine',1,0,0,BWsection,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                                    if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
                                        Tau_factor_result.Tau_W1.Tau = NaN;
                                    end
                                    Section(current_sectiondirection).transport(current_transportdirection).tortuosity(ksection,current_phase_todo+2)=Tau_factor_result.Tau_W1.Tau;
                                elseif current_transportdirection==2
                                    Tau_factor_result = TauFactor('InLine',1,0,0,BWsection,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                                    if Tau_factor_result.Tau_W2.Tau == Inf
                                        Tau_factor_result.Tau_W2.Tau = NaN;
                                    end
                                    Section(current_sectiondirection).transport(current_transportdirection).tortuosity(ksection,current_phase_todo+2)=Tau_factor_result.Tau_W2.Tau;
                                elseif current_transportdirection==3
                                    Tau_factor_result = TauFactor('InLine',1,0,0,BWsection,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                                    if Tau_factor_result.Tau_W3.Tau == Inf
                                        Tau_factor_result.Tau_W3.Tau = NaN;
                                    end
                                    Section(current_sectiondirection).transport(current_transportdirection).tortuosity(ksection,current_phase_todo+2)=Tau_factor_result.Tau_W3.Tau;
                                end

                                % Effective diffusion coefficient and Bruggeman exponent
                                Section(current_sectiondirection).transport(current_transportdirection).Deff(ksection,current_phase_todo+2) = Section(current_sectiondirection).transport(1).porosity(ksection,current_phase_todo+2) / Section(current_sectiondirection).transport(current_transportdirection).tortuosity(ksection,current_phase_todo+2);
                                Section(current_sectiondirection).transport(current_transportdirection).Bruggeman(ksection,current_phase_todo+2) = 1 - log(Section(current_sectiondirection).transport(current_transportdirection).tortuosity(ksection,current_phase_todo+2))/log(Section(current_sectiondirection).transport(1).porosity(ksection,current_phase_todo+2));
                            end
                        end
                    end
                else
                    Section(current_sectiondirection).transport(1).porosity = [];
                    Section(current_sectiondirection).transport(current_transportdirection).tortuosity = [];
                    Section(current_sectiondirection).transport(current_transportdirection).Bruggeman = [];
                    Section(current_sectiondirection).transport(current_transportdirection).Deff = [];                    
                end
            end
        else
            Section(current_sectiondirection).transport = [];
        end
    end

    % % Figures
    strunit = voxel_unit;
    if strcmp(strunit,'um') || strcmp(strunit,'micrometer') || strcmp(strunit,'Micrometer') || strcmp(strunit,'micrometers') || strcmp(strunit,'Micrometers')
        axisunit = '(\mum)';
    else
        axisunit = ['(' strunit ')'];
    end

    % Porosity
    Fig = figure; % Create figure
    Fig.Name= 'Volume fractions, section per section';
    Fig.Color='white'; % Background colour
    set(Fig,'position',round([scrsz(1) scrsz(2) 4/5*scrsz(3)  1/2*scrsz(4)]) );
    id_axe=0;
    for current_sectiondirection=1:1:number_dimension
        if p.direction_section(current_sectiondirection)>1
            id_axe=id_axe+1;
            sub_axes=subplot(1,section_direction_todo,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            h_title=title('Volume fraction');
            t_ = xlabel(' ');
            t_1 = sprintf('Position along %s ',char(infovol.directionname(current_sectiondirection)));
            t_2 = axisunit;
            t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
            ylabel('Volume fraction \epsilon');
            
            tmp = Section(current_sectiondirection).transport(1).porosity;
            xa = tmp(:,1); xb = tmp(:,2);
            current_phase_todo=0;
            for current_phase=1:1:number_phase % Loop over phases
                if p.todo(current_phase)
                    current_phase_todo=current_phase_todo+1;
                    y = tmp(:,current_phase_todo+2);
 
                    xall = []; yall = [];
                    for k=1:1:length(y)
                        xall = [xall xa(k) xb(k)]; yall = [yall y(k) y(k)];
                        if k<length(y)
                            xall = [xall xb(k) xb(k)];
                            yall = [yall y(k) y(k+1)];
                        end
                    end
                    plot(xall*voxel_size,yall,'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));
                    
                    vf_whole = Tortuositytaufactor.phase(current_phase_todo).table(1,2).Variables;
                    plot([min(xall) max(xall)]*voxel_size,[vf_whole vf_whole],'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'LineStyle','--','DisplayName',[char(infovol.phasename(current_phase,1)) ' (whole volume)' ]);

                end
            end
            legend(sub_axes,'Location','best');
            set(sub_axes,'YGrid',opt.format.grid); % Display grid
            set(sub_axes,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
    end
    sgtitle(Fig,'Volume fraction \epsilon per volume section','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
    if opt.save.savefig % Save figure
        filename= 'Transport_properties_persection_volumefractions';
        function_savefig(Fig, Current_folder, filename, opt.save); % Call function
    end
    if opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end

    % Tortuosity factor, Bruggeman exponent, and normalized effective diffusion coefficient
    strfig(1).figname = 'Tortuosity factor'; strfig(2).name = 'Bruggeman exponent'; strfig(3).name = 'Normalized effective diffusion coefficient';
    strfig(1).titlename = 'Tortuosity factor'; strfig(2).titlename = 'Bruggeman exponent'; strfig(3).titlename = 'D_{eff}/D_{bulk}';
    strfig(1).symbol = 'Tortuosity factor \tau'; strfig(2).symbol = 'Bruggeman exponent p'; strfig(3).symbol = 'D_{eff}/D_{bulk}';
    strfig(1).filename = 'Tortuosity'; strfig(2).filename = 'Bruggeman'; strfig(3).filename = 'Deff';
    for kfig=1:1:3
        Fig = figure; % Create figure
        Fig.Name= [strfig(kfig).figname ', section per section'];
        Fig.Color='white'; % Background colour
        set(Fig,'position',round([scrsz(1) scrsz(2) 4/5*scrsz(3)  4/5*scrsz(4)]) );
        id_axe=0;
        transportdirection_todo=0;
        for current_transportdirection=1:1:number_dimension
            if p.direction_todo(current_transportdirection)
                transportdirection_todo = transportdirection_todo+1;
                for current_sectiondirection=1:1:number_dimension
                    if p.direction_section(current_sectiondirection)>1
                        id_axe=id_axe+1;
                        sub_axes=subplot(number_direction_todo,section_direction_todo,id_axe,'Parent',Fig);
                        hold(sub_axes,'on'); % Active subplot
                        h_title=title([strfig(kfig).titlename ' along ' char(infovol.directionname(current_transportdirection))]);
                        t_ = xlabel(' ');
                        t_1 = sprintf('Position along %s ',char(infovol.directionname(current_sectiondirection)));
                        t_2 = axisunit;
                        t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
                        strtr = sprintf('_{%s}',char(infovol.directionname(current_transportdirection)));
                        ylabel([strfig(kfig).symbol strtr ]);
                        if kfig==1
                            tmp = Section(current_sectiondirection).transport(current_transportdirection).tortuosity;
                        elseif kfig==2
                            tmp = Section(current_sectiondirection).transport(current_transportdirection).Bruggeman;
                        else
                            tmp = Section(current_sectiondirection).transport(current_transportdirection).Deff;
                        end
                        xa = tmp(:,1); xb = tmp(:,2);
                        current_phase_todo=0;
                        for current_phase=1:1:number_phase % Loop over phases
                            if p.todo(current_phase)
                                current_phase_todo=current_phase_todo+1;
                                y = tmp(:,current_phase_todo+2);

                                xall = []; yall = [];
                                for k=1:1:length(y)
                                    xall = [xall xa(k) xb(k)]; yall = [yall y(k) y(k)];
                                    if k<length(y)
                                        xall = [xall xb(k) xb(k)];
                                        yall = [yall y(k) y(k+1)];
                                    end
                                end
                                plot(xall*voxel_size,yall,'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'DisplayName',char(infovol.phasename(current_phase,1)));

                                p_whole = Tortuositytaufactor.phase(current_phase_todo).table(transportdirection_todo,kfig+2).Variables;
                                plot([min(xall) max(xall)]*voxel_size,[p_whole p_whole],'Color', infovol.phasecolor(current_phase,:),'LineWidth',opt.format.linewidth,'LineStyle','--','DisplayName',[char(infovol.phasename(current_phase,1)) ' (whole volume)' ]);
                            end
                        end
                        legend(sub_axes,'Location','best');
                        set(sub_axes,'YGrid',opt.format.grid); % Display grid
                        set(sub_axes,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
                        set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
                        h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
                        hold(sub_axes,'off'); % Relase figure
                    end
                end
            end
        end
        sgtitle(Fig,[strfig(kfig).symbol ' per volume section'],'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
        if opt.save.savefig % Save figure
            filename= ['Transport_properties_persection_' strfig(kfig).filename];
            function_savefig(Fig, Current_folder, filename, opt.save); % Call function
        end
        if opt.format.autoclosefig
            close(Fig); % Do not keep open figures
        end
    end

    % % Tables
    clear Variable_name_table;
    Variable_name_table={['Position start ' voxel_unit] ['Position end ' voxel_unit]};
    for current_phase_todo=1:1:number_phase_todo
        Variable_name_table(current_phase_todo+2)=phasename_todo(current_phase_todo);
    end
    % Create table
    for current_sectiondirection=1:1:number_dimension
        if p.direction_section(current_sectiondirection)>1
            tmp = Section(current_sectiondirection).transport(1).porosity; tmp(:,1:2) = tmp(:,1:2)*voxel_size;
            Evolution.sectiondirection(current_sectiondirection).volumefractions_table = array2table(tmp,'VariableNames',Variable_name_table);
            for current_transportdirection=1:1:number_dimension
                if p.direction_todo(current_transportdirection)
                    tmp = Section(current_sectiondirection).transport(current_transportdirection).tortuosity; tmp(:,1:2) = tmp(:,1:2)*voxel_size;
                    Evolution.sectiondirection(current_sectiondirection).transportdirection(current_transportdirection).tortuosity_table = array2table(tmp,'VariableNames',Variable_name_table);
                    tmp = Section(current_sectiondirection).transport(current_transportdirection).Bruggeman; tmp(:,1:2) = tmp(:,1:2)*voxel_size;
                    Evolution.sectiondirection(current_sectiondirection).transportdirection(current_transportdirection).Bruggeman_table = array2table(tmp,'VariableNames',Variable_name_table);
                    tmp = Section(current_sectiondirection).transport(current_transportdirection).Deff; tmp(:,1:2) = tmp(:,1:2)*voxel_size;
                    Evolution.sectiondirection(current_sectiondirection).transportdirection(current_transportdirection).Deff_table = array2table(tmp,'VariableNames',Variable_name_table);
                end
            end
        end
    end

    % % Save tables
    if opt.save.xls
        filename = 'Transport_properties_persection'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        sheet_=0;
        for current_sectiondirection=1:1:number_dimension
            if p.direction_section(current_sectiondirection)>1
                sheet_ = sheet_+1;
                DATA_writetable.sheet(sheet_).name=['Vf, section ' num2str(current_sectiondirection)];
                DATA_writetable.sheet(sheet_).table= Evolution.sectiondirection(current_sectiondirection).volumefractions_table;
                for current_transportdirection=1:1:number_dimension
                    if p.direction_todo(current_transportdirection)
                        sheet_ = sheet_+1;
                        DATA_writetable.sheet(sheet_).name=['Tau, section ' num2str(current_sectiondirection) ', transport ' num2str(current_transportdirection)];
                        DATA_writetable.sheet(sheet_).table= Evolution.sectiondirection(current_sectiondirection).transportdirection(current_transportdirection).tortuosity_table;
                        sheet_ = sheet_+1;
                        DATA_writetable.sheet(sheet_).name=['P, section ' num2str(current_sectiondirection) ', transport ' num2str(current_transportdirection)];
                        DATA_writetable.sheet(sheet_).table= Evolution.sectiondirection(current_sectiondirection).transportdirection(current_transportdirection).Bruggeman_table;
                        sheet_ = sheet_+1;
                        DATA_writetable.sheet(sheet_).name=['Deff, section ' num2str(current_sectiondirection) ', transport ' num2str(current_transportdirection)];
                        DATA_writetable.sheet(sheet_).table= Evolution.sectiondirection(current_sectiondirection).transportdirection(current_transportdirection).Deff_table;
                    end
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
    property_voxelsizedependence = zeros(number_resize+1,number_phase_todo+1,number_direction_todo+1);
    property_voxelsizedependence(1,1,:)=voxel_size;
    for current_direction_todo=1:1:number_direction_todo
        property_voxelsizedependence(1,2:end,current_direction_todo)=Tortuosity_factor(:,current_direction_todo)';
    end
    property_voxelsizedependence(1,2:end,number_direction_todo+1)=Volumefraction';

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
                property_voxelsizedependence(current_iteration+1,current_phase_todo+1,end) = Numbervoxel_phase_tmp/numel(Phase_microstructure_resized); % Defintion of volume fraction
                current_direction_todo = 0;
                for current_direction=1:1:number_dimension
                    if p.direction_todo(current_direction)
                        current_direction_todo=current_direction_todo+1;
                        % % Algorithm
                        % Call tortuosity algorithm (cf. Tau factor manual)
                        if current_direction==1
                            Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                            if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
                                Tau_factor_result.Tau_W1.Tau = NaN;
                            end
                            property_voxelsizedependence(current_iteration+1,current_phase_todo+1,current_direction_todo)=Tau_factor_result.Tau_W1.Tau;
                        elseif current_direction==2
                            Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                            if Tau_factor_result.Tau_W2.Tau == Inf
                                Tau_factor_result.Tau_W2.Tau = NaN;
                            end
                            property_voxelsizedependence(current_iteration+1,current_phase_todo+1,current_direction_todo)=Tau_factor_result.Tau_W2.Tau;
                        elseif current_direction==3
                            Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                            if Tau_factor_result.Tau_W3.Tau == Inf
                                Tau_factor_result.Tau_W3.Tau = NaN;
                            end
                            property_voxelsizedependence(current_iteration+1,current_phase_todo+1,current_direction_todo)=Tau_factor_result.Tau_W3.Tau;
                        end
                    end
                end
                % % Time
                timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase)/number_direction_todo toc(time_clock_start_phase)/number_direction_todo]];
            end
        end

        % CPU and stopwatch time - end
        timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume)/number_direction_todo toc(time_clock_start_volume)/number_direction_todo]];

    end
    clear Phase_microstructure_resized;

    % Sort per voxel size
    for current_direction_todo=1:1:number_direction_todo
        property_voxelsizedependence(:,:,current_direction_todo) = sortrows(property_voxelsizedependence(:,:,current_direction_todo),1);  
    end
    property_voxelsizedependence(:,:,current_direction_todo+1) = sortrows(property_voxelsizedependence(:,:,current_direction_todo+1),1);       
    
    %% DEDUCE BRUGGEMAN AND NORMALIZED EFFECTIVE DIFFUSION COEFFICIENT
    eps_voxelsizedependence = property_voxelsizedependence(:,:,end);
    tau_voxelsizedependence = property_voxelsizedependence(:,:,1:end-1);
    Bruggeman_voxelsizedependence = 1 - log(tau_voxelsizedependence)./log(eps_voxelsizedependence);
    NDeff_voxelsizedependence = eps_voxelsizedependence./tau_voxelsizedependence;
    Bruggeman_voxelsizedependence(:,1,:) = tau_voxelsizedependence(:,1,:);
    NDeff_voxelsizedependence(:,1,:) = tau_voxelsizedependence(:,1,:);

    %% EXTRAPOLATION TO 0 nm
    x=property_voxelsizedependence(:,1,1);

    % Tortuosity factor
    tmp = zeros(number_resize+2,number_phase_todo+1,number_direction_todo); % + 0 nm and + initial voxel size
    if strcmp(p.scaling_extrapolation,'Linear')
        interpolation_voxelsize_order_tau=1;
    elseif strcmp(p.scaling_extrapolation,'Quadratic')
        interpolation_voxelsize_order_tau=2;
    elseif strcmp(p.scaling_extrapolation,'Cubic')
        interpolation_voxelsize_order_tau=3;
    end
    max_order = length(p.scaling)-1;
    interpolation_voxelsize_order_tau = min(interpolation_voxelsize_order_tau,max_order);
    for current_phase_todo=1:1:number_phase_todo
        current_direction_todo = 0;
        for current_direction=1:1:number_dimension
            if p.direction_todo(current_direction)
                current_direction_todo=current_direction_todo+1;
                y=tau_voxelsizedependence(:,current_phase_todo+1,current_direction_todo);
                pi = polyfit(x,y,interpolation_voxelsize_order_tau);
                vq = polyval(pi,0);
                tmp(1,current_phase_todo+1,current_direction_todo)=vq;
                interpolation_voxelsize_tau(current_phase_todo,current_direction_todo).pi=pi;
                % For correlation
                [str_direction] = function_remove_emptyandspecialcharacter_string(char(directionname_todo(current_direction_todo)));
                results_correlation(current_phase_todo).(['Tortuosity_' num2str(current_direction) '_' str_direction '_extrapolated']) = vq;
            end
        end
    end
    tmp(2:end,:,:) = tau_voxelsizedependence(:,:,:);
    tau_voxelsizedependence = tmp; clear tmp;

    % Bruggeman exponent
    tmp = zeros(number_resize+2,number_phase_todo+1,number_direction_todo); % + 0 nm and + initial voxel size
    if strcmp(p.Bruggeman_scaling_extrapolation,'Linear')
        interpolation_voxelsize_order_p=1;
    elseif strcmp(p.Bruggeman_scaling_extrapolation,'Quadratic')
        interpolation_voxelsize_order_p=2;
    elseif strcmp(p.Bruggeman_scaling_extrapolation,'Cubic')
        interpolation_voxelsize_order_p=3;
    end
    max_order = length(p.scaling)-1;
    interpolation_voxelsize_order_p = min(interpolation_voxelsize_order_p,max_order);
    for current_phase_todo=1:1:number_phase_todo
        current_direction_todo = 0;
        for current_direction=1:1:number_dimension
            if p.direction_todo(current_direction)
                current_direction_todo=current_direction_todo+1;
                y=Bruggeman_voxelsizedependence(:,current_phase_todo+1,current_direction_todo);
                pi = polyfit(x,y,interpolation_voxelsize_order_p);
                vq = polyval(pi,0);
                tmp(1,current_phase_todo+1,current_direction_todo)=vq;
                interpolation_voxelsize_p(current_phase_todo,current_direction_todo).pi=pi;
                % For correlation
                [str_direction] = function_remove_emptyandspecialcharacter_string(char(directionname_todo(current_direction_todo)));
                results_correlation(current_phase_todo).(['Bruggeman_exponent_' num2str(current_direction) '_' str_direction '_extrapolated']) = vq;
            end
        end
    end
    tmp(2:end,:,:) = Bruggeman_voxelsizedependence(:,:,:);
    Bruggeman_voxelsizedependence = tmp; clear tmp;

    % Normalized diffusivity
    tmp = zeros(number_resize+2,number_phase_todo+1,number_direction_todo); % + 0 nm and + initial voxel size
    if strcmp(p.Deff_scaling_extrapolation,'Linear')
        interpolation_voxelsize_order_deff=1;
    elseif strcmp(p.Deff_scaling_extrapolation,'Quadratic')
        interpolation_voxelsize_order_deff=2;
    elseif strcmp(p.Deff_scaling_extrapolation,'Cubic')
        interpolation_voxelsize_order_deff=3;
    end
    max_order = length(p.scaling)-1;
    interpolation_voxelsize_order_deff = min(interpolation_voxelsize_order_deff,max_order);
    for current_phase_todo=1:1:number_phase_todo
        current_direction_todo = 0;
        for current_direction=1:1:number_dimension
            if p.direction_todo(current_direction)
                current_direction_todo=current_direction_todo+1;
                y=NDeff_voxelsizedependence(:,current_phase_todo+1,current_direction_todo);
                pi = polyfit(x,y,interpolation_voxelsize_order_deff);
                vq = polyval(pi,0);
                tmp(1,current_phase_todo+1,current_direction_todo)=vq;
                interpolation_voxelsize_deff(current_phase_todo,current_direction_todo).pi=pi;
                % For correlation
                [str_direction] = function_remove_emptyandspecialcharacter_string(char(directionname_todo(current_direction_todo)));
                results_correlation(current_phase_todo).(['Normalized_diffusivity_' num2str(current_direction) '_' str_direction '_extrapolated']) = vq;
            end
        end
    end
    tmp(2:end,:,:) = NDeff_voxelsizedependence(:,:,:);
    NDeff_voxelsizedependence = tmp; clear tmp;    
     

    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={['Voxel size ' voxel_unit]}; % Columns name
    for current_direction_todo=1:1:number_direction_todo
        Variable_name_table(1+current_direction_todo)=directionname_todo(current_direction_todo);
    end
    % Table    
    for current_phase_todo = 1:1:number_phase_todo
        t=tau_voxelsizedependence(:,1,1);
        for current_direction_todo = 1:1:number_direction_todo
            t = [t tau_voxelsizedependence(:,current_phase_todo+1,current_direction_todo)];
        end
        Tau_voxelsizedependence_phase(current_phase_todo).table = array2table(t,'VariableNames',Variable_name_table);
        t=Bruggeman_voxelsizedependence(:,1,1);
        for current_direction_todo = 1:1:number_direction_todo
            t = [t Bruggeman_voxelsizedependence(:,current_phase_todo+1,current_direction_todo)];
        end
        Burggeman_voxelsizedependence_phase(current_phase_todo).table = array2table(t,'VariableNames',Variable_name_table);
        t=NDeff_voxelsizedependence(:,1,1);
        for current_direction_todo = 1:1:number_direction_todo
            t = [t NDeff_voxelsizedependence(:,current_phase_todo+1,current_direction_todo)];
        end
        NDeff_voxelsizedependence_phase(current_phase_todo).table = array2table(t,'VariableNames',Variable_name_table);      
    end    

    %% DISPLAY TEXT RESULTS
    for current_phase_todo = 1:1:number_phase_todo
        fprintf('> For %s, tortuosity factor, extrapolation to zero voxel size: polynomial of order %i\n\n',char(phasename_todo(current_phase_todo)),interpolation_voxelsize_order_tau);
        disp(Tau_voxelsizedependence_phase(current_phase_todo).table)
        fprintf('> For %s, Bruggeman exponent, extrapolation to zero voxel size: polynomial of order %i\n\n',char(phasename_todo(current_phase_todo)),interpolation_voxelsize_order_p);
        disp(Burggeman_voxelsizedependence_phase(current_phase_todo).table)
        fprintf('> For %s, normalized effective diffusivity, extrapolation to zero voxel size: polynomial of order %i\n\n',char(phasename_todo(current_phase_todo)),interpolation_voxelsize_order_deff);
        disp(NDeff_voxelsizedependence_phase(current_phase_todo).table)  
    end
    
    %% SAVE RESULTS
    Results_Tortuosity.Tau_voxelsizedependence_phase = Tau_voxelsizedependence_phase; % Save in main table result
    Results_Tortuosity.Burggeman_voxelsizedependence_phase = Burggeman_voxelsizedependence_phase;
    Results_Tortuosity.NDeff_voxelsizedependence_phase = NDeff_voxelsizedependence_phase;
    if opt.save.xls
        filename = 'Transport_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        sheet=1;
        for current_phase_todo = 1:1:number_phase_todo
            DATA_writetable.sheet(sheet).name=['Tortuosity ' char(phasename_todo(current_phase_todo))];
            DATA_writetable.sheet(sheet).table=Tau_voxelsizedependence_phase(current_phase_todo).table;
            DATA_writetable.sheet(sheet+1).name=['Bruggeman ' char(phasename_todo(current_phase_todo))];
            DATA_writetable.sheet(sheet+1).table=Burggeman_voxelsizedependence_phase(current_phase_todo).table;
            DATA_writetable.sheet(sheet+2).name=['Diffusivity ' char(phasename_todo(current_phase_todo))];
            DATA_writetable.sheet(sheet+2).table=NDeff_voxelsizedependence_phase(current_phase_todo).table;      
            sheet=sheet+3;
        end   
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    %% FIGURES
    infovolbis = infovol;
    infovolbis.phasecolor = colororder;
    infovolbis.phasename = infovol.directionname;
    for current_phase_todo=1:1:number_phase_todo
        parameters_figure.plotlog = false;
        parameters_figure.propertyname = 'Tortuosity factor';
        parameters_figure.propertyname2 = 'Bruggeman exponent';
        parameters_figure.propertyname3 = 'Normalized effective diffusion coefficient';
        parameters_figure.method = 'TauFactor (finite difference, Dirichelt BCs)';
        parameters_figure.figname = ['Transport properties, ' char(phasename_todo(current_phase_todo))];
        parameters_figure.property_voxelsizedependence = Tau_voxelsizedependence_phase(current_phase_todo).table.Variables;
        parameters_figure.property_voxelsizedependence2 = Burggeman_voxelsizedependence_phase(current_phase_todo).table.Variables;
        parameters_figure.property_voxelsizedependence3 = NDeff_voxelsizedependence_phase(current_phase_todo).table.Variables;
        parameters_figure.number_phase = number_direction_todo;
        parameters_figure.str_ylabel = '\tau';
        parameters_figure.str_ylabel2 = 'Bruggeman exponent p';
        parameters_figure.str_ylabel3 = 'D_{eff}/D_{bulk}';
        parameters_figure.propertynameunit = [];
        parameters_figure.propertynameunit2 = [];
        parameters_figure.propertynameunit3 = [];
        parameters_figure.interpolation_voxelsize.pi = []; parameters_figure.interpolation_voxelsize2.pi=[]; parameters_figure.interpolation_voxelsize3.pi=[];
        for current_direction_todo=1:1:number_direction_todo
            parameters_figure.interpolation_voxelsize.pi = [parameters_figure.interpolation_voxelsize.pi; interpolation_voxelsize_tau(current_phase_todo,current_direction_todo).pi];
            parameters_figure.interpolation_voxelsize2.pi = [parameters_figure.interpolation_voxelsize2.pi; interpolation_voxelsize_p(current_phase_todo,current_direction_todo).pi];
            parameters_figure.interpolation_voxelsize3.pi = [parameters_figure.interpolation_voxelsize3.pi; interpolation_voxelsize_deff(current_phase_todo,current_direction_todo).pi];
        end       
        parameters_figure.Current_folder = Current_folder;
        parameters_figure.filename = ['Transport_properties_' char(phasename_todo(current_phase_todo)) '_voxel_size_dependence'];
        parameters_figure.infovol = infovolbis;
        parameters_figure.opt = opt;
        parameters_figure.todo = p.direction_todo;
        Function_create_figure_voxelsizedependence(parameters_figure) % Figures
    end 
end

%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%

if p.RVE.number_RVE>0
    n_property_RVE=0;
    current_direction_todo = 0;
    for current_direction=1:1:number_dimension
        if p.direction_todo(current_direction)
            current_direction_todo=current_direction_todo+1;
            str_property(n_property_RVE+1).propertyname = ['Tortuosity factor, ' char(infovol.directionname(current_direction))];
            str_property(n_property_RVE+2).propertyname = ['Bruggeman exponent, ' char(infovol.directionname(current_direction))];
            str_property(n_property_RVE+3).propertyname = ['Normalized diffusivity, ' char(infovol.directionname(current_direction))];
            str_direction = function_remove_emptyandspecialcharacter_string(char(directionname_todo(current_direction_todo)));
            str_property(n_property_RVE+1).corrname = ['Tau_' num2str(current_direction) '_' str_direction];
            str_property(n_property_RVE+2).corrname = ['P_' num2str(current_direction) '_' str_direction];
            str_property(n_property_RVE+3).corrname = ['NDeff' num2str(current_direction) '_' str_direction];

            res(n_property_RVE+1).Wholevolume_results = Tortuosity_factor(:,current_direction_todo);
            res(n_property_RVE+2).Wholevolume_results = Bruggeman_exponent(:,current_direction_todo);
            res(n_property_RVE+3).Wholevolume_results = Normalized_effective_diffusion_coefficient(:,current_direction_todo);

            n_property_RVE=n_property_RVE+3;
        end
    end

    for k_RVE = 1:1:p.RVE.number_RVE % Loop over all RVE asked
        RVEparameters = p.RVE.RVE(k_RVE);
        for property_RVE=1:1:n_property_RVE
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
            Result_nestedRVE = zeros(n_nestedRVE+1,n_threshold+1,number_phase_todo,3,3, n_property_RVE); % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdoamin) in cubic root (=1), in square root (=2), or lenght (=3)
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
            for property_RVE=1:1:n_property_RVE
                Rk(property_RVE).RVE(k_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),All_subdomain(:,12),...
                    'VariableNames',{'Subdomain Id' 'Group Id' 'Number subdomain' 'Equivalent cubic length' 'Equivalent section length' 'Length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
            end
            [number_subdomain,~] = size(All_subdomain); % The number of subdomain
            number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size

            %% ALGORITHM
            % Colunm 1 is the subdomain id
            % Colunm 2 and 3 are the sizes of the subdomain.
            Property_eachsubdomain = zeros(number_subdomain,number_phase_todo+4, n_property_RVE);
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
                        binary_phase(current_subdomain == code_tmp) = 1; % Binary phase
                        Volumefraction_subdomain = Numbervoxel_phase_tmp/voxel_number_tmp;

                        current_direction_todo = 0;
                        pr=0;
                        for current_direction=1:1:number_dimension
                            if p.direction_todo(current_direction)
                                current_direction_todo=current_direction_todo+1;
                                % % Algorithm
                                % Call tortuosity algorithm (cf. Tau factor manual)
                                if current_direction==1
                                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                                    if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
                                        Tau_factor_result.Tau_W1.Tau = NaN;
                                    end
                                    Property_eachsubdomain(subdomain_id,current_phase_todo+4,pr+1)=Tau_factor_result.Tau_W1.Tau;                                   
                                elseif current_direction==2
                                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                                    if Tau_factor_result.Tau_W2.Tau == Inf
                                        Tau_factor_result.Tau_W2.Tau = NaN;
                                    end
                                    Property_eachsubdomain(subdomain_id,current_phase_todo+4,pr+1)=Tau_factor_result.Tau_W2.Tau;
                                elseif current_direction==3
                                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                                    if Tau_factor_result.Tau_W3.Tau == Inf
                                        Tau_factor_result.Tau_W3.Tau = NaN;
                                    end
                                    Property_eachsubdomain(subdomain_id,current_phase_todo+4,pr+1)=Tau_factor_result.Tau_W3.Tau;                                    
                                end
                                Property_eachsubdomain(subdomain_id,current_phase_todo+4,pr+2)= 1 - log(Property_eachsubdomain(subdomain_id,current_phase_todo+4,pr+1))/log(Volumefraction_subdomain);
                                Property_eachsubdomain(subdomain_id,current_phase_todo+4,pr+3)= Volumefraction_subdomain/Property_eachsubdomain(subdomain_id,current_phase_todo+4,pr+1);
                                pr = pr + 3;
                            end
                        end
                        % % Time
                        timedata_perphase = [timedata_perphase; [Numbervoxel_phase_tmp (cputime-time_cpu_start_phase)/number_direction_todo toc(time_clock_start_phase)/number_direction_todo]];

                    end
                end
                % CPU and stopwatch time - end
                timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume)/number_direction_todo toc(time_clock_start_volume)/number_direction_todo]];
            end

            for property_RVE=1:1:n_property_RVE
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
                    parameters_figure.propertynameunit = [];
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
            for property_RVE=1:1:n_property_RVE
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

    for property_RVE=1:1:n_property_RVE
        Results_Tortuosity.RVE.(str_property(property_RVE).corrname) = Rk(property_RVE).RVE; % Save in main table result
    end

end


%%
%% ENDING FUNCTION
%%

%% TIME
Table_time_pervolume = table(timedata_pervolume(:,1),timedata_pervolume(:,2),timedata_pervolume(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Tortuosity.Table_time_pervolume = Table_time_pervolume; % Save in main table result
Table_time_perphase = table(timedata_perphase(:,1),timedata_perphase(:,2),timedata_perphase(:,3),...
    'VariableNames',{'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Tortuosity.Table_time_perphase = Table_time_perphase; % Save in main table result

date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
lasted_time = date_end-date_start;
Table_date = table({char(date_start)},{char(date_end)},{char(lasted_time)},...
    'VariableNames',{'Start date' 'End date' 'Lasted time'});
Results_Tortuosity.Table_date = Table_date;

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
    Function_Writetable(Current_folder,'Tortuosity_calculation_time',DATA_writetable)
end
% Display
fprintf ('Finished the %s\n\n',date_end);
fprintf ('Lasted: %s\n\n',lasted_time);
function_time_figure(timedata_pervolume, timedata_perphase, Current_folder, 'Tortuosity_calculation_time', 'Tortuosity factor', opt);

%% SAVE CORRELATION
Current_folder = [infovol.volpath 'Correlation' separator];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end
save([Current_folder 'Correlation_tortuosity'],'results_correlation');

%% SAVE RESULTS
if opt.save.mat
    Current_folder = [infovol.volpath 'Summary' separator];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_tortuosity'],'Results_Tortuosity')
end

end