function [] = Function_Tortuosity_factor_taufactor(Phase_microstructure, PROPERTY, OPTIONS, INFO)
%Calculate tortuosity factor using tau factor
% Function_Tortuosity_factor_taufactor(array, PROPERTY, OPTIONS, INFO) - when use with the toolbox
% or
% Function_Tortuosity_factor_taufactor(array, voxelsize) - when use as a standalone function

% The discrete computation of the diffusion equation comes from:
% S.J. Cooper, A. Bertei, P.R. Shearing, J.A. Kilner, M.P. Brandon
% TauFactor: an Open-source application for calculating tortuosity factors from tomographic data
% SoftwareX 5 (2016) 203-2010

%% DEFAULT VALUES
expected_number_argument = 4;
nargin; % Number of input variable when the function is call

if nargin ~= expected_number_argument % Unexpected number of argument
    
    if nargin == 2 % Case for function called as: Function_particle_size_CPSD(Phase_microstructure, voxel_size)
        voxel_size = PROPERTY;
        clear PROPERTY OPTIONS
        
        % Set default folder
        t = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss'); % Set unique folder based on time, with second precision
        INFO.resultfoldername = ['Volume_characterization_' char(t)];
        desktop=winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); % Find desktop folder of windows user
        OPTIONS.mainsavefolder = [desktop '\' INFO.resultfoldername];
        allcolor=[get(0, 'DefaultAxesColorOrder'); rand(100,3);];
        
        % Set default phase information
        unique_code = unique(Phase_microstructure);
        INFO.number_phase = length(unique_code);
        INFO.voxel_number = numel(Phase_microstructure);
        for k=1:1:INFO.number_phase % Loop over all unique values
            INFO.phase(k).code = unique_code(k); % Assing code
            INFO.phasename(k,1) = {['Phase ' num2str(INFO.phase(k).code)]}; % Assign phase name based on code value
            INFO.phase(k).filename = ['Phase_' num2str(INFO.phase(k).code)]; % Assign phase name based on code value
            INFO.phase(k).name = INFO.phase(k).filename;
            INFO.phase(k).color = allcolor(k,:);
        end
        INFO.voxel_size = voxel_size; % nanometers
        
        % Set default direction information
        INFO.direction(1).name = 'in-plane direction 1';
        INFO.direction(2).name = 'in-plane direction 2';
        INFO.direction(3).name = 'through-plane direction';
        INFO.direction(1).filename = 'inplane_direction_1';
        INFO.direction(2).filename = 'inplane_direction_2';
        INFO.direction(3).filename = 'throughplane_direction';
        
        % Set default options
        OPTIONS.save_resultsmat = true;
        OPTIONS.save_xls = true;
        OPTIONS.save_fig = true;
        OPTIONS.savefig_infig = true;
        OPTIONS.savefig_informat = {'png'};
        
        % Set display options
        OPTIONS.fontname = 'Times New Roman';
        OPTIONS.displaytext = true;
        OPTIONS.closefigureaftercreation = false;
        OPTIONS.grid = 'on'; OPTIONS.minorgrid = 'on';
        OPTIONS.Linewidth = 2;
        OPTIONS.Fontsize_axe =  12;
        OPTIONS.Fontsize_legend =  12;
        OPTIONS.Fontsize_title =  14;
        
        % No Voxel size dependence analysis and RVE analysis
        PROPERTY.tortuosity_taufactor.voxel_size_dependence.todo = false;
        PROPERTY.tortuosity_taufactor.number_RVE = 0;
        PROPERTY.tortuosity_taufactor.todo_slice = [ {false}; {false}; {false}];
        PROPERTY.tortuosity_taufactor.sliceparameter = [{5};{5};{5}];
        PROPERTY.tortuosity_taufactor.slicechoice = 1;
        
    else % Incorrect number of argument
        disp 'Error calling Function_Tortuosity_factor_taufactor. Wrong number of argument.'
        help Function_Tortuosity_factor_taufactor
    end
    
end


%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
if ispc
    main_folder = [OPTIONS.mainsavefolder '\'];
    Sub_folder = 'Tortuosity_factor_taufactor\'; % Result are saved in this subfolder
else
    main_folder = [OPTIONS.mainsavefolder '/'];
    Sub_folder = 'Tortuosity_factor_taufactor/'; % Result are saved in this subfolder
end
Current_folder = [main_folder Sub_folder];
if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% INFO
Domain_size = size(Phase_microstructure);
if length(Domain_size)==2 % Add third column with unit value
    Domain_size = [Domain_size 1];
end
if min(Domain_size) == 1 % 2D case
    number_dimension = 2;
else
    number_dimension =3; % 3D case
end
number_phase = INFO.number_phase; % Number of phase
voxel_number = INFO.voxel_number; % Number of voxel
voxel_size = INFO.voxel_size; % Voxel size nm

%% INITIALIZE RESULTS (USE FOR CORRELATION and VISUALIZATION)
for current_phase=1:1:number_phase
    results_correlation(current_phase).name = INFO.phase(current_phase).name;
end

%% DISPLAY
if OPTIONS.displaytext==true
    disp '    TORTUOSITY FACTOR - TAUFACTOR (FINITE DIFFERENCE, Dirichlet-Dirichelt BC) METHOD';
    disp '    --------------------------------------------------------------------------------';
    disp ' ';
end

%% PARAMETERS
%Analysed_phase = PROPERTY.tortuosityfactor.taufactor.phase_to_analyse;
Analysed_phase=0:1:50;

%%
%% ALGORITHM ON WHOLE VOLUME
%%

%% CALCULATION
% Initialisation
Effective_diffusion_coefficient = zeros(number_phase,number_dimension);
Bruggeman_exponent = zeros(number_phase,number_dimension);
Tortuosity_factor = zeros(number_phase,number_dimension);
Effective_diffusion_coefficient_anisotropy = zeros(number_phase,number_dimension);
Bruggeman_exponent_anisotropy = zeros(number_phase,number_dimension);
Tortuosity_factor_anisotropy = zeros(number_phase,number_dimension);
volume_fraction_total = zeros(number_phase,number_dimension);

n_calculation = 0;
for current_phase=1:1:number_phase % Loop on every phase
    code_tmp = INFO.phase(current_phase).code; % The code of the phase
    binary_phase = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
    binary_phase( Phase_microstructure==code_tmp )=1;
    voxel_number = sum(sum(sum(binary_phase==1)));
    
    % Volume fraction
    volume_fraction_total(current_phase,:)=voxel_number/prod(Domain_size);
    
    for current_direction=1:1:number_dimension % Loop on every analysed direction
        n_calculation = n_calculation+1;
        time_cpu_start = cputime; % CPU time start
        tic; % Stopwatch start
        % Call tortuosity algorithm (cf. Tau factor manual)
        if current_direction==1
            Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
            if Tau_factor_result.Tau_W1.Tau == Inf % No percolation path. 65535 is default value for inf tortuosity factor
                Tau_factor_result.Tau_W1.Tau = NaN;
            end
            Tortuosity_factor(current_phase,current_direction)=Tau_factor_result.Tau_W1.Tau;
        elseif current_direction==2
            Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
            if Tau_factor_result.Tau_W2.Tau == Inf
                Tau_factor_result.Tau_W2.Tau = NaN;
            end            
            Tortuosity_factor(current_phase,current_direction)=Tau_factor_result.Tau_W2.Tau;
        elseif current_direction==3
            Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
            if Tau_factor_result.Tau_W3.Tau == Inf
                Tau_factor_result.Tau_W3.Tau = NaN;
            end            
            Tortuosity_factor(current_phase,current_direction)=Tau_factor_result.Tau_W3.Tau;
        end
        % Effective diffusion coefficient
        Effective_diffusion_coefficient(current_phase,current_direction) = 1*volume_fraction_total(current_phase,1)/Tortuosity_factor(current_phase,current_direction);
        tau = Tortuosity_factor(current_phase,current_direction);
        eps = volume_fraction_total(current_phase,1);
        Bruggeman_exponent(current_phase,current_direction) = 1 - log(tau)/log(eps);
        % CPU and stopwatch time - end
        time_cpu_elapsed = cputime-time_cpu_start;
        time_stopwatch_elapsed = toc;
        if n_calculation==1
            Time_measure = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
        else
            Time_tmp = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
            Time_measure = [Time_measure;Time_tmp];
        end
        
        % Save (correlation)
        [str_direction] = function_remove_emptyandspecialcharacter_string(INFO.direction(current_direction).name);
        results_correlation(current_phase).(['Tortuosity_' num2str(current_direction) '_' str_direction]) = Tortuosity_factor(current_phase,current_direction);
        results_correlation(current_phase).(['Effective_diffusion_coefficient_' num2str(current_direction) '_' str_direction]) = Effective_diffusion_coefficient(current_phase,current_direction);
        results_correlation(current_phase).(['Bruggeman_exponent_' num2str(current_direction) '_' str_direction]) = Bruggeman_exponent(current_phase,current_direction);
        
    end
    Effective_diffusion_coefficient_anisotropy(current_phase,:) = Effective_diffusion_coefficient(current_phase,:) ./  Effective_diffusion_coefficient(current_phase,end);
    Bruggeman_exponent_anisotropy(current_phase,:) = Bruggeman_exponent(current_phase,:) ./  Bruggeman_exponent(current_phase,end);
    Tortuosity_factor_anisotropy(current_phase,:) = Tortuosity_factor(current_phase,:) ./  Tortuosity_factor(current_phase,end);
    
    % Save (correlation)
    results_correlation(current_phase).Effective_diffusion_coefficient_anisotropy_1_over_3 = Effective_diffusion_coefficient_anisotropy(current_phase,1);
    results_correlation(current_phase).Effective_diffusion_coefficient_anisotropy_2_over_3 = Effective_diffusion_coefficient_anisotropy(current_phase,2);
    results_correlation(current_phase).Bruggeman_exponent_anisotropy_1_over_3 = Bruggeman_exponent_anisotropy(current_phase,1);
    results_correlation(current_phase).Bruggeman_exponent_anisotropy_2_over_3 = Bruggeman_exponent_anisotropy(current_phase,2);
    results_correlation(current_phase).Tortuosity_anisotropy_1_over_3 = Tortuosity_factor_anisotropy(current_phase,1);
    results_correlation(current_phase).Tortuosity_anisotropy_2_over_3 = Tortuosity_factor_anisotropy(current_phase,2);
end

%% TABLES
% Time
Table_time = table(Time_measure(:,1)*1e-6,Time_measure(:,2),Time_measure(:,3),...
    'VariableNames',{'Voxel_number_millions','CPU_time_s' 'Stopwatch_s'});
Results_tortuosityfactor.Table_time = Table_time; % Save in main table result

% Result calculated on whole volume
for current_phase=1:1:number_phase
    Tortuositytaufactor.phase(current_phase).table = table([{INFO.direction(1).name}; {INFO.direction(2).name}; {INFO.direction(3).name}],volume_fraction_total(current_phase,:)',Tortuosity_factor(current_phase,:)',Bruggeman_exponent(current_phase,:)',Effective_diffusion_coefficient(current_phase,:)',...
        'VariableNames',{'Direction','Volume_fraction','Tortuosity_factor','Bruggeman_exponent','Deff_Dbulk',});
    str_anisotropy_13 = {[INFO.direction(1).name ' / ' INFO.direction(3).name]};
    str_anisotropy_23 = {[INFO.direction(2).name ' / ' INFO.direction(3).name]};
    Tortuositytaufactor_anisotropy.phase(current_phase).table = table([str_anisotropy_13; str_anisotropy_23],Tortuosity_factor_anisotropy(current_phase,1:end-1)',Bruggeman_exponent_anisotropy(current_phase,1:end-1)',Effective_diffusion_coefficient_anisotropy(current_phase,1:end-1)',...
        'VariableNames',{'Direction','Tortuosity_factor_anisotropy','Bruggeman_exponent_anisotropy','Deff_Dbulk_anisotropy',});
end
Results_tortuosityfactor.Table_TortuosityFactor_taufactor = Tortuositytaufactor; % Save in main table result
Results_tortuosityfactor.Table_TortuosityFactorAnisotropy_taufactor = Tortuositytaufactor_anisotropy;

%% SAVE TABLES
if OPTIONS.save_xls==true
    filename = 'Tortuosity_factor'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    sheet_=0;
    for current_phase=1:1:number_phase
        sheet_=sheet_+1;
        DATA_writetable.sheet(sheet_).name=INFO.phase(current_phase).filename;
        DATA_writetable.sheet(sheet_).table=Tortuositytaufactor.phase(current_phase).table;
        sheet_=sheet_+1;
        DATA_writetable.sheet(sheet_).name=[INFO.phase(current_phase).filename ' anisotropy'];
        DATA_writetable.sheet(sheet_).table=Tortuositytaufactor_anisotropy.phase(current_phase).table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
if OPTIONS.displaytext==true
    fprintf('> Calculated on the whole domain:\n\n');
    disp 'Tortuosity factor';
    disp ' '
    for current_phase=1:1:number_phase
        fprintf('For the phase %s\n',INFO.phase(current_phase).name);
        disp(Tortuositytaufactor.phase(current_phase).table)
        disp(Tortuositytaufactor_anisotropy.phase(current_phase).table)
    end
    fprintf('Computation time, in seconds:\n\n');
    disp(Table_time)
end

%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME
%%
scrsz = get(0,'ScreenSize'); % Screen resolution

%% FIGURES
Fig = figure; % Create figure
Fig.Name= 'Transport properties';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',scrsz); % Full screen figure
% - Create axes as a subplot
for id_axe=1:1:8
    if id_axe ~=5
        sub_axes=subplot(2,4,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            h_title=title ('Volume fraction \epsilon'); % Title
            h_=bar(volume_fraction_total,1); % Width=1
            ylabel('Volume fraction \epsilon'); % Axis label
        elseif id_axe==2
            h_title=title ('Tortuosity factor \tau'); % Title
            h_=bar(Tortuosity_factor,1); % Width=1
            ylabel('Tortuosity factor \tau'); % Axis label
        elseif id_axe==3
            h_title=title ('Bruggeman exponent p, \tau=\epsilon^{1-p}'); % Title
            h_=bar(Bruggeman_exponent,1); % Width=1
            ylabel('Bruggeman exponent p'); % Axis label
        elseif id_axe==4
            h_title=title ('Diffusion coefficient ratio D_{eff}/D_{bulk}'); % Title
            h_=bar(Effective_diffusion_coefficient,1); % Width=1
            ylabel('D_{eff}/D_{bulk}'); % Axis label
        elseif id_axe==6
            h_title=title ('Tortuosity factor anisotropy'); % Title
            h_=bar(Tortuosity_factor_anisotropy(:,1:end-1),1); % Width=1
            ylabel('Tortuosity factor anisotropy'); % Axis label
        elseif id_axe==7
            h_title=title ('Bruggeman exponent anisotropy'); % Title
            h_=bar(Bruggeman_exponent_anisotropy(:,1:end-1),1); % Width=1
            ylabel('Bruggeman exponent anisotropy'); % Axis label
        elseif id_axe==8
            h_title=title ('Diffusion coefficient ratio anisotropy'); % Title
            h_=bar(Effective_diffusion_coefficient_anisotropy(:,1:end-1),1); % Width=1
            ylabel('D_{eff}/D_{bulk} anisotropy'); % Axis label
        end
        clear x_text
        for current_phase=1:1:number_phase
            if current_phase==1
                phase_text={'',INFO.phase(current_phase).name,''};
                x_text = phase_text;
            else
                phase_text={INFO.phase(current_phase).name,''};
                x_text=[x_text,phase_text];
            end
        end
        set(sub_axes,'xticklabel',x_text);
        % Set individual color and legend
        if id_axe<=4
            h_(1).FaceColor = [0 114 189]/255;
            h_(2).FaceColor = [217 83 25]/255;
            h_(3).FaceColor = [237 177 32]/255;
            h_legend = legend(sub_axes,INFO.direction.name,'Location','best'); % Legend
        else
            h_(1).FaceColor = [0 114 189]/255;
            h_(2).FaceColor = [217 83 25]/255;
            str_legend(1).name = char(str_anisotropy_13);
            str_legend(2).name = char(str_anisotropy_23);
            h_legend = legend(sub_axes,str_legend.name,'Location','best');
        end
        set(sub_axes,'YGrid',OPTIONS.grid); % Grid
        set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
        h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
        h_legend.FontSize = OPTIONS.Fontsize_legend; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
end
sgtitle(Fig,'Transport properties','FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
if OPTIONS.save_fig == true % Save figure
    function_savefig(Fig, Current_folder, 'Transport_properties', OPTIONS); % Call function
end
if OPTIONS.closefigureaftercreation == true
    close(Fig); % Do not keep open figures
end

%% ALONG DIRECTIONS
for current_direction=1:1:number_dimension
    if logical(cell2mat( PROPERTY.tortuosity_taufactor.todo_slice(current_direction) ))
        if PROPERTY.tortuosity_taufactor.slicechoice==1
            number_thickslice = round(cell2mat(PROPERTY.tortuosity_taufactor.sliceparameter(current_direction)));
        else
            thickness_slice = cell2mat(PROPERTY.tortuosity_taufactor.sliceparameter(current_direction));
            number_thickslice = round(Domain_size(current_direction)*voxel_size/1000 / thickness_slice);
        end        
        number_thickslice = max([1 number_thickslice]);
        position_slice = round(linspace(1,Domain_size(current_direction),number_thickslice+1));
        for current_phase=1:1:number_phase
            code_tmp = INFO.phase(current_phase).code; % The code of the phase
            binary_phase = zeros(Domain_size(1),Domain_size(2),Domain_size(3));
            binary_phase( Phase_microstructure==code_tmp )=1;
            
            val_ = zeros(number_thickslice,6); % Start position, end position, porosity,  tortuosity factor, Bruggeman exponent, Deff_Dbulk
            for current_slice=1:1:number_thickslice
                if current_direction==1
                    x0 = position_slice(current_slice); x1 = position_slice(current_slice+1);
                    y0 = 1; y1 = Domain_size(2);
                    z0 = 1; z1 = Domain_size(3);
                elseif current_direction==2
                    x0 = 1; x1 = Domain_size(1);
                    y0 = position_slice(current_slice); y1 = position_slice(current_slice+1);
                    z0 = 1; z1 = Domain_size(3);
                elseif current_direction==3
                    x0 = 1; x1 = Domain_size(1);
                    y0 = 1; y1 = Domain_size(2);
                    z0 = position_slice(current_slice); z1 = position_slice(current_slice+1);
                end
                val_(current_slice,1) = position_slice(current_slice);
                val_(current_slice,2) = position_slice(current_slice+1);
                thick_slice = binary_phase(x0:x1,y0:y1,z0:z1);
                
                voxel_number = sum(sum(sum(thick_slice==1)));
                val_(current_slice,3) = voxel_number/numel(thick_slice);
                
                time_cpu_start = cputime; % CPU time start
                tic; % Stopwatch start
                if current_direction==1
                    Tau_factor_result = TauFactor('InLine',1,0,0,thick_slice,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                    if Tau_factor_result.Tau_W1.Tau == Inf
                        Tau_factor_result.Tau_W1.Tau = NaN;
                    end
                    val_(current_slice,4)=Tau_factor_result.Tau_W1.Tau;
                elseif current_direction==2
                    Tau_factor_result = TauFactor('InLine',1,0,0,thick_slice,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                    if Tau_factor_result.Tau_W2.Tau == Inf
                        Tau_factor_result.Tau_W2.Tau = NaN;
                    end
                    val_(current_slice,4)=Tau_factor_result.Tau_W2.Tau;
                elseif current_direction==3
                    Tau_factor_result = TauFactor('InLine',1,0,0,thick_slice,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                    if Tau_factor_result.Tau_W3.Tau == Inf
                        Tau_factor_result.Tau_W3.Tau = NaN;
                    end
                    val_(current_slice,4)=Tau_factor_result.Tau_W3.Tau;
                end
                % Effective diffusion coefficient
                val_(current_slice,6) = 1*val_(current_slice,3)/val_(current_slice,4);
                tau = val_(current_slice,4);
                eps = val_(current_slice,3);
                val_(current_slice,5) = 1 - log(tau)/log(eps);
                % CPU and stopwatch time - end
                time_cpu_elapsed = cputime-time_cpu_start;
                time_stopwatch_elapsed = toc;
                Time_tmp = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
                Time_measure = [Time_measure;Time_tmp];
            end
            % Convert to micrometers
            val_(:,1) = val_(:,1)*voxel_size/1000;
            val_(:,2) = val_(:,2)*voxel_size/1000;
            val_(1,1)=0;
            % Table
            tmp_ = table(val_(:,1),val_(:,2),val_(:,3),val_(:,4),val_(:,5),val_(:,6),...
                'VariableNames',{'Start' 'End' 'Porosity' 'Tortuosity_factor' 'Bruggeman_exponent' 'Deff_Dbulk'});
            Transportslice.direction(current_direction).phase(current_phase).table = tmp_;
        end
    end
end

if logical(cell2mat(PROPERTY.tortuosity_taufactor.todo_slice(1))) || logical(cell2mat(PROPERTY.tortuosity_taufactor.todo_slice(2))) || logical(cell2mat(PROPERTY.tortuosity_taufactor.todo_slice(3)))
    Results_tortuosityfactor.Table_TortuosityFactor_slice = Transportslice; % Save in main table result
    
    % SAVE TABLES
    if OPTIONS.save_xls==true
        filename = 'Transport_slice'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        sheet_=0;
        for current_phase=1:1:number_phase
            for current_direction=1:1:number_dimension
                if logical(cell2mat( PROPERTY.tortuosity_taufactor.todo_slice(current_direction) ))
                    sheet_=sheet_+1;
                    DATA_writetable.sheet(sheet_).name=[INFO.phase(current_phase).filename ' ' INFO.direction(current_direction).filename];
                    DATA_writetable.sheet(sheet_).table=Results_tortuosityfactor.Table_TortuosityFactor_slice.direction(current_direction).phase(current_phase).table;
                end
            end
        end
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
    
    % FIGURES
    Color_order = get(0, 'DefaultAxesColorOrder');
    for current_phase=1:1:number_phase
        clear str_legend
        Fig = figure; % Create figure
        Fig.Name= ['Transport properties per thick slice, ' INFO.phase(current_phase).name];
        Fig.Color='white'; % Background colour
        scrsz = get(0,'ScreenSize'); % Screen resolution
        set(Fig,'position',scrsz); % Full screen figure
        % - Create axes as a subplot
        for id_axe=1:1:5
            k_legend=0;
            sub_axes=subplot(2,3,id_axe,'Parent',Fig);
            hold(sub_axes,'on'); % Active subplot
            if id_axe==1
                h_title=title ('Volume fraction \epsilon'); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.tortuosity_taufactor.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_tortuosityfactor.Table_TortuosityFactor_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,3);
                            y_(2*k,1)   = t_(k,3);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [volume_fraction_total(current_phase,1) volume_fraction_total(current_phase,1)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(volume_fraction_total(current_phase,1),'%1.3f')];
                    end
                end
                ylabel('Volume fraction \epsilon');
            elseif id_axe==2
                h_title=title ('Tortuosity factor \tau'); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.tortuosity_taufactor.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_tortuosityfactor.Table_TortuosityFactor_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,4);
                            y_(2*k,1)   = t_(k,4);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [Tortuosity_factor(current_phase,current_direction) Tortuosity_factor(current_phase,current_direction)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(Tortuosity_factor(current_phase,current_direction),'%1.3f')];
                    end
                    ylabel('Tortuosity factor \tau');
                end
            elseif id_axe==3
                h_title=title ('Tortuosity - volume fraction relationship'); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.tortuosity_taufactor.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_tortuosityfactor.Table_TortuosityFactor_slice.direction(current_direction).phase(current_phase).table.Variables;
                        eps = t_(:,3);
                        tau = t_(:,4);
                        p_ = polyfit(log(eps),log(tau),1);
                        x_fit = linspace(min(eps),max(eps),1000);
                        y_fit = exp(p_(2)) * x_fit.^p_(1);
                        h_dot=plot(eps,tau);
                        h_fit=plot(x_fit,y_fit);
                        set(h_dot, 'Color',Color_order(current_direction,:),'MarkerSize',OPTIONS.Fontsize_axe,'Marker','o','LineStyle','none','LineWidth',OPTIONS.Linewidth);
                        set(h_fit, 'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',OPTIONS.Linewidth);
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['\tau=' num2str(exp(p_(2)),'%1.3f') '\epsilon^{1-' num2str(1-p_(1),'%1.3f') '}'];
                    end
                    ylabel('Tortuosity factor \tau');
                end
            elseif id_axe==4
                h_title=title ('Bruggeman exponent p'); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.tortuosity_taufactor.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_tortuosityfactor.Table_TortuosityFactor_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,5);
                            y_(2*k,1)   = t_(k,5);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [Bruggeman_exponent(current_phase,current_direction) Bruggeman_exponent(current_phase,current_direction)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(Bruggeman_exponent(current_phase,current_direction),'%1.3f')];
                    end
                    ylabel('Bruggeman exponent p');
                end
            elseif id_axe==5
                h_title=title ('D_{eff}/D_{bulk}'); % Title
                for current_direction=1:1:number_dimension
                    if logical(cell2mat( PROPERTY.tortuosity_taufactor.todo_slice(current_direction) ))
                        k_legend=k_legend+1;
                        t_= Results_tortuosityfactor.Table_TortuosityFactor_slice.direction(current_direction).phase(current_phase).table.Variables;
                        [n_,~]=size(t_);
                        x_=zeros(2*n_,1); y_=zeros(2*n_,1);
                        for k=1:1:n_
                            x_(2*k-1,1) = t_(k,1);
                            x_(2*k,1)   = t_(k,2);
                            y_(2*k-1,1) = t_(k,6);
                            y_(2*k,1)   = t_(k,6);
                        end
                        h_slice=plot(x_,y_);
                        h_wholevolume = plot([0 Domain_size(current_direction)*voxel_size/1000], [Effective_diffusion_coefficient(current_phase,current_direction) Effective_diffusion_coefficient(current_phase,current_direction)]);
                        set(h_slice, 'Color',Color_order(current_direction,:),'LineWidth',(number_dimension-current_direction+1));
                        set(h_wholevolume,'Color',Color_order(current_direction,:),'LineStyle','--','LineWidth',(number_dimension-current_direction+1));
                        str_legend(2*k_legend-1).name = INFO.direction(current_direction).name;
                        str_legend(2*k_legend).name = ['Whole volume, ' num2str(Effective_diffusion_coefficient(current_phase,current_direction),'%1.3f')];
                    end
                    ylabel('D_{eff}/D_{bulk}');
                end
            end
            if id_axe~=3
                xlabel('Position along domain''s direction (\mum)');
            else
                xlabel('Volume fraction');
            end
            h_legend = legend(sub_axes,str_legend.name,'Location','best');
            % - Grid
            if strcmp(OPTIONS.grid,'on')
                grid(sub_axes,'on'); % Display grid
                set(sub_axes,'XMinorGrid',OPTIONS.minorgrid,'YMinorGrid',OPTIONS.minorgrid); % Display grid for minor thicks
            end
            set(sub_axes,'FontName',OPTIONS.fontname,'FontSize',OPTIONS.Fontsize_axe); % Fontname and fontsize
            h_title.FontSize = OPTIONS.Fontsize_title; % Set title fontsize
            h_legend.FontSize = OPTIONS.Fontsize_legend; % Set title fontsize
            hold(sub_axes,'off'); % Relase figure
        end
        sgtitle(Fig,['Transport properties along directions, ' INFO.phase(current_phase).name] ,'FontWeight','bold','FontSize',OPTIONS.Fontsize_title+2,'FontName',OPTIONS.fontname);
        if OPTIONS.save_fig == true % Save figure
            function_savefig(Fig, Current_folder, ['Transport_slice_' INFO.phase(current_phase).filename], OPTIONS); % Call function
        end
        if OPTIONS.closefigureaftercreation == true
            close(Fig); % Do not keep open figures
        end
    end
end

%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%
interpolation_voxelsize_order=1;

if PROPERTY.tortuosity_taufactor.voxel_size_dependence.todo % Check if voxel size analysis is asked
    size_choice = PROPERTY.tortuosity_taufactor.voxel_size_dependence.voxel;
    size_choice = sort(size_choice);
    size_choice(size_choice==1)=[];
    number_resize=length(size_choice); % Number of different voxel size that will be analyzed

    %% CALCULATION
    % Initialization
    for current_direction=1:1:number_dimension
        direction(current_direction).property_voxelsizedependence = zeros(number_resize+1,number_phase+1,4);
        direction(current_direction).property_voxelsizedependence(1,1,:)=voxel_size;
        direction(current_direction).property_voxelsizedependence(1,2:end,1)=Tortuosity_factor(:,current_direction)'; % Tortuosity factor
        direction(current_direction).property_voxelsizedependence(1,2:end,2)=Bruggeman_exponent(:,current_direction)'; % Bruggmeman exponent
        direction(current_direction).property_voxelsizedependence(1,2:end,3)=Effective_diffusion_coefficient(:,current_direction)'; % Deff_Dulk
        direction(current_direction).property_voxelsizedependence(1,2:end,4)=volume_fraction_total(:,current_direction)'; % Volume fraction
    end
        
    % Loop on each voxel size
    for current_iteration=1:1:number_resize
        % New voxel size
        current_voxel_size = size_choice(current_iteration)*voxel_size;
        % Microstructure resized
        [Phase_microstructure_resized] = function_scale_array(Phase_microstructure, voxel_size, current_voxel_size, INFO.phaseinfo);        
        Current_domain_size=size(Phase_microstructure_resized);
        for current_phase=1:1:number_phase
            code_tmp = INFO.phase(current_phase).code; % Code of the phase
            binary_phase=zeros(Current_domain_size(1),Current_domain_size(2),Current_domain_size(3)); % Binary phase
            binary_phase(Phase_microstructure_resized == code_tmp) = 1;
            
            voxel_number = sum(sum(sum(binary_phase==1)));
            eps=voxel_number/numel(binary_phase);
            for current_direction=1:1:number_dimension % Loop on every analysed direction
                time_cpu_start = cputime; % CPU time start
                tic; % Stopwatch start
                direction(current_direction).property_voxelsizedependence(current_iteration+1,1,:)=current_voxel_size;
                % Call tortuosity algorithm (cf. Tau factor manual)
                if current_direction==1
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                    if Tau_factor_result.Tau_W1.Tau == Inf
                        Tau_factor_result.Tau_W1.Tau = NaN;
                    end                    
                    direction(current_direction).property_voxelsizedependence(current_iteration+1,current_phase+1,1)=Tau_factor_result.Tau_W1.Tau;
                elseif current_direction==2
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                    if Tau_factor_result.Tau_W2.Tau == Inf
                        Tau_factor_result.Tau_W2.Tau = NaN;
                    end                    
                    direction(current_direction).property_voxelsizedependence(current_iteration+1,current_phase+1,1)=Tau_factor_result.Tau_W2.Tau;
                elseif current_direction==3
                    Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                    if Tau_factor_result.Tau_W3.Tau == Inf
                        Tau_factor_result.Tau_W3.Tau = NaN;
                    end                    
                    direction(current_direction).property_voxelsizedependence(current_iteration+1,current_phase+1,1)=Tau_factor_result.Tau_W3.Tau;
                end
                tau = direction(current_direction).property_voxelsizedependence(current_iteration+1,current_phase+1,1);
                direction(current_direction).property_voxelsizedependence(current_iteration+1,current_phase+1,2) = 1 - log(tau)/log(eps); 
                direction(current_direction).property_voxelsizedependence(current_iteration+1,current_phase+1,3) = 1*eps/tau;
                direction(current_direction).property_voxelsizedependence(current_iteration+1,current_phase+1,4) = eps;
                % CPU and stopwatch time - end
                time_cpu_elapsed = cputime-time_cpu_start;
                time_stopwatch_elapsed = toc;
                Time_tmp = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
                Time_measure = [Time_measure;Time_tmp];
            end
        end        
    end
    clear Phase_microstructure_resized;
    
    %% EXTRAPOLATION TO 0 nm
    str_correlation(1).name = 'Tortuosity_factor_1_taufactor';
    str_correlation(2).name = 'Tortuosity_factor_2_taufactor';
    str_correlation(3).name = 'Tortuosity_factor_3_taufactor';
    str_correlation(4).name = 'Bruggmeman_exponent_1_taufactor';
    str_correlation(5).name = 'Bruggmeman_exponent_2_taufactor';
    str_correlation(6).name = 'Bruggmeman_exponent_3_taufactor';
    str_correlation(7).name = 'Effective_diffusion_coefficient_1_taufactor';
    str_correlation(8).name = 'Effective_diffusion_coefficient_2_taufactor';
    str_correlation(9).name = 'Effective_diffusion_coefficient_3_taufactor';
    
    direction_tmp(1).property_voxelsizedependence = zeros(number_resize+2,number_phase+1,4);
    direction_tmp(2).property_voxelsizedependence = zeros(number_resize+2,number_phase+1,4);
    direction_tmp(3).property_voxelsizedependence = zeros(number_resize+2,number_phase+1,4);

    for current_phase=1:1:number_phase
        k_=0;
        for property_resized = 1:1:3
            for current_direction=1:1:number_dimension
                k_=k_+1;
                x=direction(current_direction).property_voxelsizedependence(:,1,property_resized);
                y=direction(current_direction).property_voxelsizedependence(:,current_phase+1,property_resized);
                p = polyfit(x,y,interpolation_voxelsize_order);
                vq = polyval(p,0);
                direction_tmp(current_direction).property_voxelsizedependence(1,current_phase+1,property_resized)=vq;
                direction_tmp(current_direction).interpolation_voxelsize(current_phase,property_resized).p=p;
                results_correlation(current_phase).([str_correlation(k_).name '_extrapolated']) = vq; % For correlation
            end
        end
    end
    direction_tmp(1).property_voxelsizedependence(2:end,:,:) = direction(1).property_voxelsizedependence(:,:,:);    
    direction_tmp(2).property_voxelsizedependence(2:end,:,:) = direction(2).property_voxelsizedependence(:,:,:);  
    direction_tmp(3).property_voxelsizedependence(2:end,:,:) = direction(3).property_voxelsizedependence(:,:,:);      
    direction = direction_tmp; clear direction_tmp;
            
    %% MANAGING RESULTS
    % Results are saved in a table
    Variable_name_table={'Voxel_size_nm'}; % Columns name
    k_=0;
    for current_direction=1:1:number_dimension
        for current_phase=1:1:number_phase
            k_=k_+1;
            Variable_name_table(1+k_)={[INFO.phase(current_phase).filename '_' INFO.direction(current_direction).filename]};
        end
    end
    % Table
    Table_tortuosity_taufactor_voxelsizedependence = array2table([direction(1).property_voxelsizedependence(:,:,1) direction(2).property_voxelsizedependence(:,2:end,1) direction(3).property_voxelsizedependence(:,2:end,1)],...
        'VariableNames',Variable_name_table);
    Table_Bruggeman_taufactor_voxelsizedependence = array2table([direction(1).property_voxelsizedependence(:,:,2) direction(2).property_voxelsizedependence(:,2:end,2) direction(3).property_voxelsizedependence(:,2:end,2)],...
        'VariableNames',Variable_name_table);
    Table_DeffDbulk_taufactor_voxelsizedependence = array2table([direction(1).property_voxelsizedependence(:,:,3) direction(2).property_voxelsizedependence(:,2:end,3) direction(3).property_voxelsizedependence(:,2:end,3)],...
        'VariableNames',Variable_name_table);
    
    %% DISPLAY TEXT RESULTS
    if (OPTIONS.displaytext==1)
        fprintf('> Tortuosity factor dependence with the voxel size:\n\n');
        disp(Table_tortuosity_taufactor_voxelsizedependence)
        fprintf('> Bruggeman exponent dependence with the voxel size:\n\n');
        disp(Table_Bruggeman_taufactor_voxelsizedependence)        
        fprintf('> Normalized effective diffusion coefficient dependence with the voxel size:\n\n');
        disp(Table_DeffDbulk_taufactor_voxelsizedependence)          
    end
    
    %% SAVE RESULTS
    Results_tortuosityfactor.voxelsizedependence_tau = Table_tortuosity_taufactor_voxelsizedependence; % Save in main table result
    Results_tortuosityfactor.voxelsizedependence_Bruggeman = Table_Bruggeman_taufactor_voxelsizedependence;
    Results_tortuosityfactor.voxelsizedependence_DeffDbulk = Table_DeffDbulk_taufactor_voxelsizedependence;
    if OPTIONS.save_xls==true
        filename = 'Transport_voxel_size_dependence'; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name='Tortuosity_factor';
        DATA_writetable.sheet(1).table=Table_tortuosity_taufactor_voxelsizedependence;
        DATA_writetable.sheet(2).name='Bruggeman_exponent';
        DATA_writetable.sheet(2).table=Table_Bruggeman_taufactor_voxelsizedependence;    
        DATA_writetable.sheet(3).name='Deff_Dbulk';
        DATA_writetable.sheet(3).table=Table_DeffDbulk_taufactor_voxelsizedependence;           
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end       
    
    %% FIGURES
    parameters_figure.number_phase = number_phase;
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.INFO = INFO;
    parameters_figure.OPTIONS = OPTIONS;
    
    parameters_figure.propertyname = 'Tortuosity factor';
    parameters_figure.method = 'Tau factor';
    parameters_figure.interpolation_voxelsize = direction(1).interpolation_voxelsize(:,1);
    parameters_figure.property_voxelsizedependence = direction(1).property_voxelsizedependence(:,:,1);
    parameters_figure.interpolation_voxelsize2 = direction(2).interpolation_voxelsize(:,1);
    parameters_figure.property_voxelsizedependence2 = direction(2).property_voxelsizedependence(:,:,1);
    parameters_figure.interpolation_voxelsize3 = direction(3).interpolation_voxelsize(:,1);
    parameters_figure.property_voxelsizedependence3 = direction(3).property_voxelsizedependence(:,:,1);
    parameters_figure.str_ylabel = 'Tortuosity factor \tau';
    parameters_figure.propertynameunit = '';
    parameters_figure.filename = 'Tortuosityfactor_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures
    
    parameters_figure.propertyname = 'Bruggeman exponent';
    parameters_figure.method = 'Tau factor';
    parameters_figure.interpolation_voxelsize = direction(1).interpolation_voxelsize(:,2);
    parameters_figure.property_voxelsizedependence = direction(1).property_voxelsizedependence(:,:,2);
    parameters_figure.interpolation_voxelsize2 = direction(2).interpolation_voxelsize(:,2);
    parameters_figure.property_voxelsizedependence2 = direction(2).property_voxelsizedependence(:,:,2);
    parameters_figure.interpolation_voxelsize3 = direction(3).interpolation_voxelsize(:,2);
    parameters_figure.property_voxelsizedependence3 = direction(3).property_voxelsizedependence(:,:,2);
    parameters_figure.str_ylabel = 'Bruggeman exponent p';
    parameters_figure.propertynameunit = '';
    parameters_figure.filename = 'Bruggeman_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures
    
    parameters_figure.propertyname = 'Normalized effective diffusion coefficient';
    parameters_figure.method = 'Tau factor';
    parameters_figure.interpolation_voxelsize = direction(1).interpolation_voxelsize(:,3);
    parameters_figure.property_voxelsizedependence = direction(1).property_voxelsizedependence(:,:,3);
    parameters_figure.interpolation_voxelsize2 = direction(2).interpolation_voxelsize(:,3);
    parameters_figure.property_voxelsizedependence2 = direction(2).property_voxelsizedependence(:,:,3);
    parameters_figure.interpolation_voxelsize3 = direction(3).interpolation_voxelsize(:,3);
    parameters_figure.property_voxelsizedependence3 = direction(3).property_voxelsizedependence(:,:,3);
    parameters_figure.str_ylabel = 'D_{eff}/D_{bulk}';
    parameters_figure.propertynameunit = '';
    parameters_figure.filename = 'DeffDbulk_voxel_size_dependence';
    Function_create_figure_voxelsizedependence(parameters_figure) % Figures
    
    parameters_figure.direction = direction;
    parameters_figure.filename = 'TauEpsCorrelation_voxel_size_dependence';
    Function_create_figure_TauPorosity(parameters_figure,'Voxel_size_dependence') % Figures
    
    
end

%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) ANALYSIS
%%

if PROPERTY.tortuosity_taufactor.number_RVE>0
    for n_RVE=1:1:PROPERTY.tortuosity_taufactor.number_RVE % Loop over all RVE asked
        RVEparameters.name = PROPERTY.tortuosity_taufactor.RVE(n_RVE).name;
        RVEparameters.savename = PROPERTY.tortuosity_taufactor.RVE(n_RVE).savename;
        RVEparameters.type = PROPERTY.tortuosity_taufactor.RVE(n_RVE).type;
        RVEparameters.divisions = PROPERTY.tortuosity_taufactor.RVE(n_RVE).divisions;
        RVEparameters.subs2 = PROPERTY.tortuosity_taufactor.RVE(n_RVE).subs2;
        RVEparameters.subs4 = PROPERTY.tortuosity_taufactor.RVE(n_RVE).subs4;
        RVEparameters.Aspectratio = PROPERTY.tortuosity_taufactor.RVE(n_RVE).Aspectratio;
        if  strcmp(PROPERTY.tortuosity_taufactor.RVE(n_RVE).type,'A')
            RVEparameters.Aspectratio_name = [num2str(Domain_size(1)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(2)/Domain_size(3),'%1.3f\t') ' ' num2str(Domain_size(3)/Domain_size(3),'%1.3f\t')];
        elseif strcmp(PROPERTY.tortuosity_taufactor.RVE(n_RVE).type,'B') || strcmp(PROPERTY.tortuosity_taufactor.RVE(n_RVE).type,'D')
            RVEparameters.Aspectratio_name = [num2str(RVEparameters.Aspectratio(1)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(2)/RVEparameters.Aspectratio(3),'%1.3f\t') ' ' num2str(RVEparameters.Aspectratio(3)/RVEparameters.Aspectratio(3),'%1.3f\t')];
        end
        RVEparameters.Constantdirection = PROPERTY.tortuosity_taufactor.RVE(n_RVE).Constantdirection;
        RVEparameters.Growthdirection = PROPERTY.tortuosity_taufactor.RVE(n_RVE).Growthdirection;
        RVEparameters.Growthperstep = PROPERTY.tortuosity_taufactor.RVE(n_RVE).Growthperstep;
        RVEparameters.Growthrelativeto = PROPERTY.tortuosity_taufactor.RVE(n_RVE).Growthrelativeto;
        RVEparameters.threshold_std = PROPERTY.tortuosity_taufactor.RVE(n_RVE).threshold_std;
        RVEparameters.threshold_numbersubvolumes = PROPERTY.tortuosity_taufactor.RVE(n_RVE).threshold_numbersubvolumes;
        RVEparameters.firstuniquevolume_size = PROPERTY.tortuosity_taufactor.RVE(n_RVE).firstuniquevolume_size;
        RVEparameters.firstuniquevolume_unit = PROPERTY.tortuosity_taufactor.RVE(n_RVE).firstuniquevolume_unit;
        
        if OPTIONS.save_xls==true || OPTIONS.save_fig == true
            if ispc
                Sub_folder_RVE = [Current_folder RVEparameters.savename '\'];
            else
                Sub_folder_RVE = [Current_folder RVEparameters.savename '/'];
            end
            if exist(Sub_folder_RVE,'dir')==0 % Folder existence is checked, and created if necessary
                mkdir(Sub_folder_RVE);
            end
        end
        
        [All_subdomain,GROUP_SUBDOMAIN, Wholevolume_size] = Function_get_Subdomains(RVEparameters, Domain_size, voxel_size); % Location of all subdomains
        Wholevolume_size(1:2) = Wholevolume_size(1:2)*voxel_size/1000; % Size are set in micrometer
        
        % Information about subdomains
        tmp = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),...
            'VariableNames',{'Subdomain_Id' 'Group_Id' 'Number_subdomain' 'Equivalent_cubic_length' 'Section_length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
         for property_RVE=1:1:3
            for current_direction=1:1:number_dimension
                res(property_RVE).direction(current_direction).RVE(n_RVE).RVEparameters = RVEparameters;
                res(property_RVE).direction(current_direction).RVE(n_RVE).info = tmp;
            end
        end
        
        [number_subdomain,~] = size(All_subdomain); % The number of subdomain
        number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size
        
        %% ALGORITHM
        % Initialisation
        % Colunm 1 is the subdomain id
        % Colunm 2 and 3 are the sizes of the subdomain.
        Property_eachsubdomain = zeros(number_subdomain,number_phase+3,3,number_dimension); %1 Tortuosity factor, 2 Bruggeman exponent, 3 Deff_Dbulk
        % Property calculated for each subdomain
        for subdomain_id = 1:1:number_subdomain
            % Boundary of the subdomain
            x0 = All_subdomain(subdomain_id,6); x1 = All_subdomain(subdomain_id,7);
            y0 = All_subdomain(subdomain_id,8); y1 = All_subdomain(subdomain_id,9);
            z0 = All_subdomain(subdomain_id,10); z1 = All_subdomain(subdomain_id,11);
            clear current_subdomain;
            current_subdomain = Phase_microstructure(x0:x1,y0:y1,z0:z1);
            Current_domain_size = size(current_subdomain);
            current_voxel_number = numel(current_subdomain);
            Property_eachsubdomain(subdomain_id,1,:,:)=subdomain_id;
            % Equivalent size of the subdomain
            Property_eachsubdomain(subdomain_id,2,:,:)=All_subdomain(subdomain_id,4)*voxel_size/1000; % Size are set in micrometer
            Property_eachsubdomain(subdomain_id,3,:,:)=All_subdomain(subdomain_id,5)*voxel_size/1000;            
            % CPU and stopwatch time - start
            time_cpu_start = cputime;
            tic;
            for current_phase=1:1:number_phase % Loop over all phases
                code_tmp = INFO.phase(current_phase).code; % Code of the phase
                binary_phase=zeros(Current_domain_size(1),Current_domain_size(2),Current_domain_size(3)); % Initialization
                binary_phase(current_subdomain == code_tmp) = 1;
                
                voxel_number = sum(sum(sum(binary_phase==1)));
                eps=voxel_number/numel(binary_phase);
                for current_direction=1:1:number_dimension % Loop on every analysed direction
                    time_cpu_start = cputime; % CPU time start
                    tic; % Stopwatch start
                    % Call tortuosity algorithm (cf. Tau factor manual)
                    if current_direction==1
                        Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;1 0 0],[1 1 1]);
                        if Tau_factor_result.Tau_W1.Tau == Inf
                            Tau_factor_result.Tau_W1.Tau = NaN;
                        end
                        Property_eachsubdomain(subdomain_id,current_phase+3,1,current_direction)=Tau_factor_result.Tau_W1.Tau;
                    elseif current_direction==2
                        Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 1 0],[1 1 1]);
                        if Tau_factor_result.Tau_W2.Tau == Inf
                            Tau_factor_result.Tau_W2.Tau = NaN;
                        end
                        Property_eachsubdomain(subdomain_id,current_phase+3,1,current_direction)=Tau_factor_result.Tau_W2.Tau;
                    elseif current_direction==3
                        Tau_factor_result = TauFactor('InLine',1,0,0,binary_phase,[0 0 0;0 0 0;0 0 1],[1 1 1]);
                        if Tau_factor_result.Tau_W3.Tau == Inf
                            Tau_factor_result.Tau_W3.Tau = NaN;
                        end
                        Property_eachsubdomain(subdomain_id,current_phase+3,1,current_direction)=Tau_factor_result.Tau_W3.Tau;
                    end
                    tau = Property_eachsubdomain(subdomain_id,current_phase+3,1,current_direction);
                    Property_eachsubdomain(subdomain_id,current_phase+3,2,current_direction) = 1 - log(tau)/log(eps);
                    Property_eachsubdomain(subdomain_id,current_phase+3,3,current_direction) = 1*eps/tau;
                    % CPU and stopwatch time - end
                    time_cpu_elapsed = cputime-time_cpu_start;
                    time_stopwatch_elapsed = toc;
                    Time_tmp = [voxel_number time_cpu_elapsed time_stopwatch_elapsed];
                    Time_measure = [Time_measure;Time_tmp];
                end
            end
        end
        
        %% STATISTICAL ANALYSIS and RVE SIZE
        str_RVE(1).name = 'Tortuosityfactor_RVE_';
        str_RVE(2).name = 'Bruggeman_RVE_';
        str_RVE(3).name = 'DeffDbulk_RVE_';
        str_RVE(1).propertyname = 'Tortuosity factor';
        str_RVE(2).propertyname = 'Bruggeman exponent';  
        str_RVE(3).propertyname = 'Normalized eff. diff. coeff.';      
        str_RVE(1).propertynameunit = '';
        str_RVE(2).propertynameunit = '';
        str_RVE(3).propertynameunit = '';
        
        res(1).Wholevolume_results = Tortuosity_factor(:,:);
        res(2).Wholevolume_results = Bruggeman_exponent(:,:);
        res(3).Wholevolume_results = Effective_diffusion_coefficient(:,:);  
        
        for property_RVE=1:1:3
            for current_direction=1:1:number_dimension
                [res(property_RVE).direction(current_direction).Property_subdomains_statistics, res(property_RVE).direction(current_direction).Size_RVE] = Function_subdomains_statistical_analysis(number_group_size,number_phase,GROUP_SUBDOMAIN,Property_eachsubdomain(:,:,property_RVE,current_direction),voxel_size,RVEparameters);
                % For correlation
                for current_phase=1:1:number_phase
                    if max(res(property_RVE).direction(current_direction).Size_RVE(2,:,1))~=0
                        results_correlation(current_phase).([str_RVE(property_RVE).name INFO.direction(current_direction).filename '_'  RVEparameters.savename]) = res(property_RVE).direction(current_direction).Size_RVE(2,current_phase,1);
                    end
                    if strcmp(RVEparameters.type,'C')
                        if max(res(property_RVE).direction(current_direction).Size_RVE(2,:,2))~=0
                            results_correlation(current_phase).([str_RVE(property_RVE).name INFO.direction(current_direction).filename '_length_' RVEparameters.savename]) = res(property_RVE).direction(current_direction).Size_RVE(2,current_phase,2);
                        end
                    end
                end
                
                %% MANAGING RESULTS
                [res(property_RVE).direction(current_direction).RVE] = Function_subdomains_manage_results(Property_eachsubdomain(:,:,property_RVE,current_direction),...
                                                                                                          res(property_RVE).direction(current_direction).Property_subdomains_statistics,...
                                                                                                          RVEparameters,...
                                                                                                          res(property_RVE).direction(current_direction).RVE,...
                                                                                                          res(property_RVE).direction(current_direction).Size_RVE,...
                                                                                                          n_RVE,number_phase,INFO);
                
                %% TEXT DISPLAY AND SAVE RESULTS
                propertyname=[str_RVE(property_RVE).propertyname ' ' INFO.direction(current_direction).name];
                if property_RVE>1 || current_direction>1
                    INFO.showrveparameters=false;
                else
                    INFO.showrveparameters=true;
                end
                Function_subdomains_display_and_save(OPTIONS,INFO,res(property_RVE).direction(current_direction).RVE,n_RVE,RVEparameters,number_phase,propertyname,Sub_folder_RVE)
                
                %% FIGURES
                parameters_figure.propertyname = propertyname;
                parameters_figure.propertynameunit = str_RVE(property_RVE).propertynameunit;
                parameters_figure.RVE = RVEparameters;
                parameters_figure.Criterion=[RVEparameters.threshold_std RVEparameters.threshold_numbersubvolumes];
                parameters_figure.savefolder = Sub_folder_RVE;
                parameters_figure.OPTIONS = OPTIONS;
                parameters_figure.INFO = INFO;
                parameters_figure.number_phase = number_phase;
                parameters_figure.Property_subdomains_statistics = res(property_RVE).direction(current_direction).Property_subdomains_statistics;
                parameters_figure.Property_eachsubdomain = Property_eachsubdomain(:,:,property_RVE,current_direction);
                parameters_figure.Size_RVE = res(property_RVE).direction(current_direction).Size_RVE;
                parameters_figure.Wholevolume_results = res(property_RVE).Wholevolume_results(:,current_direction);
                parameters_figure.Wholevolume_size =  Wholevolume_size;
                Function_create_figures_RVE(parameters_figure) % Figures
            end
        end
        
    end
    Results_tortuosityfactor.RVE.tortuosityfactor = res(1).direction; % Save in main table result
    Results_tortuosityfactor.RVE.Bruggeman = res(2).direction;
    Results_tortuosityfactor.RVE.DeffDbulk = res(3).direction;
end


%%
%% ENDING FUNCTION
%%

%% TIME
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
filename = 'Taufactor_calculation_time';
[Results_tortuosityfactor] = function_time_figure(Time_measure,date_start, date_end, Results_tortuosityfactor, Current_folder, filename, 'Tortuosity factor (taufactor)', OPTIONS);
 
%% SAVE RESULTS
if OPTIONS.save_resultsmat == true
    Sub_folder = 'Summary\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Results_tortuosityfactor_taufactor'],'Results_tortuosityfactor')
    Sub_folder = 'Correlation\'; % Result are saved in this subfolder
    Current_folder = [main_folder Sub_folder];
    if exist(Current_folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(Current_folder);
    end
    save([Current_folder 'Correlation_tortuosityfactor_taufactor'],'results_correlation')
end

end