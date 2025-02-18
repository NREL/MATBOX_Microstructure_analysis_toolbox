function [] = function_microstructure_characterization_segmentedvolume(infovol, optc, opts)
% Calculate microstructure properties useful for energy storage macroscale models.
% For instance, for Pseudo-2D battery models.
% Intended use of this function is through the GUI of src\Microstructure_characterization\Microstructure_characterization.mlapp
% infovol : volume information
% optc    : options common to all volumes
% opts    : options volume specific 

% Tables below are re-created as they are saved in an excel file for each volume
date_start = datetime('now','TimeZone','local','Format','eeee, HH:mm:ss Z, MMMM d, y');
Table_date = table({char(date_start)},...
    'VariableNames',{'Date'});
Table_toolbox = table({'Toolbox';'Version';'GitHub repository';'Official webpage';'Journal article'},{optc.toolbox.toolboxname;optc.toolbox.version;optc.toolbox.githubrepo;optc.toolbox.officialwebpage;optc.toolbox.article},...
    'VariableNames',{'Item','Value'});
Table_author = table(optc.toolbox.authorlist(:,1),optc.toolbox.authorlist(:,2),...
    'VariableNames',{'Author','Affiliation'});
username = getenv('USERNAME');
computername = getenv('COMPUTERNAME');
os = getenv('OS');
Table_User = table({username},{computername},{os},...
    'VariableNames',{'User_name','Computer_name','Operating_system'});

%% IMPORT
if max(max(infovol.initial_phaselabel)) <= 255
    fileformat = 'uint8';
else
    fileformat = 'uint16';
end
disp(['Loading file: ' infovol.loadingpath ', format: ' fileformat]);
disp 'Please wait...';
[Phase_microstructure,outcome] = function_load_tif(infovol.loadingpath,fileformat);
if ~outcome % fail to import
   return
end
disp 'Import successful !';

%% FOLDER
if ~exist(infovol.volpath,'dir') % Check existence of the main folder
    mkdir(infovol.volpath); % Create folder if not exist
end
if ispc
    separator = '\';
else
    separator = '/';
end

%% SAVE LOADED VOLUME
volume_folder = [infovol.volpath 'Volume' separator];
if ~exist(volume_folder,'dir') % Check existence of the main folder
    mkdir(volume_folder); % Create folder if not exist
end
function_save_tif(Phase_microstructure, [volume_folder 'loaded_volume.tif']);

%% GENERAL INFORMATION ABOUT THE VOLUME
Domain_size = size(Phase_microstructure); % Size of the loaded domain (number of voxel)
Domain_size_unit = Domain_size*infovol.initial_voxelsize; % Size of the loaded domain
number_direction = length(Domain_size); % Number of direction (2D or 3D)
initial_numbervoxel = prod(Domain_size);
% Create table
if number_direction==2
    equivalent_square_size_unit = prod(Domain_size_unit)^(1/2);
    Table_domainsize = table([1;2],[infovol.directionname(1); infovol.directionname(2)], [Domain_size(1); Domain_size(2)],[Domain_size_unit(1); Domain_size_unit(2)], [equivalent_square_size_unit;equivalent_square_size_unit],...
        'VariableNames',{'Direction', 'Name', 'Number of voxel' ['Length ' infovol.unit] ['Equivalent squareroot length ' infovol.unit]});
else
    equivalent_cubic_size_unit = prod(Domain_size_unit)^(1/3);
    Table_domainsize = table([1;2;3],[infovol.directionname(1); infovol.directionname(2); infovol.directionname(3)], [Domain_size(1); Domain_size(2); Domain_size(3)],[Domain_size_unit(1); Domain_size_unit(2); Domain_size_unit(3)], [equivalent_cubic_size_unit; equivalent_cubic_size_unit; equivalent_cubic_size_unit],...
        'VariableNames',{'Direction', 'Name', 'Number of voxel' ['Length ' infovol.unit] ['Equivalent cubicroot length ' infovol.unit]});
end
Table_voxel = table([{['Voxel size (' infovol.unit ')']};{'Number of voxel'}], [infovol.initial_voxelsize; prod(Domain_size)],...
    'VariableNames',{'Information', 'Value'});

% Size of the loaded data
tmp=whos('Phase_microstructure'); data_MB=tmp.bytes*9.53674e-7;
Table_inputoutput = table([{'File loaded'};{'Array size (MB)'};{'Save folder'}], [{infovol.loadingpath}; num2str(data_MB,'%1.1f'); {infovol.volpath}],...
    'VariableNames',{'Information', 'Text'});

Table_phase = table(infovol.initial_phasename,infovol.initial_phaselabel(:,1),...
    'VariableNames',{'Phase name' 'Phase label'});

Table_volumeinformation = table(infovol.volumeinformation(:,1),infovol.volumeinformation(:,2),...
    'VariableNames',{'Information' 'Text'});

disp ' ';
disp(Table_inputoutput);disp ' ';
disp(Table_volumeinformation);disp ' ';
disp(Table_phase);disp ' ';
disp(Table_voxel);disp ' ';
disp(Table_domainsize);disp ' ';

if optc.save.xls
    % Intialization
    sheet_number=0;
    clear DATA_writetable
    % Input/output
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='IO';
    DATA_writetable.sheet(sheet_number).table=Table_inputoutput;    
    % Phase name
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Phase';
    DATA_writetable.sheet(sheet_number).table=Table_phase;        
    % Volume information    
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Volume_information';
    DATA_writetable.sheet(sheet_number).table=Table_volumeinformation;    
    % Domain size
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Domain_size';
    DATA_writetable.sheet(sheet_number).table=Table_domainsize;    
    % Voxel
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Voxel';
    DATA_writetable.sheet(sheet_number).table=Table_voxel;          
    % User information
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='User_information';
    DATA_writetable.sheet(sheet_number).table=Table_User;    
    % Date
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Start_date';
    DATA_writetable.sheet(sheet_number).table=Table_date;
    % Version information
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Toolbox_version';
    DATA_writetable.sheet(sheet_number).table=Table_toolbox;
    % Author information
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Toolbox_authors';
    DATA_writetable.sheet(sheet_number).table=Table_author;        

    % Filename without extension
    filename = 'General_information';
    % Save function
    Function_Writetable(infovol.volpath,filename,DATA_writetable)
end

%% TIME
time_cpu_start = cputime; % CPU start
tStart = tic; % Stopwatch start
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % Date start

%%
%% SETUP VOLUME
%%
disp '********************';
disp '>>> SETUP VOLUME <<<';
disp '********************';
disp ' ';

%% REGION OF INTEREST
% Update tables
ROI_size = infovol.ROI(:,2)-infovol.ROI(:,1)+1;
ROI_size_unit = ROI_size*infovol.initial_voxelsize;
if number_direction==2
    equivalent_squareROI_size_unit = prod(ROI_size_unit)^(1/2);
    Table_domainsize_ROI = table([1;2],[infovol.directionname(1); infovol.directionname(2)], infovol.ROI(:,1),infovol.ROI(:,2), ROI_size, ROI_size_unit, [equivalent_squareROI_size_unit; equivalent_squareROI_size_unit],...
        'VariableNames',{'Direction', 'Name', 'ROI Start', 'ROI End', 'Number of voxel', ['Length ' infovol.unit], ['Equivalent squareroot length ' infovol.unit]});
else
    equivalent_cubicROI_size_unit = prod(ROI_size_unit)^(1/3);
    Table_domainsize_ROI = table([1;2;3],[infovol.directionname(1); infovol.directionname(2); infovol.directionname(3)], infovol.ROI(:,1),infovol.ROI(:,2), ROI_size, ROI_size_unit, [equivalent_cubicROI_size_unit; equivalent_cubicROI_size_unit; equivalent_cubicROI_size_unit],...
        'VariableNames',{'Direction', 'Name', 'ROI Start', 'ROI End', 'Number of voxel', ['Length ' infovol.unit], ['Equivalent cubicroot length ' infovol.unit]});
end
Table_voxel_ROI = table([{['Voxel size (' infovol.unit ')']};{'Number of voxel'}], [infovol.initial_voxelsize; prod(Domain_size)],...
    'VariableNames',{'Information', 'Value'});

% Check if crop required
if number_direction==2
    check_ROI_start = min(infovol.ROI(:,1)==[1;1]);
    check_ROI_end = min(infovol.ROI(:,2)==[Domain_size(1);Domain_size(2)]);
else
    check_ROI_start = min(infovol.ROI(:,1)==[1;1;1]);
    check_ROI_end = min(infovol.ROI(:,2)==[Domain_size(1);Domain_size(2);Domain_size(3)]);
end

if check_ROI_start*check_ROI_end==0 % Compare with region of interest
    % Crop
    Phase_microstructure = Phase_microstructure( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2) , infovol.ROI(3,1):infovol.ROI(3,2));    
    Domain_size = size(Phase_microstructure);
    disp 'Volume has been cropped'; disp ' ';
    disp(Table_domainsize_ROI); disp ' ';
    disp(Table_voxel_ROI);disp ' ';
else
    disp 'Volume has NOT been crooped';  disp ' ';
end
tmp=whos('Phase_microstructure'); data_cropMB=tmp.bytes*9.53674e-7; % Keep track of data size    
cropped_numbervoxel = prod(Domain_size);

%% RE-LABEL PHASES
Table_phase_relabel = table(infovol.phasename,infovol.phaselabel,...
    'VariableNames',{'Phase name' 'Phase label'});
if sum(infovol.initial_phaselabel(:,1)~=infovol.initial_phaselabel(:,2)) % At least one phase is re-assigned
    number_initialphase = length(infovol.initial_phasename); % Number of assigned phase
    tmp = Phase_microstructure; % Temporary variable
    for current_phase=1:1:number_initialphase
        old_label = infovol.initial_phaselabel(current_phase,1);
        new_label = infovol.initial_phaselabel(current_phase,2);
        index_ = find(Phase_microstructure==old_label);
        tmp(index_)=new_label;
    end
    Phase_microstructure = tmp; % Assign
    clear tmp % Clean temporary variable
    disp 'Phase label has been re-assigned';  disp ' ';
    disp(Table_phase_relabel);disp ' ';
else
    disp 'Phase label has NOT been re-assigned';  disp ' ';
end
Phase_microstructure=uint8(Phase_microstructure);

%% SCALING (IMAGE RESOLUTION)
if infovol.scaling_factor ~= 1.0
    % Set parameters
    parameters_scaling.scaling_factor = infovol.scaling_factor;
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = min(infovol.phaselabel);
    % Scale
    Phase_microstructure = function_scaling(Phase_microstructure,parameters_scaling);
    Domain_size = size(Phase_microstructure);
    Domain_size_unit = Domain_size*infovol.voxelsize;
    % Create table
    if number_direction==2
        equivalent_square_size_unit = prod(Domain_size_unit)^(1/2);
        Table_domainsize_ROI_scaling = table([1;2],[infovol.directionname(1); infovol.directionname(2)], [Domain_size(1); Domain_size(2)],[Domain_size_unit(1); Domain_size_unit(2)], [equivalent_square_size_unit;equivalent_square_size_unit],...
            'VariableNames',{'Direction', 'Name', 'Number of voxel' ['Length ' infovol.unit] ['Equivalent squareroot length ' infovol.unit]});
    else
        equivalent_cubic_size_unit = prod(Domain_size_unit)^(1/3);
        Table_domainsize_ROI_scaling = table([1;2;3],[infovol.directionname(1); infovol.directionname(2); infovol.directionname(3)], [Domain_size(1); Domain_size(2); Domain_size(3)],[Domain_size_unit(1); Domain_size_unit(2); Domain_size_unit(3)], [equivalent_cubic_size_unit; equivalent_cubic_size_unit; equivalent_cubic_size_unit],...
            'VariableNames',{'Direction', 'Name', 'Number of voxel' ['Length ' infovol.unit] ['Equivalent cubicroot length ' infovol.unit]});
    end
    Table_voxel_ROI_scaling = table([{['Voxel size (' infovol.unit ')']};{'Number of voxel'}], [infovol.voxelsize; prod(Domain_size)],...
        'VariableNames',{'Information', 'Value'});

    disp 'Voxel size has been changed';  disp ' ';
    disp(Table_domainsize_ROI_scaling);disp ' ';
    disp(Table_voxel_ROI_scaling);disp ' ';
else
    Table_domainsize_ROI_scaling = Table_domainsize_ROI;
    Table_voxel_ROI_scaling = Table_voxel_ROI;
    disp 'Voxel size has NOT been changed';  disp ' ';
end
tmp=whos('Phase_microstructure'); data_voxelresizeMB=tmp.bytes*9.53674e-7; % Keep track of data size
voxelresize_numbervoxel = prod(Domain_size);

%% MEMORY
Table_memory = table([{'Loaded data'};{'Region of interest'};{'Voxel resize'}],[initial_numbervoxel;cropped_numbervoxel;voxelresize_numbervoxel], [data_MB;data_cropMB;data_voxelresizeMB],...
    'VariableNames',{'Microstructure array','Number of voxel', 'Memory MB'});
disp 'Memory';  disp ' ';
disp(Table_memory);disp ' ';

%% SAVE INVESTIGATED VOLUME
function_save_tif(Phase_microstructure, [volume_folder 'Investigated_volume.tif']);

%% SAVE INFORMATION IN EXCEL SHEET
if optc.save.xls
    % Intialization
    sheet_number=0;
    clear DATA_writetable
    % Region of interest: domain size
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='ROI_domainsize';
    DATA_writetable.sheet(sheet_number).table=Table_domainsize_ROI;
    % Region of interest: Voxel
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='ROI_voxelsize';
    DATA_writetable.sheet(sheet_number).table=Table_voxel_ROI;    
    % Phase name
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Relabel_phase';
    DATA_writetable.sheet(sheet_number).table=Table_phase_relabel;
    % Voxel resize: domain size
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Scaling_domainsize';
    DATA_writetable.sheet(sheet_number).table=Table_domainsize_ROI_scaling;
    % Voxel resize: Voxel
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Scaling_voxelsize';
    DATA_writetable.sheet(sheet_number).table=Table_voxel_ROI_scaling;       
    % Size
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Memory';
    DATA_writetable.sheet(sheet_number).table=Table_memory;  
    
    % Filename without extension
    filename = 'Volume_setup';
    % Save function
    Function_Writetable(infovol.volpath,filename,DATA_writetable)
end

%%
%% CHARACTERIZATION
%%

disp '********************';
disp '>>> CALCULATIONS <<<';
disp '********************';
disp ' ';

foo=1; % Un-used inputs to have a different number of argument between GUI and standalone use for some functions.

%% VOLUME FRACTIONS
if sum(opts.volumefraction.todo)
    Function_Volume_fractions(Phase_microstructure, infovol, optc, opts.volumefraction) % Call function
end

%% CONNECTIVITY
if opts.connectivity.todo
     Function_connectivity(Phase_microstructure, infovol, optc, opts.connectivity) % Call function
end

%% TORTUOSITY FACTORS
if sum(opts.tortuosity.todo)
    Function_Tortuosity_factor_taufactor(Phase_microstructure, infovol, optc, opts.tortuosity) % Call function
end

%% SPECIFIC SURFACE AREA, DIRECT
if sum(opts.Sp_direct.todo)
    Function_Specificsurface_direct(Phase_microstructure, infovol, optc, opts.Sp_direct, foo, foo) % Call function
end

%% SPECIFIC INTERFACE AREA, DIRECT
if opts.Int_direct.todo
    Function_Specificinterface_direct(Phase_microstructure, infovol, optc, opts.Int_direct, foo, foo) % Call function
end

%% TRIPLE PHASE BOUNDARY LENGTH, DIRECT
if opts.TPBL_direct.todo
    Function_TPBL_direct(Phase_microstructure, infovol, optc, opts.TPBL_direct, foo, foo) % Call function
end

%% DIAMETER, C-PSD
if sum(opts.D_cpsd.todo)
    Function_particle_size_CPSD(Phase_microstructure, infovol, optc, opts.D_cpsd, foo) % Call function
end

%% DIAMETER, EDMF
if sum(opts.D_edmf.todo)
    Function_particle_size_EDMF(Phase_microstructure, infovol, optc, opts.D_edmf) % Call function
end

%% DIAMETER, CHORD
if sum(opts.D_chord.todo)
    Function_particle_size_CHORD(Phase_microstructure, infovol, optc, opts.D_chord) % Call function
end

%% DIAMETER, WATERSHED
if sum(opts.D_watershed.todo)
    Function_particle_size_Watershed(Phase_microstructure, infovol, optc, opts.D_watershed) % Call function
end

% Chord: harmonize with EDMF

% Tau

% Convergence(RVE), and RVE(RVE)

% Article on RVE

% Article on fractal/voxel size
% I know what fractal dimension will bring: we can compare Sp or TPBL from two different materials, only if they
% share same voxel size, and same fractal dimension (indicator of the surface roughness)

% Publish on Github (legacy and new version in choice menu)

% Connectivity

% Publish on Github

% Watershed
% PCRF

% Publish on Github

% Tutorial video

% Correlation module

%%
%% SUMMARY
%%

%% TIME
summary_folder = [infovol.volpath 'Summary' separator]; 
if ~exist(summary_folder,'dir') % Check existence of the main folder
    mkdir(summary_folder); % Create folder if not exist
end

time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
time_stopwatch_elapsed = toc(tStart); % Stopwatch elapsed time
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
% time table
Table_time = table({'Date start';'Date end';'Elasped time (hh:mm:ss)';'Elasped time (s), tic-toc';'CPU time (s)'} ,{char(date_start);char(date_end);char(date_end-date_start);num2str(time_stopwatch_elapsed,'%1.1f');num2str(time_cpu_elapsed,'%1.1f')},...
    'VariableNames',{'Information' 'text'});
disp ' ';
disp '***********************';
disp '>>> RESULTS SUMMARY <<<';
disp '***********************';
disp ' ';
disp(Table_time)
if optc.save.xls
    filename = 'Time'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Time';
    DATA_writetable.sheet(1).table=Table_time;
    % Save function
    Function_Writetable(summary_folder,filename,DATA_writetable)
end

%% SUMMARIZE RESULTS
if optc.save.mat              
    function_summary_onevolume(summary_folder, infovol, optc, opts);
end

end