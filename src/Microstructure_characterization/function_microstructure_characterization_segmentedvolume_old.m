function [] = function_microstructure_characterization_segmentedvolume(OPTIONS, INFO, PROPERTY)
%Calculate classic microstructure properties for battery and fuel cells macromodels:
% - volume fraction
% - connectivity
% - specific surface area
% - triple phase boundary
% - particle diameter
% - tortuosity factor
%You can choose among a set of different numerical methods to finely evaluate properties value.
%Perform Representative Volume Element (RVE) analysis and image resolution sensitivity analysis.
%Calculate morphology parameters (particle shpericity and particle elongation).
%Calculate geometric tortuosity with the graph representation of the microstructure.
%Generate python files to be used with the FEniCS Finite Element Solver for homogeneisation calculation.'

% Version
MAT_charac_version = '0.9';
Table_version = table({char(MAT_charac_version)},...
    'VariableNames',{'Characterization_toolbox_version'});

%% FOLDER
main_folder = [OPTIONS.mainsavefolder '\']; % Main folder that will contain all the results for this microstructure
if ~exist(main_folder,'dir') % Check existence of the main folder
    mkdir(main_folder); % Create folder if not exist
end

%% DATE
date_start = datetime('now','TimeZone','local','Format','eeee, HH:mm:ss Z, MMMM d, y');
Table_date = table({char(date_start)},...
    'VariableNames',{'Date'});

%% USERNAME
username = getenv('USERNAME');
computername = getenv('COMPUTERNAME');
os = getenv('OS');
Table_User = table({username},{computername},{os},...
     'VariableNames',{'User_name','Computer_name','Operating_system'});
  
%% CHECK REQURIED TOOLBOX
is_installed_imagetoolbox = license('test','Image_Toolbox');
is_installed_reportgenerator = license('test','MATLAB_Report_Gen');

%% INITIALIZE ERROR REPORT (BETA)
% if OPTIONS.report.todo == 1 && is_installed_reportgenerator == 1
%     import mlreportgen.dom.*; % Import the DOM API
%     errorreport_path = [main_folder 'Error_Report']; % Save path of the error report
%     OPTIONS.errorreportdoc = Document(errorreport_path, OPTIONS.report.rpt_type); % Create an empty document with the specified type
%     open(OPTIONS.errorreportdoc); % Open file    
%    
%     % Title
%     title = append(OPTIONS.errorreportdoc, Paragraph('Microstructure characterization: Error report'));
%     title.Style = OPTIONS.report.style.title;
%     % Title 2
%     str = sprintf ('%s',INFO.filename_input);
%     title = append(OPTIONS.errorreportdoc, Paragraph(str));
%     title.Style = OPTIONS.report.style.title2;
%     % Empty text
%     standard_text = append(OPTIONS.errorreportdoc, Paragraph('-')); standard_text.Style = OPTIONS.report.style.standardtext; standard_text.HAlign='Center';
%      % Volume basic information
%     str1 = sprintf ('File loaded from: %s', OPTIONS.loadingpath);
%     str2 = sprintf ('Results saved in: %s', main_folder);
%     cellstr = {str1, str2};
%     for k=1:1:length(cellstr) % Loop over all strings contained in the cell
%         text = append(OPTIONS.errorreportdoc, Paragraph( char(cellstr(k))) ); % Append string
%         text.Style = OPTIONS.report.style.standardtext; % Set style
%     end        
%     % Empty text
%     standard_text = append(OPTIONS.errorreportdoc, Paragraph('-')); standard_text.Style = OPTIONS.report.style.standardtext; standard_text.HAlign='Center';
%     % Date, username, machine and version
%     str1 = sprintf ('Started %s',date_start);
%     str2 = sprintf ('By user %s, on machine %s, with the OS %s',username,computername,os);    
%     str3 = sprintf ('Matlab version: %s',version);
%     str4 = sprintf ('Microstructure analysis tool, characterization module version: %s',MAT_charac_version);
%     str5 = 'Developed by F. Usseglio-Viretta, NREL';
%     cellstr = {str1, str2, str3, str4, str5};
%     for k=1:1:length(cellstr) % Loop over all strings contained in the cell
%         text = append(OPTIONS.errorreportdoc, Paragraph( char(cellstr(k))) ); % Append string
%         text.Style = OPTIONS.report.style.standardtext; % Set style
%     end        
%     % Close report
%     close(OPTIONS.errorreportdoc)
% end

%% LOAD FILE
% Call function
[Phase_microstructure, outcome] = function_loadvolume(OPTIONS.loadingpath, OPTIONS.arrayformat, OPTIONS.structurefieldname);
if outcome.success==false
    str1 = sprintf ('Error during loading file');
    str2 = sprintf ('Error message: %s', outcome.reason);
    cellstr = {str1, str2};
%     if is_installed_reportgenerator == 1 % Write error report
%         open(OPTIONS.errorreportdoc); % Open file
%         for k=1:1:length(cellstr) % Loop over all strings contained in the cell
%             text = append(OPTIONS.errorreportdoc, Paragraph( char(cellstr(k))) ); % Append string
%             text.Style = OPTIONS.report.style.standardtext; % Set style
%         end
%     elseif OPTIONS.displaytext
%         fprintf('%s. %s\n',str1, str2);
%         pause(5); % Time to read the error
%     end
    if OPTIONS.displaytext
        fprintf('%s. %s\n',str1, str2);
        pause(5); % Time to read the error
    end
    return % Return to main function
end

%% SAVE LOADED VOLUME

volume_folder = [main_folder 'Volume\'];
if ~exist(volume_folder,'dir') % Check existence of the main folder
    mkdir(volume_folder); % Create folder if not exist
end
function_save_tif(Phase_microstructure, [volume_folder 'loaded_volume.tif'])

%% DIMENSIONS
Domain_size = size(Phase_microstructure); % Size of the loaded domain (number of voxel)
initial_numbervoxel = prod(Domain_size);
Domain_size_um = Domain_size*INFO.initial_voxelsize/1000; % Size of the loaded domain (in micrometers)
number_direction = length(INFO.direction); % Number of direction (2D or 3D)
ROI_size = INFO.RegionOfInterest(:,2)-INFO.RegionOfInterest(:,1)+1;
ROI_size_um = ROI_size*INFO.initial_voxelsize/1000;
ratio_voxelsize = INFO.asked_voxelsize/INFO.initial_voxelsize;
% Create table
if number_direction==2
    equivalent_square_size_um = prod(Domain_size_um)^(1/2);
    numbervoxel_newvoxelsize = initial_numbervoxel/(ratio_voxelsize^2);
    Table_domainsize = table([1;2],[{INFO.direction(1).consistentname}; {INFO.direction(2).consistentname}], [Domain_size(1); Domain_size(2)],[Domain_size_um(1); Domain_size_um(2)], [equivalent_square_size_um;equivalent_square_size_um],...
        'VariableNames',{'Direction', 'Name', 'Number_of_voxel' 'Length_micrometers' 'Equivalent_square_length_micrometers'});
    equivalent_squareROI_size_um = prod(ROI_size_um)^(1/2);
    numbervoxelROI_newvoxelsize = prod(ROI_size)/(ratio_voxelsize^2);
    Table_Region_of_Interest = table([1;2;3],[{INFO.direction(1).consistentname}; {INFO.direction(2).consistentname}], INFO.RegionOfInterest(:,1),INFO.RegionOfInterest(:,2), ROI_size, ROI_size_um, [equivalent_squareROI_size_um; equivalent_squareROI_size_um],...
        'VariableNames',{'Direction', 'Name', 'ROI_Start', 'ROI_End', 'Number_of_voxel', 'Length_micrometers', 'Equivalent_square_length_micrometers'});
else
    equivalent_cubic_size_um = prod(Domain_size_um)^(1/3);
    numbervoxel_newvoxelsize = initial_numbervoxel/(ratio_voxelsize^3);
    Table_domainsize = table([1;2;3],[{INFO.direction(1).consistentname}; {INFO.direction(2).consistentname}; {INFO.direction(3).consistentname}], [Domain_size(1); Domain_size(2); Domain_size(3)],[Domain_size_um(1); Domain_size_um(2); Domain_size_um(3)], [equivalent_cubic_size_um; equivalent_cubic_size_um; equivalent_cubic_size_um],...
        'VariableNames',{'Direction', 'Name', 'Number_of_voxel' 'Length_micrometers' 'Equivalent_square_length_micrometers'});
    equivalent_cubicROI_size_um = prod(ROI_size_um)^(1/3);
    numbervoxelROI_newvoxelsize = prod(ROI_size)/(ratio_voxelsize^3);
    Table_Region_of_Interest = table([1;2;3],[{INFO.direction(1).consistentname}; {INFO.direction(2).consistentname}; {INFO.direction(3).consistentname}], INFO.RegionOfInterest(:,1),INFO.RegionOfInterest(:,2), ROI_size, ROI_size_um, [equivalent_cubicROI_size_um; equivalent_cubicROI_size_um; equivalent_cubicROI_size_um],...
        'VariableNames',{'Direction', 'Name', 'ROI_Start', 'ROI_End', 'Number_of_voxel', 'Length_micrometers', 'Equivalent_square_length_micrometers'});
end
Table_voxel = table([{'Voxel size (micrometers)'};{'Number of voxel'}], [INFO.initial_voxelsize; initial_numbervoxel],...
    'VariableNames',{'Information', 'Text'});

%% GENERAL INFORMATION

% Size of the loaded data
tmp=whos('Phase_microstructure'); data_MB=tmp.bytes*9.53674e-7;
Table_inputoutput = table([{'File loaded'};{'Array size (MB)'};{'Save folder'}], [{OPTIONS.loadingpath}; num2str(data_MB,'%1.1f'); {main_folder}],...
    'VariableNames',{'Information', 'Text'});
Table_phase = table(INFO.initial_phasename,INFO.initial_phasecode,...
    'VariableNames',{'Phase_name' 'Phase_code'});
Table_volumeinformation = table(INFO.volumeinformation(:,1),INFO.volumeinformation(:,2),...
    'VariableNames',{'Information' 'Text'});

if OPTIONS.displaytext
    disp '****************************';
    disp '>>> GENERAL INFORMATIONS <<<';
    disp '****************************';
    disp ' ';
    disp(Table_date);disp ' ';
    disp(Table_User);disp ' ';
    disp(Table_version);disp ' ';
    disp(Table_inputoutput);disp ' ';
    disp(Table_volumeinformation);disp ' ';
    disp(Table_phase);disp ' ';
    disp(Table_voxel);disp ' ';
    disp(Table_domainsize);disp ' ';
end

if OPTIONS.save_xls
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
    % Version information
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Version';
    DATA_writetable.sheet(sheet_number).table=Table_version;
    % User information
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='User_information';
    DATA_writetable.sheet(sheet_number).table=Table_User;    
    % Date
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Start_date';
    DATA_writetable.sheet(sheet_number).table=Table_date;

    % Filename without extension
    filename = 'General_information';
    % Save function
    Function_Writetable(main_folder,filename,DATA_writetable)
end

%% INITIALIZE REPORT (BETA)
% if OPTIONS.report.todo == 1 && is_installed_reportgenerator == 1
%     report_path = [main_folder 'Report']; % Save path of the report
%     OPTIONS.reportdoc = Document(report_path, OPTIONS.report.rpt_type); % Create an empty document with the specified type
%     open(OPTIONS.reportdoc); % Open file
%     
%     % Title
%     title = append(OPTIONS.reportdoc, Paragraph('Microstructure characterization'));
%     title.Style = OPTIONS.report.style.title;
%     % Title 2
%     str = sprintf ('%s',INFO.filename_input);
%     title = append(OPTIONS.reportdoc, Paragraph(str));
%     title.Style = OPTIONS.report.style.title2;
%     
%     % Decription
%     str1 = sprintf ('Source  : %s', INFO.volume_source);
%     str2 = sprintf ('Material: %s', INFO.volume_material);
%     str3 = sprintf ('Material group: %s', INFO.volume_groupmaterial);
%     str4 = sprintf ('Obversation method: %s', INFO.volume_experimental);
%     str5 = sprintf ('%s', INFO.volume_description);
%     str6 = sprintf ('%s', INFO.analysis_description);   
%     cellstr = {str1, str2, str3, str4, str5, str6};
%     for k=1:1:length(cellstr) % Loop over all strings contained in the cell
%         text = append(OPTIONS.reportdoc, Paragraph( char(cellstr(k))) ); % Append string
%         text.Style = OPTIONS.report.style.standardtext; % Set style
%     end    
%     
%     % Empty text
%     standard_text = append(OPTIONS.reportdoc, Paragraph('-')); standard_text.Style = OPTIONS.report.style.standardtext; standard_text.HAlign='Center';
%      % Volume basic information
%     str1 = sprintf ('File loaded from: %s', INFO.loadingpath);
%     str2 = sprintf ('Results saved in: %s', main_folder);
%     cellstr = {str1, str2};
%     for k=1:1:length(cellstr) % Loop over all strings contained in the cell
%         text = append(OPTIONS.reportdoc, Paragraph( char(cellstr(k))) ); % Append string
%         text.Style = OPTIONS.report.style.standardtext; % Set style
%     end        
%      
%     % Empty text
%     standard_text = append(OPTIONS.reportdoc, Paragraph('-')); standard_text.Style = OPTIONS.report.style.standardtext; standard_text.HAlign='Center';
%     % Date, username
%     str1 = sprintf ('Started %s',date_start);
%     str2 = sprintf ('By user %s, on machine %s, with the OS %s',username,computername,os);    
%     cellstr = {str1, str2};
%     for k=1:1:length(cellstr) % Loop over all strings contained in the cell
%         text = append(OPTIONS.reportdoc, Paragraph( char(cellstr(k))) ); % Append string
%         text.Style = OPTIONS.report.style.standardtext; % Set style
%     end    
%     
%     % Empty text
%     standard_text = append(OPTIONS.reportdoc, Paragraph('-')); standard_text.Style = OPTIONS.report.style.standardtext; standard_text.HAlign='Center';
%     % machine and version
%     str1 = sprintf ('Matlab version: %s',version);
%     str2 = sprintf ('Microstructure analysis tool, characterization module version: %s',MAT_charac_version);
%     str3 = '   Developed by F. Usseglio-Viretta, NREL';
%     str4 = 'This report has been programmatically generated';
%     cellstr = {str1, str2, str3, str4};
%     for k=1:1:length(cellstr) % Loop over all strings contained in the cell
%         text = append(OPTIONS.reportdoc, Paragraph( char(cellstr(k))) ); % Append string
%         text.Style = OPTIONS.report.style.standardtext; % Set style
%         text.WhiteSpace='preserve';
%     end
%     
%     % Close report
%     close(OPTIONS.reportdoc)    
% end

%% TIME
time_cpu_start = cputime; % CPU start
tStart = tic; % Stopwatch start
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % Date stard

%%
%% SETUP VOLUME
%%
if OPTIONS.displaytext
    disp '********************';
    disp '>>> SETUP VOLUME <<<';
    disp '********************';
    disp ' ';
end

%% REGION OF INTEREST
% Check if crop required
Domain_size = size(Phase_microstructure); % Size of the loaded domain
check_ROI_start = min(INFO.RegionOfInterest(:,1)==[1;1;1]);
check_ROI_end = min(INFO.RegionOfInterest(:,2)==[Domain_size(1);Domain_size(2);Domain_size(3)]);
Table_voxel_ROI = Table_voxel; % Initialize
if check_ROI_start*check_ROI_end==0 % Compare with region of interest
    % Resize
    tmp_ = Phase_microstructure( INFO.RegionOfInterest(1,1):INFO.RegionOfInterest(1,2) , INFO.RegionOfInterest(2,1):INFO.RegionOfInterest(2,2) , INFO.RegionOfInterest(3,1):INFO.RegionOfInterest(3,2));
    Phase_microstructure = tmp_; clear tmp_;
    % Update informations
    Domain_size = size(Phase_microstructure);
    Domain_size_um = Domain_size*INFO.initial_voxelsize/1000; % in micrometers
    Table_voxel_ROI = table([{'Voxel size (micrometers)'};{'Number of voxel'}], [INFO.initial_voxelsize; prod(Domain_size)],...
        'VariableNames',{'Information', 'Text'});
    if OPTIONS.displaytext
        disp 'Volume has been cropped'; disp ' ';
        disp(Table_Region_of_Interest);disp ' ';
        disp(Table_voxel_ROI);disp ' ';
    end
else
    if OPTIONS.displaytext
        disp 'Volume has NOT been crooped';  disp ' ';
    end
end
tmp=whos('Phase_microstructure'); data_cropMB=tmp.bytes*9.53674e-7; % Keep track of data size    
cropped_numbervoxel = prod(Domain_size);

%% RE-ASSIGN PHASES
Table_phase = table(INFO.phasename,[INFO.phase(:).code]',...
    'VariableNames',{'Phase_name' 'Phase_code'});
% Check if re-assign is required
if min(INFO.initial_phasecode==INFO.assigned_phasecode)==0
    number_initialphase = length(INFO.initial_phasecode); % Number of assigned phase    
    tmp = Phase_microstructure; % Temporary variable
    for current_phase=1:1:number_initialphase % Loop over all initial phase
       old_code = INFO.initial_phasecode(current_phase);
       new_code = INFO.assigned_phasecode(current_phase);
       index_ = Phase_microstructure==old_code;
       tmp(index_)=new_code;
    end
    Phase_microstructure = tmp; % Assign
    clear tmp % Clean temporary variable
    if OPTIONS.displaytext
        disp 'Phase code has been re-assigned';  disp ' ';
        disp(Table_phase);disp ' ';
    end
else
    if OPTIONS.displaytext
        disp 'Phase code has NOT been re-assigned';  disp ' ';
    end
end

%% IMAGE RESOLUTION
Table_voxel_RESIZE = Table_voxel_ROI; % Initialize
Table_domainsize = Table_Region_of_Interest; % Initialize
if INFO.initial_voxelsize ~= INFO.asked_voxelsize
    [Phase_microstructure] = function_scale_array(Phase_microstructure, INFO.initial_voxelsize, INFO.asked_voxelsize, INFO.phaseinfo);
    % Update informations
    Domain_size = size(Phase_microstructure);
    Domain_size_um = Domain_size*INFO.asked_voxelsize/1000; % in micrometers
    Table_voxel_RESIZE = table([{'Voxel size (micrometers)'};{'Number of voxel'}], [INFO.asked_voxelsize; prod(Domain_size)],...
        'VariableNames',{'Information', 'Text'});
    if number_direction==2
        equivalent_square_size_um = prod(Domain_size_um)^(1/2);
        Table_domainsize = table([1;2],[{INFO.direction(1).consistentname}; {INFO.direction(2).consistentname}], [Domain_size(1); Domain_size(2)],[Domain_size_um(1); Domain_size_um(2)], [equivalent_square_size_um;equivalent_square_size_um],...
            'VariableNames',{'Direction', 'Name', 'Number_of_voxel' 'Length_micrometers' 'Equivalent_square_length_micrometers'});
    else
        equivalent_cubic_size_um = prod(Domain_size_um)^(1/3);
        Table_domainsize = table([1;2;3],[{INFO.direction(1).consistentname}; {INFO.direction(2).consistentname}; {INFO.direction(3).consistentname}], [Domain_size(1); Domain_size(2); Domain_size(3)],[Domain_size_um(1); Domain_size_um(2); Domain_size_um(3)], [equivalent_cubic_size_um; equivalent_cubic_size_um; equivalent_cubic_size_um],...
            'VariableNames',{'Direction', 'Name', 'Number_of_voxel' 'Length_micrometers' 'Equivalent_square_length_micrometers'});
    end
    if OPTIONS.displaytext
        disp 'Voxel size has been changed';  disp ' ';
        disp(Table_domainsize);disp ' ';
        disp(Table_voxel_RESIZE);disp ' ';
    end
else
    if OPTIONS.displaytext
        disp 'Voxel size has NOT been changed';  disp ' ';
    end
end
tmp=whos('Phase_microstructure'); data_voxelresizeMB=tmp.bytes*9.53674e-7; % Keep track of data size
voxelresize_numbervoxel = prod(Domain_size);

%% MEMORY
Table_memory = table([{'Loaded data'};{'Region of interest'};{'Voxel resize'}],[initial_numbervoxel;cropped_numbervoxel;voxelresize_numbervoxel], [data_MB;data_cropMB;data_voxelresizeMB],...
    'VariableNames',{'Microstructure_array','Number_of_voxel', 'Memory_MB'});
if OPTIONS.displaytext
    disp 'Memory';  disp ' ';
    disp(Table_memory);disp ' ';
end

%% SAVE INVESTIGATED VOLUME
function_save_tif(Phase_microstructure, [volume_folder 'Investigated_volume.tif'])

%% SAVE INFORMATION IN EXCEL SHEET
if OPTIONS.save_xls
    % Intialization
    sheet_number=0;
    clear DATA_writetable
    % Region of interest: domain size
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='ROI_domainsize';
    DATA_writetable.sheet(sheet_number).table=Table_Region_of_Interest;
    % Region of interest: Voxel
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='ROI_voxel';
    DATA_writetable.sheet(sheet_number).table=Table_voxel_ROI;    
    % Phase name
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Reassign_phase';
    DATA_writetable.sheet(sheet_number).table=Table_phase;
    % Voxel resize: domain size
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Voxelresize_domainsize';
    DATA_writetable.sheet(sheet_number).table=Table_domainsize;
    % Voxel resize: Voxel
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Voxelresize_voxel';
    DATA_writetable.sheet(sheet_number).table=Table_voxel_RESIZE;       
    % Size
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Memory_size';
    DATA_writetable.sheet(sheet_number).table=Table_memory;  
    
    % Filename without extension
    filename = 'Volume_setup';
    % Save function
    Function_Writetable(main_folder,filename,DATA_writetable)
end

%%
%% CALCULATIONS
%%

% Update common info
Domain_size = size(Phase_microstructure);
INFO.number_phase = length(INFO.phase); % Number of phase
INFO.Domain_size = Domain_size; % Domain size
INFO.voxel_number = prod(Domain_size); % Number of voxel
INFO.voxel_size = INFO.asked_voxelsize; % Voxel size nm

if OPTIONS.displaytext
    disp '********************';
    disp '>>> CALCULATIONS <<<';
    disp '********************';
    disp ' ';
end

%% VOLUME FRACTIONS
if PROPERTY.volumefractions.todo
    Function_Volume_fractions(Phase_microstructure, PROPERTY, OPTIONS, INFO); % Call function
end

%% TORTUOSITY (tau factor)
if PROPERTY.tortuosity_taufactor.todo
    Function_Tortuosity_factor_taufactor(Phase_microstructure, PROPERTY, OPTIONS, INFO); % Call function
end

%% SPECIFIC SURFACE AREA (counting surfaces)
if PROPERTY.specificsurfacearea_directmethod.todo
    Function_Specificsurface_direct(Phase_microstructure, PROPERTY, OPTIONS, INFO); % Call function
end

%% SPECIFIC INTERFACE AREA (counting surfaces)
if INFO.number_phase>2
    if PROPERTY.specificinterfacearea_directmethod.todo
        Function_Specificinterface_direct(Phase_microstructure, PROPERTY, OPTIONS, INFO); % Call function
    end
end

%% PARTICLE SIZE (Continuum Particle-size distribution C-PSD)
if PROPERTY.particlesize_cpsd.todo
    Function_particle_size_CPSD(Phase_microstructure, PROPERTY, OPTIONS, INFO); % Call function
end

%% PARTICLE SIZE (distance map cumulative function fitting)
if PROPERTY.particlesize_dmap.todo
    Function_particle_size_distancemap(Phase_microstructure, PROPERTY, OPTIONS, INFO); % Call function
end

%% CONNECTIVITY
if PROPERTY.connectivity.todo
    Function_connectivity(Phase_microstructure, PROPERTY, OPTIONS, INFO); % Call function
end 

%% PARTICLE SIZE (discrete particle size distribution, based on watershed)
if PROPERTY.particlesize_watershed.todo
    Function_Discrete_particle_size_watershed_immersion(Phase_microstructure, PROPERTY, OPTIONS, INFO)
end

%% PARTICLE SIZE (discrete particle size distribution, based on Pseudo Coulomb Repulsive Field, PCRF)
if PROPERTY.particlesize_PCRF.todo
    Function_Discrete_particle_size_PCRF(Phase_microstructure, PROPERTY, OPTIONS, INFO)
end


%%
%% SUMMARY
%%

%% TIME
summary_folder = [main_folder 'Summary\']; 
if ~exist(summary_folder,'dir') % Check existence of the main folder
    mkdir(summary_folder); % Create folder if not exist
end

time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
time_stopwatch_elapsed = toc(tStart); % Stopwatch elapsed time
date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
% time table
Table_time = table({'Date start';'Date end';'Elasped time (hh:mm:ss)';'Elasped time (s), tic-toc';'CPU time (s)'} ,{char(date_start);char(date_end);char(date_end-date_start);num2str(time_stopwatch_elapsed,'%1.1f');num2str(time_cpu_elapsed,'%1.1f')},...
    'VariableNames',{'Information' 'text'});
if OPTIONS.displaytext==true
    disp ' ';
    disp '***********************';
    disp '>>> RESULTS SUMMARY <<<';
    disp '***********************';
    disp ' ';
    disp(Table_time)
end
if OPTIONS.save_xls==true
    filename = 'Time'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    DATA_writetable.sheet(1).name='Time';
    DATA_writetable.sheet(1).table=Table_time;
    % Save function
    Function_Writetable(summary_folder,filename,DATA_writetable)
end


%% SUMMARIZE RESULTS
% if OPTIONS.save_resultsmat                
%     function_summary_onevolume(summary_folder, INFO, OPTIONS);
% end


%% CLOSE REPORT
% if OPTIONS.report.todo == 1 && is_installed_reportgenerator==1
%       close(OPTIONS.reportdoc)
% end

end



