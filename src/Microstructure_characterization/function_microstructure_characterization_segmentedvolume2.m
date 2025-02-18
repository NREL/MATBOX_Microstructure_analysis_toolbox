function [] = function_microstructure_characterization_segmentedvolume2(mode, segmentation_type, infovol, optc, opts_sem,opts_ins,tag,tag_unit)
% Calculate microstructure properties useful for energy storage macroscale models.
% For instance, for Pseudo-2D battery models (i.e., Doyle–Fuller–Newman models).
% Intended use of this function is through the GUI of src\Microstructure_characterization\Microstructure_characterization.mlapp
% Inputs:
%   infovol  : volume information
%   optc     : options common to all volumes
%   opts_sem : options volume specific for semantic segmentation file
%   opts_ins : options volume specific for instance segmentation file

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
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    disp(['Loading semantic segmentation file: ' infovol.loadingpath_semantic]);
    disp '   Please wait...';
    [M_semantic,outcome] = function_load_tif(infovol.loadingpath_semantic,'uint8');
    if ~outcome % fail to import
        return
    end
    disp(['   Import successful ! Format: ' class(M_semantic)]);
    Domain_size = size(M_semantic); % Size of the loaded domain (number of voxel)

    M_nanoporosity = [];
    M_bulkdiffusivity = [];
    if strcmp(infovol.scale,'Dual scale (heterogenous)')
        if ~isempty(infovol.loadingpath_nanoporosity)
            disp(['Loading nanoporosity file: ' infovol.loadingpath_nanoporosity]);
            disp '   Please wait...';
            [M_nanoporosity,outcome] = function_load_tif(infovol.loadingpath_nanoporosity,'single');
            if ~outcome % fail to import
                return
            end
            if max(max(max(M_nanoporosity)))>1
                M_nanoporosity = M_nanoporosity / max(max(max(M_nanoporosity)));
            end
            disp(['   Import successful ! Format: ' class(M_nanoporosity)]);
        end

        if ~isempty(infovol.loadingpath_bulkdiffusivity)
            disp(['Loading bulk diffusivity file: ' infovol.loadingpath_bulkdiffusivity]);
            disp '   Please wait...';
            [M_bulkdiffusivity,outcome] = function_load_tif(infovol.loadingpath_bulkdiffusivity,'single');
            if ~outcome % fail to import
                return
            end
            if max(max(max(M_bulkdiffusivity)))>1
                M_bulkdiffusivity = M_bulkdiffusivity / max(max(max(M_bulkdiffusivity)));
            end
            disp(['   Import successful ! Format: ' class(M_bulkdiffusivity)]);
        end
    end
end

if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    disp(['Loading instance segmentation file: ' infovol.loadingpath_instance]);
    disp '   Please wait...';
    [M_instance,outcome] = function_load_tif(infovol.loadingpath_instance,'uint32');
    if ~outcome % fail to import
        return
    end
    [M_instance] = fct_intconvert(M_instance);
    disp(['   Import successful ! Format: ' class(M_instance)]);
    Domain_size = size(M_instance); % Size of the loaded domain (number of voxel)
end

%% OVERWRITTE IF infovol NOT RELEVANT
if ~strcmp(mode,'Set parameters for each volume independently')
    % Dimension
    number_direction = length(Domain_size);
    if number_direction<length(infovol.directionname)
        infovol.directionname = infovol.directionname(1:number_direction);
    elseif number_direction>length(infovol.directionname)
        for k = length(infovol.directionname)+1:1:number_direction
            infovol.directionname(k) = {['Direction ' num2str(k)]};
        end
    end

    if (strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')) && (isempty(infovol.nphase) || ~isfield(infovol,'phasename_semantic'))
        unis = unique(M_semantic);
        infoval.nphase = length(unis);

        if ~isfield(infovol,'phasename_semantic')
            infovol.phasename_semantic = infovol.initial_phasename_semantic;
            infovol.phaselabel_semantic = infovol.initial_phaselabel_semantic(:,1);
        end

        if infoval.nphase < length(infovol.phasename_semantic)
            infovol.initial_phasename_semantic = infovol.initial_phasename_semantic(1:infoval.nphase);
            infovol.phasename_semantic = infovol.phasename_semantic(1:infoval.nphase);
            infovol.initial_phaselabel_semantic = infovol.initial_phaselabel_semantic(1:infoval.nphase,:);
            infovol.phaselabel_semantic = infovol.phaselabel_semantic(1:infoval.nphase);
            infovol.phasecolor_semantic = infovol.phasecolor_semantic(1:infoval.nphase,:);
        elseif infoval.nphase > length(infovol.phasename_semantic)
            for k = length(infovol.phasename_semantic)+1:1:infoval.nphase
                infovol.initial_phasename_semantic(k) = {['Label ' num2str(unis(k))]};
                infovol.phasename_semantic(k) = {['Label ' num2str(unis(k))]};
                infovol.initial_phaselabel_semantic(k,:) = [unis(k) unis(k)];
                infovol.phaselabel_semantic(k) = unis(k);
                if k<7
                    tmp = colororder;
                    infovol.phasecolor_semantic(k,:) = tmp(k,:);
                else
                    infovol.phasecolor_semantic(k,:) = rand(1,3);
                end
            end
        end

    elseif (strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')) && isempty(infovol.ninstance)
        % No issue here: it is always a binary image
        foo = 1;
    end
end

%% FOLDER
if ~exist(infovol.volpath,'dir') % Check existence of the main folder
    mkdir(infovol.volpath); % Create folder if not exist
end

%% SAVE LOADED VOLUME
volume_folder = fullfile(infovol.volpath,'Volume');
if ~exist(volume_folder,'dir')
    mkdir(volume_folder);
end
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    function_save_tif(M_semantic, fullfile(volume_folder,'loaded_semanticvolume.tif'));
    if strcmp(infovol.scale,'Dual scale (heterogenous)')
        if ~isempty(infovol.loadingpath_nanoporosity)
            function_save_tif(M_nanoporosity, fullfile(volume_folder,'loaded_nanoporosity.tif'));
            save(fullfile(volume_folder,'loaded_nanoporosity.mat'),'M_nanoporosity');
        end
        if ~isempty(infovol.loadingpath_bulkdiffusivity)
            function_save_tif(M_bulkdiffusivity, fullfile(volume_folder,'loaded_bulkdiffusivity.tif'));
            save(fullfile(volume_folder,'loaded_bulkdiffusivity.mat'),'M_bulkdiffusivity');
        end
    end
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    function_save_tif(M_instance, fullfile(volume_folder,'loaded_instancevolume.tif'));
end

%% GENERAL INFORMATION ABOUT THE VOLUME
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
if strcmp(segmentation_type,'semantic')
    tmp=whos('M_semantic'); data_MB_sem=tmp.bytes*9.53674e-7;
    c1 = [{'File loaded'};{'Array size (MB)'}];
    c2 = [{infovol.loadingpath_semantic}; num2str(data_MB_sem,'%1.1f')];
    if strcmp(infovol.scale,'Dual scale (heterogenous)')
        if ~isempty(infovol.loadingpath_nanoporosity)
            c1(end+1) = {'Nanoporosity'};
            c2(end+1) = {infovol.loadingpath_nanoporosity};
        end
        if ~isempty(infovol.loadingpath_bulkdiffusivity)
            c1(end+1) = {'Bulk diffusivity'};
            c2(end+1) = {infovol.loadingpath_bulkdiffusivity};
        end
    end
    c1(end+1) = {'Save folder'};
    c2(end+1) = {infovol.volpath};
    Table_inputoutput = table(c1, c2, 'VariableNames',{'Information', 'Text'});

elseif strcmp(segmentation_type,'instance')
    tmp=whos('M_instance'); data_MB_ins=tmp.bytes*9.53674e-7;
    Table_inputoutput = table([{'File loaded'};{'Array size (MB)'};{'Save folder'}], [{infovol.loadingpath_instance}; num2str(data_MB_ins,'%1.1f'); {infovol.volpath}],...
        'VariableNames',{'Information', 'Text'});    

elseif strcmp(segmentation_type,'semantic and instance')
    tmp=whos('M_semantic'); data_MB_sem=tmp.bytes*9.53674e-7;
    tmp=whos('M_instance'); data_MB_ins=tmp.bytes*9.53674e-7;
    c1 = [{'File loaded (semantic)'};{'Array size (MB)'}];
    c2 = [{infovol.loadingpath_semantic}; num2str(data_MB_sem,'%1.1f')];
    if strcmp(infovol.scale,'Dual scale (heterogenous)')
        if ~isempty(infovol.loadingpath_nanoporosity)
            c1(end+1) = {'Nanoporosity'};
            c2(end+1) = {infovol.loadingpath_nanoporosity};
        end
        if ~isempty(infovol.loadingpath_bulkdiffusivity)
            c1(end+1) = {'Bulk diffusivity'};
            c2(end+1) = {infovol.loadingpath_bulkdiffusivity};
        end
    end
    c1(end+1) = {'File loaded (instance)'};
    c2(end+1) = {infovol.loadingpath_instance};
    c1(end+1) = {'Array size (MB)'};
    c2(end+1) = num2str(data_MB_ins,'%1.1f');    
    c1(end+1) = {'Save folder'};
    c2(end+1) = {infovol.volpath};
    Table_inputoutput = table(c1, c2, 'VariableNames',{'Information', 'Text'});
end

Table_nanoporosity = [];
Table_bulkdiffusivity = [];
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    if strcmp(infovol.scale,'Dual scale (heterogenous)')
        labels = infovol.initial_phaselabel_semantic;
        if ~isempty(infovol.loadingpath_nanoporosity)
            for k = 1:length(labels)
                idx = find(M_semantic == infovol.initial_phaselabel_semantic(k));
                avg_nano(k) = double(mean( M_nanoporosity(idx) ));
                std_nano(k) = double(std( M_nanoporosity(idx) ));
                min_nano(k) = double(min( M_nanoporosity(idx) ));
                max_nano(k) = double(max( M_nanoporosity(idx) ));
            end
            Table_nanoporosity = table(infovol.initial_phasename_semantic,infovol.initial_phaselabel_semantic,min_nano',max_nano',std_nano',avg_nano',...
                'VariableNames',{'Phase' 'Label' 'Minimum' 'Maximum' 'Standard deviation' 'Average'});%
        end
        if ~isempty(infovol.loadingpath_bulkdiffusivity)
            for k = 1:length(labels)
                idx = find(M_semantic == infovol.initial_phaselabel_semantic(k));
                avg_bulk(k) = double(mean( M_nanoporosity(idx) ));
                std_bulk(k) = double(std( M_nanoporosity(idx) ));
                min_bulk(k) = double(min( M_nanoporosity(idx) ));
                max_bulk(k) = double(max( M_nanoporosity(idx) ));
            end
            Table_bulkdiffusivity = table(infovol.initial_phasename_semantic,infovol.initial_phaselabel_semantic,min_bulk',max_bulk',std_bulk',avg_bulk',...
                'VariableNames',{'Phase' 'Label' 'Minimum' 'Maximum' 'Standard deviation' 'Average'});%
        end

    elseif strcmp(infovol.scale,'Dual scale (uniform per phase)')
            Table_nanoporosity = table(infovol.initial_phasename_semantic,infovol.initial_phaselabel_semantic,infovol.nanoporosity,...
                'VariableNames',{'Phase' 'Label' 'Nanoporosity'});%        
    end
end

if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    Table_phase_semantic = table(infovol.initial_phasename_semantic,infovol.initial_phaselabel_semantic(:,1),...
        'VariableNames',{'Phase name (semantic)' 'Phase label'});
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    Table_phase_BWinstance = table(infovol.initial_phasename_instance,infovol.initial_phaselabel_instance(:,1),...
        'VariableNames',{'Phase name (binary instance)' 'Phase label'});
    Table_phase_instance = table([{'Min label'};{'Max label'};{'Number of labels'}],[min(min(min(M_instance)));max(max(max(M_instance)));infovol.ninstance],...
        'VariableNames',{'Phase info (instance)' 'Value'});    
end

Table_volumeinformation = table(infovol.volumeinformation(:,1),infovol.volumeinformation(:,2),...
    'VariableNames',{'Information' 'Text'});

disp ' ';
disp(Table_inputoutput);disp ' ';

disp(Table_volumeinformation);disp ' ';

if infovol.isbackground
    disp('First label is background and will be excluded from the analysis');
end

if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    disp(Table_phase_semantic);disp ' ';
    if strcmp(infovol.scale,'Dual scale (heterogenous)') || strcmp(infovol.scale,'Dual scale (uniform per phase)')
        disp(['    Semantic volume lower scale representation: ' infovol.scale])
        if ~isempty(Table_nanoporosity)
            disp('       - nanoporosity');
            disp(Table_nanoporosity);disp ' ';
        end
        if ~isempty(Table_bulkdiffusivity)
            disp('       - bulk diffusivity');
            disp(Table_bulkdiffusivity);disp ' ';
        end
    end
end

if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    disp(Table_phase_instance);disp ' ';
    disp(Table_phase_BWinstance);disp ' ';
end
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
    % Volume information
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='Volume_information';
    DATA_writetable.sheet(sheet_number).table=Table_volumeinformation;
    % Phase name
    if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
        sheet_number=sheet_number+1;
        DATA_writetable.sheet(sheet_number).name='Phase (semantic)';
        DATA_writetable.sheet(sheet_number).table=Table_phase_semantic;
    end
    if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
        sheet_number=sheet_number+1;
        DATA_writetable.sheet(sheet_number).name='Phase (instance)';
        DATA_writetable.sheet(sheet_number).table=Table_phase_instance;
        sheet_number=sheet_number+1;
        DATA_writetable.sheet(sheet_number).name='Phase (binary instance)';
        DATA_writetable.sheet(sheet_number).table=Table_phase_BWinstance;
    end
    if ~isempty(Table_nanoporosity)
        sheet_number=sheet_number+1;
        DATA_writetable.sheet(sheet_number).name='Nanoporosity';
        DATA_writetable.sheet(sheet_number).table=Table_nanoporosity;
    end
    if ~isempty(Table_bulkdiffusivity)
        sheet_number=sheet_number+1;
        DATA_writetable.sheet(sheet_number).name='Bulk diffusivity';
        DATA_writetable.sheet(sheet_number).table=Table_bulkdiffusivity;
    end
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
    DATA_writetable.sheet(sheet_number).name='MATBOX version';
    DATA_writetable.sheet(sheet_number).table=Table_toolbox;
    % Author information
    sheet_number=sheet_number+1;
    DATA_writetable.sheet(sheet_number).name='MATBOX authors';
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
if strcmp(mode,'Set parameters for each volume independently')
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
        if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
            if number_direction==2
                M_semantic = M_semantic( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2));
                if strcmp(infovol.scale,'Dual scale (heterogenous)') || strcmp(infovol.scale,'Dual scale (uniform per phase)')
                    if ~isempty(M_nanoporosity)
                        M_nanoporosity = M_nanoporosity( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2));
                    end
                    if ~isempty(M_bulkdiffusivity)
                        M_bulkdiffusivity = M_bulkdiffusivity( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2));
                    end
                end
            else
                M_semantic = M_semantic( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2) , infovol.ROI(3,1):infovol.ROI(3,2));
                if strcmp(infovol.scale,'Dual scale (heterogenous)') || strcmp(infovol.scale,'Dual scale (uniform per phase)')
                    if ~isempty(M_nanoporosity)
                        M_nanoporosity = M_nanoporosity( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2) , infovol.ROI(3,1):infovol.ROI(3,2));
                    end
                    if ~isempty(M_bulkdiffusivity)
                        M_bulkdiffusivity = M_bulkdiffusivity( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2) , infovol.ROI(3,1):infovol.ROI(3,2));
                    end
                end                
            end
            Domain_size = size(M_semantic);
        end
        if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
            if number_direction==2
                M_instance = M_instance( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2));
            else
                M_instance = M_instance( infovol.ROI(1,1):infovol.ROI(1,2) , infovol.ROI(2,1):infovol.ROI(2,2) , infovol.ROI(3,1):infovol.ROI(3,2));
            end
            Domain_size = size(M_instance);
        end        

        Table_voxel_ROI = table([{['Voxel size (' infovol.unit ')']};{'Number of voxel'}], [infovol.initial_voxelsize; prod(Domain_size)],...
            'VariableNames',{'Information', 'Value'});

        disp 'Volume has been cropped'; disp ' ';
        disp(Table_domainsize_ROI); disp ' ';
        disp(Table_voxel_ROI);disp ' ';
    else
        disp 'Volume has NOT been crooped';  disp ' ';
        Table_voxel_ROI = Table_voxel;
        Table_domainsize_ROI = Table_domainsize;
    end
else
    disp 'Volume has NOT been crooped';  disp ' ';
    Table_voxel_ROI = Table_voxel;
    Table_domainsize_ROI = Table_domainsize;
end

if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    tmp=whos('M_semantic'); data_cropMB_sem=tmp.bytes*9.53674e-7;
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    tmp=whos('M_instance'); data_cropMB_ins=tmp.bytes*9.53674e-7;
end
cropped_numbervoxel = prod(Domain_size);


%% RE-LABEL PHASES
Table_phase_relabel = [];
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    if ~isempty(infovol.nphase_reassign)
        Table_phase_relabel = table(infovol.phasename_semantic,infovol.phaselabel_semantic,...
            'VariableNames',{'Phase name' 'Phase label'});
        if sum(infovol.initial_phaselabel_semantic(:,1)~=infovol.initial_phaselabel_semantic(:,2)) % At least one phase is re-assigned
            number_initialphase = length(infovol.initial_phasename_semantic); % Number of assigned phase
            tmp = M_semantic; % Temporary variable
            for current_phase=1:1:number_initialphase
                old_label = infovol.initial_phaselabel_semantic(current_phase,1);
                new_label = infovol.initial_phaselabel_semantic(current_phase,2);
                index_ = find(M_semantic==old_label);
                tmp(index_)=new_label;
            end
            M_semantic = tmp; % Assign
            clear tmp % Clean temporary variable
            disp 'Phase label has been re-assigned';  disp ' ';
            disp(Table_phase_relabel);disp ' ';
        else
            disp 'Phase label has NOT been re-assigned';  disp ' ';
        end
        M_semantic=uint8(M_semantic);
    end
end

%% SCALING (IMAGE RESOLUTION)
if infovol.scaling_factor ~= 1.0
    % Set parameters
    parameters_scaling.scaling_factor = infovol.scaling_factor;
    parameters_scaling.label_or_greylevel = 'Label';
    % Scale
    if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
        parameters_scaling.background = min(infovol.phaselabel_semantic);
        M_semantic = function_scaling(M_semantic,parameters_scaling);
        Domain_size = size(M_semantic);
        if ~isempty(M_nanoporosity)
            parameters_scaling.label_or_greylevel = 'Grey level';
            M_nanoporosity = function_scaling(M_nanoporosity,parameters_scaling);
        end
        if ~isempty(M_bulkdiffusivity)
            parameters_scaling.label_or_greylevel = 'Grey level';
            M_bulkdiffusivity = function_scaling(M_bulkdiffusivity,parameters_scaling);
        end
    end
    if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
        parameters_scaling.label_or_greylevel = 'Label';
        parameters_scaling.background = min(infovol.phaselabel_instance);
        M_instance = function_scaling(M_instance,parameters_scaling);
        Domain_size = size(M_semantic);
    end
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

if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    tmp=whos('M_semantic'); data_voxelresizeMB_sem=tmp.bytes*9.53674e-7;
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    tmp=whos('M_instance'); data_voxelresizeMB_ins=tmp.bytes*9.53674e-7;
end
voxelresize_numbervoxel = prod(Domain_size);

%% MEMORY
if strcmp(segmentation_type,'semantic')
    Table_memory = table([{'Loaded data'};{'Region of interest'};{'Voxel resize'}],[initial_numbervoxel;cropped_numbervoxel;voxelresize_numbervoxel], [data_MB_sem;data_cropMB_sem;data_voxelresizeMB_sem],...
        'VariableNames',{'Microstructure array','Number of voxel', 'Memory MB'});
elseif strcmp(segmentation_type,'instance')
    Table_memory = table([{'Loaded data'};{'Region of interest'};{'Voxel resize'}],[initial_numbervoxel;cropped_numbervoxel;voxelresize_numbervoxel], [data_MB_ins;data_cropMB_ins;data_voxelresizeMB_ins],...
        'VariableNames',{'Microstructure array','Number of voxel', 'Memory MB'});
elseif strcmp(segmentation_type,'semantic and instance')
    Table_memory = table([{'Loaded data'};{'Region of interest'};{'Voxel resize'}],[initial_numbervoxel;cropped_numbervoxel;voxelresize_numbervoxel], [{[num2str(data_MB_sem) '/' num2str(data_MB_ins)]};{[num2str(data_cropMB_sem) '/' num2str(data_cropMB_ins)]};{[num2str(data_voxelresizeMB_sem) '/' num2str(data_voxelresizeMB_ins)]}],...
        'VariableNames',{'Microstructure array','Number of voxel', 'Memory MB (semantic/instance)'});
end
disp 'Memory';  disp ' ';
disp(Table_memory);disp ' ';


%% SAVE INVESTIGATED VOLUME
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    function_save_tif(M_semantic, fullfile(volume_folder,'Investigated_semanticvolume.tif'));
    if strcmp(infovol.scale,'Dual scale (heterogenous)')
        if ~isempty(M_nanoporosity)
            function_save_tif(M_nanoporosity, fullfile(volume_folder,'Investigated_nanoporosity.tif'));
            save(fullfile(volume_folder,'Investigated_nanoporosity.mat'),'M_nanoporosity');
        end
        if ~isempty(M_bulkdiffusivity)
            function_save_tif(M_bulkdiffusivity, fullfile(volume_folder,'Investigated_bulkdiffusivity.tif'));
            save(fullfile(volume_folder,'Investigated_bulkdiffusivity.mat'),'M_bulkdiffusivity');
        end
    end 
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    function_save_tif(M_instance, fullfile(volume_folder,'Investigated_instancevolume.tif'));
end

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
    if ~isempty(Table_phase_relabel)
        sheet_number=sheet_number+1;
        DATA_writetable.sheet(sheet_number).name='Relabel_phase';
        DATA_writetable.sheet(sheet_number).table=Table_phase_relabel;
    end
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

if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    infovol_sem = infovol;
    if strcmp(segmentation_type,'semantic and instance')
        infovol_sem.sub = 'semantic';
    else
        infovol_sem.sub = [];
    end        
    infovol_sem.phaselabel = infovol.phaselabel_semantic;
    infovol_sem.phasename = infovol.phasename_semantic;
    infovol_sem.phasecolor = infovol.phasecolor_semantic;
end

if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    infovol_ins = infovol;
    if strcmp(segmentation_type,'semantic and instance')
        infovol_ins.sub = 'Instance';
    else
        infovol_ins.sub = [];
    end    
    infovol_ins.phaselabel = infovol.phaselabel_instance;
    infovol_ins.phasename = infovol.phasename_instance;
    infovol_ins.phasecolor = infovol.phasecolor_instance;

    M_binaryinstance = M_instance;
    if infovol.isbackground
        M_binaryinstance(M_binaryinstance>1)=2;
    else        
        M_binaryinstance(M_binaryinstance~=0)=1;
    end
end


%% VOLUME FRACTIONS
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    if sum(opts_sem.volumefraction.todo)
        disp '    VOLUME FRACTIONS: semantic file';
        disp '    -------------------------------';
        disp ' ';
        Charact_Volumefractions(M_semantic, M_nanoporosity, infovol_sem, optc, opts_sem.volumefraction) % Call function
    end
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    if sum(opts_ins.volumefraction.todo)
        disp '    VOLUME FRACTIONS: Binary instance file';
        disp '    --------------------------------------';
        disp ' ';        
        Charact_Volumefractions(M_binaryinstance, [], infovol_ins, optc, opts_ins.volumefraction) % Call function
    end
end

%% TORTUOSITY
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    if sum(opts_sem.tortuosity.todo)
        disp '    TORTUOSITY FACTOR: semantic file';
        disp '    --------------------------------';
        disp ' ';
        Charact_Tortuosity(M_semantic, infovol_sem, optc, opts_sem.tortuosity) % Call function
    end
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    if sum(opts_ins.tortuosity.todo)
        disp '    TORTUOSITY FACTOR: Binary instance file';
        disp '    ---------------------------------------';
        disp ' ';        
        Charact_Tortuosity(M_binaryinstance, infovol_ins, optc, opts_ins.tortuosity) % Call function
    end
end

%% Continuous particle size distribution CPSD
if strcmp(segmentation_type,'semantic') || strcmp(segmentation_type,'semantic and instance')
    if sum(opts_sem.cpsd.todo)
        disp '    Continuous particle size distribution (CPSD): semantic file';
        disp '    -----------------------------------------------------------';
        disp ' ';
        Charact_Cpsd(M_semantic, infovol_sem, optc, opts_sem.cpsd, foo) % Call function
    end
end
if strcmp(segmentation_type,'instance') || strcmp(segmentation_type,'semantic and instance')
    if sum(opts_ins.cpsd.todo)
        disp '    Continuous particle size distribution (CPSD): Binary instance file';
        disp '    ------------------------------------------------------------------';
        disp ' ';        
        Charact_Cpsd(M_binaryinstance, infovol_ins, optc, opts_ins.cpsd, M_instance) % Call function
    end
end

%% METRICS PER PARTICLE
if ~isempty(opts_ins) && opts_ins.metricsperparticle.todo
    if strcmp(segmentation_type,'instance')
        disp '    Metrics per particle (instance only)';
        disp '    ------------------------------------';
        disp ' ';        
        Charact_Metricsperparticle(M_instance, [], infovol_ins, optc, opts_ins.metricsperparticle) % Call function
    elseif strcmp(segmentation_type,'semantic and instance')
        disp '    Metrics per particle (semantic and instance)';
        disp '    --------------------------------------------';
        disp ' ';           
        Charact_Metricsperparticle(M_instance, M_semantic, infovol_ins, optc, opts_ins.metricsperparticle) % Call function
    end
end


%% TAG
if ~isempty(tag)
    results_correlation(1).tag = tag;
    results_correlation(1).tag_unit = tag_unit;
    folder = fullfile(infovol.volpath, 'Correlation');
    if exist(folder,'dir')==0 % Folder existence is checked, and created if necessary
        mkdir(folder);
    end
    save(fullfile(folder, 'Correlation_tag.mat'),'results_correlation');
end



% 
% 
% %% CONNECTIVITY
% if opts.connectivity.todo
%      Function_connectivity(M_semantic, infovol, optc, opts.connectivity) % Call function
% end
% 
% %% TORTUOSITY FACTORS
% if sum(opts.tortuosity.todo)
%     Function_Tortuosity_factor_taufactor(M_semantic, infovol, optc, opts.tortuosity) % Call function
% end
% 
% %% SPECIFIC SURFACE AREA, DIRECT
% if sum(opts.Sp_direct.todo)
%     Function_Specificsurface_direct(M_semantic, infovol, optc, opts.Sp_direct, foo, foo) % Call function
% end
% 
% %% SPECIFIC INTERFACE AREA, DIRECT
% if opts.Int_direct.todo
%     Function_Specificinterface_direct(M_semantic, infovol, optc, opts.Int_direct, foo, foo) % Call function
% end
% 
% %% TRIPLE PHASE BOUNDARY LENGTH, DIRECT
% if opts.TPBL_direct.todo
%     Function_TPBL_direct(M_semantic, infovol, optc, opts.TPBL_direct, foo, foo) % Call function
% end
% 
% %% DIAMETER, C-PSD
% if sum(opts.D_cpsd.todo)
%     Function_particle_size_CPSD(M_semantic, infovol, optc, opts.D_cpsd, foo) % Call function
% end
% 
% %% DIAMETER, EDMF
% if sum(opts.D_edmf.todo)
%     Function_particle_size_EDMF(M_semantic, infovol, optc, opts.D_edmf) % Call function
% end
% 
% %% DIAMETER, CHORD
% if sum(opts.D_chord.todo)
%     Function_particle_size_CHORD(M_semantic, infovol, optc, opts.D_chord) % Call function
% end
% 
% %% DIAMETER, WATERSHED
% if sum(opts.D_watershed.todo)
%     Function_particle_size_Watershed(M_semantic, infovol, optc, opts.D_watershed) % Call function
% end
% 
% % Chord: harmonize with EDMF
% 
% % Tau
% 
% % Convergence(RVE), and RVE(RVE)
% 
% % Article on RVE
% 
% % Article on fractal/voxel size
% % I know what fractal dimension will bring: we can compare Sp or TPBL from two different materials, only if they
% % share same voxel size, and same fractal dimension (indicator of the surface roughness)
% 
% % Publish on Github (legacy and new version in choice menu)
% 
% % Connectivity
% 
% % Publish on Github
% 
% % Watershed
% % PCRF
% 
% % Publish on Github
% 
% % Tutorial video
% 
% % Correlation module
% 
% %%
% %% SUMMARY
% %%
% 
% %% TIME
% summary_folder = [infovol.volpath 'Summary' separator]; 
% if ~exist(summary_folder,'dir') % Check existence of the main folder
%     mkdir(summary_folder); % Create folder if not exist
% end
% 
% time_cpu_elapsed = cputime-time_cpu_start; % CPU elapsed time
% time_stopwatch_elapsed = toc(tStart); % Stopwatch elapsed time
% date_end = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
% % time table
% Table_time = table({'Date start';'Date end';'Elasped time (hh:mm:ss)';'Elasped time (s), tic-toc';'CPU time (s)'} ,{char(date_start);char(date_end);char(date_end-date_start);num2str(time_stopwatch_elapsed,'%1.1f');num2str(time_cpu_elapsed,'%1.1f')},...
%     'VariableNames',{'Information' 'text'});
% disp ' ';
% disp '***********************';
% disp '>>> RESULTS SUMMARY <<<';
% disp '***********************';
% disp ' ';
% disp(Table_time)
% if optc.save.xls
%     filename = 'Time'; % Filename without extension
%     % Prepare the data
%     clear DATA_writetable
%     DATA_writetable.sheet(1).name='Time';
%     DATA_writetable.sheet(1).table=Table_time;
%     % Save function
%     Function_Writetable(summary_folder,filename,DATA_writetable)
% end
% 
% %% SUMMARIZE RESULTS
% if optc.save.mat              
%     function_summary_onevolume(summary_folder, infovol, optc, opts);
% end
% 
% end