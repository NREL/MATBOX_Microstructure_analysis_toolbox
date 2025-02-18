function [results_main, results_correlation, timedata_perphase, timedata_pervolume] = ImageResolution_erroranalysis(pMET,Msem,Mins,voxel_size, number_domain, number_domain_todo, todo, scaling, scaling_extrapolation, Current_folder, opt, infovol, results_main, results_correlation, timedata_perphase, timedata_pervolume)

fprintf('> Image resolution analysis:\n\n');

size_choice = scaling;
size_choice = sort(size_choice);
number_resize=length(size_choice); % Number of different voxel size that will be analyzed

nMetric = length(pMET.metric);

if isempty(Mins)
    sz = size(Msem);
else
    sz = size(Mins);
end
dimension = length(sz);
if dimension == 2
    sz(3)=1;
end

%% CALCULATION
% Initialization
property_voxelsizedependence = zeros(number_resize+1,number_domain_todo+1,nMetric);
property_voxelsizedependence(1,1)=voxel_size;
for kMetric = 1:nMetric
    property_voxelsizedependence(1,2:end,kMetric)=pMET.metric(1).result_initialres;
end

% Loop on each voxel size
for current_iteration=1:1:number_resize

    % % Microstructure scaling
    % New voxel size
    property_voxelsizedependence(current_iteration+1,1,:)=size_choice(current_iteration)*voxel_size;
    % Set parameters
    parameters_scaling.scaling_factor = size_choice(current_iteration);
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = min(infovol.phaselabel);
    % Scale
    if ~isempty(Msem)
        Msem_resized = function_scaling(Msem,parameters_scaling);
        voxel_number_tmp=numel(Msem_resized);
    end
    if ~isempty(Mins)
        Mins_resized = function_scaling(Mins,parameters_scaling);
        voxel_number_tmp=numel(Mins_resized);
    end

    % CPU and stopwatch time - start
    time_cpu_start_volume = cputime; % CPU start
    time_clock_start_volume = tic; % Stopwatch start
 
    current_domain_todo = 0;
    for current_domain=1:1:number_domain % Loop over all phases
        if todo(current_domain)
            time_cpu_start_phase = cputime; % CPU start
            time_clock_start_phase = tic; % Stopwatch start
            current_domain_todo=current_domain_todo+1;

            % % Algorithm: SPECIFIC FOR EACH FILE
            if strcmp(pMET.fct_name,'Volume fraction')
                if infovol.isbackground
                    nvoxel_background = sum(sum(sum(Msem_resized==0)));
                    vf_background = nvoxel_background/voxel_number_tmp;
                end
                code_tmp = infovol.phaselabel(current_domain);
                Numbervoxel_domain_tmp= sum(sum(sum(Msem_resized==code_tmp )));
                property_voxelsizedependence(current_iteration+1,current_domain_todo+1,1)=Numbervoxel_domain_tmp/voxel_number_tmp;
                if infovol.isbackground
                    property_voxelsizedependence(current_iteration+1,current_domain_todo+1,1) = property_voxelsizedependence(current_iteration+1,current_domain_todo+1,1) / (1-vf_background);
                end
            end

            % Time
            timedata_perphase = [timedata_perphase; [Numbervoxel_domain_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];

        end
    end
    % CPU and stopwatch time - end
    timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];

end

% Sort per voxel size
for kMetric = 1:nMetric
    property_voxelsizedependence(:,:,kMetric) = sortrows(property_voxelsizedependence(:,:,kMetric),1);
end

%% EXTRAPOLATION TO 0 nm
if ~strcmp(scaling_extrapolation,'None')
    property_voxelsizedependence_new = zeros(number_resize+2,number_domain_todo+1,nMetric);
    interpolation_voxelsize_order = zeros(nMetric,1);
    for kMetric = 1:nMetric
        tmp = zeros(number_resize+2,number_domain_todo+1); % + 0 nm and + initial voxel size
        x=property_voxelsizedependence(:,1,kMetric);
        max_order = length(scaling)-1;

        if strcmp(scaling_extrapolation,'Best')

        else
            if strcmp(scaling_extrapolation,'Linear')
                interpolation_voxelsize_order(kMetric)=1;
            elseif strcmp(scaling_extrapolation,'Quadratic')
                interpolation_voxelsize_order(kMetric)=2;
            elseif strcmp(scaling_extrapolation,'Cubic')
                interpolation_voxelsize_order(kMetric)=3;
            end
            interpolation_voxelsize_order(kMetric) = min(interpolation_voxelsize_order,max_order);
            % fprintf('  %s: extrapolation to zero voxel size: polynomial of order %i\n\n',pMET.metric(kMetric).name,interpolation_voxelsize_order);
            for current_domain_todo=1:1:number_domain_todo
                y=property_voxelsizedependence(:,current_domain_todo+1,kMetric);
                pi = polyfit(x,y,interpolation_voxelsize_order);
                vq = polyval(pi,0);
                tmp(1,current_domain_todo+1)=vq;
                results_correlation(current_domain_todo).([pMET.metric(kMetric).shortname_correlation '_extrapolated']) = vq;  
                interpolation_voxelsize(kMetric).domain(current_domain_todo).pi=pi;
            end
        end
        tmp(2:end,:) = property_voxelsizedependence(:,:,kMetric);
        property_voxelsizedependence_new(:,:,kMetric) = tmp; clear tmp;
    end
    property_voxelsizedependence = property_voxelsizedependence_new;
    clear property_voxelsizedependence_new;
end

%% FRACTAL DIMENSION
if pMET.fractal_bertei
    if strcmp(pMET.fct_name,'Volume fraction')
        topology_dimension = ones(number_domain_todo,nMetric) * dimension;
    end
    fractal_dimension = zeros(number_domain_todo,nMetric);
    fractal_propensity = zeros(number_domain_todo,nMetric);

    resolution = property_voxelsizedependence(:,1,1);
    for kMetric = 1:nMetric
        for current_domain_todo=1:1:number_domain_todo
            values = property_voxelsizedependence(:,current_domain_todo+1,kMetric);
            [fractal_dimension(current_domain_todo,kMetric), fractal_propensity(current_domain_todo,kMetric), logs(current_domain_todo,kMetric).vals] = RichardsonMandelbrot_fractaldimension(topology_dimension(1),resolution,values);
            results_correlation(current_domain_todo).([pMET.metric(kMetric).shortname_correlation '_fractaldimension']) = fractal_dimension(current_domain_todo,kMetric);             
            results_correlation(current_domain_todo).([pMET.metric(kMetric).shortname_correlation '_fractalpropensity']) = fractal_propensity(current_domain_todo,kMetric);             
        end
    end
end

%% MANAGING RESULTS

current_domain_todo = 0;
for current_domain=1:1:number_domain % Loop over all phases
    if todo(current_domain)
        current_domain_todo = current_domain_todo + 1;
        domainname_todo(current_domain_todo,1) = infovol.phasename(current_domain,1);
    end
end

% Results are saved in a table
Variable_name_table={['Voxel size ' infovol.unit]}; % Columns name
for current_domain_todo=1:1:number_domain_todo
    Variable_name_table(1+current_domain_todo)=domainname_todo(current_domain_todo);
end
% Table
for kMetric = 1:nMetric
    voxelsizedependence(kMetric).table = array2table(property_voxelsizedependence(:,:,kMetric),...
        'VariableNames',Variable_name_table);
    if pMET.fractal_bertei
        Variable_name_table={'Voxel size log'}; % Columns name
        tmp = [];
        for current_domain_todo=1:1:number_domain_todo
            Variable_name_table(1+current_domain_todo)={[char(domainname_todo(current_domain_todo)) ' log']};
            tmp = [tmp logs(current_domain_todo,kMetric).vals(:,2)];
        end
        voxelsizedependence(kMetric).table_log = array2table([logs(current_domain_todo,kMetric).vals(:,1) tmp],'VariableNames',Variable_name_table);
        fractaldimension(kMetric).Table = table(domainname_todo(:,1),fractal_dimension(:,kMetric),topology_dimension(:,kMetric),fractal_propensity(:,kMetric),...
            'VariableNames',{'Phase','Fractal dimension','Topology dimension','Fractal propensity'});
    end
end

%% DISPLAY TEXT RESULTS
for kMetric = 1:nMetric
    if strcmp(scaling_extrapolation,'None')
        fprintf('  %s\n\n',pMET.metric(kMetric).name);
    else
        fprintf('  %s: extrapolation to zero voxel size: polynomial of order %i\n\n',pMET.metric(kMetric).name,interpolation_voxelsize_order);
    end
    disp(voxelsizedependence(kMetric).table)
    if pMET.fractal_bertei
        disp(voxelsizedependence(kMetric).table_log)
        str = ['    Richardson/Mandelbrot formula: Log(property) = m + (' num2str(topology_dimension(1,kMetric)) '-fractal dimension)*log(voxel size))'];
        disp(str);
        disp(fractaldimension(kMetric).Table)
    end
end

%% SAVE RESULTS
results_main.voxelsizedependence = voxelsizedependence; % Save in main table result
if pMET.fractal_bertei
    results_main.fractaldimension = fractaldimension;
end

if opt.save.xls
    for kMetric = 1:nMetric
        filename = [pMET.metric(kMetric).shortname_correlation '_voxelsize_dependence']; % Filename without extension
        % Prepare the data
        clear DATA_writetable
        DATA_writetable.sheet(1).name=pMET.metric(kMetric).name;
        DATA_writetable.sheet(1).table=voxelsizedependence(kMetric).table;
        if pMET.fractal_bertei
            DATA_writetable.sheet(2).name='Log log';
            DATA_writetable.sheet(2).table=voxelsizedependence(kMetric).table_log;
            DATA_writetable.sheet(3).name='Fractal dimension Mandelbrot';
            DATA_writetable.sheet(3).table=fractaldimension(kMetric).Table;
        end
        % Save function
        Function_Writetable(Current_folder,filename,DATA_writetable)
    end
end

%% FIGURES
for kMetric = 1:nMetric
    parameters_figure.propertyname = pMET.metric(kMetric).name;
    parameters_figure.method = [];
    parameters_figure.property_voxelsizedependence = property_voxelsizedependence(:,:,kMetric);
    parameters_figure.number_domain = number_domain_todo;
    parameters_figure.str_ylabel = pMET.metric(kMetric).name;
    parameters_figure.propertynameunit = pMET.metric(kMetric).unit;
    if ~strcmp(scaling_extrapolation,'None')
        parameters_figure.interpolation_voxelsize = interpolation_voxelsize(kMetric);
    else
        parameters_figure.interpolation_voxelsize = [];
    end
    parameters_figure.Current_folder = Current_folder;
    parameters_figure.filename = [pMET.metric(kMetric).shortname_correlation '_voxelsize_dependence'];
    parameters_figure.infovol = infovol;
    parameters_figure.opt = opt;
    parameters_figure.todo = todo;
    ImageResolution_erroranalysis_figures(parameters_figure) % Figures
end

end