function [results_main, results_correlation, timedata_perphase, timedata_pervolume] = RVE_main(pRVEs,pMET,Msem,Mins,voxel_size,number_domain,number_domain_todo,Current_folder,opt,infovol,results_main, results_correlation, timedata_perphase, timedata_pervolume)

% To do:
% Add: RVE_figures           : add violon plot
% Fix: RVE_convergence_Table : (not critical) Correct FOV_length and RVE_length (currently tables are not saved)
% Fix: Result_RVEconv        : (not critial) remove x axis length as already put inFOV_length instead 

if isempty(Mins)
    sz = size(Msem);
else
    sz = size(Mins);
end
dimension = length(sz);
if dimension == 2
    sz(3)=1;
end

nMetric = length(pMET.metric);

for k_RVE = 1:1:length(pRVEs) % Loop over all RVE analysis
    pRVE = pRVEs(k_RVE); % Select RVE parameters
    RVE(k_RVE).RVEparameters = pRVE; % For result structure
    RVE(k_RVE).attempted = false;
    
    %% CHECK IF ANALYSIS POSSIBLE
    if dimension == 2
        if strcmp(pRVE.type,'C') && strcmp(pRVE.Constantdirection,'Direction 3')
            warning('Impossible RVE analysis: domain is 2D but constant direction is set along axe 3!');
            continue
        elseif strcmp(pRVE.type,'D')
            warning('Impossible RVE analysis: domain is 2D but with two constant directions!');
            continue
        elseif strcmp(pRVE.type,'G') && (strcmp(pRVE.Growthdirection,'Direction 3, from z min to z max') || strcmp(pRVE.Growthdirection,'Direction 3, from z max to z min'))
            warning('Impossible convergence analysis: domain is 2D but growing direction is set along axe 3!');
            continue
        end
        if strcmp(pRVE.analysis,'Independent subvolumes') && pRVE.RVEconvergence_with_FOV
            if strcmp(pRVE.RVEconvergence_Crop,'1 direction') && (strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from z min to z max') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from z max to z min') || strcmp(pRVE.RVEconvergence_thatis,'Shrink along direction 3, from both ends'))
                warning('Impossible RVE=f(FOV) analysis: domain is 2D but shrinking direction is set along axe 3!');
                continue
            elseif strcmp(pRVE.RVEconvergence_Crop,'2 direction') && (strcmp(pRVE.RVEconvergence_thatis,'And keep direction 1 uncropped') || strcmp(pRVE.RVEconvergence_thatis,'And keep direction 2 uncropped'))
                warning('Impossible RVE=f(FOV) analysis: domain is 2D but shrinking direction is set along axe 3!');
                continue
            elseif strcmp(pRVE.RVEconvergence_Crop,'3 direction')
                warning('Impossible RVE=f(FOV) analysis: domain is 2D but shrinking direction is set along axe 3!');
                continue
            end
        end
    end
    RVE(k_RVE).attempted = true;

    %% CREATE FOLDER
    if opt.save.xls || opt.save.savefig
        Sub_folder_RVE = fullfile(Current_folder, pRVE.savename);
        v=0;
        while exist(Sub_folder_RVE,'dir')
            v=v+1;
            Sub_folder_RVE = fullfile(Current_folder, [pRVE.savename '_bis' num2str(v)]);
        end
        mkdir(Sub_folder_RVE);
    end

    %% SELECT THRESHOLDS
    if ischar(pRVE.threshold_std) % 'n/a'
        thresholds = pRVE.threshold_reldiff;
    else
        thresholds = pRVE.threshold_std;
    end
    n_threshold = length(thresholds);

    %% RVE CONVERGENCE ANALYSIS
    if strcmp(pRVE.analysis,'Independent subvolumes') && pRVE.RVEconvergence_with_FOV
        FOVbounds = RVE_convergence_FOVbounds(pRVE,sz);
        if ~isempty(FOVbounds)
            [nFOV,~] = size(FOVbounds); % Number of FOV
            length_FOV = zeros(nFOV+1,2);
            Result_RVEconv = zeros(nFOV+1,n_threshold+1,number_domain_todo,3,2,nMetric);
            % FOV size / number of threshold / phase / subdomain RVE or convergence size <, = , > /  size (both FOV and subdomain) in cubic or square root (=1), square root or length (=2)
        end
    end
    if exist('FOVbounds','var') && ~isempty(FOVbounds)
        pRVE.RVEconvergence_with_FOV = true; % Overwritte
    else
        pRVE.RVEconvergence_with_FOV = false;
        nFOV = 0;
    end

    %% LOOP OVER FOVs
    for k_FOV = 0:1:nFOV

        %% SELECT FOV
        if k_FOV>0
            x0 = FOVbounds(k_FOV,1); x1 = FOVbounds(k_FOV,2);
            y0 = FOVbounds(k_FOV,3); y1 = FOVbounds(k_FOV,4);
            z0 = FOVbounds(k_FOV,5); z1 = FOVbounds(k_FOV,6);
            if ~isempty(Msem)
                Msem_FOV = Msem(x0:x1,y0:y1,z0:z1);
                sz_FOV = size(Msem_FOV);
            end
            if ~isempty(Mins)
                Mins_FOV = Mins(x0:x1,y0:y1,z0:z1);
                sz_FOV = size(Mins_FOV);
            end
        else
            sz_FOV = sz;
        end
        FOVsize = (prod(sz_FOV))^(1/dimension) * voxel_size;
        if dimension == 3
            FOVlength_str = 'cubic root';
            if strcmp(pRVE.analysis,'Independent subvolumes')
                RVE(k_RVE).Representativity_analysis = 'Representative volume element (RVE)';
            else
                RVE(k_RVE).Representativity_analysis = 'Volume convergence';
            end
        else
            FOVlength_str = 'square root';
            if strcmp(pRVE.analysis,'Independent subvolumes')
                RVE(k_RVE).Representativity_analysis = 'Representative Section Area (RSA)';
            else
                RVE(k_RVE).Representativity_analysis = 'Section area convergence';
            end
        end
        RVE(k_RVE).FOVlength_str = FOVlength_str;

        %% SELECT SUBDOMAINS FOR THIS FOV
        [All_subdomain,GROUP_SUBDOMAIN, REP2ndsize, Representativity_2ndlength_str, Representativity_analysis_size2] = RVE_Subdomains(pRVE, sz_FOV, voxel_size); % Location of all subdomains
        if isempty(All_subdomain)
            warning('FOV is too small to satisfy the minimum subomain condition.')
            continue
        end

        RVE(k_RVE).Representativity_2ndlength_str = Representativity_2ndlength_str;
        RVE(k_RVE).Representativity_analysis_size2 = Representativity_analysis_size2;

        if pRVE.RVEconvergence_with_FOV
            [FOV2ndsize, FOV2ndsize_lengthstr] = RVE_2ndsize(pRVE,sz_FOV,voxel_size); % Get FOV second size

            length_FOV(k_FOV+1,1) = FOVsize;
            %Result_RVEconv(k_FOV+1,1,:,:,1,:) = FOVsize; % FOV size
            if k_FOV ==1
                fprintf('       Representativity convergence analysis\n');
            end
            if k_FOV > 0
                if isempty(FOV2ndsize)
                    fprintf('          Cropping iteration #%i/%i: FOV %s =%1.2f %s\n',k_FOV,nFOV,FOVlength_str,FOVsize,infovol.unit);
                else
                    fprintf('          Cropping iteration #%i/%i: FOV %s =%1.2f %s, FOV 2nd size %s =%1.2f %s\n',k_FOV,nFOV,FOVlength_str,FOVsize,infovol.unit,FOV2ndsize_lengthstr,FOV2ndsize,infovol.unit);
                end
            end
            if ~isempty(FOV2ndsize)
                length_FOV(k_FOV+1,2) = FOV2ndsize;
                %Result_RVEconv(k_FOV+1,1,:,:,2,:) = FOV2ndsize; % FOV second size
                nsize = 2;
            else
                nsize = 1;
            end
        end

        % Information about subdomains
        RVE(k_RVE).info = table(All_subdomain(:,1),All_subdomain(:,2),All_subdomain(:,3),All_subdomain(:,4),All_subdomain(:,5),All_subdomain(:,6),All_subdomain(:,7),All_subdomain(:,8),All_subdomain(:,9),All_subdomain(:,10),All_subdomain(:,11),All_subdomain(:,12),...
            'VariableNames',{'Subdomain Id' 'Group Id' 'Number subdomain' 'Cubic root' 'Square root' 'Length' 'x0' 'x1' 'y0' 'y1' 'z0' 'z1'});
        [number_subdomain,~] = size(All_subdomain); % The number of subdomain
        number_group_size = length(GROUP_SUBDOMAIN.id); % the number of group of subdomains sharing the same size

        %% APPLY ALGORITHM ON SUBDOMAINS
        % Colunm 1 is the subdomain id
        % Colunm 2-4 are the sizes of the subdomain.
        Property_eachsubdomain = zeros(number_subdomain, number_domain_todo+4, nMetric);

        % Property calculated for each subdomain
        for subdomain_id = 1:1:number_subdomain
            Property_eachsubdomain(subdomain_id,1,:)=subdomain_id;
            % Equivalent size of the subdomain
            Property_eachsubdomain(subdomain_id,2,:)=All_subdomain(subdomain_id,4)*voxel_size; % Cubic root length
            Property_eachsubdomain(subdomain_id,3,:)=All_subdomain(subdomain_id,5)*voxel_size; % Square root length
            Property_eachsubdomain(subdomain_id,4,:)=All_subdomain(subdomain_id,6)*voxel_size; % Length
            % Boundary of the subdomain
            x0 = All_subdomain(subdomain_id,7); x1 = All_subdomain(subdomain_id,8);
            y0 = All_subdomain(subdomain_id,9); y1 = All_subdomain(subdomain_id,10);
            z0 = All_subdomain(subdomain_id,11); z1 = All_subdomain(subdomain_id,12);
            % Crop FOV
            if k_FOV == 0
                if ~isempty(Msem)
                    current_subdomain_sem = Msem(x0:x1,y0:y1,z0:z1);
                    voxel_number_tmp=numel(current_subdomain_sem);
                end
                if ~isempty(Mins)
                    current_subdomain_ins = Mins(x0:x1,y0:y1,z0:z1);
                    voxel_number_tmp=numel(current_subdomain_ins);
                end
            else
                if ~isempty(Msem)
                    current_subdomain_sem = Msem_FOV(x0:x1,y0:y1,z0:z1);
                    voxel_number_tmp=numel(current_subdomain_sem);
                end
                if ~isempty(Mins)
                    current_subdomain_ins = Mins_FOV(x0:x1,y0:y1,z0:z1);
                    voxel_number_tmp=numel(current_subdomain_ins);
                end
            end

            % CPU and stopwatch time - start
            time_cpu_start_volume = cputime; % CPU start
            time_clock_start_volume = tic; % Stopwatch start

            current_domain_todo = 0;
            for current_domain=1:1:number_domain % Loop over all phases
                if pMET.p.todo(current_domain)
                    time_cpu_start_phase = cputime; % CPU start
                    time_clock_start_phase = tic; % Stopwatch start
                    current_domain_todo=current_domain_todo+1;

                    % % Algorithm: SPECIFIC FOR EACH FILE
                    if strcmp(pMET.fct_name,'Volume fraction')
                        if infovol.isbackground
                            nvoxel_background = sum(sum(sum(current_subdomain_sem==0)));
                            vf_background = nvoxel_background/voxel_number_tmp;
                        end
                        code_tmp = infovol.phaselabel(current_domain);
                        Numbervoxel_domain_tmp= sum(sum(sum(current_subdomain_sem==code_tmp )));
                        Property_eachsubdomain(subdomain_id,current_domain_todo+4)=Numbervoxel_domain_tmp/voxel_number_tmp;
                        if infovol.isbackground
                            if vf_background==1
                                Property_eachsubdomain(subdomain_id,current_domain_todo+4) = NaN;
                            else
                                Property_eachsubdomain(subdomain_id,current_domain_todo+4) = Property_eachsubdomain(subdomain_id,current_domain_todo+4) / (1-vf_background);
                            end
                        end
                    end

                    % % Time
                    timedata_perphase = [timedata_perphase; [Numbervoxel_domain_tmp (cputime-time_cpu_start_phase) toc(time_clock_start_phase)]];
                end
            end

            % CPU and stopwatch time - end
            timedata_pervolume = [timedata_pervolume; [voxel_number_tmp (cputime-time_cpu_start_volume) toc(time_clock_start_volume)]];
        end

        %% STATISTICAL ANALYSIS and RVE SIZE
        for kMetric=1:1:nMetric
            RVE(k_RVE).Metric(kMetric).name = pMET.metric(kMetric).name;
            [Property_subdomains_statistics, Size_RVE, relativedifference_convergence, Size_convergence] = RVE_statisticalanalysis(number_group_size,number_domain_todo,GROUP_SUBDOMAIN,Property_eachsubdomain(:,:,kMetric),voxel_size,pRVE,dimension);

            %% SAVE FOR CORRELATION
            if k_FOV == 0
                for k_threshold=1:1:n_threshold
                    current_domain_todo = 0;
                    for current_domain=1:1:number_domain
                        if pMET.p.todo(current_domain)
                            current_domain_todo=current_domain_todo+1;
                            if strcmp(pRVE.analysis,'Independent subvolumes')
                                if Size_RVE(k_threshold,current_domain_todo,2,1)~=0
                                    s1 = strrep(FOVlength_str, ' ', '');
                                    str = [pMET.metric(kMetric).shortname_correlation '_RVE_' pRVE.savename '_' s1 '_thresh' num2str(thresholds(k_threshold),'%1.2f')  ];
                                    str = strrep(str, '.', 'p');
                                    results_correlation(current_domain_todo).(str) = Size_RVE(k_threshold,current_domain_todo,2,1);
                                end
                                if ischar(Representativity_2ndlength_str) && Size_RVE(k_threshold,current_domain_todo,2,2)~=0
                                    s1 = strrep(Representativity_2ndlength_str, ' ', '');
                                    str = [pMET.metric(kMetric).shortname_correlation '_RVE_' pRVE.savename '_' s1 '_thresh' num2str(thresholds(k_threshold),'%1.2f')  ];
                                    str = strrep(str, '.', 'p');
                                    results_correlation(current_domain_todo).(str) = Size_RVE(k_threshold,current_domain_todo,2,2);
                                end

                            else
                                if Size_convergence(k_threshold,current_domain_todo,2,1)~=0
                                    s1 = strrep(FOVlength_str, ' ', '');
                                    str = [pMET.metric(kMetric).shortname_correlation '_conv_' pRVE.savename '_' s1 '_thresh' num2str(thresholds(k_threshold),'%1.2f')  ];
                                    str = strrep(str, '.', 'p');
                                    results_correlation(current_domain_todo).(str) = Size_convergence(k_threshold,current_domain_todo,2,1);
                                end
                                if ischar(Representativity_2ndlength_str) && Size_convergence(k_threshold,current_domain_todo,2,2)~=0
                                    s1 = strrep(Representativity_2ndlength_str, ' ', '');
                                    str = [pMET.metric(kMetric).shortname_correlation '_conv_' pRVE.savename '_' s1 '_thresh' num2str(thresholds(k_threshold),'%1.2f')  ];
                                    str = strrep(str, '.', 'p');
                                    results_correlation(current_domain_todo).(str) = Size_convergence(k_threshold,current_domain_todo,2,2);
                                end
                            end

                        end
                    end
                end
            end

            %% MANAGING RESULTS
            [Res] = RVE_tableresult(Property_eachsubdomain(:,:,kMetric), Property_subdomains_statistics, Size_RVE, relativedifference_convergence, Size_convergence, pRVE, number_domain,number_domain_todo,FOVlength_str,Representativity_2ndlength_str,infovol,pMET.p);
            RVE(k_RVE).Metric(kMetric).res = Res;

            if pRVE.RVEconvergence_with_FOV
                for k_threshold=1:1:n_threshold
                    current_domain_todo = 0;
                    for current_domain=1:1:number_domain
                        if pMET.p.todo(current_domain)
                            current_domain_todo=current_domain_todo+1;

                            for ksize = 1:nsize
                                if strcmp(pRVE.analysis,'Independent subvolumes')
                                    Result_RVEconv(k_FOV+1,k_threshold+1,current_domain_todo,1,ksize,kMetric) = Size_RVE(k_threshold,current_domain_todo,1,ksize);
                                    Result_RVEconv(k_FOV+1,k_threshold+1,current_domain_todo,2,ksize,kMetric) = Size_RVE(k_threshold,current_domain_todo,2,ksize);
                                    Result_RVEconv(k_FOV+1,k_threshold+1,current_domain_todo,3,ksize,kMetric) = Size_RVE(k_threshold,current_domain_todo,3,ksize);
                                else
                                    Result_RVEconv(k_FOV+1,k_threshold+1,current_domain_todo,1,ksize,kMetric) = Size_convergence(k_threshold,current_domain_todo,1,ksize);
                                    Result_RVEconv(k_FOV+1,k_threshold+1,current_domain_todo,2,ksize,kMetric) = Size_convergence(k_threshold,current_domain_todo,2,ksize);
                                    Result_RVEconv(k_FOV+1,k_threshold+1,current_domain_todo,3,ksize,kMetric) = Size_convergence(k_threshold,current_domain_todo,3,ksize);
                                end
                            end
                        end
                    end
                end
            end

            %% TEXT DISPLAY AND SAVE RESULTS
            if k_FOV == 0
                pRVE.disp_parRVE = true;
                RVE_displayandsave(RVE,k_RVE,kMetric,pRVE,number_domain,number_domain_todo,pMET.metric(kMetric).name,Sub_folder_RVE,opt,infovol,pMET.p);
            end

            %% FIGURES
            if k_FOV == 0
                parameters_figure.propertyname = pMET.metric(kMetric).name;
                parameters_figure.propertynameunit = pMET.metric(kMetric).unit;
                parameters_figure.RVE = pRVE;
                parameters_figure.Criterion=[pRVE.threshold_std 0];
                parameters_figure.savefolder = Sub_folder_RVE;
                parameters_figure.number_domain = number_domain;
                parameters_figure.number_domain_todo = number_domain_todo;
                parameters_figure.Property_subdomains_statistics = Property_subdomains_statistics;
                parameters_figure.Property_eachsubdomain = Property_eachsubdomain(:,:,kMetric);
                parameters_figure.relativedifference_convergence = relativedifference_convergence;
                parameters_figure.Size_RVE = Size_RVE;
                parameters_figure.convergence_criterion = pRVE.threshold_reldiff;
                parameters_figure.Size_convergence = Size_convergence;
                parameters_figure.Analysis_name = {RVE(k_RVE).Representativity_analysis, Representativity_analysis_size2};
                parameters_figure.Length_name = {FOVlength_str, Representativity_2ndlength_str};
                parameters_figure.Wholevolume_size = [FOVsize REP2ndsize];
                parameters_figure.Wholevolume_results = pMET.metric(kMetric).result_wholevolume;
                parameters_figure.infovol = infovol;
                parameters_figure.todo = pMET.p.todo;
                parameters_figure.opt = opt;
                RVE_figures(parameters_figure) % Figures
            end
        end       
    end

    %% CONVERGENCE ANALYSIS RESULT
    if pRVE.RVEconvergence_with_FOV
        for kMetric=1:1:nMetric
            % Table
            parameters_figure.nsize = nsize;
            parameters_figure.RVE_length_name = {FOVlength_str, Representativity_2ndlength_str};
            parameters_figure.FOV_length_name = {FOVlength_str, FOV2ndsize_lengthstr};
            parameters_figure.infovol = infovol;
            parameters_figure.thresholds = thresholds;
            parameters_figure.todo = pMET.p.todo;
            parameters_figure.propertyname = pMET.metric(kMetric).name;
            parameters_figure.Result_RVEconv = Result_RVEconv(:,:,:,:,:,kMetric);
            parameters_figure.length_FOV = length_FOV;
            parameters_figure.number_domain = number_domain;
            parameters_figure.opt = opt;
            parameters_figure.Sub_folder_RVE = Sub_folder_RVE;
            %[RVE(k_RVE).Metric(kMetric).res] = RVE_convergence_Table(parameters_figure ,RVE(k_RVE).Metric(kMetric).res);

            % Figure
            parameters_figure.Result_RVEconv = Result_RVEconv(:,:,:,:,:,kMetric);
            RVE_convergence_figures(parameters_figure) % Figures
            % Save   
            RVE(k_RVE).Metric(kMetric).res.RVEconvergence = Result_RVEconv(:,:,:,:,:,kMetric);
            RVE(k_RVE).Metric(kMetric).res.length_FOV = length_FOV;

        end
    end

end

results_main.RVE = RVE; % Save in main table result

end