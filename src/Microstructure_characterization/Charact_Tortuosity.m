function [] = Charact_Tortuosity(M,infovol,opt,p)

%% DEFAULT VALUES
% (M,label,direction)
% (M,lower_scale_diffusivity,direction)


%% DATE
date_start = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');

%% FOLDER
Current_folder = fullfile(infovol.volpath,'Tortuosity',infovol.sub);
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

%% VOLUME INFORMATION
sz = size(M);
dimension = length(sz);
if length(sz)==2
    sz = [sz 1];
end
number_phase = length(infovol.phaselabel); % Number of phase
voxel_number = prod(sz); % Number of voxel
voxel_size = infovol.voxelsize;
voxel_unit = infovol.unit;

%% INITIALIZE RESULTS (USE FOR CORRELATION)
if strcmp(p.scale,'One scale')
    current_phase_todo = 0;

    L = length(p.todo);
    if L < number_phase
        p.todo = [p.todo; ones(number_phase-L)];
    elseif L > number_phase
        p.todo = p.todo(1:number_phase);
    end

    if infovol.isbackground
        p.todo(1) = 0;
    end

    for current_phase=1:1:number_phase
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;
            results_correlation(current_phase_todo).name = infovol.phasename(current_phase,1);
        end
    end
else
    results_correlation.name = 'Full volume';
end


%%
%% ALGORITHM ON WHOLE VOLUME
%%

% Algorithm background parameters
p.isbackground = infovol.isbackground;
p.background_label = min(min(min(M)));
p.iterate_figure = true;
p.folder = Current_folder;

%% CALCULATION
if strcmp(p.scale,'One scale')
    number_phase_todo = sum(p.todo); % How many phase are we going to analyse ?
    directionname_todo = cell(number_phase_todo,dimension);    
    
    phaselabel = zeros(number_phase_todo,1);
    phasename_todo = cell(number_phase_todo,1);
    timedata = zeros(number_phase_todo,3);
    timedata_domain = cell(number_phase_todo,1);
   
    Deff = zeros(number_phase_todo,dimension);
    Mac = zeros(number_phase_todo,dimension);
    Tau = zeros(number_phase_todo,dimension);
    Bruggeman = zeros(number_phase_todo,dimension);
    vf = zeros(number_phase_todo,1);

    current_phase_todo = 0;
    for current_phase=1:1:number_phase % Loop over all phases
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;
            phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);
            number_direction_todo = sum(p.directions(current_phase).vals);


            idx = find(p.directions(current_phase).vals);
            direction_todo(current_phase_todo).name = infovol.directionname(idx);
            direction_todo(current_phase_todo).number = length(idx);
            direction_todo(current_phase_todo).idx = idx;

            current_direction_todo = 0;
            for dir = 1:dimension
                if p.directions(current_phase).vals(dir)
                    current_direction_todo = current_direction_todo + 1;
                    directionname_todo(current_phase_todo,current_direction_todo) = infovol.directionname(dir,1);

                    time_cpu_start_phase = cputime; % CPU start
                    time_clock_start_phase = tic; % Stopwatch start

                    % % Algorithm
                    p.direction = dir;
                    phaselabel(current_phase_todo,1) = infovol.phaselabel(current_phase);
                    p.onescale_label = phaselabel(current_phase_todo,1);
                    p.name = [char(phasename_todo(current_phase_todo,1)) '_direction' num2str(dir)];
                    [Deff(current_phase_todo,dir), Tau(current_phase_todo,dir), Bruggeman(current_phase_todo,dir), vf(current_phase_todo)] = Tortuosity_algorithm(M,p);
                    Mac(current_phase_todo,dir) = 1/Deff(current_phase_todo,dir);

                    % % Correlation
                    [str_direction] = function_remove_emptyandspecialcharacter_string(char(infovol.directionname(dir)));
                    results_correlation(current_phase_todo).('Vf') = vf(current_phase_todo);
                    results_correlation(current_phase_todo).(['Tortuosity_ax' num2str(dir) '_' str_direction]) = Tau(current_phase_todo,dir);
                    results_correlation(current_phase_todo).(['Normalized_diffusivity_ax' num2str(dir) '_' str_direction]) = Deff(current_phase_todo,dir);
                    results_correlation(current_phase_todo).(['MacMullin_number_ax' num2str(dir) '_' str_direction]) = Mac(current_phase_todo,dir);
                    results_correlation(current_phase_todo).(['Bruggeman_exponent_ax' num2str(dir) '_' str_direction]) = Bruggeman(current_phase_todo,dir);
                end
            end

            if number_direction_todo==3
                results_correlation(current_phase_todo).Tortuosity_ax1over2 = Tau(current_phase_todo,1)/Tau(current_phase_todo,2);
                results_correlation(current_phase_todo).Tortuosity_ax1over3 = Tau(current_phase_todo,1)/Tau(current_phase_todo,3);
                results_correlation(current_phase_todo).Tortuosity_ax2over3 = Tau(current_phase_todo,2)/Tau(current_phase_todo,3);
            end

            % % Time
            timedata_domain(current_phase_todo,1) = infovol.phasename(current_phase,1);
            timedata(current_phase_todo,1) = sum(sum(sum(M==phaselabel(current_phase_todo,1) )));
            timedata(current_phase_todo,2) = (cputime-time_cpu_start_phase)/number_direction_todo; % CPU elapsed time
            timedata(current_phase_todo,3) = (toc(time_clock_start_phase))/number_direction_todo; % CPU elapsed time            
        end
    end

elseif strcmp(p.scale,'Dual scale')
    time_cpu_start_volume = cputime; % CPU start
    time_clock_start_volume = tic; % Stopwatch start

    timedata = zeros(1,3); timedata(1,1) = voxel_number;
    timedata_domain = cell(1,1);
   
    labels = reshape(infovol.phaselabel,[1 length(infovol.phaselabel)]);
    if strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
        vals = reshape(p.lowerscale_porosity,[1 length(p.lowerscale_porosity)]);
        p.lowerscale_porosity = [double(labels); double(vals)];
        vals = reshape(p.lowerscale_Bruggeman,[1 length(p.lowerscale_Bruggeman)]);
        p.lowerscale_Bruggeman = [double(labels); double(vals)];
        p.lowerscale_diffusivity = p.lowerscale_porosity;
        p.lowerscale_diffusivity(2,:) = p.lowerscale_porosity(2,:) .^ p.lowerscale_Bruggeman(2,:);
    elseif strcmp(p.dualscale_input,'Lower scale diffusivity')
        vals = reshape(p.lowerscale_diffusivity,[1 length(p.lowerscale_diffusivity)]);
        p.lowerscale_diffusivity = [double(labels); double(vals)];
    end

    number_direction_todo = sum(p.directions);

    Deff = zeros(number_direction_todo,1);
    Mac = zeros(number_direction_todo,1);
    Tau = zeros(number_direction_todo,1);
    Bruggeman = zeros(number_direction_todo,1);
    vf = zeros(number_direction_todo,1);

    % Volume fractions (for table)
    volfraction = zeros(length(p.lowerscale_diffusivity(1,:)),1);
    for k=1:length(p.lowerscale_diffusivity(1,:))
        volfraction(k) = sum(sum(sum( M == p.lowerscale_diffusivity(1,k))));
    end
    volfraction = volfraction / numel(M);
    if p.isbackground
        volfraction = volfraction/(1-volfraction(1));
        volfraction(1) = []; % With GUI, first phase is always the background
    end

    directionname_todo = cell(number_direction_todo,1);
    current_direction_todo = 0;
    for dir = 1:dimension
        if p.directions(dir)
            current_direction_todo = current_direction_todo + 1;
            directionname_todo(current_direction_todo,1) = infovol.directionname(dir,1);

            time_cpu_start_phase = cputime; % CPU start
            time_clock_start_phase = tic; % Stopwatch start

            % % Algorithm
            p.direction = dir;
            p.name = ['Direction' num2str(dir)];
            [Deff(current_direction_todo), Tau(current_direction_todo), Bruggeman(current_direction_todo), vf(current_direction_todo)] = Tortuosity_algorithm(M,p);
            Mac(current_direction_todo) = 1/Deff(current_direction_todo);

            % % Correlation
            [str_direction] = function_remove_emptyandspecialcharacter_string(char(infovol.directionname(dir)));
            results_correlation.('Vf') = vf(1);
            results_correlation.(['Tortuosity_ax' num2str(dir) '_' str_direction]) = Tau(current_direction_todo);
            results_correlation.(['Normalized_diffusivity_ax' num2str(dir) '_' str_direction]) = Deff(current_direction_todo);
            results_correlation.(['MacMullin_number_ax' num2str(dir) '_' str_direction]) = Mac(current_direction_todo);
            results_correlation.(['Bruggeman_exponent_ax' num2str(dir) '_' str_direction]) = Bruggeman(current_direction_todo);
        end
    end

    if number_direction_todo==3
        results_correlation.Tortuosity_ax1over2 = Tau(1)/Tau(2);
        results_correlation.Tortuosity_ax1over3 = Tau(1)/Tau(3);
        results_correlation.Tortuosity_ax2over3 = Tau(2)/Tau(3);
    end

    % % Time
    timedata_domain(1,1) = {'Full volume'};
    timedata(1,2) = (cputime-time_cpu_start_volume) /number_direction_todo; % CPU elapsed time
    timedata(1,3) = toc(time_clock_start_volume) / number_direction_todo; % Stopwatch elapsed time
end

%% TABLES
% Time
Table_time = table(timedata_domain(:,1), timedata(:,1),timedata(:,2),timedata(:,3),...
    'VariableNames',{'Domain', 'Number of voxel','CPU time s' 'Stopwatch s'});
Results_Tortuosity.Table_time = Table_time; % Save in main table result

if strcmp(p.scale,'One scale')
    for current_phase_todo=1:1:number_phase_todo
        Tortuositytaufactor.phase(current_phase_todo).table = table(direction_todo(current_phase_todo).name,ones(direction_todo(current_phase_todo).number,1)*vf(current_phase_todo),Tau(current_phase_todo,direction_todo(current_phase_todo).idx)',Bruggeman(current_phase_todo,direction_todo(current_phase_todo).idx)',Deff(current_phase_todo,direction_todo(current_phase_todo).idx)',Mac(current_phase_todo,direction_todo(current_phase_todo).idx)',...
            'VariableNames',{'Direction','Volume fraction','Tortuosity factor','Bruggeman exponent','Normalized Deff','MacMullin number'});
    end
elseif strcmp(p.scale,'Dual scale')
    if strcmp(p.dualscale_input,'Lower scale diffusivity')
        if p.isbackground
            Lowerscalecoefficients_table = table(infovol.phasename(2:end),p.lowerscale_diffusivity(1,2:end)',volfraction,p.lowerscale_diffusivity(2,2:end)',...
                'VariableNames',{'Phase' 'Label' 'Volume fraction' 'Lower scale diffusivity'});%
        else
            Lowerscalecoefficients_table = table(infovol.phasename,p.lowerscale_diffusivity(1,:)',volfraction,p.lowerscale_diffusivity(2,:)',...
                'VariableNames',{'Phase' 'Label' 'Volume fraction' 'Lower scale diffusivity'});%
        end
    elseif strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
        if p.isbackground
            Lowerscalecoefficients_table = table(infovol.phasename(2:end),p.lowerscale_porosity(1,2:end)',volfraction,p.lowerscale_porosity(2,2:end)', p.lowerscale_Bruggeman(2,2:end)',p.lowerscale_diffusivity(2,2:end)',...
                'VariableNames',{'Phase' 'Label' 'Volume fraction', 'Lower scale porosity' 'Lower scale Bruggeman exponent' 'Lower scale diffusivity'});%
        else
            Lowerscalecoefficients_table = table(infovol.phasename,p.lowerscale_porosity(1,:)',volfraction,p.lowerscale_porosity(2,:)',p.lowerscale_Bruggeman(2,:)',p.lowerscale_diffusivity(2,2:end)',...
                'VariableNames',{'Phase' 'Label' 'Volume fraction', 'Lower scale Porosity' 'Lower scale Bruggeman exponent' 'Lower scale diffusivity'});%
        end
    end
    Tortuositytaufactor_table = table(directionname_todo,vf,Tau,Bruggeman,Deff,Mac,...
        'VariableNames',{'Direction','Porosity','Tortuosity factor','Bruggeman exponent','Normalized Deff','MacMullin number'});
end

%% SAVE TABLES
if opt.save.xls
    filename = 'Tortuosity'; % Filename without extension
    % Prepare the data
    clear DATA_writetable
    if strcmp(p.scale,'One scale')
        for current_phase_todo=1:1:number_phase_todo
            DATA_writetable.sheet(current_phase_todo).name=char(phasename_todo(current_phase_todo,1));
            DATA_writetable.sheet(current_phase_todo).table=Tortuositytaufactor.phase(current_phase_todo).table;
        end
    elseif strcmp(p.scale,'Dual scale')
        DATA_writetable.sheet(1).name='Lower scale coefficients';
        DATA_writetable.sheet(1).table=Lowerscalecoefficients_table;
        DATA_writetable.sheet(2).name='Full volume';
        DATA_writetable.sheet(2).table=Tortuositytaufactor_table;
    end
    % Save function
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%% DISPLAY
fprintf('> Calculated on the whole domain:\n\n');
if strcmp(p.scale,'One scale')
    for current_phase_todo=1:1:number_phase_todo
        fprintf('  For the phase %s:\n\n',char(phasename_todo(current_phase_todo,1)));
        disp(Tortuositytaufactor.phase(current_phase_todo).table)
    end
elseif strcmp(p.scale,'Dual scale')
    fprintf('  For the full volume:\n\n');
    fprintf('  Lower scale coefficients (%s):\n\n',p.dualscale_input);
    disp(Lowerscalecoefficients_table)
    fprintf('  Effective coefficients:\n\n');
    disp(Tortuositytaufactor_table)
end
fprintf('Computation time, in seconds:\n\n');
disp(Table_time)

%% FIGURE
scrsz = get(0,'ScreenSize'); % Screen resolution
Fig = figure; % Create figure
Fig.Name= 'Transport properties';
Fig.Color='white'; % Background colour

if strcmp(p.scale,'One scale')
    set(Fig,'position',round([scrsz(1) scrsz(2) 2/3*scrsz(3)  2/3*scrsz(4)]) ); % Full screen figure
    % - Create axes as a subplot
    X = categorical(phasename_todo);
    X = reordercats(X,phasename_todo);
    for id_axe=1:1:5
        sub_axes=subplot(2,3,id_axe,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        if id_axe==1
            h_title=title ('Volume fraction \epsilon'); % Title
            h_=bar(X,vf,0.5); % Width=1
            ylabel('Volume fraction \epsilon'); % Axis label
        elseif id_axe==2
            h_title=title ('Tortuosity factor \tau'); % Title
            h_=bar(X,Tau,1); % Width=1
            ylabel('Tortuosity factor \tau'); % Axis label
            ylim([1 Inf]);
        elseif id_axe==3
            h_title=title ('Bruggeman exponent p, \tau=\epsilon^{1-p}'); % Title
            h_=bar(X,Bruggeman,1); % Width=1
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
            h_=bar(X,Deff,1); % Width=1
            ylabel('D_{eff}/D_{bulk}'); % Axis label
        elseif id_axe==5
            h_title=title ('MacMullin number'); % Title
            h_=bar(X,Mac,1); % Width=1
            ylabel('D_{bulk}/D_{eff}'); % Axis label            
        end
        if id_axe>1
            h_legend = legend(sub_axes,infovol.directionname,'Location','best'); % Legend
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

elseif strcmp(p.scale,'Dual scale')

    if p.isbackground
        X = categorical(infovol.phasename(2:end));
        X = reordercats(X,infovol.phasename(2:end));
    else
        X = categorical(infovol.phasename);
        X = reordercats(X,infovol.phasename);
    end
    Y = categorical({'Full volume'});
    if strcmp(p.dualscale_input,'Lower scale diffusivity')
        set(Fig,'position',round([scrsz(1) scrsz(2) 2/3*scrsz(3)  2/3*scrsz(4)]) ); % Full screen figure
        for id_axe=1:1:6
            if id_axe~=3
                sub_axes=subplot(2,3,id_axe,'Parent',Fig);
                hold(sub_axes,'on'); % Active subplot
            end
            if id_axe==1
                h_title=title ('Volume fraction \epsilon'); % Title
                h_=bar(X,volfraction,0.5); % Width=1
                ylabel('\epsilon'); % Axis label
            elseif id_axe==2
                h_title=title ('Lower scale diffusivity'); % Title
                h_=bar(X,Lowerscalecoefficients_table.("Lower scale diffusivity"),0.5); % Width=1
                ylabel('D_{lower}'); % Axis label                   
            elseif id_axe==4
                h_title=title ('Tortuosity factor \tau'); % Title
                h_=bar(Y,Tau,1); % Width=1
                ylabel('\tau'); % Axis label
            elseif id_axe==5
                h_title=title ('Normalized effective diffusion coefficient'); % Title
                h_=bar(Y,Deff,1); % Width=1
                ylabel('D_{eff}/D_{bulk}'); % Axis label                
            elseif id_axe==6
                h_title=title ('MacMullin number'); % Title
                h_=bar(Y,Mac,1); % Width=1
                ylabel('D_{bulk}/D_{eff}'); % Axis label                
            end
            if id_axe<=3
                h_.FaceColor = [0.5 0.5 0.5];
            end
            if id_axe>3
                h_legend = legend(sub_axes,infovol.directionname,'Location','best'); % Legend
            end
            set(sub_axes,'YGrid',opt.format.grid); % Display grid
            set(sub_axes,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
            if id_axe>3
                h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
            end
            hold(sub_axes,'off'); % Relase figure
        end   

    elseif strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
        set(Fig,'position',round([scrsz(1) scrsz(2) 1*scrsz(3)  2/3*scrsz(4)]) ); % Full screen figure
        for id_axe=1:1:10
            if id_axe~=5
                sub_axes=subplot(2,5,id_axe,'Parent',Fig);
                hold(sub_axes,'on'); % Active subplot
            end
            if id_axe==1
                h_title=title ('Volume fraction \epsilon'); % Title
                h_=bar(X,volfraction,0.5); % Width=1
                ylabel('\epsilon'); % Axis label
            elseif id_axe==2
                h_title=title ('Lower scale porosity \epsilon_{lower}'); % Title
                h_=bar(X,Lowerscalecoefficients_table.("Lower scale porosity"),0.5); % Width=1
                ylabel('\epsilon_{lower}'); % Axis label
            elseif id_axe==3
                h_title=title ('Lower scale Bruggeman exponent'); % Title
                h_=bar(X,Lowerscalecoefficients_table.("Lower scale Bruggeman exponent"),0.5); % Width=1
                ylabel('p_{lower}'); % Axis label   
                xls = xlim;
                plot([xls(1) xls(2)],[1 1],'Color','k','LineStyle','--','LineWidth',2);
                plot([xls(1) xls(2)],[1.5 1.5],'Color','k','LineStyle','--','LineWidth',2);
                plot([xls(1) xls(2)],[2.0 2.0],'Color','k','LineStyle','--','LineWidth',2);
                xt = [xls(2) xls(2) xls(2)];
                yt = [1 1.5 2.0];
                str = {'Rule of mixture','Spheres','Cylinders'};
                text(xt,yt,str,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize);                
            elseif id_axe==4
                h_title=title ('Lower scale diffusivity'); % Title
                h_=bar(X,Lowerscalecoefficients_table.("Lower scale diffusivity"),0.5); % Width=1
                ylabel('D_{lower}'); % Axis label                   
            elseif id_axe==6
                h_title=title ('Porosity \epsilon'); % Title
                h_=bar(Y,vf,1); % Width=1
                ylabel('\epsilon'); % Axis label               
            elseif id_axe==7
                h_title=title ('Tortuosity factor \tau'); % Title
                h_=bar(Y,Tau,1); % Width=1
                ylabel('\tau'); % Axis label
            elseif id_axe==8
                h_title=title ('Bruggeman exponent p, \tau=\epsilon^{1-p}'); % Title
                h_=bar(Y,Bruggeman,1); % Width=1
                ylabel('p'); % Axis label      
                xls = xlim;
                plot([xls(1) xls(2)],[1 1],'Color','k','LineStyle','--','LineWidth',2);
                plot([xls(1) xls(2)],[1.5 1.5],'Color','k','LineStyle','--','LineWidth',2);
                plot([xls(1) xls(2)],[2.0 2.0],'Color','k','LineStyle','--','LineWidth',2);
                xt = [xls(2) xls(2) xls(2)];
                yt = [1 1.5 2.0];
                str = {'Rule of mixture','Spheres','Cylinders'};
                text(xt,yt,str,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize);                   
            elseif id_axe==9
                h_title=title ('Normalized effective diffusion coefficient'); % Title
                h_=bar(Y,Deff,1); % Width=1
                ylabel('D_{eff}/D_{bulk}'); % Axis label                
            elseif id_axe==10
                h_title=title ('MacMullin number'); % Title
                h_=bar(Y,Mac,1); % Width=1
                ylabel('D_{bulk}/D_{eff}'); % Axis label                
            end
            if id_axe<=4
                h_.FaceColor = [0.5 0.5 0.5];
            end
            if id_axe>5
                h_legend = legend(sub_axes,infovol.directionname,'Location','best'); % Legend
            end
            set(sub_axes,'YGrid',opt.format.grid); % Display grid
            set(sub_axes,'YMinorGrid',opt.format.minorgrid); % Display grid for minor thicks
            set(sub_axes,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
            h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
            if id_axe>5
                h_legend.FontSize = opt.format.legendfontsize; % Set title fontsize
            end
            hold(sub_axes,'off'); % Relase figure
        end
    end
end
sgtitle(Fig,'Transport properties','FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
if opt.save.savefig % Save figure
    filename= 'Transport_properties';
    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
end
if opt.format.autoclosefig
    %close(Fig); % Do not keep open figures
end


%%
%% ADDITIONAL RESULTS ON THE WHOLE VOLUME 
%%

p.iterate_figure = false;

%% ALONG DIRECTIONS

sheet = 0;
clear DATA_writetable

if strcmp(p.scale,'One scale')
    current_phase_todo = 0;
    for current_phase=1:1:number_phase % Loop over all phases
        DeffS = zeros(3,3) - 1;
        if p.todo(current_phase)
            current_phase_todo=current_phase_todo+1;
            phasename_todo(current_phase_todo,1) = infovol.phasename(current_phase,1);
            for section_dir = 1:dimension
                nsection = p.sections(current_phase).vals(section_dir);
                bounds = round(linspace(0,sz(section_dir),nsection+1));
                if nsection>1
                    for tau_dir = 1:1:dimension
                        if p.directions(current_phase).vals(tau_dir)
                            fprintf('> Calculated section per section for phase %s (cut along direction %i with %i sections) for transport along direction %i:\n\n',char(phasename_todo(current_phase_todo,1)),section_dir,nsection,tau_dir);
                            for ksection = 1:nsection
                                if section_dir == 1
                                    x0 = bounds(ksection)+1; x1 = bounds(ksection+1);
                                    y0 = 1; y1 = sz(2);
                                    z0 = 1; z1 = sz(3);
                                elseif section_dir == 2
                                    x0 = 1; x1 = sz(1);
                                    y0 = bounds(ksection)+1; y1 = bounds(ksection+1);
                                    z0 = 1; z1 = sz(3);
                                else
                                    x0 = 1; x1 = sz(1);
                                    y0 = 1; y1 = sz(2);
                                    z0 = bounds(ksection)+1; z1 = bounds(ksection+1);
                                end
                                MM = M(x0:x1,y0:y1,z0:z1);

                                % % Algorithm
                                p.direction = tau_dir;
                                p.onescale_label = infovol.phaselabel(current_phase);
                                [DeffS(section_dir,tau_dir,ksection), TauS(section_dir,tau_dir,ksection), BruggemanS(section_dir,tau_dir,ksection), vfS(section_dir,tau_dir,ksection)] = Tortuosity_algorithm(MM,p);
                                MacS(section_dir,tau_dir,ksection) = 1/DeffS(section_dir,tau_dir,ksection);
                            end
                        end
                    end
                end
            end
        

            % Table
            for section_dir = 1:dimension
                nsection = p.sections(current_phase).vals(section_dir);
                bounds = round(linspace(0,sz(section_dir),nsection+1));
                for tau_dir = 1:1:dimension
                    if DeffS(section_dir,tau_dir,1) ~= -1
                        sheet = sheet+1;

                        c1 = (bounds(1:end-1)+1) * voxel_size;
                        c2 = bounds(2:end) * voxel_size;
                        vfvals = reshape( vfS(section_dir,tau_dir,1:nsection),[nsection,1]);
                        Deffvals = reshape( DeffS(section_dir,tau_dir,1:nsection),[nsection,1]);
                        Macvals = reshape( MacS(section_dir,tau_dir,1:nsection),[nsection,1]);
                        Tauvals = reshape( TauS(section_dir,tau_dir,1:nsection),[nsection,1]);
                        Brugvals = reshape( BruggemanS(section_dir,tau_dir,1:nsection),[nsection,1]);

                        Variable_name_table={['Section start ' voxel_unit],['Section end ' voxel_unit],'Volume fraction','Effective diffusivity','MacMullin number','Tortuosity factor','Bruggeman exponent'};
                        tmp = array2table([c1',c2',vfvals,Deffvals,Macvals,Tauvals,Brugvals],'VariableNames',Variable_name_table);

                        DATA_writetable.sheet(sheet).name = [char(phasename_todo(current_phase_todo,1)) '_CutalongAxe' num2str(section_dir) '_TransportalongAxe' num2str(tau_dir)];
                        DATA_writetable.sheet(sheet).table = tmp;
                    end
                end
            end

            % Figure
            m = sum( DeffS(:,1,1)~=-1 );
            if m>0
                Fig = figure;
                Fig.Name= ['Transport properties calculated per section for ' char(phasename_todo(current_phase_todo,1))];
                Fig.Color='white'; % Background colour
                set(Fig,'position',round([scrsz(1) scrsz(2) scrsz(3)  round(m/3*scrsz(4))]) ); % Full screen figure
                t = tiledlayout(m,5,'Parent',Fig);
                %t.TileSpacing = 'compact';
                %t.Padding = 'compact';
                for section_dir = 1:dimension
                    if DeffS(section_dir,1,1)~=-1
                        nsection = p.sections(current_phase).vals(section_dir);
                        ntaudir = sum( DeffS(section_dir,:,1)~=-1 );

                        sectioname = [];
                        for s=1:nsection
                            sectioname = [sectioname {['section #' num2str(s)]}];
                        end
                        X = categorical(sectioname);
                        X = reordercats(X,sectioname);

                        for k=1:5
                            ax = nexttile;
                            hold(ax,'on'); % Active subplot
                            if k==1
                                y = reshape( vfS(section_dir,1,1:nsection),[nsection,1]);
                                bar(X,y)
                                strylabel = 'Volume fraction \epsilon';
                                xls = xlim;
                                plot([xls(1) xls(2)],[vf(current_phase_todo) vf(current_phase_todo)],'Color','k','LineStyle','--','LineWidth',2);
                                text(xls(2),vf(current_phase_todo),{'Full volume'},'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize);
                            else
                                if k==2
                                    tmp = DeffS;
                                    strylabel = 'Effective diffusivity D_{eff}/D_{bulk}';

                                    % xls = xlim;
                                    % for dir = 1:dimension
                                    %     if p.directions(dir)
                                    %         plot([xls(1) xls(2)],[Deff(current_phase_todo,dir) Deff(current_phase_todo,dir)],'Color','k','LineStyle','--','LineWidth',2);
                                    %     end
                                    % end
                                    % text(xls(2),vf(current_phase_todo),{'Full volume'},'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize);

                                elseif k==3
                                    tmp = MacS;
                                    strylabel = 'MacMullin number D_{bulk}/D_{eff}';
                                elseif k==4
                                    tmp = TauS;
                                    strylabel = 'Tortuosity factor \tau';
                                elseif k==5
                                    tmp = BruggemanS;
                                    strylabel = 'Bruggeman exponent p';
                                end
                                y = zeros(nsection,ntaudir);
                                for s=1:nsection
                                    y(s,:) = reshape( tmp(section_dir,1:ntaudir,s),[ntaudir,1]);
                                end
                                bar(X,y,1.0)
                                h_legend = legend(infovol.directionname,'Location','best'); % Legend
                                h_legend.FontSize = opt.format.legendfontsize;

                            end
                            xlabel(['Cut along ' char(infovol.directionname(section_dir))])
                            ylabel(strylabel)
                            set(ax,'YGrid',opt.format.grid)
                            set(ax,'YMinorGrid',opt.format.minorgrid)
                            set(ax,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
                        end

                    end
                end
                sgtitle(Fig,Fig.Name,'FontWeight','bold','FontSize',opt.format.sgtitlefontsize,'FontName',opt.format.fontname);
                if opt.save.savefig % Save figure
                    filename= Fig.Name;
                    function_savefig(Fig, Current_folder, filename, opt.save); % Call function
                end
                if opt.format.autoclosefig
                    %    close(Fig); % Do not keep open figures
                end
            end
        end




    end

elseif strcmp(p.scale,'Dual scale')


end

if opt.save.xls && sheet>0
    filename = 'Transport_along_directions'; % Filename without extension
    Function_Writetable(Current_folder,filename,DATA_writetable)
end

%%
%% IMAGE RESOLUTION SENSITIVITY ANALYSIS
%%


%%
%% REPRESENTATITVE VOLUME ELEMENT (RVE) AND CONVERGENCE ANALYSIS
%%


%%
%% ENDING FUNCTION
%%

end