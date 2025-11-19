function [] = function_summary_onevolume(folder, infovol, optc)
%function_summary_onevolume summarize results in table

keyboard

TO BE REDONE

% Find all m files
MyFolderInfo=dir(folder);
n_files=length(MyFolderInfo);
kk=0;
for k=1:1:n_files
    [~,name,ext] = fileparts(MyFolderInfo(k).name);
    if strcmp(ext,'.mat')
        sub_name = name(1:7);
        if strcmp(sub_name,'Results')
            kk=kk+1;
            results_filenames(kk)={[name '.mat']};
        end
    end
end

if kk>0

    string_precision = 4;
    number_results = length(results_filenames); % Number of results to summarize
    number_phase = length(infovol.phasename);
    scrsz = get(0,'ScreenSize'); % Screen resolution

    if number_results>=1
        line_table=0;
        for k=1:1:number_results % Loop over all saved results
            % Load result
            current_result=char(results_filenames(k));

            % Case volume fraction
            if strcmp(current_result,'Results_volume_fraction.mat')
                pathname = fullfile(folder, current_result);
                datamat = load(pathname);
                line_table=line_table+1;
                Propertyname(line_table) = {'Volume fraction'}; % Name
                Propertymethod(line_table) = {'Voxel summation'}; % Name
                T = datamat.Results_Volumefraction.Table_Volumefraction; % Table
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.volumefraction.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        array = T{:,3};
                        mean_value = array(current_phase_todo);
                        PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                        PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                        PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                        PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'[]'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
            end

            % Case Tortuosity
            if strcmp(current_result,'Results_tortuosity.mat')
                pathname = [folder current_result];
                datamat = load(pathname);
                for kproperty=1:1:3
                    line_table=line_table+1;
                    if kproperty==1
                        Propertyname(line_table) = {'Tortuosity factor'}; % Name
                    elseif kproperty==2
                        Propertyname(line_table) = {'Bruggeman exponent'}; 
                    else
                        Propertyname(line_table) = {'Normalized effective diffusivity'}; 
                    end
                    Propertymethod(line_table) = {'TauFactor (FD, Dirichlet BCs)'}; % Name
                    current_phase_todo = 0;
                    for current_phase=1:1:number_phase
                        if opts.tortuosity.todo(current_phase)
                            current_phase_todo=current_phase_todo+1;
                            T = datamat.Results_Tortuosity.Tortuositytaufactor.phase(current_phase_todo).table; % Table
                            array = T{:,kproperty+2};
                            current_direction_todo = 0;
                            for current_direction=1:1:length(opts.tortuosity.direction_todo)
                                if opts.tortuosity.direction_todo(current_direction)
                                    current_direction_todo = current_direction_todo+1;
                                    if current_direction==1
                                        str = num2str(array(current_direction_todo),string_precision);
                                    else
                                        str = [str ' / ' num2str(array(current_direction_todo),string_precision)];
                                    end
                                else
                                    if current_direction==1
                                        str = '-';
                                    else
                                        str = [str ' / ' '-'];
                                    end
                                end
                            end
                            PropertyMean(line_table,current_phase) = {str}; % Mean
                            PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                            PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                            PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                            PropertyUnit(line_table,current_phase) = {'[]'}; % Unit
                            PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                        else % Not calculated
                            PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                            PropertyMin(line_table,current_phase) = {'-'}; % Min
                            PropertyMax(line_table,current_phase) = {'-'}; % Max
                            PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                            PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                            PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                        end
                    end
                end  
            end

            % Case specific surface area (direct)
            if strcmp(current_result,'Results_Sp_direct.mat')
                Spunit = [infovol.unit '^-1'];
                pathname = [folder current_result];
                datamat = load(pathname);
                line_table=line_table+1;
                Propertyname(line_table) = {'Phase surface / phase volume'}; % Name
                Propertymethod(line_table) = {'Face summation'}; % Name
                if isfield(datamat.Results_specificsurfacearea,'Table_Sp_phasevolume_direct_removedcompvol')
                    T = datamat.Results_specificsurfacearea.Table_Sp_phasevolume_direct_removedcompvol; % Table
                else
                    T = datamat.Results_specificsurfacearea.Table_Sp_phasevolume_direct; % Table
                end
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.Sp_direct.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        array = T{:,7};
                        mean_value = array(current_phase_todo);
                        PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                        PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                        PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                        PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {Spunit}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
                line_table=line_table+1;
                Propertyname(line_table) = {'Phase surface / domain volume'}; % Name
                Propertymethod(line_table) = {'Face summation'}; % Name
                T = datamat.Results_specificsurfacearea.Table_Sp_domainvolume_direct; % Table
                if isfield(datamat.Results_specificsurfacearea,'Table_Sp_phasevolume_direct_removedcompvol')
                    T = datamat.Results_specificsurfacearea.Table_Sp_domainvolume_direct_removedcompvol; % Table
                else
                    T = datamat.Results_specificsurfacearea.Table_Sp_domainvolume_direct; % Table
                end
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.Sp_direct.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        array = T{:,7};
                        mean_value = array(current_phase_todo);
                        PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                        PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                        PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                        PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {Spunit}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
            end

            % Case specific interface area (direct)
            if strcmp(current_result,'Results_Int_direct.mat')
                Spunit = [infovol.unit '^-1'];
                pathname = [folder current_result];
                datamat = load(pathname);
                T = datamat.Results_specificinterfacearea.Table_Int_direct; % Table
                interface_names = T(:,1).Variables; % Name
                [number_interface,~] = size(interface_names);
                for current_interface=1:1:number_interface
                    line_table=line_table+1;
                    Propertyname(line_table) = {[char(interface_names(current_interface,:)) ' area / domain volume']}; % Name
                    Propertymethod(line_table) = {'Face summation'}; % Name
                    for current_phase=1:1:number_phase
                        if contains(char(interface_names(current_interface)), char(infovol.phasename(current_phase)))
                            PropertyMean(line_table,current_phase) = {num2str(T{current_interface,8},string_precision)}; % Mean
                            PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                            PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                            PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                            PropertyUnit(line_table,current_phase) = {Spunit}; % Unit
                            PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                        else % Not calculated
                            PropertyMean(line_table,current_phase) = {'n/a'}; % Mean
                            PropertyMin(line_table,current_phase) = {'-'}; % Min
                            PropertyMax(line_table,current_phase) = {'-'}; % Max
                            PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                            PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                            PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                        end
                    end
                end
            end

            % Case Triple phase boundary length (direct)
            if strcmp(current_result,'Results_TPBL_direct.mat')
                TPBLunit = [infovol.unit '^-2'];
                pathname = [folder current_result];
                datamat = load(pathname);
                T = datamat.Results_TPBL.Table_TPBL_direct; % Table
                line_names = T(:,1).Variables; % Name
                [number_line,~] = size(line_names);                
                for current_line=1:1:number_line
                    line_table=line_table+1;
                    Propertyname(line_table) = {[char(line_names(current_line,:)) ' line lenght / domain volume']}; % Name
                    Propertymethod(line_table) = {'Line summation'}; % Name
                    for current_phase=1:1:number_phase
                        if contains(char(line_names(current_line)), char(infovol.phasename(current_phase)))
                            PropertyMean(line_table,current_phase) = {num2str(T{current_line,8},string_precision)}; % Mean
                            PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                            PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                            PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                            PropertyUnit(line_table,current_phase) = {TPBLunit}; % Unit
                            PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                        else % Not calculated
                            PropertyMean(line_table,current_phase) = {'n/a'}; % Mean
                            PropertyMin(line_table,current_phase) = {'-'}; % Min
                            PropertyMax(line_table,current_phase) = {'-'}; % Max
                            PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                            PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                            PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                        end
                    end
                end
            end

            % Case mean diameter c-psd
            if strcmp(current_result,'Results_D50_Cpsd.mat')
                Dpunit = infovol.unit;
                pathname = [folder current_result];
                datamat = load(pathname);
                line_table=line_table+1;
                Propertyname(line_table) = {'Equivalent diameter'}; % Name
                Propertymethod(line_table) = {'C-psd'}; % Name
                T = datamat.Results_cpsd.Table_cpsd; % Table
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.D_cpsd.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        array = T{current_phase_todo,3:7};
                        PropertyMin(line_table,current_phase) = {num2str(array(1),string_precision)}; % Min
                        PropertyMean(line_table,current_phase) = {num2str(array(2),string_precision)}; % Mean
                        PropertyMax(line_table,current_phase) = {num2str(array(3),string_precision)}; % Max
                        PropertyStd(line_table,current_phase) = {num2str(array(4),string_precision)}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {Dpunit}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {num2str(array(5),string_precision)}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
                line_table=line_table+1;
                Propertyname(line_table) = {'Particle level of details'}; % Name
                Propertymethod(line_table) = {'C-psd'}; % Name
                T = datamat.Results_cpsd.Table_cpsd_lvldetail; % Table
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.D_cpsd.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        array = T{current_phase_todo,3};
                        PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                        PropertyMean(line_table,current_phase) = {num2str(array,string_precision)}; % Mean
                        PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                        PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'[0,1]'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
            end

            % Case mean diameter edmf
            if strcmp(current_result,'Results_D50_edmf.mat')
                Dpunit = infovol.unit;
                pathname = [folder current_result];
                datamat = load(pathname);
                line_table=line_table+1;
                Propertyname(line_table) = {'Equivalent diameter'}; % Name
                Propertymethod(line_table) = {'EDMF'}; % Name
                T = datamat.Results_edmf.Table_edmf_d50; % Table
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.D_edmf.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        array = T{current_phase_todo,3};
                        PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                        PropertyMean(line_table,current_phase) = {num2str(array,string_precision)}; % Mean
                        PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                        PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {Dpunit}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
                line_table=line_table+1;
                Propertyname(line_table) = {'Distance to surface'}; % Name
                Propertymethod(line_table) = {'Distance map'}; % Name
                T = datamat.Results_edmf.Table_dmap; % Table
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.D_edmf.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        array = T{current_phase_todo,3:7};
                        PropertyMin(line_table,current_phase) = {num2str(array(1),string_precision)}; % Min
                        PropertyMean(line_table,current_phase) = {num2str(array(2),string_precision)}; % Mean
                        PropertyMax(line_table,current_phase) = {num2str(array(3),string_precision)}; % Max
                        PropertyStd(line_table,current_phase) = {num2str(array(4),string_precision)}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {Dpunit}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {num2str(array(5),string_precision)}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
            end

            % Case mean diameter chord
            if strcmp(current_result,'Results_D50_chord.mat')
                Dpunit = infovol.unit;
                pathname = [folder current_result];
                datamat = load(pathname);
                line_table=line_table+1;
                Propertyname(line_table) = {'Equivalent diameter'}; % Name
                Propertymethod(line_table) = {'CLMF'}; % Name
                current_phase_todo = 0;
                for current_phase=1:1:number_phase
                    if opts.D_chord.todo(current_phase)
                        current_phase_todo=current_phase_todo+1;
                        T = datamat.Results_chord.Table_chord_d50; % Table
                        array = T{current_phase_todo,2:end};
                        str = [num2str(array(1),string_precision) ' / ' num2str(array(2),string_precision) ' / ' num2str(array(3),string_precision)];
                        PropertyMean(line_table,current_phase) = {str}; % Mean
                        PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                        PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                        PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {Dpunit}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                    else % Not calculated
                        PropertyMean(line_table,current_phase) = {'not calculated'}; % Mean
                        PropertyMin(line_table,current_phase) = {'-'}; % Min
                        PropertyMax(line_table,current_phase) = {'-'}; % Max
                        PropertyStd(line_table,current_phase) = {'-'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'-'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'-'}; % Standard deviation in percents of the mean
                    end
                end
            end



        end

        if exist('Propertyname','var')
            Propertyname_sav = Propertyname; Propertymethod_sav = Propertymethod; PropertyMean_sav = PropertyMean; PropertyMin_sav = PropertyMin; PropertyMax_sav = PropertyMax; PropertyStd_std = PropertyStd; PropertyUnit_sav=PropertyUnit; PropertyStdpercent_sav=PropertyStdpercent;
            for current_phase=1:1:number_phase % Loop over all phases

                Propertyname = Propertyname_sav; Propertymethod = Propertymethod_sav; PropertyMean = PropertyMean_sav; PropertyMin = PropertyMin_sav; PropertyMax = PropertyMax_sav; PropertyStd = PropertyStd_std; PropertyUnit=PropertyUnit_sav; PropertyStdpercent=PropertyStdpercent_sav;
                has_been_modified=true;
                while has_been_modified
                    n=length(PropertyMean(:,current_phase));
                    for k=1:1:n
                        if isempty(cell2mat( PropertyMean(k,current_phase) ))
                            Propertyname(k)=[]; Propertymethod(k)=[]; PropertyMean(k,:)=[]; PropertyMin(k,:)=[]; PropertyMax(k,:)=[]; PropertyStd(k,:)=[]; PropertyUnit(k,:)=[]; PropertyStdpercent(k,:)=[];
                            has_been_modified=true; break;
                        end
                        has_been_modified=false;
                    end
                end

                Table_summary_results.phase(current_phase).table = table(Propertyname',Propertymethod',PropertyMean(:,current_phase),PropertyMin(:,current_phase),PropertyMax(:,current_phase),PropertyStd(:,current_phase),PropertyUnit(:,current_phase),PropertyStdpercent(:,current_phase),...
                    'VariableNames',{'Property' 'Method' 'Mean' 'Min' 'Max' 'Std' 'Unit' 'Std_percents'});
                fprintf('       For the phase: %s:\n\n',char(infovol.phasename(current_phase)));
                disp(Table_summary_results.phase(current_phase).table); disp ' ';

                Fig = figure; % Create figure
                Fig.Name= sprintf('Results of phase %s',char(infovol.phasename(current_phase))); % Figure name
                Fig.Color='white'; % Background colour
                % Remove toolbar
                Fig.ToolBar='none';
                Fig.MenuBar='none';
                Fig.NumberTitle='off';
                % Position unit
                Fig.Units='normalized';
                % Position
                set(Fig, 'Position', [0.1 0.1 0.8 0.375]);
                length_fig_x = scrsz(3)*0.8;
                % Create uitable
                % Cell data
                Cell_data = [Propertyname',Propertymethod',PropertyMean(:,current_phase),PropertyMin(:,current_phase),PropertyMax(:,current_phase),PropertyStd(:,current_phase),PropertyUnit(:,current_phase),PropertyStdpercent(:,current_phase)];
                % Coloum name
                hs='<html><font size=5><font face="Times">';
                he='</h1></html>';
                Column_name_table={[hs 'Property' he],[hs 'Method' he],[hs 'Mean' he],[hs 'Min' he],[hs 'Max' he],[hs 'Std' he],[hs 'Unit' he],[hs 'Std_percents' he]};
                % Uitable
                t_basic = uitable('Parent', Fig,...
                    'RowName',[],...
                    'ColumnName',Column_name_table,...
                    'ColumnFormat',({'char' 'char' 'char' 'char' 'char' 'char' 'char' 'char'}),...
                    'ColumnEditable', true,...
                    'RowStriping','on',...
                    'Data', Cell_data);

                % Format uitable
                t_basic.Units='normalized';
                t_basic.Position=[0.025 0.025 0.95 0.95];
                t_basic.FontName='Times New Roman';
                t_basic.FontSize=14;
                t_basic.FontSize=14;
                length_tab_x = 0.95*length_fig_x;
                t_basic.ColumnWidth ={length_tab_x*0.25 length_tab_x*0.125 length_tab_x*0.2 length_tab_x*0.075 length_tab_x*0.075 length_tab_x*0.075 length_tab_x*0.1 length_tab_x*0.1};
                if optc.save.savefig % Save figure
                    filename= sprintf('Results_of_%s',char(infovol.phasename(current_phase)));
                    function_savefig(Fig, folder, filename, optc.save); % Call function
                end
                if optc.format.autoclosefig
                    close(Fig); % Do not keep open figures
                end
            end

            % Save
            if optc.save.mat
                save([folder 'Table_summary_results'],'Table_summary_results')
            end
            if optc.save.xls
                filename = 'Table_summary_results'; % Filename without extension
                % Prepare the data
                clear DATA_writetable
                sheetnumber=0;
                for current_phase=1:1:number_phase
                    sheetnumber=sheetnumber+1;
                    DATA_writetable.sheet(sheetnumber).name=char(infovol.phasename(current_phase));
                    DATA_writetable.sheet(sheetnumber).table=Table_summary_results.phase(current_phase).table;
                end
                % Save function
                Function_Writetable(folder,filename,DATA_writetable)
            end
        end

    end

end

end

