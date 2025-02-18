function [] = function_summary_onevolume(folder, INFO, OPTIONS)
%function_summary_onevolume summarize results in table

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

string_precision = 4;
number_phase = length(INFO.phasename); % Number of phase
number_results = length(results_filenames); % Number of results to summarize
scrsz = get(0,'ScreenSize'); % Screen resolution

if number_results>=1
    line_table=0;
    for k=1:1:number_results % Loop over all saved results
        % Load result
        current_result=char(results_filenames(k));
        
        % Case volume fraction
        if strcmp(current_result,'Results_volume_fraction.mat')
            pathname = [folder current_result];
            datamat = load(pathname);
            line_table=line_table+1;
            Propertyname(line_table) = {'Volume fraction'}; % Name
            Propertymethod(line_table) = {'Voxel summation'}; % Name
            T = datamat.Results_Volumefraction.Table_Volumefraction; % Table
            for current_phase=1:1:number_phase
                array = T{:,3};
                mean_value = array(current_phase);
                PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
        end
        
        % Case specific surface area (direct)
        if strcmp(current_result,'Results_specific_surface_direct.mat')
            pathname = [folder current_result];
            datamat = load(pathname);
            line_table=line_table+1;            
            Propertyname(line_table) = {'Phase surface / phase volume'}; % Name
            Propertymethod(line_table) = {'Face summation'}; % Name
            T = datamat.Results_specificsurfacearea.Table_Specificsurface_phasevolume_direct; % Table
            for current_phase=1:1:number_phase
                array = T{:,7};
                mean_value = array(current_phase);
                PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[um^-1]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
            line_table=line_table+1;
            Propertyname(line_table) = {'Phase surface / domain volume'}; % Name
            Propertymethod(line_table) = {'Face summation'}; % Name
            T = datamat.Results_specificsurfacearea.Table_Specificsurface_domainvolume_direct; % Table
            for current_phase=1:1:number_phase
                array = T{:,7};
                mean_value = array(current_phase);
                PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[um^-1]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
        end
        
        % Case specific interface area (direct)
        if strcmp(current_result,'Results_specific_interface_direct.mat')
            pathname = [folder current_result];
            datamat = load(pathname);
            T = datamat.Results_specificinterfacearea.Table_Specificinterface_direct; % Table
            interface_names = T(:,1).Variables;
            [number_interface,~] = size(interface_names);
            for current_interface=1:1:number_interface
                for current_phase=1:1:number_phase
                    if contains(interface_names(current_interface,:), INFO.phase(current_phase).name)
                        line_table=line_table+1;
                        Propertyname(line_table) = {[interface_names(current_interface,:) ' area / domain volume']}; % Name
                        Propertymethod(line_table) = {'Face summation'}; % Name
                        mean_value = T{current_interface,6};
                        PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                        PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                        PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                        PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                        PropertyUnit(line_table,current_phase) = {'[um^-1]'}; % Unit
                        PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
                    end
                end
            end
        end        

        % Case C-PSD
        if strcmp(current_result,'Results_particlediameter_CPSD.mat')
            pathname = [folder current_result];
            datamat = load(pathname);
            line_table=line_table+1;
            Propertyname(line_table) = {'Equivalent diameter'}; % Name
            Propertymethod(line_table) = {'C-PSD'}; % Name            
            T = datamat.Results_cpsd.Table_cpsd; % Table
            for current_phase=1:1:number_phase
                array = T{current_phase,2:end};
                PropertyMin(line_table,current_phase) = {num2str(array(1),string_precision)}; % Min
                PropertyMean(line_table,current_phase) = {num2str(array(2),string_precision)}; % Mean
                PropertyMax(line_table,current_phase) = {num2str(array(3),string_precision)}; % Max
                PropertyStd(line_table,current_phase) = {num2str(array(4),string_precision)}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[um]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {num2str(array(5),string_precision)}; % Standard deviation in percents of the mean                
            end            
            line_table=line_table+1;
            Propertyname(line_table) = {'Particle level of details'}; % Name
            Propertymethod(line_table) = {'C-PSD'}; % Name
            T = datamat.Results_cpsd.Table_cpsd_lvldetail; % Table
            for current_phase=1:1:number_phase
                array = T{:,2};
                mean_value = array(current_phase);
                PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
%             line_table=line_table+1; % Example how to add new summary
%             Propertyname(line_table) = {'Equivalent radius'}; % Name
%             Propertymethod(line_table) = {'C-PSD'}; % Name
%             T = datamat.Results_cpsd.Table_radius; % Table
%             for current_phase=1:1:number_phase
%                 array = T{current_phase,2:end};
%                 PropertyMin(line_table,current_phase) = {num2str(array(1),string_precision)}; % Min
%                 PropertyMean(line_table,current_phase) = {num2str(array(2),string_precision)}; % Mean
%                 PropertyMax(line_table,current_phase) = {num2str(array(3),string_precision)}; % Max
%                 PropertyStd(line_table,current_phase) = {num2str(array(4),string_precision)}; % Standard deviation
%                 PropertyUnit(line_table,current_phase) = {'[um]'}; % Unit
%                 PropertyStdpercent(line_table,current_phase) = {num2str(array(5),string_precision)}; % Standard deviation in percents of the mean
%             end
        end        
        
        % Case Distance map
        if strcmp(current_result,'Results_particlediameter_Dmap.mat')
            pathname = [folder current_result];
            datamat = load(pathname);
            line_table=line_table+1;
            Propertyname(line_table) = {'Distance to surface'}; % Name
            Propertymethod(line_table) = {'Distance map'}; % Name            
            T = datamat.Results_dmap.Table_dmap; % Table
            for current_phase=1:1:number_phase
                array = T{current_phase,2:end};
                PropertyMin(line_table,current_phase) = {num2str(array(1),string_precision)}; % Min
                PropertyMean(line_table,current_phase) = {num2str(array(2),string_precision)}; % Mean
                PropertyMax(line_table,current_phase) = {num2str(array(3),string_precision)}; % Max
                PropertyStd(line_table,current_phase) = {num2str(array(4),string_precision)}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[um]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {num2str(array(5),string_precision)}; % Standard deviation in percents of the mean                
            end            
            line_table=line_table+1;
            Propertyname(line_table) = {'Equivalent diameter'}; % Name
            Propertymethod(line_table) = {'Distance map'}; % Name
            T = datamat.Results_dmap.Table_dmap_d50; % Table
            for current_phase=1:1:number_phase
                array = T{:,2};
                mean_value = array(current_phase);
                PropertyMean(line_table,current_phase) = {num2str(mean_value,string_precision)}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[um]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end            
        end                     
                
        % Case Tortuosity factor
        if strcmp(current_result,'Results_tortuosityfactor_taufactor.mat')
            pathname = [folder current_result];
            datamat = load(pathname);
            line_table=line_table+1;
            Propertyname(line_table) = {'Tortuosity factor 1/2/3'}; % Name
            Propertymethod(line_table) = {'Tau factor'}; % Name
            for current_phase=1:1:number_phase
                T = datamat.Results_tortuosityfactor.Table_TortuosityFactor_taufactor.phase(current_phase).table; % Table
                array = T{:,3};
                str_mean=[num2str(array(1),string_precision) ' / ' num2str(array(2),string_precision) ' / ' num2str(array(3),string_precision)];
                PropertyMean(line_table,current_phase) = {str_mean}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
            line_table=line_table+1;
            Propertyname(line_table) = {'Bruggeman exponent 1/2/3'}; % Name
            Propertymethod(line_table) = {'Tau factor'}; % Name
            for current_phase=1:1:number_phase
                T = datamat.Results_tortuosityfactor.Table_TortuosityFactor_taufactor.phase(current_phase).table; % Table
                array = T{:,4};
                str_mean=[num2str(array(1),string_precision) ' / ' num2str(array(2),string_precision) ' / ' num2str(array(3),string_precision)];
                PropertyMean(line_table,current_phase) = {str_mean}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
            line_table=line_table+1;
            Propertyname(line_table) = {'Normalized Eff. Diff. Coeff. 1/2/3'}; % Name
            Propertymethod(line_table) = {'Tau factor'}; % Name
            for current_phase=1:1:number_phase
                T = datamat.Results_tortuosityfactor.Table_TortuosityFactor_taufactor.phase(current_phase).table; % Table
                array = T{:,5};
                str_mean=[num2str(array(1),string_precision) ' / ' num2str(array(2),string_precision) ' / ' num2str(array(3),string_precision)];
                PropertyMean(line_table,current_phase) = {str_mean}; % Mean
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
        end
        
        % Case Connectivity
        if strcmp(current_result,'Results_connectivity.mat')
            pathname = [folder current_result];
            datamat = load(pathname);
            line_table=line_table+1;
            Propertyname(line_table) = {'Main cluster connectivity'}; % Name
            Propertymethod(line_table) = {'6-connected'}; % Name
            T = datamat.Results_connectivity.Table_connectivity; % Table
            for current_phase=1:1:number_phase
                array = T{current_phase,2:end};
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMean(line_table,current_phase) = {num2str(array(1),string_precision)}; % Mean
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[%]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'[%]'}; % Standard deviation in percents of the mean
            end
            line_table=line_table+1;
            Propertyname(line_table) = {'Face2Face connectivity 1/2/3'}; % Name
            Propertymethod(line_table) = {'6-connected'}; % Name
            T = datamat.Results_connectivity.Table_face2face; % Table
            for current_phase=1:1:number_phase
                array = T{current_phase,2:end};
                str_mean=[num2str(array(1),string_precision) ' / ' num2str(array(4),string_precision) ' / ' num2str(array(7),string_precision)];
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMean(line_table,current_phase) = {str_mean}; % Mean
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[%]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
            line_table=line_table+1;
            Propertyname(line_table) = {'From 1st face connectivity 1/2/3'}; % Name
            Propertymethod(line_table) = {'6-connected'}; % Name
            T = datamat.Results_connectivity.Table_fromface1; % Table
            for current_phase=1:1:number_phase
                array = T{current_phase,2:end};
                str_mean=[num2str(array(1),string_precision) ' / ' num2str(array(4),string_precision) ' / ' num2str(array(7),string_precision)];
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMean(line_table,current_phase) = {str_mean}; % Mean
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[%]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
            line_table=line_table+1;
            Propertyname(line_table) = {'From last face connectivity 1/2/3'}; % Name
            Propertymethod(line_table) = {'6-connected'}; % Name
            T = datamat.Results_connectivity.Table_fromface2; % Table
            for current_phase=1:1:number_phase
                array = T{current_phase,2:end};
                str_mean=[num2str(array(1),string_precision) ' / ' num2str(array(4),string_precision) ' / ' num2str(array(7),string_precision)];
                PropertyMin(line_table,current_phase) = {'n/a'}; % Min
                PropertyMean(line_table,current_phase) = {str_mean}; % Mean
                PropertyMax(line_table,current_phase) = {'n/a'}; % Max
                PropertyStd(line_table,current_phase) = {'n/a'}; % Standard deviation
                PropertyUnit(line_table,current_phase) = {'[%]'}; % Unit
                PropertyStdpercent(line_table,current_phase) = {'n/a'}; % Standard deviation in percents of the mean
            end
        end

    end
    
    if exist('Propertyname','var') == 1
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
            if (OPTIONS.displaytext==1)
                fprintf('       For the phase: %s:\n\n',INFO.phase(current_phase).name);
                disp(Table_summary_results.phase(current_phase).table); disp ' ';
            end
            
            Fig = figure; % Create figure
            Fig.Name= sprintf('Results of phase %s',INFO.phase(current_phase).name); % Figure name
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
            if OPTIONS.save_fig == true % Save figure
                filename= sprintf('Results_of_%s',INFO.phase(current_phase).filename);
                function_savefig(Fig, folder, filename, OPTIONS); % Call function
            end
            if OPTIONS.closefigureaftercreation == true
                close(Fig); % Do not keep open figures
            end
        end
        
        % Save
        if OPTIONS.save_resultsmat == true
            save([folder 'Table_summary_results'],'Table_summary_results')
        end
        if OPTIONS.save_xls==true
            filename = 'Table_summary_results'; % Filename without extension
            % Prepare the data
            clear DATA_writetable
            sheetnumber=0;
            for current_phase=1:1:number_phase
                sheetnumber=sheetnumber+1;
                DATA_writetable.sheet(sheetnumber).name=INFO.phase(current_phase).filename;
                DATA_writetable.sheet(sheetnumber).table=Table_summary_results.phase(current_phase).table;
            end
            % Save function
            Function_Writetable(folder,filename,DATA_writetable)
        end
    end
    
end

end

