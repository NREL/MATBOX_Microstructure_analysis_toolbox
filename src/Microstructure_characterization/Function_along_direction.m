function [Table_evolution] = Function_along_direction(p)

%% Calculation
for current_phase=1:1:p.number_phase
    data_ = p.data(current_phase).(p.field);
    for current_direction=1:1:p.number_dimension
        % Initialization
        % :,1,: position
        % :,2 : min
        % :,3,: mean
        % :,4,: max
        % :,5,: std
        Alongdirection.phase(current_phase).direction(current_direction).val = zeros(p.Domain_size(current_direction),4);
        Alongdirection.phase(current_phase).direction(current_direction).val(:,1) = (1:1:p.Domain_size(current_direction))*p.voxel_size;
        for position=1:1:p.Domain_size(current_direction)
            % Current slice of data
            if current_direction==1
                slice_ = data_(position,:,:);
            elseif current_direction==2
                slice_ = data_(:,position,:);
            elseif current_direction==3
                slice_ = data_(:,:,position);
            end
            different_values =  unique(slice_); % Get all valuies
            
            % Remove specific values from the analysis
            for k=1:1:length(p.ignore_value)
                different_values(different_values==p.ignore_value(k)) = [];
            end
            
            if ~isempty(different_values)
                % Minimum
                Alongdirection.phase(current_phase).direction(current_direction).val(position,2) = min(different_values);
                
                % Maximum
                Alongdirection.phase(current_phase).direction(current_direction).val(position,4) = max(different_values);
                
                % Weighed values
                % Initialisation
                values=zeros(length(different_values),2);
                for current_val=1:1:length(different_values)
                    values(current_val,1)=different_values(current_val);
                    values(current_val,2)=sum(sum( slice_== different_values(current_val)));
                end
                % Mean
                Alongdirection.phase(current_phase).direction(current_direction).val(position,3) = sum(values(:,1).*values(:,2))/sum(values(:,2));
                
                % Standard deviation
                % The weighted starndard deviation formula is
                % sqrt( sum(wi((xi-<x>)^2)  / ( (n-1)/n * sum wi ) )
                % With wi the weight of the xi, and <x> the weighted mean (mean_size)
                wi = values(:,2);
                xi = values(:,1);
                n = length(xi);
                mean_ = Alongdirection.phase(current_phase).direction(current_direction).val(position,3);
                Alongdirection.phase(current_phase).direction(current_direction).val(position,5) = sqrt( sum( wi.*((xi-mean_).^2)) / ( (n-1)/n * sum(wi)  ));
            else
                Alongdirection.phase(current_phase).direction(current_direction).val(position,2) = NaN;
                Alongdirection.phase(current_phase).direction(current_direction).val(position,3) = NaN;
                Alongdirection.phase(current_phase).direction(current_direction).val(position,4) = NaN;
                Alongdirection.phase(current_phase).direction(current_direction).val(position,5) = NaN;
            end
            
        end
    end
end


%% Table
for current_phase=1:1:p.number_phase
    for current_direction=1:1:p.number_dimension
        data_ = Alongdirection.phase(current_phase).direction(current_direction).val;
        if p.ignore_min
            Table_evolution.direction(current_direction).phase(current_phase).table = array2table([data_(:,1) data_(:,3) data_(:,4) data_(:,5)],...
                'VariableNames',p.Variable_name_table);
        else
            Table_evolution.direction(current_direction).phase(current_phase).table = array2table([data_(:,1) data_(:,2) data_(:,3) data_(:,4) data_(:,5)],...
                'VariableNames',p.Variable_name_table);            
        end
    end
end


%% Save tables
if p.OPTIONS.save_xls==true
    for current_phase = 1:1:p.number_phase
        % Prepare the data
        clear DATA_writetable
        filename = [p.Table_filename '_along_directions_' p.INFO.phase(current_phase).filename]; % Filename without extension
        sheet_=0;
        for direction=1:1:p.number_dimension
            sheet_=sheet_+1;
            % Data : Evolution of volume fraction along direction 1
            DATA_writetable.sheet(sheet_).name =p.INFO.direction(direction).filename;
            DATA_writetable.sheet(sheet_).table=Table_evolution.direction(direction).phase(current_phase).table;
        end
        % Save function
        Function_Writetable(p.Current_folder,filename,DATA_writetable)
    end
end

%% Figures
scrsz = get(0,'ScreenSize'); % Screen resolution
for current_phase=1:1:p.number_phase
    Fig = figure; % Create figure
    Fig.Name= [p.figure_name p.INFO.phase(current_phase).name]; % Figure name
    Fig.Color='white'; % Background colour
    set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)*p.number_dimension/3 scrsz(4)*1/2]); % Full screen figure
    for current_direction=1:1:p.number_dimension % Iterate over axe
        sub_axes=subplot(1,p.number_dimension,current_direction,'Parent',Fig);
        hold(sub_axes,'on'); % Active subplot
        h_title=title (p.axe_title); % Set title font
        % Plot graphs
        data_ = Alongdirection.phase(current_phase).direction(current_direction).val;
        x_=data_(:,1);
        y_min = data_(:,2);
        y_mean = data_(:,3);
        y_max = data_(:,4);
        y_std = data_(:,5);
        
        % Mean
        h_mean=plot(x_,y_mean); % For the legend order
        % Extremums
        if ~p.ignore_min
            h_min=plot(x_,y_min);
        end
        h_max=plot(x_,y_max);
        
        % Colors, thickness, markers
        set(h_mean, 'Color', 'k','LineWidth',p.OPTIONS.Linewidth);
        if ~p.ignore_min
            set(h_min, 'Color', p.INFO.phase(current_phase).color,'LineStyle','--','LineWidth',p.OPTIONS.Linewidth);
        end
        set(h_max, 'Color', p.INFO.phase(current_phase).color,'LineStyle','--','LineWidth',p.OPTIONS.Linewidth);
        % Mean with error bar (+- standard deviation)
        h_mean_witherrorbar = errorbar(x_,y_mean,y_std);
        set(h_mean_witherrorbar, 'Color', p.INFO.phase(current_phase).color,'LineWidth',p.OPTIONS.Linewidth);
        h_mean=plot(x_,y_mean); % Plot over the other
        set(h_mean, 'Color', 'k','LineWidth',p.OPTIONS.Linewidth);
        
        % Axis label
        t_ = xlabel(' ');
        t_1 = sprintf('Position along %s ',p.INFO.direction(current_direction).name);
        t_2 = '(\mum)';
        t_.String= [t_1 t_2]; % Sprintf does not accept greek characters
        t_ = ylabel(p.ylabel);
        % Legend
        if p.ignore_min
            str_legend(1).name = ['Mean ' p.legendname ', ' num2str(p.mean_val(current_phase,1),'%1.3f') p.ylabel_unit];
            str_legend(2).name = ['Maximun ' p.legendname];
        else
            str_legend(1).name = ['Mean ' p.legendname ', ' num2str(p.mean_val(current_phase,1),'%1.3f') p.ylabel_unit];
            str_legend(2).name = ['Extremums ' p.legendname];
        end
        h_legend = legend(sub_axes,str_legend.name,'Location','best');
        % - Grid
        if strcmp(p.OPTIONS.grid,'on')
            grid(sub_axes,'on'); % Display grid
            set(sub_axes,'XMinorGrid',p.OPTIONS.minorgrid,'YMinorGrid',p.OPTIONS.minorgrid); % Display grid for minor thicks
        end
        set(sub_axes,'FontName',p.OPTIONS.fontname,'FontSize',p.OPTIONS.Fontsize_axe); % Fontname and fontsize
        h_title.FontSize = p.OPTIONS.Fontsize_title; % Set title fontsize
        h_legend.FontSize = p.OPTIONS.Fontsize_legend; % Set title fontsize
        hold(sub_axes,'off'); % Relase figure
    end
    sgtitle(Fig,[p.figure_title ' ' p.INFO.phase(current_phase).name] ,'FontWeight','bold','FontSize',p.OPTIONS.Fontsize_title+2,'FontName',p.OPTIONS.fontname);
    if p.OPTIONS.save_fig == true % Save figure
        function_savefig(Fig, p.Current_folder, [p.figure_filename p.INFO.phase(current_phase).name], p.OPTIONS); % Call function
    end
    if p.OPTIONS.closefigureaftercreation == true
        close(Fig); % Do not keep open figures
    end
end


end

