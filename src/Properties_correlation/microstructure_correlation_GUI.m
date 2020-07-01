function microstructure_correlation_GUI
%Graphic user interface of the microstructure correlation toolbox module

close all
clear all

%%
%% GUI DISPLAY OPTIONS
%%

font_name_GUI ='Times New Roman'; % Font
% Font size
font_size_small_GUI =10;
font_size_large_GUI =14;
% Background and text color
background_description_tab = [0 0.5 0];
ForegroundColor_description_tab = [1 1 1];
error_message_red = [1 0 0];
error_message_orange = [1 0.64 0];


%%
%% PRESELECTED OPTIONS AND VARIABLES
%%

scrsz = get(0,'ScreenSize');
if ispc
    folder_separation = '\';
else
    folder_separation = '/';
end

% Initialize variables here is similar to define them as global variables
allresults={};
number_results = [];
number_parameter_id = 100;
number_maxvolumes = 100;
volume_name=[];
number_volume_loaded=0;
number_volume_loaded_tmp=0;
color_phase = [get(0, 'DefaultAxesColorOrder'); rand(number_maxvolumes,3);]; % Color attributed to each phase. One row per phase.
marker_order = {'o','s','d','^','v','>','<','p','h','x','+','*'}';
while length(marker_order)<number_maxvolumes
    marker_order=[marker_order;marker_order];
end

plot_per_volume = true;
plot_per_group = false;
Save_folder = 'None';
Save_filename = 'None';

axe_fontsize=10;
legend_fontsize=8;
diagonale_fontsize=14;
title_fontsize=25;
figure_Linewidth=2;
figure_Markersize=10;
figure_grid='on';
figure_Minorgrid='off';


%%
%% CREATE GUI: OBJECT
%%

%% MAIN FIGURE
main_figure = figure; % Create figure
main_figure.Name= 'Microstructure correlation'; % Set name
main_figure.NumberTitle='off'; % Remove number from title name
main_figure.Color='white'; % Background colour
main_figure.MenuBar='none'; % Remove menubar and toolbar
main_figure.ToolBar='none';
main_figure.Units='normalized'; % Set unit
lfigure=[0.1 0.1 0.7 0.8];
main_figure.Position=lfigure; % Set Position

volume_figure = figure; % Create figure
volume_figure.NumberTitle='off'; % Remove number from title name
volume_figure.Color='white'; % Background colour
volume_figure.MenuBar='none'; % Remove menubar and toolbar
volume_figure.ToolBar='none';
volume_figure.Units='normalized'; % Set unit
lvolfigure=[0.2 0.2 0.4 0.7];
volume_figure.Position=lvolfigure; % Set Position
volume_figure.Visible = 'off';


Text_module_name = uicontrol('Parent', main_figure, 'Style', 'text','Units','normalized','Position',[0 0.97 1 0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',background_description_tab,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Microstructure correlation');

Text_error = uicontrol('Parent', main_figure, 'Style', 'text','Units','normalized','Position',[0.01 0.01 0.98 0.03],...
    'FontName',font_name_GUI,'FontSize',font_size_large_GUI,'BackgroundColor',error_message_red,'ForegroundColor',ForegroundColor_description_tab,...
    'String','Error message','Visible','off','Enable','off');

l=[0.5-0.025-0.1 0.9 0.1 0.04];
Button_add_volume = uicontrol('Parent', main_figure, 'Style', 'pushbutton', 'String', 'Load volume',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', l,'Callback',{@add_volume_Callback});

l(3)=l(1)-0.01-0.025; l(1)=0.025;
str='Assign parameters with a name in the table below and choose which ones you want to plot. Then load volume main folders one after one.';
Text_parameter_Id = uicontrol('Parent', main_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String',str,...
    'HorizontalAlignment','left','Units','normalized','Position', l);

ltable_parameter_id = [0.025 0.3 0.45 0.585];
table_parameter_Id = uitable('Parent', main_figure,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',ltable_parameter_id);
table_parameter_Id.ColumnName = {'Parameter Id','Parameter name in figure','Unit name','Number of values','Plot'}; % Column name
table_parameter_Id.ColumnEditable = [false true true false true]; % Select column editable
ltable=scrsz(3)*lfigure(3)*ltable_parameter_id(3);
table_parameter_Id.ColumnWidth = {ltable*0.2, ltable*0.35, ltable*0.125, ltable*0.2, ltable*0.1}; % Auto width
table_parameter_Id.RowName = []; % Remove row name
% Initialize data
for k_id=1:1:number_parameter_id
    column_1(k_id,1) = {num2str(k_id)};
    column_2(k_id,1) = {['Parameter ' num2str(k_id)]};
    column_3(k_id,1) = {['']};
    column_4(k_id,1) = {'0'};
    column_5(k_id,1) = {false};
    parameters_input(k_id).name = column_2(k_id,1);
    parameters_input(k_id).to_plot = 0;
    parameters_input(k_id).value = ones(1,number_maxvolumes)*NaN;
end
table_parameter_Id.Data=[column_1 column_2 column_3 column_4 column_5];

parameters_input(1).volumes_name = cell(1,number_maxvolumes);
for k_vol=1:1:number_maxvolumes
    parameters_input(1).volumes_name(k_vol)={'None'};
end
parameters_input_tmp = parameters_input;

str='Assign volume to a group id and choose which ones you want to plot. Similarly, chooses which group you want to plot.';
Text_volumegroup = uicontrol('Parent', main_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String',str,...
    'HorizontalAlignment','left','Units','normalized','Position', [0.525 0.9 0.45 0.04],'visible','off');

ltable_volumes = [0.525 0.505 0.45 0.38];
table_volumes = uitable('Parent', main_figure,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',ltable_volumes,'visible','off','enable','off');
table_volumes.ColumnName = {'Volume name','Legend in figure','Group Id','Parameters','Plot'}; % Column name
table_volumes.ColumnEditable = [false true true false true]; % Select column editable
ltable=scrsz(3)*lfigure(3)*ltable_volumes(3);
table_volumes.ColumnWidth = {ltable*0.29, ltable*0.29, ltable*0.1, ltable*0.2, ltable*0.1 }; % Auto width
table_volumes.RowName = []; % Remove row name
% Initialize data
for k_vol=1:1:number_maxvolumes
    column_vol_1(k_vol,1) = {'None'};
    column_vol_2(k_vol,1) = {'None'};
    column_vol_3(k_vol,1) = {[k_vol]};
    column_vol_4(k_vol,1) = {[0]};
    column_vol_5(k_vol,1) = {false};
end
table_volumes.Data=[column_vol_1 column_vol_2 column_vol_3 column_vol_4 column_vol_5];

ltable_group = [0.525 0.3 0.45 0.2];
table_group = uitable('Parent', main_figure,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',ltable_group,'visible','off','enable','off');
table_group.ColumnName = {'Group id','Legend in figure','Marker','Filled','R','G','B','Plot'}; % Column name
table_group.ColumnEditable = [false true true true true true true true]; % Select column editable
ltable=scrsz(3)*lfigure(3)*ltable_group(3);
table_group.ColumnWidth = {ltable*0.19, ltable*0.19, ltable*0.1, ltable*0.1, ltable*0.1, ltable*0.1, ltable*0.1, ltable*0.1 }; % Auto width
table_group.RowName = []; % Remove row name
% Initialize data
for k_group=1:1:number_maxvolumes
    column_group_1(k_group,1) = {[k_group]};
    column_group_2(k_group,1) = {['Group ' num2str(k_group)]};
    column_group_3(k_group,1) = marker_order(k_group);
    column_group_4(k_group,1) = {false};
    column_group_5(k_group,1) = {[color_phase(k_group,1)]};
    column_group_6(k_group,1) = {[color_phase(k_group,2)]};
    column_group_7(k_group,1) = {[color_phase(k_group,3)]};
    column_group_8(k_group,1) = {false};
end
column_group_8(1,1) = {true};
table_group.Data=[column_group_1 column_group_2 column_group_3 column_group_4 column_group_5 column_group_6 column_group_7 column_group_8];


Text_volumefigure_instructions = uicontrol('Parent', volume_figure, 'Style', 'text','Units','normalized','Position',[0.05 0.915 0.9 0.085],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'HorizontalAlignment','left',...
    'String','Select the phase you want to correlate. Then assign a parameter Id to the parameter you want to correlate (>=1). An Id equals to 0 means the parameter will not be correlated. Then, click on the button ''Save'' if you want to repeat with another phase of this volume. Otherwise click on the button ''Close'' to keep your changes or ''Cancel and close'' to discard them.');

Popup_selectphase = uicontrol('Parent', volume_figure,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', 'None','Units','normalized','Position', [0.05 0.865 0.9 0.05],'Callback',{@popup_selectphase_Callback});

Button_save_parameterphase = uicontrol('Parent', volume_figure, 'Style', 'pushbutton', 'String', 'Save','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.05 0.02 0.08 0.05],'Callback',{@Button_save_parameterphase_Callback});
Button_close_parameterphase = uicontrol('Parent', volume_figure, 'Style', 'pushbutton', 'String', 'Close','enable','off',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.14 0.02 0.18 0.05],'Callback',{@Button_close_parameterphase_Callback});
Button_cancelclose_parameterphase = uicontrol('Parent', volume_figure, 'Style', 'pushbutton', 'String', 'Cancel and close',...
    'FontSize',font_size_small_GUI,'FontName',font_name_GUI,'Units','normalized','Position', [0.33 0.02 0.18 0.05],'Callback',{@Button_cancelclose_parameterphase_Callback});

Text_error_parameterphase = uicontrol('Parent', volume_figure, 'Style', 'text','Units','normalized','Position',[0.53 0.02 0.42 0.05],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'BackgroundColor',error_message_red,'ForegroundColor',ForegroundColor_description_tab,'String','Error message','Visible','off');

ltable_parameter_choice = [0.05 0.345 0.9 0.52];
table_parameter_choice = uitable('Parent', volume_figure,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',ltable_parameter_choice,'CellEditCallback',@cellsection_tableparameterschoice_Callback, 'CellSelectionCallback',@cellsection_tableparameterschoice_Callback );
table_parameter_choice.ColumnName = {'Parameter name','Value','Parameter Id'}; % Column name
table_parameter_choice.ColumnEditable = [false true true]; % Select column editable
ltable=scrsz(3)*lvolfigure(3)*ltable_parameter_choice(3);
table_parameter_choice.ColumnWidth = {ltable*0.56, ltable*0.2, ltable*0.2}; % Auto width
table_parameter_choice.RowName = []; % Remove row name

ltable_parameter_additionalchoice = [0.05 0.1 0.9 0.225];
table_parameter_additionalchoice = uitable('Parent', volume_figure,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',ltable_parameter_additionalchoice, 'CellEditCallback',@cellsection_tableparameterschoice_Callback, 'CellSelectionCallback',@cellsection_tableparameterschoice_Callback );
table_parameter_additionalchoice.ColumnName = {'Additional parameter you may manually add','Value','Parameter Id'}; % Column name
table_parameter_additionalchoice.ColumnEditable = [true true true]; % Select column editable
ltable=scrsz(3)*lvolfigure(3)*ltable_parameter_additionalchoice(3);
table_parameter_additionalchoice.ColumnWidth = {ltable*0.56, ltable*0.2, ltable*0.2}; % Auto width
table_parameter_additionalchoice.RowName = []; % Remove row name


y_ = 0.3;
annotation(main_figure,'line',[0.025 0.975],[y_-0.025 y_-0.025]);

Button_createfigure = uicontrol('Parent', main_figure, 'Style', 'pushbutton', 'String', 'Create figure',...
    'FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','Units','normalized','BackgroundColor',[0.1961    0.8039    0.1961],'Position', [0.975-0.2 y_-0.1 0.2 0.06],'Callback',{@button_createfigure_Callback},'enable','off');

Button_visualizematrix = uicontrol('Parent', main_figure, 'Style', 'pushbutton', 'String', 'Visualize matrix',...
    'FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','Units','normalized','BackgroundColor',[0.6784    1.0000    0.1843],'Position', [0.975-0.2 y_-0.16 0.2 0.06],'Callback',{@button_visualize_matrix_Callback},'enable','off');

Button_exportmatrix = uicontrol('Parent', main_figure, 'Style', 'pushbutton', 'String', 'Export matrix',...
    'FontSize',font_size_large_GUI,'FontName',font_name_GUI,'FontWeight','bold','Units','normalized','Position', [0.975-0.2 y_-0.22 0.2 0.06],'Callback',{@button_export_matrix_Callback},'enable','off');

Button_select_savefolder = uicontrol('Parent', main_figure,'Style', 'pushbutton','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Select save folder',...
    'HorizontalAlignment','left','Units','normalized','Position', [0.025 y_-0.08 0.1 0.04],'Callback',{@button_select_savefolder_Callback});
Text_savefolder = uicontrol('Parent', main_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',Save_folder,...
    'HorizontalAlignment','left','BackgroundColor','w','Units','normalized','Position', [0.13 y_-0.08 0.345 0.04]);

Text_filename_instructions = uicontrol('Parent', main_figure,'Style', 'text','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'FontWeight','bold','String','Write filename',...
    'HorizontalAlignment','center','Units','normalized','Position', [0.025 y_-0.12 0.1 0.04]);
Text_filename = uicontrol('Parent', main_figure,'Style', 'Edit','FontSize',font_size_small_GUI,'FontName',font_name_GUI,'String',Save_filename,...
    'HorizontalAlignment','left','BackgroundColor','w','Units','normalized','Position', [0.13 y_-0.12 0.345 0.04]);

checkbox_plot_pervolume = uicontrol('Parent', main_figure, 'Style', 'checkbox','Units','normalized','Position',[0.025 y_-0.16 0.1 0.04],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Plot per volume','Value',plot_per_volume,'Callback',{@checkbox_plotpervolume_Callback});
checkbox_plot_pergroup = uicontrol('Parent', main_figure, 'Style', 'checkbox','Units','normalized','Position',[0.125 y_-0.16 0.1 0.04],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold','BackgroundColor','w',...
    'String','Plot per group','Value',plot_per_group,'Callback',{@checkbox_plotpergroup_Callback});

text_displayoption_fontname = uicontrol('Parent', main_figure, 'Style', 'text','Units','normalized','Position',[0.225 y_-0.16 0.1 0.04],...
    'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'FontWeight','bold',...
    'String','Select fontname','HorizontalAlignment','left');
Popup_displayoptions_fontname = uicontrol('Parent', main_figure,'Style', 'popup','FontSize',font_size_small_GUI,'FontName',font_name_GUI,...
    'String', listfonts,'Value',203,'Units','normalized','Position', [0.325 y_-0.16 0.15 0.04]);

ltable_figureoption = [0.025 y_-0.24 0.6 0.07];
table_figureoption = uitable('Parent', main_figure,'FontName',font_name_GUI,'FontSize',font_size_small_GUI,'Units','normalized','Position',ltable_figureoption);
table_figureoption.ColumnName = {'Axe fontsize','Legend fontsize','Diagonale fontsize','Title fontsize','Linewidth','MarkerSize','Grid','Minorgird'}; % Column name
ltable=scrsz(3)*lfigure(3)*ltable_figureoption(3);
table_figureoption.ColumnWidth = {ltable*0.125, ltable*0.125, ltable*0.125, ltable*0.125, ltable*0.125,ltable*0.125, ltable*0.12, ltable*0.12}; % Auto width
table_figureoption.ColumnEditable = [true true true true true true true true]; % Select column editable

table_figureoption.RowName = []; % Remove row name
% Initialize data
if strcmp(figure_grid,'on')
    figure_grid_boolean = true;
    if strcmp(figure_Minorgrid,'on')
        figure_minorgrid_boolean = true;
    else
        figure_minorgrid_boolean = false;
    end
else
    figure_grid_boolean = false;
    figure_minorgrid_boolean = false;
end
table_figureoption.Data=[{axe_fontsize} {legend_fontsize} {diagonale_fontsize} {title_fontsize} {figure_Linewidth} {figure_Markersize} {figure_grid_boolean} {figure_minorgrid_boolean}];


%%
%% CALLBACK FUNCTIONS
%%

    function checkbox_plotpervolume_Callback(~,~)
        plot_per_volume = logical(checkbox_plot_pervolume.Value);
        plot_per_group = ~plot_per_volume;
        checkbox_plot_pergroup.Value = plot_per_group;
    end
    function checkbox_plotpergroup_Callback(~,~)
        plot_per_group = logical(checkbox_plot_pergroup.Value);
        plot_per_volume = ~plot_per_group;
        checkbox_plot_pervolume.Value = plot_per_volume;
    end

    function button_select_savefolder_Callback(~,~)
        % Set string of the dialog box
        str_dialogbox = 'Select the save folder where figure and results will be saved';
        % Open dialog box to choose file path
        Save_folder_tmp = uigetdir(str_dialogbox);
        if Save_folder_tmp==0
            Save_folder = 'None';
        else
            Save_folder = [Save_folder_tmp folder_separation];
        end
        set(Text_savefolder,'String',Save_folder);
    end

    function button_createfigure_Callback(~,~)
        r =function_prepare_inputdata();
        
        parameters_figure.Save_folder = Save_folder;
        parameters_figure.save_fig = true;
        parameters_figure.plot_per_volume = plot_per_volume;
        parameters_figure.plot_per_group = plot_per_group;
        parameters_figure.figurename = Text_filename.String;
        parameters_figure.matrix = r.matrix;
        parameters_figure.parameternames = r.parameternames;
        parameters_figure.parameterunitnames = r.parameterunitnames;
        parameters_figure.volumenames = r.volumenames;
        parameters_figure.volumegroup = r.volumegroup;
        parameters_figure.group = r.group;
        parameters_figure.correlation = r.Kendalltauranking;
        parameters_figure.correlation_group = r.Kendalltauranking_averagedpergroup;
        
        parameters_figure.marker_order = marker_order;
        parameters_figure.fontname = char(Popup_displayoptions_fontname.String(Popup_displayoptions_fontname.Value));
        parameters_figure.Fontsize_axe = cell2mat(table_figureoption.Data(1));
        parameters_figure.Fontsize_legend = cell2mat(table_figureoption.Data(2));
        parameters_figure.font_diagonale = cell2mat(table_figureoption.Data(3));
        parameters_figure.Fontsize_title = cell2mat(table_figureoption.Data(4));
        parameters_figure.Linewidth = cell2mat(table_figureoption.Data(5));
        parameters_figure.MarkerSize = cell2mat(table_figureoption.Data(6));
        if cell2mat(table_figureoption.Data(7))
            parameters_figure.grid = 'on';
            if cell2mat(table_figureoption.Data(8))
                parameters_figure.minorgrid = 'on';
            else
                parameters_figure.minorgrid = 'off';
            end
        else
            parameters_figure.grid = 'off';
            parameters_figure.minorgrid = 'off';
        end
        
        function_figure_correlation_matrix(parameters_figure)
    end

    function button_visualize_matrix_Callback(~,~)
        r = function_prepare_inputdata();
        disp '     PARAMETERS';
        disp '     ----------';
        disp(r.T)
        disp '     CORRELATION: Kendall rank correlation coefficient';
        disp '     -------------------------------------------------';
        disp(r.T_Kendalltauranking)
        disp '     CORRELATION: Kendall rank correlation coefficient, averaged per group';
        disp '     ----------------------------------------------------------------------';
        disp(r.T_Kendalltauranking_averagedpergroup)        
    end

    function button_export_matrix_Callback(~,~)
        r = function_prepare_inputdata();
        % .mat
        table_parameter_Id_Data = table_parameter_Id.Data;
        table_volumes_Data = table_volumes.Data;
        table_group_Data = table_group.Data;
        
        Save_filename = function_remove_emptyandspecialcharacter_string(Text_filename.String);
        
        save([Save_folder Save_filename] ,'table_parameter_Id_Data','table_volumes_Data','table_group_Data');
        % Excel data sheet
        clear DATA_writetable
        DATA_writetable.sheet(1).name = 'Parameters';
        DATA_writetable.sheet(1).table = r.Texport;
        DATA_writetable.sheet(2).name = 'Kendall_correlation';
        DATA_writetable.sheet(2).table = r.T_Kendalltauranking_export;
        DATA_writetable.sheet(3).name = 'Kendall_correlation_group_avg';
        DATA_writetable.sheet(3).table = r.T_Kendalltauranking_averagedpergroup_export;        
        % Save function
        Function_Writetable(Save_folder,Save_filename,DATA_writetable)
    end

    function [r] = function_prepare_inputdata
        matrix_correlation_tmp = zeros(number_maxvolumes,number_parameter_id);
        for k_vol=1:1:number_maxvolumes
            for k_id=1:1:number_parameter_id
                matrix_correlation_tmp(k_vol,k_id) = parameters_input(k_id).value(k_vol);
            end
        end
        % Remove parameters with 0 values and volumes with 0 parameters
        %n_values = str2num(cell2mat(table_parameter_Id.Data(:,4)));
        n_values = str2double(table_parameter_Id.Data(:,4));
        zeros_location = find(n_values==0);
        matrix_correlation_tmp(:,zeros_location)=NaN;
        n_parameters = cell2mat(table_volumes.Data(:,4));
        zeros_location = find(n_parameters==0);
        matrix_correlation_tmp(zeros_location,:)=NaN;
        
        % Find parameters and volumes for which the plot checkbox is true
        plot_parameter = cell2mat(table_parameter_Id.Data(:,5)) ;
        zeros_location = find(plot_parameter==0);
        matrix_correlation_tmp(:,zeros_location)=NaN;
        plot_volumes = cell2mat(table_volumes.Data(:,5)) ;
        zeros_location = find(plot_volumes==0);
        matrix_correlation_tmp(zeros_location,:)=NaN;
        
        idx = find(~isnan(matrix_correlation_tmp));
        [I,J] = ind2sub(size(matrix_correlation_tmp),idx);
        Vol = unique(I); Par = unique(J);
        r.matrix = matrix_correlation_tmp(Vol , Par );
        r.parameternames = table_parameter_Id.Data(Par,2);
        r.parameterunitnames = table_parameter_Id.Data(Par,3);
        r.volumenames = table_volumes.Data(Vol,2);
        r.volumegroup = table_volumes.Data(Vol,3);
        uniquegroup = unique(cell2mat(r.volumegroup));
        r.group = table_group.Data(uniquegroup,:);
        
        % Create Table
        [n_vol, n_par] = size(r.matrix);
        varNames = cell(1,n_par);
        rowNames = cell(n_vol,1);
        for k_par = 1:1:n_par
            varNames(k_par)={ function_remove_emptyandspecialcharacter_string(char(r.parameternames(k_par))) };
        end
        for k_vol = 1:1:n_vol
            rowNames(k_vol,1)={ function_remove_emptyandspecialcharacter_string(char(r.volumenames(k_vol))) };
        end
        rowNames(k_vol+1,1) = {'Unit'};
        
        %T = array2table(r.PLOT.matrix, 'VariableNames',varNames,'RowNames',rowNames )
        tmp = cell(n_vol+1,n_par);
        tmp(1:n_vol,:) = num2cell(r.matrix);
        tmp(end,:) =  r.parameterunitnames;
        
        tmp2 = cell(n_vol+1,n_par+1);
        tmp2(1:n_vol,2:n_par+1) = num2cell(r.matrix);
        tmp2(end,2:end) =  r.parameterunitnames;
        tmp2(1:end-1,1) =  r.volumenames;
        tmp2(end,1) =  {'Unit'};
        varNames2 = cell(n_par+1,1);
        varNames2(1,1)= {'Volume'};
        varNames2(2:end,1) = varNames;
        r.T = cell2table(tmp, 'VariableNames',varNames,'RowNames',rowNames );
        r.Texport = cell2table(tmp2, 'VariableNames',varNames2);
        
        % Correlation matrix (Kendall tau ranking)
        r.Kendalltauranking = zeros(n_par,n_par);
        for row=1:1:n_par-1
            for column=row+1:1:n_par
                par1 = r.matrix(:,row);
                par2 = r.matrix(:,column);
                r.Kendalltauranking(row,column) = Function_Kendall_tau_ranking(par1,par2);
                r.Kendalltauranking(column,row) = r.Kendalltauranking(row,column); % Symmetric
            end
        end
        r.T_Kendalltauranking = cell2table(num2cell(r.Kendalltauranking), 'VariableNames',varNames,'RowNames',varNames );
        tmp2 = cell(n_par,n_par+1);
        tmp2(:,2:n_par+1) = num2cell(r.Kendalltauranking);
        tmp2(:,1) =  varNames;
        varNames2 = cell(n_par+1,1);
        varNames2(1,1)= {'Parameters'};
        varNames2(2:end,1) = varNames;
        r.T_Kendalltauranking_export = cell2table(tmp2, 'VariableNames',varNames2);
        
        % Correlation maxtrix (Kendall tau ranking), per group, averaged
        plot_group = find(cell2mat(r.group(:,end)));
        number_group = length(plot_group);
        for k_group=1:1:number_group
            current_group = cell2mat (r.group(k_group,1));
            group(k_group).idx = find( cell2mat(r.volumegroup) == current_group );
        end
        r.Kendalltauranking_averagedpergroup = zeros(n_par,n_par);
        for row=1:1:n_par-1
            for column=row+1:1:n_par
                tmp_kendall=zeros(number_group,1);
                for k_group=1:1:number_group
                    x=r.matrix(group(k_group).idx,column);
                    y=r.matrix(group(k_group).idx,row);
                    tmp_kendall(k_group,1) = Function_Kendall_tau_ranking(x,y);
                end
                avg_ = nanmean([tmp_kendall]); % Kandall = NaN if only evaluated with 1 point.
                if isnan(avg_)
                    avg_=0;
                end
                r.Kendalltauranking_averagedpergroup(row,column) = avg_;
                r.Kendalltauranking_averagedpergroup(column,row) = r.Kendalltauranking_averagedpergroup(row,column); % Symmetric
            end
        end
        r.T_Kendalltauranking_averagedpergroup = cell2table(num2cell(r.Kendalltauranking_averagedpergroup), 'VariableNames',varNames,'RowNames',varNames );
        tmp2 = cell(n_par,n_par+1);
        tmp2(:,2:n_par+1) = num2cell(r.Kendalltauranking_averagedpergroup);
        tmp2(:,1) =  varNames;
        varNames2 = cell(n_par+1,1);
        varNames2(1,1)= {'Parameters'};
        varNames2(2:end,1) = varNames;
        r.T_Kendalltauranking_averagedpergroup_export = cell2table(tmp2, 'VariableNames',varNames2);       
    end


    function add_volume_Callback(~,~)
        % Set string of the dialog box
        str_dialogbox = 'Select the main folder of the volume you want to investigate';
        % Open dialog box to choose file path
        volume_main_folder = uigetdir(str_dialogbox);
        if volume_main_folder==0
            % User clicked cancel button or closed the dialog box
        else
            current_folder = [volume_main_folder folder_separation 'Correlation' folder_separation];
            % Find all mat files that start with 'Correlation'
            MyFolderInfo=dir(current_folder);
            n_files=length(MyFolderInfo);
            detected_files=0;
            for k_file=1:1:n_files
                [~,filename,ext] = fileparts(MyFolderInfo(k_file).name);
                if strcmp(ext,'.mat')
                    sub_name = filename(1:11);
                    if strcmp(sub_name,'Correlation')
                        detected_files=detected_files+1;
                        results_filenames(detected_files)={[filename '.mat']};
                    end
                end
            end
            if detected_files
                set(Text_error,'String','Error message','BackgroundColor',error_message_red,'Visible','off','Enable','off');
                set(volume_figure,'Name',['Volume selected: ' volume_main_folder],'visible','on'); % Set name
                
                tmp = find(volume_main_folder==folder_separation);
                volume_name = volume_main_folder(tmp(end)+1:length(volume_main_folder));
                
                number_results = length(results_filenames); % Number of results to summarize
                allresults={};
                k_data=0;
                for current_matfile=1:1:number_results % Loop over all saved results
                    current_result=char(results_filenames(current_matfile)); % Load result
                    pathname = [current_folder current_result];
                    datamat = load(pathname);
                    if isfield(datamat,'results_correlation')
                        data = datamat.results_correlation;
                        number_phase = length(data);
                        for current_phase=1:1:number_phase
                            phase_choice(current_phase,1)={data(current_phase).name};
                        end
                        k_data=k_data+1;
                        allresults(k_data).phase=data;
                    end
                end
                set(Popup_selectphase,'String', phase_choice);
                input_ = [{'-'} {[0]} {[0]}];
                table_parameter_choice.Data=input_;
                
                parameters_input_tmp=parameters_input;
                for k_vol=1:1:number_maxvolumes
                    if strcmp( char(parameters_input(1).volumes_name(k_vol)),volume_name) % User has selected a volume that has been already loaded before
                        number_volume_loaded_tmp = k_vol;
                        for k_id =1:1:number_parameter_id
                            parameters_input_tmp(k_id).value(number_volume_loaded_tmp) = NaN; % Reset data
                        end
                        break % Exit for loop
                    else
                        number_volume_loaded_tmp = number_volume_loaded+1; % Increment volume number
                    end
                end
                
            else
                set(Text_error,'String','Wrong folder! ...\Main_folder_to_select\Correlation\Correlation*.mat','BackgroundColor',error_message_orange,'Visible','on','Enable','on');
            end
        end
    end

    function Button_save_parameterphase_Callback(~,~)
        table_choice_data = table_parameter_choice.Data;
        table_additionalchoice_data = table_parameter_additionalchoice.Data;
        table_data = [cell2mat(table_choice_data(:,2:3)); cell2mat(table_additionalchoice_data(:,2:3))];
        idx=find( table_data(:,2)~=0 ); % Find all parameters to be correlated
        if ~isempty(idx)
            for k_idx=1:1:length(idx)
                parameter_id = table_data(idx(k_idx),2);
                parameter_value = table_data(idx(k_idx),1);
                parameters_input_tmp(parameter_id).value(number_volume_loaded_tmp) = parameter_value;
            end
            parameters_input_tmp(1).volumes_name(number_volume_loaded_tmp) = {volume_name};
        end
    end

    function Button_close_parameterphase_Callback(~,~)
        % Update
        parameters_input = parameters_input_tmp;
        number_volume_loaded = max([number_volume_loaded_tmp number_volume_loaded]);
        set(volume_figure,'visible','off')
        
        current_volume = number_volume_loaded_tmp;
        
        % Update table_parameter_Id
        counter_forthisvolume=0;
        allcolumns = table_parameter_Id.Data;
        for k_id= 1:1:number_parameter_id
            counter=0;
            for k_volume=1:1:number_maxvolumes
                if ~isnan( parameters_input(k_id).value(k_volume))
                    counter=counter+1;
                    if k_volume==current_volume
                        counter_forthisvolume = counter_forthisvolume+1;
                    end
                end
            end
            allcolumns(k_id,4) = {num2str(counter)};
            if counter>=1
                allcolumns(k_id,5) = {true};
            end
        end
        table_parameter_Id.Data=[allcolumns(:,1) allcolumns(:,2) allcolumns(:,3) allcolumns(:,4) allcolumns(:,5)];
        
        % Update table_volumes and table_group
        set(Text_volumegroup,'visible','on')
        
        allcolumns = table_volumes.Data;
        allcolumns(current_volume,1) = {volume_name};
        allcolumns(current_volume,2) = {volume_name};
        allcolumns(current_volume,3) = {[1]};
        allcolumns(current_volume,4) = {[counter_forthisvolume]};
        allcolumns(current_volume,5) = {true};
        table_volumes.Data=[allcolumns(:,1) allcolumns(:,2) allcolumns(:,3) allcolumns(:,4) allcolumns(:,5)];
        set(table_volumes,'visible','on','enable','on');
        
        % Update table group
        set(table_group,'visible','on','enable','on');
        
        if number_volume_loaded>=2
            set(Button_createfigure,'enable','on');
            set(Button_visualizematrix,'enable','on');
            set(Button_exportmatrix,'enable','on');
        end
        
    end

    function Button_cancelclose_parameterphase_Callback(~,~)
        set(volume_figure,'visible','off')
    end

    function popup_selectphase_Callback(~,~)
        % Initialize data
        current_phase = Popup_selectphase.Value;
        k_parameter=0;
        for k_result=1:1:number_results
            parameter_name = fieldnames(allresults(k_result).phase);
            for k_name=1:1:length(parameter_name)
                current_name = char(parameter_name(k_name));
                if ~strcmp(current_name,'name') % name is keyword for the phase name, it is then ignored.
                    k_parameter=k_parameter+1;
                    column1(k_parameter,1) = {current_name};
                    allvalues = struct2cell(allresults(k_result).phase(current_phase));
                    column2(k_parameter,1) = allvalues(k_name);
                    column3(k_parameter,1) = {[0]}; % By default
                end
            end
        end
        input_ = [column1 column2 column3];
        table_parameter_choice.Data=input_;
        set(Text_error_parameterphase,'String','No parameter set to be correlated. Only parameters with an Id >0 will be correlated.','BackgroundColor',error_message_orange,'Visible','on');
        set(Button_save_parameterphase,'Enable','off');
        set(Button_close_parameterphase,'Enable','off');
        
        % Initialize data
        column__1 = [{'-'} {'-'} {'-'} {'-'} {'-'} {'-'} {'-'} {'-'} {'-'} {'-'}]';
        column__23 = [{[0]} {[0]} {[0]} {[0]} {[0]} {[0]} {[0]} {[0]} {[0]} {[0]}]';
        table_parameter_additionalchoice.Data=[column__1 column__23 column__23];
    end

    function cellsection_tableparameterschoice_Callback(~,~)
        % Get back data of the table
        table_choice_data = table_parameter_choice.Data;
        table_additionalchoice_data = table_parameter_additionalchoice.Data;
        table_data = [cell2mat(table_choice_data(:,2:3)); cell2mat(table_additionalchoice_data(:,2:3))];
        stopuser=false;
        if max(max( isnan(table_data) )) == 1
            set(Text_error_parameterphase,'String','Error: do not use NaN!','BackgroundColor',error_message_red,'Visible','on');
            stopuser=true;
        else
            if max( table_data(:,end)<0 ) ==1 || min( table_data(:,end)==uint8(table_data(:,end)) )==0
                set(Text_error_parameterphase,'String','Error: Parameter Id must be positive integer!','BackgroundColor',error_message_red,'Visible','on');
                stopuser=true;
            else
                nonzero = table_data(:,end); nonzero(nonzero==0)=[];
                if isempty(nonzero)
                    stopuser=true;
                    set(Text_error_parameterphase,'String','No parameter set to be correlated. Only parameters with an Id >0 will be correlated.','BackgroundColor',error_message_orange,'Visible','on');
                else
                    if length(nonzero)>length(unique(nonzero))
                        stopuser=true;
                        set(Text_error_parameterphase,'String','Error: Parameters Id must be different.','BackgroundColor',error_message_red,'Visible','on');
                    else
                        if max(nonzero)>number_parameter_id
                            stopuser=true;
                            set(Text_error_parameterphase,'String','Error: Parameters Id must be lower or equal with the maximum parameter Id of the main figure table.','BackgroundColor',error_message_red,'Visible','on');
                        end
                    end
                end
            end
        end
        if stopuser
            set(Button_save_parameterphase,'Enable','off');
            set(Button_close_parameterphase,'Enable','off');
        else
            set(Text_error_parameterphase,'String','Error message','BackgroundColor',error_message_red,'Visible','off');
            set(Button_save_parameterphase,'Enable','on');
            set(Button_close_parameterphase,'Enable','on');
        end
    end


end