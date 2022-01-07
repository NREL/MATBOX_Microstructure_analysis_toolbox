function [] = Microstructure_basic_visualization_interface(array, p)

nargin; % Number of input variable when the function is call

if nargin == 1 % Set default name for the axis
    direction_name = {'Normal to axe 1','Normal to axe 2','Normal to axe 3'};
    customcmap = [];
else
    if isfield(p,'direction_name')
        direction_name=p.direction_name;
    else
        direction_name = {'Normal to axe 1','Normal to axe 2','Normal to axe 3'};
    end
    if isfield(p,'custom_colormap')
        customcmap=p.custom_colormap;
    else
        customcmap = [];
    end
end


%% PARAMETERS
font_name_GUI = 'Times new roman';
font_size_GUI = 12;

% Default colormap
color_phase_default = round(colororder*255);
tmp=randi(255,10000,3);
color_phase_default=[color_phase_default;tmp];
color_phase_default=color_phase_default/255; % Normalized for Matlab
choice_colormap = 'MATLAB default';


%% PHASE INFORMATION
Phase_code=unique(array); % Get phase code
number_phase=length(Phase_code); % Get number of phase
Domain_size = size(array);
% Default Colors
if nargin==2 && isfield(p,'custom_colormap')
    for current_phase=1:1:number_phase
        RGB_phase.index(current_phase).rgb = [customcmap(current_phase,1) customcmap(current_phase,2) customcmap(current_phase,3)];
    end
else
    for current_phase=1:1:number_phase
        RGB_phase.index(current_phase).rgb = [color_phase_default(current_phase,1) color_phase_default(current_phase,2) color_phase_default(current_phase,3)];
    end
end

%% FIGURE
Fig = figure; % Create figure
Fig.Name= 'Visualize microstructure';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]); % Full screen figure
% Figure for volume view
axes_2dview = axes('Parent', Fig,'Visible','off','FontName',font_name_GUI,'Units','normalized','Position', [0 0.125 1 0.85]);
% Remove tick and label
set(axes_2dview,'xtick',[],'ytick',[]);
set(axes_2dview,'xticklabel',[],'yticklabel',[]);
% Box on
box(axes_2dview,'on');
% Fit the axes box
axis(axes_2dview,'tight');
% Aspect ratio is 1:1
axis(axes_2dview,'equal');

% Slider
Slider_axes_2dview = uicontrol('Parent', Fig,'Style', 'slider','Min',1,'Max',100,'Value',1,'Units','normalized','Position', [0.3 0.05 0.4 0.04],'Callback', @slider_axes2dview_Callback,...
    'Visible','on','enable','on');
% Text (slider position)
Text_slider = uicontrol('Parent', Fig,'Style', 'text','FontSize',font_size_GUI,'FontName',font_name_GUI,'String','Position: -/-','Visible','on','enable','on',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.3 0.0 0.4 0.04]);
% Popup menu (slider direction)
Popup_slider = uicontrol('Parent', Fig,'Style', 'popup','FontSize',font_size_GUI,'FontName',font_name_GUI,...
    'String', direction_name,'Value',3,'Units','normalized','Position', [0.05 0.05 0.215 0.04],'enable','on','Visible','on','Callback', @popup_axes2dview_Callback);

direction = 3; pos_=1;
set(Text_slider,'String',['Slice: ' num2str(pos_) '/' num2str(Domain_size(direction))]);
minor_step = 1/(Domain_size(direction)-1);
major_step = max([minor_step 0.1]);
set(Slider_axes_2dview,'Min',1,'Max',Domain_size(direction),'SliderStep', [minor_step, major_step],'Value',1);

% Colormap
list_cmap = {'MATLAB default','gray','bone','copper','turbo','jet','parula','Random'};
if ~isempty(customcmap)    
    list_cmap = {'Custom','MATLAB default','gray','bone','copper','turbo','jet','parula','Random'};
end
Popup_colormap = uicontrol('Parent', Fig,'Style', 'popup','FontSize',font_size_GUI,'FontName',font_name_GUI,...
    'String',list_cmap ,'Value',1,'Units','normalized','Position', [0.725 0.05 0.215 0.04],'enable','on','Visible','on','Callback', @popup_colormap_Callback);

% Slider
    function slider_axes2dview_Callback(source,~)
        % Get position value
        pos_=round(source.Value);
        direction = Popup_slider.Value;
        % Update text
        set(Text_slider,'String',['Slice: ' num2str(pos_) '/' num2str(Domain_size(direction))]);
        % Update position array
        % Position_slice(direction)=pos_;
        % Update figure
        update_figure
    end

% Select direction
    function popup_colormap_Callback(source,~)
        % Get colormap
        color_map = char(Popup_colormap.String(source.Value));
        if strcmp(color_map,'MATLAB default')
            for current_phase=1:1:number_phase
                RGB_phase.index(current_phase).rgb = [color_phase_default(current_phase,1) color_phase_default(current_phase,2) color_phase_default(current_phase,3)];
            end
        elseif strcmp(color_map,'Random')
            color_phase_random=randi(255,10000,3);
            color_phase_random=color_phase_random/255; % Normalized for Matlab
            for current_phase=1:1:number_phase
                RGB_phase.index(current_phase).rgb = [color_phase_random(current_phase,1) color_phase_random(current_phase,2) color_phase_random(current_phase,3)];
            end            
        else
            if strcmp(color_map,'Custom')
                custom_colormap = customcmap;
            else
                custom_colormap = eval([color_map '(' num2str(number_phase) ')']);
            end
            for current_phase=1:1:number_phase
                RGB_phase.index(current_phase).rgb = [custom_colormap(current_phase,1) custom_colormap(current_phase,2) custom_colormap(current_phase,3)];
            end
        end
        % Update figure
        update_figure
    end

    function popup_axes2dview_Callback(source,~)
        % Get direction
        direction=source.Value;
        % Set slider min, max
        minor_step = 1/(Domain_size(direction)-1);
        major_step = max([minor_step 0.1]);
        set(Slider_axes_2dview,'Min',1,'Max',Domain_size(direction),'SliderStep', [minor_step, major_step],'Value',1);
        % Update text
        set(Text_slider,'String',['Slice: ' num2str(1) '/' num2str(Domain_size(direction))]);
        % Update figure
        update_figure
    end

    function update_figure
        direction = Popup_slider.Value;
        pos_ = round(Slider_axes_2dview.Value);
        Position_slice(direction)=pos_;
        if direction==1
            % Initializaion
            slice_color = zeros(Domain_size(2),Domain_size(3),3); % RGB color map
            slice_r = zeros(Domain_size(2),Domain_size(3)); % Red color map
            slice_g = zeros(Domain_size(2),Domain_size(3)); % Green color map
            slice_b = zeros(Domain_size(2),Domain_size(3)); % Blue color map
            % Attribute RGB colors for each voxel
            for current_phase=1:1:number_phase
                code_tmp =Phase_code(current_phase); % Current phase code
                slice_r(array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                slice_g(array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                slice_b(array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
            end
        elseif direction==2
            % Initializaion
            slice_color = zeros(Domain_size(1),Domain_size(3),3); % RGB color map
            slice_r = zeros(Domain_size(1),Domain_size(3)); % Red color map
            slice_g = zeros(Domain_size(1),Domain_size(3)); % Green color map
            slice_b = zeros(Domain_size(1),Domain_size(3)); % Blue color map
            % Attribute RGB colors for each voxel
            for current_phase=1:1:number_phase
                code_tmp =Phase_code(current_phase); % Current phase code
                slice_r(array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                slice_g(array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                slice_b(array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
            end
        elseif direction==3
            % Initializaion
            slice_color = zeros(Domain_size(1),Domain_size(2),3); % RGB color map
            slice_r = zeros(Domain_size(1),Domain_size(2)); % Red color map
            slice_g = zeros(Domain_size(1),Domain_size(2)); % Green color map
            slice_b = zeros(Domain_size(1),Domain_size(2)); % Blue color map
            % Attribute RGB colors for each voxel
            for current_phase=1:1:number_phase
                code_tmp =Phase_code(current_phase); % Current phase code
                slice_r(array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                slice_g(array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                slice_b(array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
            end
        end
        slice_color(:,:,1)=slice_r; % Attribute RGB color
        slice_color(:,:,2)=slice_g;
        slice_color(:,:,3)=slice_b;
        % Display the slice
        image(slice_color,'parent',axes_2dview);
        % Get axis position
        axe_position =axes_2dview.Position;
        set(axes_2dview ,'Position', axe_position);
        % Remove tick and label
        set(axes_2dview,'xtick',[],'ytick',[]);
        set(axes_2dview,'xticklabel',[],'yticklabel',[]);
        % Fit the axes box
        axis(axes_2dview,'tight');
        % Aspect ratio is 1:1
        axis(axes_2dview,'equal');
    end

% Initialize
update_figure
end

