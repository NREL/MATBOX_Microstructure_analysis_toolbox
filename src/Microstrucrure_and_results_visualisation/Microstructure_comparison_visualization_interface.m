function [] = Microstructure_comparison_visualization_interface(varargin)
% Microstructure_comparison_visualization_interface(a,b,c,...)
% or
% Microstructure_comparison_visualization_interface(a,b,c,...,{'Axe 1 name','Axe 2 name','Axe 3 name'})
% Wtih a,b,c,... 3D arrays

nargin; % Number of input variable when the function is call
tmp = varargin(nargin); last_variable=tmp{1};
if iscell(last_variable)
    number_volume=nargin-1;
    direction_name = varargin(nargin); direction_name=direction_name{1};
else
    number_volume=nargin;
    direction_name = {'Normal to axe 1','Normal to axe 2','Normal to axe 3'}; % Defaut string
end

%% PARAMETERS
font_name_GUI = 'Times new roman';
font_size_GUI = 12;

% Default colormap
color_phase_default = round(colororder*255);
color_phase_default(2,:) = [127 127 127];
tmp=randi(255,1e6,3);
color_phase_default=[color_phase_default;tmp];
color_phase_default=color_phase_default/255; % Normalized for Matlab
choice_colormap = 'MATLAB default';

% Spacing
dw = 0.1; % Width spacing as fraction of axe width
dh = 0.1; % Height spacing as fraction of axe height
sh = 0.125; % Bottom space reserved for interface

%% VOLUME AND PHASE INFORMATION
Phase_code = [];
max_domainsize = [0 0 0];
for k=1:1:number_volume
    tmp = varargin(k); Arrays(k).array = tmp{1}; % More convenient that using varargin
    unis = unique(Arrays(k).array);
    Phase_code=unique( [unique(Phase_code); unis] ); % Get phase code
    if min(Phase_code)>=0 && max(Phase_code)<255
        if sum(Phase_code==round(Phase_code))==length(Phase_code)
            Arrays(k).array = uint8(Arrays(k).array);
        end
    end
    Arrays(k).domainsize = size(Arrays(k).array);
    max_domainsize(1) = max([max_domainsize(1) Arrays(k).domainsize(1)]);
    max_domainsize(2) = max([max_domainsize(2) Arrays(k).domainsize(2)]);
    if length(Arrays(k).domainsize)==3
        max_domainsize(3) = max([max_domainsize(3) Arrays(k).domainsize(3)]);
    else
        max_domainsize(3) = max([max_domainsize(3) 1]);
    end

    % Choose of colormap
    if length(unis)>10
        data_type(k) = {'Grey level'};
    else
        data_type(k) = {'Segmented'};
    end


end

clear varargin
 
number_phase=length(Phase_code); % Get number of phase
% Default Colors
for current_phase=1:1:number_phase
    RGB_phase.index(current_phase).rgb = [color_phase_default(current_phase,1) color_phase_default(current_phase,2) color_phase_default(current_phase,3)];
end

%% FIGURE
Fig = figure; % Create figure
Fig.Name= 'Compare microstructure';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) scrsz(3)/2 scrsz(4)/2]); % Full screen figure
% Figure for volume view
n_row = floor(sqrt(number_volume));
if floor(sqrt(number_volume)) < sqrt(number_volume)
    n_col = ceil(number_volume/n_row);
else
    n_col = n_row;
end

% Control axes location
w = 1 / (n_col+(n_col+1)*dw);
dw=dw*w;
h = (1-sh) / (n_row+n_row*dh);
dh=dh*h;

idx=0;
for k_row=1:1:n_row
    k=n_row-k_row+1;
    y=((k-1)*dh + (k-1)*h)+sh;
    for k_col=1:1:n_col
        x=k_col*dw + (k_col-1)*w;
        idx=idx+1;
        axes_2dview(idx).ax_ = axes('Parent', Fig,'Visible','off','FontName',font_name_GUI,'Units','normalized','Position', [x y w h]);
        % Remove tick and label
        set(axes_2dview(idx).ax_,'xtick',[],'ytick',[]);
        set(axes_2dview(idx).ax_,'xticklabel',[],'yticklabel',[]);
        % Box on
        box(axes_2dview(idx).ax_,'on');
        % Fit the axes box
        axis(axes_2dview(idx).ax_,'tight');
        % Aspect ratio is 1:1
        axis(axes_2dview(idx).ax_,'equal');
    end
end

% Slider
Slider_axes_2dview = uicontrol('Parent', Fig,'Style', 'slider','Min',1,'Max',100,'Value',1,'Units','normalized','Position', [0.3 0.06 0.4 0.04],'Callback', @slider_axes2dview_Callback,...
    'Visible','on','enable','on');
% Text (slider position)
Text_slider = uicontrol('Parent', Fig,'Style', 'text','FontSize',font_size_GUI,'FontName',font_name_GUI,'String','Position: -/-','Visible','on','enable','on',...
    'BackgroundColor','w','HorizontalAlignment','left','Units','normalized','Position', [0.3 0.0 0.4 0.04]);
% Popup menu (slider direction)
Popup_slider = uicontrol('Parent', Fig,'Style', 'popup','FontSize',font_size_GUI,'FontName',font_name_GUI,...
    'String', direction_name,'Value',3,'Units','normalized','Position', [0.06 0.05 0.215 0.04],'enable','on','Visible','on','Callback', @popup_axes2dview_Callback);

direction = 3; pos_=1;
set(Text_slider,'String',['Slice: ' num2str(pos_) '/' num2str(max_domainsize(direction))]);
minor_step = 1/(max_domainsize(direction)-1);
major_step = max([minor_step 0.1]);
if max_domainsize(direction)>1
    set(Slider_axes_2dview,'Min',1,'Max',max_domainsize(direction),'SliderStep', [minor_step, major_step],'Value',1);
else
    set(Slider_axes_2dview,'Enable','off');
end

% Colormap
Popup_colormap_seg = uicontrol('Parent', Fig,'Style', 'popup','FontSize',font_size_GUI,'FontName',font_name_GUI,...
    'String', {'MATLAB default','gray','bone','copper','jet','turbo','parula','Random'},'Value',1,'Units','normalized','Position', [0.725 0.06 0.215 0.04],'enable','on','Visible','on','Tooltip','Segmened colormap','Callback', @popup_colormap_seg_Callback);
%Popup_colormap_grey = uicontrol('Parent', Fig,'Style', 'popup','FontSize',font_size_GUI,'FontName',font_name_GUI,...
%    'String', {'gray','bone','copper','jet','turbo','parula'},'Value',1,'Units','normalized','Position', [0.8325 0.06 0.1075 0.04],'enable','on','Visible','on','Tooltip','Grey level colormap','Callback', @popup_colormap_seg_Callback);


% Pixel grid checkbox
Checkbox_pixelgrid = uicontrol('Parent', Fig,'Style', 'checkbox','FontSize',font_size_GUI,'FontName',font_name_GUI,...
    'String','Pixel grid' ,'Value',0,'Units','normalized','Position', [0.725 0.01 0.215 0.04],'enable','on','Visible','on','Callback', @checkbox_pixelgrid_Callback);

% Checkbox
    function checkbox_pixelgrid_Callback(source,~)
        % Get position value
        if source.Value
            pixelgrid
        end    
        update_figure
    end

% Slider
    function slider_axes2dview_Callback(source,~)
        % Get position value
        pos_=round(source.Value);
        direction = Popup_slider.Value;
        % Update text
        % set(Text_slider,'String',['Slice: ' num2str(pos_) '/' num2str(max_domainsize(direction))]);
        % Update position array
        % Position_slice(direction)=pos_;
        % Update figure
        update_figure
    end

% Select direction
    function popup_colormap_seg_Callback(source,~)
        % Get colormap
        color_map = char(Popup_colormap_seg.String(source.Value));
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
            custom_colormap = eval([color_map '(' num2str(number_phase) ')']);
            for current_phase=1:1:number_phase
                
                RGB_phase.index(current_phase).rgb = [custom_colormap(current_phase,1) custom_colormap(current_phase,2) custom_colormap(current_phase,3)];
            end
        end

        %RGB_phase.index(1).rgb = [1 1 1];

        % Update figure
        update_figure
    end

    function popup_axes2dview_Callback(source,~)
        % Get direction
        direction=source.Value;
        % Set slider min, max
        minor_step = 1/(max_domainsize(direction)-1);
        major_step = max([minor_step 0.1]);
        set(Slider_axes_2dview,'Min',1,'Max',max_domainsize(direction),'SliderStep', [minor_step, major_step],'Value',1);
        % Update text
        set(Text_slider,'String',['Slice: ' num2str(1) '/' num2str(max_domainsize(direction))]);
        % Update figure
        update_figure
    end

    function update_figure
        direction = Popup_slider.Value;
        pos_ = round(Slider_axes_2dview.Value);
        idx = 0;
        str_slice = 'Slice:';
        for k_row=1:1:n_row
            for k_col=1:1:n_col
                idx=idx+1;
                if idx<=number_volume
                    Domain_size = Arrays(idx).domainsize;
                    if length(Domain_size)==2
                        Domain_size = [Domain_size 1];
                    end
                    a = (Domain_size(direction) - 1) / (max_domainsize(direction) - 1);
                    b = 1-a*1;
                    pos__ = round(a*pos_+b);
                    Position_slice(direction)=pos__;
                    if Domain_size(3)==1
                        Position_slice(3)=1;
                    end
                    
                    str_slice = [str_slice '  ' num2str(pos__) '/' num2str(Domain_size(direction))];
                    
                    if direction==1
                        % Initializaion
                        slice_color = zeros(Domain_size(2),Domain_size(3),3); % RGB color map
                        slice_r = zeros(Domain_size(2),Domain_size(3)); % Red color map
                        slice_g = zeros(Domain_size(2),Domain_size(3)); % Green color map
                        slice_b = zeros(Domain_size(2),Domain_size(3)); % Blue color map
                        % Attribute RGB colors for each voxel
                        for current_phase=1:1:number_phase
                            code_tmp =Phase_code(current_phase); % Current phase code
                            slice_r(Arrays(idx).array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                            slice_g(Arrays(idx).array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                            slice_b(Arrays(idx).array(Position_slice(1),:,:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
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
                            slice_r(Arrays(idx).array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                            slice_g(Arrays(idx).array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                            slice_b(Arrays(idx).array(:,Position_slice(2),:)==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
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
                            slice_r(Arrays(idx).array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(1); % Attribute red color
                            slice_g(Arrays(idx).array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(2); % Attribute green color
                            slice_b(Arrays(idx).array(:,:,Position_slice(3))==code_tmp)= RGB_phase.index(current_phase).rgb(3); % Attribute blue color
                        end
                    end
                    slice_color(:,:,1)=slice_r; % Attribute RGB color
                    slice_color(:,:,2)=slice_g;
                    slice_color(:,:,3)=slice_b;
                    % Display the slice
                    image(slice_color,'parent',axes_2dview(idx).ax_);
                    % Get axis position
                    axe_position =axes_2dview(idx).ax_.Position;
                    set(axes_2dview(idx).ax_ ,'Position', axe_position);
                    % Remove tick and label
                    set(axes_2dview(idx).ax_,'xtick',[],'ytick',[]);
                    set(axes_2dview(idx).ax_,'xticklabel',[],'yticklabel',[]);
                    % Fit the axes box
                    axis(axes_2dview(idx).ax_,'tight');
                    % Aspect ratio is 1:1
                    axis(axes_2dview(idx).ax_,'equal');
                    if Checkbox_pixelgrid.Value
                        pixelgrid(axes_2dview(idx).ax_)
                    end
                end
            end
        end
        set(Text_slider,'String',str_slice);

    end

% Initialize
update_figure


end

