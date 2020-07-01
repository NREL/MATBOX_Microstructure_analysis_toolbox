function [array_] = Function_generate_strips(parameters)

dimension = parameters.dimension;
direction = parameters.direction;
alignement = parameters.alignement;
distancefromedge = parameters.distancefromedge;
centered = parameters.centered;
stripwidth = parameters.stripwidth;
striplength = parameters.striplength;
inplanedistance = parameters.inplanedistance;
thicknesslayers1357 = parameters.thicknesslayers1357;
thicknesslayers2468 = parameters.thicknesslayers2468;
thicknessnolayer = parameters.thicknessnolayer;
only_alllayers = parameters.only_alllayers;
centered_alllayers = parameters.centered_alllayers;

% Set all directions
if direction==1
    through_plane = 1;
    inplane_a = 2;
    inplane_b = 3;
    if alignement==2
        striplength = round(striplength*dimension(inplane_a));
        axis_x = inplane_b;
        axis_y = inplane_a;
    elseif alignement==3
        striplength = round(striplength*dimension(inplane_b));
        axis_x = inplane_a;
        axis_y = inplane_b;        
    end
    
elseif direction==2
    through_plane = 2;
    inplane_a = 1;
    inplane_b = 3;
    if alignement==1
        striplength = round(striplength*dimension(inplane_a));
        axis_x = inplane_b;
        axis_y = inplane_a;  
    elseif alignement==3
        striplength = round(striplength*dimension(inplane_b));
        axis_x = inplane_a;
        axis_y = inplane_b;           
    end 
    
else
    through_plane = 3;
    inplane_a = 1;
    inplane_b = 2;
    if alignement==1
        striplength = round(striplength*dimension(inplane_a));
        axis_x = inplane_b;
        axis_y = inplane_a;         
    elseif alignement==2
        striplength = round(striplength*dimension(inplane_b));
        axis_x = inplane_a;
        axis_y = inplane_b;          
    end 
end

% Initialization
array_ = zeros(dimension);
dimension_layer = [dimension(inplane_a) dimension(inplane_b)];
layer_1 = zeros(dimension_layer);
layer_2468 = zeros(dimension_layer);
layer_3 = zeros(dimension_layer);
layer_5 = zeros(dimension_layer);
layer_7 = zeros(dimension_layer);

% Layer 1 (from bottom)
x0 = distancefromedge+1; x1 = x0+stripwidth-1; % Initialization x
while x1<=dimension(axis_x)-distancefromedge
    y0 = 1; y1 = y0+striplength-1;
    if axis_x==inplane_a
        layer_1(x0:x1,y0:y1) = 1;
    else
        layer_1(y0:y1,x0:x1) = 1;
    end
    % Next iteration
    x0 = x1 + inplanedistance + 1; x1 = x0 + stripwidth-1;
end
if centered
    index_ = find(layer_1==1);
    [I,J] = ind2sub(dimension_layer,index_);
    x_min = min(I); x_max = max(I);
    y_min = min(J); y_max = max(J);
    dx_min = x_min-1; dx_max = dimension_layer(1)-x_max;
    dy_min = y_min-1; dy_max = dimension_layer(2)-y_max;
    if dx_min~=dx_max && axis_x==inplane_a
        layer_1_tmp = zeros(dimension_layer); % Initialize
        if dx_min<dx_max
            x0 = x_min + round((dx_max-dx_min)/2);
            x1 = x_max + round((dx_max-dx_min)/2);            
        elseif dx_min>dx_max
            x0 = x_min - round((dx_min-dx_max)/2);
            x1 = x_max - round((dx_min-dx_max)/2);
        end
        values = layer_1(x_min:x_max, :);
        layer_1_tmp(x0:x1,:) = values; % Assign
        layer_1 = layer_1_tmp; % Replace
    elseif dy_min~=dy_max && axis_x==inplane_b
        layer_1_tmp = zeros(dimension_layer); % Initialize
        if dy_min<dy_max
            y0 = y_min + round((dy_max-dy_min)/2);
            y1 = y_max + round((dy_max-dy_min)/2);
        elseif dy_min>dy_max
            y0 = y_min - round((dy_min-dy_max)/2);
            y1 = y_max - round((dy_min-dy_max)/2);
        end
        values = layer_1(:,y_min:y_max);
        layer_1_tmp(:,y0:y1) = values; % Assign
        layer_1 = layer_1_tmp; % Replace
    end
end

% Layer 5 (from bottom)
index_ = find(layer_1==1);
[I,J] = ind2sub(dimension_layer,index_);
x_min = min(I); x_max = max(I);
y_min = min(J); y_max = max(J);
if axis_x==inplane_a
    x0 = round(x_min + 2/3*stripwidth);
    y0 = y_min;
    x1 = round(x_max - 2/3*stripwidth);
    y1 = y_max;
elseif axis_x==inplane_b
    x0 = x_min;
    y0 = round(y_min + 2/3*stripwidth);
    x1 = x_max;
    y1 = round(y_max - 2/3*stripwidth);
end
dx = x1-x0;
dy = y1-y0;
values = layer_1(x_min:x_min+dx, y_min:y_min+dy);
layer_5(x0:x1,y0:y1) = values; % Assign

% Layer 3 (from top)
if axis_x==inplane_a
    x0 = round(x_min+1/3*stripwidth); x1 = x0+stripwidth-1; % Initialization x
    while x1<=dimension(axis_x)-distancefromedge
        y0 = dimension(axis_y) - striplength + 1; y1 = dimension(axis_y);
        layer_3(x0:x1,y0:y1) = 1;
        % Next iteration
        x0 = x1 + inplanedistance + 1; x1 = x0 + stripwidth-1;
    end
elseif axis_x==inplane_b
    x0 = round(y_min+1/3*stripwidth); x1 = x0+stripwidth-1; % Initialization x
    while x1<=dimension(axis_x)-distancefromedge
        y0 = dimension(axis_y) - striplength + 1; y1 = dimension(axis_y);
        layer_3(y0:y1,x0:x1) = 1;
        % Next iteration
        x0 = x1 + inplanedistance + 1; x1 = x0 + stripwidth-1;
    end    
end

% Layer 7 (from top)
if axis_x==inplane_a
    x0 = round(x_min+stripwidth); x1 = x0+stripwidth-1; % Initialization x
    while x1<=dimension(axis_x)-distancefromedge
        y0 = dimension(axis_y) - striplength + 1; y1 = dimension(axis_y);
        layer_7(x0:x1,y0:y1) = 1;
        % Next iteration
        x0 = x1 + inplanedistance + 1; x1 = x0 + stripwidth-1;
    end
elseif axis_x==inplane_b
    x0 = round(y_min+stripwidth); x1 = x0+stripwidth-1; % Initialization x
    while x1<=dimension(axis_x)-distancefromedge
        y0 = dimension(axis_y) - striplength + 1; y1 = dimension(axis_y);
        layer_7(y0:y1,x0:x1) = 1;
        % Next iteration
        x0 = x1 + inplanedistance + 1; x1 = x0 + stripwidth-1;
    end    
end



% 3D volume
current_layer_id=1;
current_thickness=thicknesslayers1357;
n_thickness=0;
for k=1:1:dimension(through_plane)
    if k<=thicknessnolayer || k>=dimension(through_plane)-thicknessnolayer+1
        foo=1; % Empty layer
    else
        n_thickness=n_thickness+1;
        if n_thickness>current_thickness
            current_layer_id=current_layer_id+1; % Next layer
            n_thickness = 1;
            if current_layer_id==9
                current_layer_id=1; % Start over form first layer
            end
        end
        if current_layer_id==1
            current_layer = layer_1;
            current_thickness = thicknesslayers1357;
            if n_thickness==current_thickness
                current_layer=current_layer*2; % Mark end of layers pattern
            end
        elseif current_layer_id==2
            current_layer = layer_2468;
            current_thickness = thicknesslayers2468;
        elseif current_layer_id==3
            current_layer = layer_3;
            current_thickness = thicknesslayers1357;
        elseif current_layer_id==4
            current_layer = layer_2468;
            current_thickness = thicknesslayers2468;
        elseif current_layer_id==5
            current_layer = layer_5;
            current_thickness = thicknesslayers1357;
        elseif current_layer_id==6
            current_layer = layer_2468;
            current_thickness = thicknesslayers2468;
         elseif current_layer_id==7
            current_layer = layer_7;
            current_thickness = thicknesslayers1357;           
          elseif current_layer_id==8
            current_layer = layer_2468;
            current_thickness = thicknesslayers2468;          
        end
        if direction==1
            array_(k,:,:) = current_layer;
        elseif direction==2
            array_(:,k,:) = current_layer;
        elseif direction==3
            array_(:,:,k) = current_layer;
        end
    end
end

if only_alllayers % Keep only pattern layers 1,2,3,4,5,6,7,8,1 
    index_ = find(array_==2); % End of layer 1
    [I,J,K] = ind2sub(dimension,index_);
    x_max = max(I); y_max = max(J); z_max = max(K);
    if direction==1
        if x_max<dimension(through_plane)
            array_(x_max+1:dimension(through_plane),:,:)=0;
        end
    elseif direction==2
        if y_max<dimension(through_plane)
            array_(:,y_max+1:dimension(through_plane),:)=0;
        end
    elseif direction==3
        if z_max<dimension(through_plane)
            array_(:,:,z_max+1:dimension(through_plane))=0;
        end
    end
    array_(array_~=0)=1; % Re-assign
    
    if centered_alllayers
        index_ = find(array_==1);
        [I,J,K] = ind2sub(dimension,index_);
        x_min = min(I); y_min = min(J); z_min = min(K);
        x_max = max(I); y_max = max(J); z_max = max(K);
        if direction==1
            dx_min = x_min-1; dx_max = dimension(through_plane)-x_max;
            if dx_min~=dx_max
                array_tmp = zeros(dimension); % Initialize
                advance_x = round((dx_min+dx_max)/2);
                if dx_min<dx_max
                    x0 = x_min + round((dx_max-dx_min)/2);
                    x1 = x_max + round((dx_max-dx_min)/2);
                elseif dx_min>dx_max
                    x0 = x_min - round((dx_min-dx_max)/2);
                    x1 = x_max - round((dx_min-dx_max)/2);
                end
                values = array_(x_min:x_max, :,:);
                array_tmp(x0:x1,:,:) = values; % Assign
                array_ = array_tmp; % Replace
            end
        elseif direction==2
            dy_min = y_min-1; dy_max = dimension(through_plane)-y_max;
            if dy_min~=dy_max
                array_tmp = zeros(dimension); % Initialize
                advance_y = round((dy_min+dy_max)/2);
                if dy_min<dy_max
                    y0 = y_min + round((dy_max-dy_min)/2);
                    y1 = y_max + round((dy_max-dy_min)/2);
                elseif dy_min>dy_max
                    y0 = y_min - round((dy_min-dy_max)/2);
                    y1 = y_max - round((dy_min-dy_max)/2);
                end
                values = array_(:,y_min:y_max,:);
                array_tmp(:,y0:y1,:) = values; % Assign
                array_ = array_tmp; % Replace
            end            
        elseif direction==3
            dz_min = z_min-1; dz_max = dimension(through_plane)-z_max;
            if dz_min~=dz_max
                array_tmp = zeros(dimension); % Initialize
                advance_z = round((dz_min+dz_max)/2);
                if dz_min<dz_max
                    z0 = z_min + round((dz_max-dz_min)/2);
                    z1 = z_max + round((dz_max-dz_min)/2);
                elseif dz_min>dz_max
                    z0 = z_min - round((dz_min-dz_max)/2);
                    z1 = z_max - round((dz_min-dz_max)/2);
                end
                values = array_(:,:,z_min:z_max);
                array_tmp(:,:,z0:z1) = values; % Assign
                array_ = array_tmp; % Replace
            end      
        end
    end
    
end
array_(array_~=0)=1; % Re-assign

end

