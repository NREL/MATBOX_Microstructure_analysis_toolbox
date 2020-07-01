function [array_] = Function_generate_connectedplates(parameters)

dimension = parameters.dimension;
direction = parameters.direction;
distancefromedge = parameters.distancefromedge;
centered = parameters.centered;
inplanelength = parameters.inplanelength;
inplanedistance = parameters.inplanedistance;
thicknesslayers135 = parameters.thicknesslayers135;
thicknesslayers246 = parameters.thicknesslayers246;
thicknessnolayer = parameters.thicknessnolayer;
only_alllayers = parameters.only_alllayers;
centered_alllayers = parameters.centered_alllayers;

% Set all directions
if direction==1
    through_plane = 1;
    inplane_a = 2;
    inplane_b = 3;
elseif direction ==2
    through_plane = 2;
    inplane_a = 1;
    inplane_b = 3;
else
    through_plane = 3;
    inplane_a = 1;
    inplane_b = 2;
end

% Initialization
array_ = zeros(dimension);
dimension_layer = [dimension(inplane_a) dimension(inplane_b)];
layer_1 = zeros(dimension_layer);
layer_2 = zeros(dimension_layer);
layer_3 = zeros(dimension_layer);
layer_5a = zeros(dimension_layer);
layer_5b = zeros(dimension_layer);

% Layer  1
x0 = distancefromedge+1; x1 = x0+inplanelength-1; % Initialization x
while x1<=dimension(inplane_a)-distancefromedge
    y0 = distancefromedge+1; y1 = y0+inplanelength-1;
    while y1<=dimension(inplane_b)-distancefromedge
        layer_1(x0:x1,y0:y1) = 1;
        % Next iteration
        y0 = y1 + inplanedistance + 1; y1 = y0 + inplanelength-1;        
    end
    % Next iteration
    x0 = x1 + inplanedistance + 1; x1 = x0 + inplanelength-1;
end
if centered
    index_ = find(layer_1==1);
    [I,J] = ind2sub(dimension_layer,index_);
    x_min = min(I); x_max = max(I);
    y_min = min(J); y_max = max(J);
    dx_min = x_min-1; dx_max = dimension_layer(1)-x_max;
    dy_min = y_min-1; dy_max = dimension_layer(2)-y_max;
    if dx_min~=dx_max || dy_min~=dy_max
        layer_1_tmp = zeros(dimension_layer); % Initialize
        %advance_x = round((dx_min+dx_max)/2);
        %advance_y = round((dy_min+dy_max)/2);
        if dx_min<dx_max
            %x0 = x_min + advance_x;
            %x1 = x_max + advance_x;  
            x0 = x_min + round((dx_max-dx_min)/2);
            x1 = x_max + round((dx_max-dx_min)/2);            
        elseif dx_min>dx_max
            %x0 = x_min - advance_x;
            %x1 = x_max - advance_x;
            x0 = x_min - round((dx_min-dx_max)/2);
            x1 = x_max - round((dx_min-dx_max)/2);
        elseif dx_min==dx_max
            x0 = x_min;
            x1 = x_max;
        end
        if dy_min<dy_max
            %y0 = y_min + advance_y;
            %y1 = y_max + advance_y;
            y0 = y_min + round((dy_max-dy_min)/2);
            y1 = y_max + round((dy_max-dy_min)/2);             
        elseif dy_min>dy_max
            %y0 = y_min - advance_y;
            %y1 = y_max - advance_y;  
            y0 = y_min - round((dy_min-dy_max)/2);
            y1 = y_max - round((dy_min-dy_max)/2);            
        elseif dy_min==dy_max
            y0 = y_min;
            y1 = y_max;
        end
        values = layer_1(x_min:x_max, y_min:y_max);
        layer_1_tmp(x0:x1,y0:y1) = values; % Assign
        layer_1 = layer_1_tmp; % Replace
    end
end

% Layer 3
index_ = find(layer_1==1);
[I,J] = ind2sub(dimension_layer,index_);
x_min = min(I); x_max = max(I);
y_min = min(J); y_max = max(J);
x0 = min(x_min+inplanelength, round(x_min + 2/3*inplanelength));
y0 = min(y_min+inplanelength, round(y_min + 2/3*inplanelength));
x1 = max(x_max-inplanelength, round(x_max - 2/3*inplanelength));
y1 = max(y_max-inplanelength, round(y_max - 2/3*inplanelength));
dx = x1-x0;
dy = y1-y0;
values = layer_1(x_min:x_min+dx, y_min:y_min+dy);
layer_3(x0:x1,y0:y1) = values; % Assign

% Layer 5a
index_ = find(layer_1==1);
[I,J] = ind2sub(dimension_layer,index_);
x_min = min(I); x_max = max(I);
y_min = min(J); y_max = max(J);
x0 = min(x_min+inplanelength, round(x_min + 2/3*inplanelength));
y0 = y_min;
x1 = max(x_max-inplanelength, round(x_max - 2/3*inplanelength));
y1 = y_max;
dx = x1-x0;
dy = y1-y0;
values = layer_1(x_min:x_min+dx, y_min:y_min+dy);
layer_5a(x0:x1,y0:y1) = values; % Assign

% Layer 5b
index_ = find(layer_1==1);
[I,J] = ind2sub(dimension_layer,index_);
x_min = min(I); x_max = max(I);
y_min = min(J); y_max = max(J);
x0 = x_min;
y0 = min(y_min+inplanelength, round(y_min + 2/3*inplanelength));
x1 = x_max;
y1 = max(y_max-inplanelength, round(y_max - 2/3*inplanelength));
dx = x1-x0;
dy = y1-y0;
values = layer_1(x_min:x_min+dx, y_min:y_min+dy);
layer_5b(x0:x1,y0:y1) = values; % Assign
layer_5 = layer_5a + layer_5b;
layer_5(layer_5~=0)=1;

% Layer 2, 4, 6
index_ = layer_1+layer_3==2;
layer_2(index_)=1;
layer_4 = layer_2;
layer_6 = layer_2;

% 3D volume
current_layer_id=1;
current_thickness=thicknesslayers135;
n_thickness=0;
for k=1:1:dimension(through_plane)
    if k<=thicknessnolayer || k>=dimension(through_plane)-thicknessnolayer+1
        foo=1; % Empty layer
    else
        n_thickness=n_thickness+1;
        if n_thickness>current_thickness
            current_layer_id=current_layer_id+1; % Next layer
            n_thickness = 1;
            if current_layer_id==7
                current_layer_id=1; % Start over form first layer
            end
        end
        if current_layer_id==1
            current_layer = layer_1;
            current_thickness = thicknesslayers135;
            if n_thickness==current_thickness
                current_layer=current_layer*2; % Mark end of layers pattern
            end
        elseif current_layer_id==2
            current_layer = layer_2;
            current_thickness = thicknesslayers246;
        elseif current_layer_id==3
            current_layer = layer_3;
            current_thickness = thicknesslayers135;
        elseif current_layer_id==4
            current_layer = layer_4;
            current_thickness = thicknesslayers246;
         elseif current_layer_id==5
            current_layer = layer_5;
            current_thickness = thicknesslayers135;           
         elseif current_layer_id==6
            current_layer = layer_6;
            current_thickness = thicknesslayers246;           
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

