function [Fig_] = function_evolution_along_direction(Microstructure)
%function_evolution_along_direction calculates and plot evolution of array values along its axes
%[Fig_] = function_evolution_along_direction(Microstructure)


%% CALCULATION
domain_size=size(Microstructure);
dimension=length(domain_size);
% direction(*).greylevel
% :,1,: position in voxel
% :,2,: mean
% :,3,: min
% :,4,: max
% :,5,: std

% % direction 1
% Position
direction(1).greylevel=zeros(domain_size(1),5);
direction(1).greylevel(:,1)=(1:1:domain_size(1));
for position=1:1:domain_size(1)
    % Current slice of data
    slice_ = Microstructure(position,:,:);
    % Convert in an 1D array
    data_ = reshape(slice_,[numel(slice_),1]);
    % Min and max
    direction(1).greylevel(position,3) = min(data_);
    direction(1).greylevel(position,4) = max(data_);
    % Mean
    direction(1).greylevel(position,2) = mean(data_);
    % Std
    direction(1).greylevel(position,5) = std(double(data_));
end

% % direction 2
direction(2).greylevel=zeros(domain_size(2),5);
direction(2).greylevel(:,1)=(1:1:domain_size(2));

for position=1:1:domain_size(2)
    slice_ = Microstructure(:,position,:);
    data_ = reshape(slice_,[numel(slice_),1]);
    direction(2).greylevel(position,3) = min(data_);
    direction(2).greylevel(position,4) = max(data_);
    direction(2).greylevel(position,2) = mean(data_);
    direction(2).greylevel(position,5) = std(double(data_));
end

% % direction 3
if dimension==3
    direction(3).greylevel=zeros(domain_size(3),5);
    direction(3).greylevel(:,1)=(1:1:domain_size(3));
    for position=1:1:domain_size(3)
        slice_ = Microstructure(:,:,position);
        data_ = reshape(slice_,[numel(slice_),1]);
        direction(3).greylevel(position,3) = min(data_);
        direction(3).greylevel(position,4) = max(data_);
        direction(3).greylevel(position,2) = mean(data_);
        direction(3).greylevel(position,5) = std(double(data_));
    end
end

%% FIGURE
Fig_ = figure;
Fig_.Name= 'Grey level along direction';
Fig_.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig_,'position',[scrsz(1) scrsz(2) scrsz(3)*dimension/3 round(3/5*scrsz(4))]); % Full screen figure
for id_axe=1:1:dimension
    sub_axes=subplot(1,3,id_axe,'Parent',Fig_); 
    hold(sub_axes,'on');
    
    title (sprintf('Along direction %i',id_axe),'FontName','Times New Roman','FontSize',16);
    
    % Data
    x_=direction(id_axe).greylevel(:,1);
    y_mean = direction(id_axe).greylevel(:,2);
    y_min = direction(id_axe).greylevel(:,3);
    y_max = direction(id_axe).greylevel(:,4);
    y_std = direction(id_axe).greylevel(:,5);
    % Mean
    h_mean=plot(x_,y_mean); % For the legend order
    % Extremums
    h_min=plot(x_,y_min);
    h_max=plot(x_,y_max);
    % Colors, thickness, markers
    set(h_mean, 'Color', 'b','LineWidth',2,'MarkerSize',12,'Marker','none');
    set(h_min, 'Color', 'r','LineStyle','--','MarkerSize',12,'Marker','none');
    set(h_max, 'Color', 'r','LineStyle','--','MarkerSize',12,'Marker','none');
    % Mean with error bar (+- standard deviation)
    h_mean_witherrorbar = errorbar(x_,y_mean,y_std);
    set(h_mean_witherrorbar, 'Color', 'k','LineWidth',1,'MarkerSize',12,'Marker','none');
    h_mean=plot(x_,y_mean); % Plot over the other
    set(h_mean, 'Color', 'b','LineWidth',2,'MarkerSize',12,'Marker','none');
    
    legend(sub_axes,'Mean grey-level','Minimum grey-level','Maximum grey-level','Location','best');
    
    xlabel('Slice');
    ylabel('Grey level value');
    
    grid(sub_axes,'on'); % Display grid
    set(sub_axes,'XMinorGrid','on','YMinorGrid','on'); % Display grid for minor thicks also
    
    set(sub_axes,'FontName','Times New Roman','FontSize',14);
    
    hold(sub_axes,'off');
end


end

