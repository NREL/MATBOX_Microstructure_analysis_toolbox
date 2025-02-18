function [Mseg,allthresholds] = fct_localthresholding_Manual(Mgrey,p)

sz = size(Mgrey);
dimension = length(sz);

if p.excludebackground
    if strcmp(p.selectbackground,'From label')
        idbackground = find(Mgrey==p.background_label);
    else
        idbackground = find(p.background_volume==1);
    end
end

dir = str2num(p.axe(end));

% Get all thresholds per slice
xq=1:1:sz(dir);
xs = p.manualsthresh(:,1)*sz(dir);
xs(1) = 1 ;
thresholds = interp1(xs,p.manualsthresh(:,2:end),xq,'linear');
thresholds=[1e-9*ones(sz(dir),1) thresholds];

% Assign voxels to phase based on threshold
Mseg = zeros(sz);
if dir==1
    for x=1:1:sz(1)
        slice_ = squeeze(Mgrey(x,:,:));
        for k=1:1:p.numberofphase
            condition_1 = slice_ > thresholds(x,k);
            condition_2 = slice_ <= thresholds(x,k+1);
            slice_( condition_1+condition_2==2 ) = p.labels(k);
        end
        Mseg(x,:,:) = slice_;
    end    

elseif dir==2
    for y=1:1:sz(2)
        slice_ = squeeze(Mgrey(:,y,:));
        for k=1:1:p.numberofphase
            condition_1 = slice_ > thresholds(y,k);
            condition_2 = slice_ <= thresholds(y,k+1);
            slice_( condition_1+condition_2==2 ) = p.labels(k);
        end
        Mseg(:,y,:) = slice_;
    end  

elseif dir==3
    for z=1:1:sz(3)
        slice_ = Mgrey(:,:,z);
        for k=1:1:p.numberofphase
            condition_1 = slice_ > thresholds(z,k);
            condition_2 = slice_ <= thresholds(z,k+1);
            slice_( condition_1+condition_2==2 ) = p.labels(k);
        end
        Mseg(:,:,z) = slice_;
    end
end

thresholds(:,1) = [];

allthresholds.x = [1:1:sz(dir)]';
allthresholds.thresholds = thresholds;


if p.excludebackground && ~isempty(idbackground)
    if min(p.labels)>0
        Mseg(idbackground)=0;
    else
        Mseg = Mseg+1;
        Mseg(idbackground)=0;
    end
end

if p.plotprofile
    col = colororder;
    if strcmp(p.voxelunit,'um')
        p.voxelunit = '\mum';
    end
    if p.excludebackground
        str_sub = [p.sub ' (background excluded)'];
    else
        str_sub = p.sub;
    end

    Fig = figure;
    Fig.Color='white'; % Background colour
    ax = axes('Parent',Fig);
    hold(ax,'on');

    x = [1:1:sz(dir)]*p.voxelsize;

    for current_phase = 1:1:p.numberofphase
        h = plot(x,thresholds(:,current_phase),'LineWidth',2,'Color',col(current_phase,:),'DisplayName',['Phase#' num2str(current_phase,'%i') ', label=' num2str(p.labels(current_phase),'%i') ]);
        if sz(dir)==1
            set(h,'Marker','o');
        end              
    end
    xlabel(['Position along axe ' num2str(dir,'%i') ' (' p.voxelunit ')' ]);
    ylabel('Threshold (grey level)');
  
    legend(ax)

    grid(ax,'on');
    axis tight;
    set(ax,'FontName','Times new roman','FontSize',12);
    title(ax,p.filename,'FontName','Times new roman','FontSize',14);
    subtitle(ax,str_sub,'FontName','Times new roman','FontSize',12);
    hold(ax,'off');
end

% Convert to uint8 or uint16
[Mseg] = fct_intconvert(Mseg);

end