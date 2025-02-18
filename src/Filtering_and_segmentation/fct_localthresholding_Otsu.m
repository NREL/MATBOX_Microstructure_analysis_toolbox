function [Mseg,allthresholds] = fct_localthresholding_Otsu(Mgrey,p)

sz = size(Mgrey);
dimension = length(sz);

idbackground = [];
n_background = 0;
Mgrey = double(Mgrey);
if p.excludebackground
    if strcmp(p.selectbackground,'From label')
        idbackground = find(Mgrey==p.background_label);
        Mgrey(idbackground) = NaN;
    else
        idbackground = find(p.background_volume==1);
        Mgrey(idbackground) = NaN;
    end
    n_background = length(idbackground);
end

min_ = double( min(min(min( Mgrey,[],'omitnan' ),[],'omitnan'),[],'omitnan') );
max_ = double( max(max(max( Mgrey,[],'omitnan' ),[],'omitnan'),[],'omitnan') );

if dimension==2
    sz = [sz 1];
end

if strcmp(p.axe,'Per slice along axis 1') || strcmp(p.axe,'Per slice along axis 2') || strcmp(p.axe,'Per slice along axis 3')
    dir = str2num(p.axe(end));
    thresholds=zeros(sz(dir),p.numberofphase);
    best_ratio=zeros(sz(dir),1);
    for k=1:1:sz(dir)
        if dir==1
            tmp = Mgrey(k,:,:);
        elseif dir==2
            tmp = Mgrey(:,k,:);
        elseif dir==3
            tmp = Mgrey(:,:,k);
        end
        tmp = reshape(tmp,[numel(tmp),1]);

        if p.excludebackground && n_background>0
            tmp(isnan(tmp)) = [];
        end

        if ~isempty(tmp) && length(unique(tmp))>=p.numberofphase
            [histogram_Mgrey] = function_calculate_histogram(tmp);
            % Otsu's algorithm applied to the whole volume
            Otsu_result_wholevolume=Function_otsu_algorithm(histogram_Mgrey,p.numberofphase);
            % Best threshold
            all_threshold=Otsu_result_wholevolume.allpermuation(:,2:p.numberofphase);
            all_ratio=Otsu_result_wholevolume.allratio;
            best_ratio(k) = max(all_ratio);
            idx = find(all_ratio==best_ratio(k));
            best_threshold = all_threshold(idx(1),:);
            thresholds(k,:) = [best_threshold'; max_];
        else
            tmp = linspace(min_,max_,p.numberofphase+1);
            thresholds(k,:) = tmp(2:end);
        end
    end

    thresholds_smoothed = thresholds;
    if p.smoothrange~=0 && sz(dir)>1
        for current_phase=1:1:p.numberofphase-1
            t_=thresholds(:,current_phase);
            for z=1:1:sz(dir)
                min_z=max(1,z-p.smoothrange);
                max_z=min(sz(dir),z+p.smoothrange);
                thresholds_smoothed(z,current_phase)= mean(t_(min_z:max_z));
            end
        end
    end

    thresholds = [ones(sz(dir),1)*(-1e9) thresholds];
    thresholds_smoothed = [ones(sz(dir),1)*(-1e9) thresholds_smoothed];


    Mseg = zeros(sz);
    if ~p.translateprofile
        translated_thresholds = [];
        for k=1:1:sz(dir)
            if dir==1
                tmp = Mgrey(k,:,:);
            elseif dir==2
                tmp = Mgrey(:,k,:);
            elseif dir==3
                tmp = Mgrey(:,:,k);
            end
            tmp_seg = zeros(size(tmp));
            
            % Assign voxels to phase based on threshold
            for current_phase=1:1:p.numberofphase
                if p.smoothrange~=0
                    condition_1 = tmp > thresholds_smoothed(k,current_phase);
                    condition_2 = tmp <= thresholds_smoothed(k,current_phase+1);
                else
                    condition_1 = tmp > thresholds(k,current_phase);
                    condition_2 = tmp <= thresholds(k,current_phase+1);
                end
                tmp_seg( condition_1+condition_2==2 ) = p.labels(current_phase);
            end

            if dir==1
                Mseg(k,:,:) = tmp_seg;
            elseif dir==2
                Mseg(:,k,:) = tmp_seg;
            elseif dir==3
                Mseg(:,:,k) = tmp_seg;
            end
        end        

    else
        % Threshold will be translated to match expected volume fraction
        if p.excludebackground && n_background>0
            background_vf = n_background / numel(Mgrey);
            p.expectedvf = (1-background_vf)*p.expectedvf;
        end
        initial_thresholds = thresholds_smoothed;
        translated_thresholds = thresholds_smoothed;
        number_voxel=numel(Mgrey);
        max_iteration = 50;
        min_error = 1e-3;
        r = max_-min_;        
        for current_phase=1:1:p.numberofphase-1 % No need to check last phase.
            % Dichotomy algorithm
            error_vf=1e9;
            n_iteration=1;
            threshold_translation = 0;
            while abs(error_vf)>min_error && n_iteration<max_iteration
                Mseg_tmp = zeros(sz);
                for k=1:1:sz(dir)
                    if dir==1
                        tmp = Mgrey(k,:,:);
                    elseif dir==2
                        tmp = Mgrey(:,k,:);
                    elseif dir==3
                        tmp = Mgrey(:,:,k);
                    end
                    
                    tmp2 = zeros(size(tmp));
                    condition_1 = tmp > translated_thresholds(k,current_phase);
                    condition_2 = tmp <= translated_thresholds(k,current_phase+1) + threshold_translation;
                    if p.excludebackground && n_background>0
                        condition_3 = ~isnan(tmp);
                        tmp2( condition_1+condition_2+condition_3==3 ) = 1;
                    else
                        tmp2( condition_1+condition_2==2 ) = 1;
                    end
                    
                    if dir==1
                        Mseg_tmp(k,:,:) = tmp2;
                    elseif dir==2
                        Mseg_tmp(:,k,:) = tmp2;
                    elseif dir==3
                        Mseg_tmp(:,:,k) = tmp2;
                    end
                end

                obtained_vf = sum(sum(sum( Mseg_tmp==1 ))) / number_voxel;
                error_vf = p.expectedvf(current_phase)-obtained_vf;
                if abs(error_vf)>min_error && n_iteration<max_iteration
                    if error_vf>0 % High threshold must increase
                        threshold_translation = threshold_translation + r/(2^n_iteration);
                    else % High threshold must decrease
                        threshold_translation = threshold_translation - r/(2^n_iteration);
                    end
                end
                n_iteration=n_iteration+1;
            end
            translated_thresholds(:,current_phase+1) = double(initial_thresholds(:,current_phase+1)) + threshold_translation;
        end

        for k=1:1:sz(dir)
            if dir==1
                tmp = Mgrey(k,:,:);
            elseif dir==2
                tmp = Mgrey(:,k,:);
            elseif dir==3
                tmp = Mgrey(:,:,k);
            end
            tmp_seg = zeros(size(tmp));

            % Assign voxels to phase based on threshold
            for current_phase=1:1:p.numberofphase
                condition_1 = tmp > translated_thresholds(k,current_phase);
                condition_2 = tmp <= translated_thresholds(k,current_phase+1);
                if p.excludebackground && n_background>0
                    condition_3 = ~isnan(tmp);
                    tmp_seg( condition_1+condition_2+condition_3==3 ) = p.labels(current_phase);
                else
                    tmp_seg( condition_1+condition_2==2 ) = p.labels(current_phase);
                end
            end

            if dir==1
                Mseg(k,:,:) = tmp_seg;
            elseif dir==2
                Mseg(:,k,:) = tmp_seg;
            elseif dir==3
                Mseg(:,:,k) = tmp_seg;
            end
        end
        translated_thresholds(:,1) = [];

    end

else


end

thresholds(:,1 ) =[];
thresholds_smoothed(:,1) = [];

allthresholds.x = [1:1:sz(dir)]';
allthresholds.thresholds = thresholds;
allthresholds.best_ratio = best_ratio;
allthresholds.thresholds_smoothed = thresholds_smoothed;
allthresholds.translated_thresholds = translated_thresholds;

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

    yyaxis left
    for current_phase = 1:1:p.numberofphase
        h = plot(x,thresholds(:,current_phase),'LineWidth',2,'Color',col(current_phase,:),'DisplayName',['Phase#' num2str(current_phase,'%i') ', label=' num2str(p.labels(current_phase),'%i') ]);
        if sz(dir)==1
            set(h,'Marker','o');
        end
        if p.smoothrange~=0 && sz(dir)>1
            plot(x,thresholds_smoothed(:,current_phase),'LineWidth',2,'LineStyle','--','Color',col(current_phase,:),'DisplayName','Smoothed');
        end
        if p.translateprofile
             plot(x,translated_thresholds(:,current_phase),'LineWidth',2,'LineStyle','-.','Color',col(current_phase,:),'DisplayName','Translated');
        end           
    end
    xlabel(['Position along axe ' num2str(dir,'%i') ' (' p.voxelunit ')' ]);
    ylabel('Threshold (grey level)');
  
    yyaxis right
    h = plot(x,best_ratio,'LineWidth',2,'Color','k','DisplayName','Separability criterion \eta');
    set(ax,'YColor',[0 0 0]);

    legend(ax)

    grid(ax,'on');
    axis tight;
    set(ax,'FontName','Times new roman','FontSize',12);
    title(ax,p.filename,'FontName','Times new roman','FontSize',14);
    subtitle(ax,str_sub,'FontName','Times new roman','FontSize',12);
    hold(ax,'off');
end

if p.excludebackground && ~isempty(idbackground)
    if min(p.labels)>0
        Mseg(idbackground)=0;
    else
        Mseg = Mseg+1;
        Mseg(idbackground)=0;
    end
end

% Convert to uint8 or uint16
[Mseg] = fct_intconvert(Mseg);

end