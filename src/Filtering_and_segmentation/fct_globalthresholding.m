function [Mseg,newtype,thresholds] = fct_globalthresholding(Mgrey,p)
sz = size(Mgrey);
dimension = length(sz);
newtype = 'Segmented (phase)';

unis = unique(Mgrey);
n_greylevel = length(unis);

% May have too modifed format if parameters were loaded from a modified macro
if isfield(p,'labels') && ischar(p.labels)
    p.labels = str2num(p.labels);
end
if isfield(p,'thresholds') && ischar(p.thresholds)
    p.thresholds = str2num(p.thresholds);
end
p.thresholds = reshape(p.thresholds,[length(p.thresholds),1]);
p.labels = reshape(p.labels,[length(p.labels),1]);
p.expectedvf = reshape(p.expectedvf,[length(p.expectedvf),1]);

if n_greylevel >= p.numberofphase

    if p.excludebackground
        if strcmp(p.selectbackground,'From label')
            idbackground = find(Mgrey==p.background_label);
        else
            idbackground = find(p.background_volume==1);
        end
    end

    if strcmp(p.method,'Manual')
        thresholds = p.thresholds;

    elseif strcmp(p.method,'Linear')
        if p.excludebackground && ~isempty(idbackground)
            tmp = Mgrey;
            tmp(idbackground)=[];
            max_=max(tmp);
            min_=min(tmp);
        else
            max_=max(max(max(Mgrey)));
            min_=min(min(min(Mgrey)));
        end
        initial_delta=max_-min_;
        final_delta=p.numberofphase-1;
        Mseg =  round( ((double((Mgrey-min_)) ./ double(initial_delta)) .* final_delta)+1 );
        thresholds = p.thresholds; % Not used

        unis = unique(Mseg);
        if sum(unis==p.labels)~=length(unis)
            tmp = Mseg;
            for k=1:1:length(unis)
                Mseg(tmp==unis(k)) = p.labels(k);
            end
            clear tmp;
        end

    elseif strcmp(p.method,'Fit to reach volume fractions')
        if p.excludebackground && ~isempty(idbackground)
            tmp = Mgrey;
            tmp(idbackground)=[];
            histogram_Mgrey = function_calculate_histogram(tmp);
            n = numel(tmp);
        else
            histogram_Mgrey = function_calculate_histogram(Mgrey);
            n = numel(Mgrey);
        end

        histogram_Mgrey(:,3) = histogram_Mgrey(:,2)/n;
        % Cumulative function
        n=length(histogram_Mgrey(:,3));
        cumulative_fct = zeros(n,1);
        cumulative_fct(end,1)=histogram_Mgrey(end,3);
        for current_value=n-1:-1:1
            cumulative_fct(current_value,1)=cumulative_fct(current_value+1,1)+histogram_Mgrey(current_value,3);
        end
        % Determine thresholds
        target = 1-p.expectedvf(1);
        thresholds = zeros(p.numberofphase,1);
        for k=1:1:p.numberofphase-1
            idx = find( abs(cumulative_fct-target) == min(abs(cumulative_fct-target)) );
            thresholds(k)=histogram_Mgrey(idx(1)-1,1);
            target = target-p.expectedvf(k+1);
        end
        thresholds(end) = max(max(max( Mgrey )));

    elseif strcmp(p.method,'Otsu')
        if p.excludebackground && ~isempty(idbackground)
            tmp = Mgrey;
            tmp(idbackground)=[];
            m = max(tmp);
            [histogram_Mgrey] = function_calculate_histogram(tmp);
        else
            m = max(max(max( Mgrey)));
            [histogram_Mgrey] = function_calculate_histogram(Mgrey);
        end
        % Otsu's algorithm applied to the whole volume
        if n_greylevel>p.numberofphase
            Otsu_result_wholevolume=Function_otsu_algorithm(histogram_Mgrey,p.numberofphase);
            % Best threshold
            all_threshold_wholevolume=Otsu_result_wholevolume.allpermuation(:,2:p.numberofphase);
            all_ratio_wholevolume=Otsu_result_wholevolume.allratio;
            best_ratio_wholevolume = max(all_ratio_wholevolume);
            idx = find(all_ratio_wholevolume==best_ratio_wholevolume);
            best_threshold_wholevolume = all_threshold_wholevolume(idx(1),:);
            thresholds = [best_threshold_wholevolume'; m];
        else
            thresholds = unis;
        end

    elseif strcmp(p.method,'Range-constrained Otsu')
        % Characteristic analysis of Otsu threshold and its applications
        % X. Xu S. Xu, L. Jin, E. Song, Pattern Recognition Letters, 32, 956-961, 2011
        if p.excludebackground && ~isempty(idbackground)
            Mgrey(idbackground)=[];
        end
        min_ = min(min(min( Mgrey )));
        max_ = max(max(max( Mgrey )));

        % Baseline Otsu
        [histogram_Mgrey] = function_calculate_histogram(Mgrey);
        Otsu_result_wholevolume=Function_otsu_algorithm(histogram_Mgrey,p.numberofphase);
        all_threshold_wholevolume=Otsu_result_wholevolume.allpermuation(:,2:p.numberofphase);
        all_ratio_wholevolume=Otsu_result_wholevolume.allratio;
        best_ratio_wholevolume = max(all_ratio_wholevolume);
        idx = find(all_ratio_wholevolume==best_ratio_wholevolume);
        best_threshold_wholevolume = all_threshold_wholevolume(idx(1),:);
        thresholds = [best_threshold_wholevolume'; max_];

        for t=1:p.numberofphase-1
            vars = Otsu_result_wholevolume.sigmaw2(t:t+1);
            if vars(1)<vars(2)
                bounds = Otsu_result_wholevolume.thresholdbounds(t,:);
            else
                bounds = Otsu_result_wholevolume.thresholdbounds(t+1,:);
            end

            tmp = Mgrey;
            tmp(tmp<bounds(1))=[];
            tmp(tmp>bounds(2))=[];

            [histogram_tmp] = function_calculate_histogram(tmp);
            Otsu_result_tmp=Function_otsu_algorithm(histogram_tmp,2);

            thresholds(t) = Otsu_result_tmp.threshold(2);
        end
    end
 
    if ~strcmp(p.method,'Linear')
        % Assign voxels to phase based on threshold
        Mseg = zeros(sz);
        thresholds=[-1e9;thresholds];
        for k=1:1:p.numberofphase
            condition_1 = Mgrey > thresholds(k);
            condition_2 = Mgrey <= thresholds(k+1);
            Mseg( condition_1+condition_2==2 ) = p.labels(k);
        end
        thresholds(1)=[];
    end

    if p.excludebackground && ~isempty(idbackground)
        if min(p.labels)>0
            Mseg(idbackground)=0;
        else
            Mseg = Mseg+1;
            Mseg(idbackground)=0;
        end
    end

else
    Mseg = Mgrey;
    warndlg('No segmentation performed: number of phases > number of grey level','Segmentation warning')
end

% Convert to uint8 or uint16
[Mseg] = fct_intconvert(Mseg);

end