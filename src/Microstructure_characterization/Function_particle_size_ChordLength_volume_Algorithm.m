function [fitted_diameter,Rs,difference,CDF,Length_along1,Length_along2,Length_along3] = Function_particle_size_ChordLength_volume_Algorithm(BW,size_filter)

% Parameter (analytical mesh)
res = 400;
nR = 1000;

sz = size(BW);

% Direction 1
Length_along1 = zeros(sz);
%L1 = [];
for ky=1:1:sz(2)
    for kz=1:1:sz(3)
        OneDarray = reshape(BW(:,ky,kz),[sz(1),1]);
        %[Length_along1] = get_lengtharray(OneDarray,Length_along1,1,ky,kz); Using function is too slow...
        % Identify starting and ending points
        d = diff(OneDarray);
        idstart = find(d==1);
        idend = find(d==-1);
        if ~isempty(idstart) && ~isempty(idend)
            % Remove truncated (edge effect)
            if idend(1) < idstart(1)
                idend(1)=[];
            end
            if length(idstart)>length(idend)
                idstart(end)=[];
            end
            % Length
            current_L1 = idend-idstart;
            %L1 = [L1; current_L1];
            % Location
            for k=1:1:length(idstart)
                Length_along1(idstart(k)+1:idend(k),ky,kz) = current_L1(k);
            end
        end
    end
end

% Direction 2
Length_along2 = zeros(sz);
%L2 = [];
for kx=1:1:sz(1)
    for kz=1:1:sz(3)
        OneDarray = reshape(BW(kx,:,kz),[sz(2),1]);
        %[Length_along2] = get_lengtharray(OneDarray,Length_along2,2,kx,kz);  Using function is too slow...
        % Identify starting and ending points
        d = diff(OneDarray);
        idstart = find(d==1);
        idend = find(d==-1);
        if ~isempty(idstart) && ~isempty(idend)
            % Remove truncated (edge effect)
            if idend(1) < idstart(1)
                idend(1)=[];
            end
            if length(idstart)>length(idend)
                idstart(end)=[];
            end
            % Length
            current_L2 = idend-idstart;
            %L2 = [L2; current_L2];
            % Location
            for k=1:1:length(idstart)
                Length_along2(kx,idstart(k)+1:idend(k),kz) = current_L2(k);
            end
        end
    end
end

% Direction 3
Length_along3 = zeros(sz);
%L3 = [];
for kx=1:1:sz(1)
    for ky=1:1:sz(2)
        OneDarray = reshape(BW(kx,ky,:),[sz(3),1]);
        %[Length_along3] = get_lengtharray(OneDarray,Length_along3,3,kx,ky);  Using function is too slow...
        % Identify starting and ending points
        d = diff(OneDarray);
        idstart = find(d==1);
        idend = find(d==-1);
        if ~isempty(idstart) && ~isempty(idend)
            % Remove truncated (edge effect)
            if idend(1) < idstart(1)
                idend(1)=[];
            end
            if length(idstart)>length(idend)
                idstart(end)=[];
            end
            % Length
            current_L3 = idend-idstart;
            %L3 = [L3; current_L3];
            % Location
            for k=1:1:length(idstart)
                Length_along3(kx,ky,idstart(k)+1:idend(k)) = current_L3(k);
            end
        end
    end
end

Length_along1_1darray = Length_along1; Length_along1_1darray(Length_along1_1darray==0)=[];
Length_along2_1darray = Length_along2; Length_along2_1darray(Length_along2_1darray==0)=[];
Length_along3_1darray = Length_along3; Length_along3_1darray(Length_along3_1darray==0)=[];

%% FILTER
if size_filter>0
    %L1(L1<=size_filter)=[];
    %L2(L2<=size_filter)=[];
    %L3(L3<=size_filter)=[];
    Length_along1(Length_along1<=size_filter)=0;
    Length_along2(Length_along2<=size_filter)=0;
    Length_along3(Length_along3<=size_filter)=0;
    Length_along1_1darray(Length_along1_1darray<=size_filter)=[];
    Length_along2_1darray(Length_along2_1darray<=size_filter)=[];
    Length_along3_1darray(Length_along3_1darray<=size_filter)=[];
end

%% CALCULATE PROBABILITY AND CUMULATIVE DISTRIBUTION FUNCTION
[~,C_L1] = get_cumulativedistfct(Length_along1_1darray);
[~,C_L2] = get_cumulativedistfct(Length_along2_1darray);
[~,C_L3] = get_cumulativedistfct(Length_along3_1darray);

%% COMPARISON WITH ANALYTCAL CUMULATIVE DISTRIBUTION FUNCTION FOR A SPHERE
unique_L1 = unique(Length_along1_1darray);
unique_L2 = unique(Length_along2_1darray);
unique_L3 = unique(Length_along3_1darray);
Rmax = max( [max(unique_L1),max(unique_L2),max(unique_L3)] )/2;
Rs = linspace(Rmax/50,Rmax*5,nR);

% Loop over all potential radius
difference = zeros(nR,3);
for k=1:1:nR
    [~, C] = function_cumulativefct_chordlength_sphere(Rs(k),res);
    x_analytical = C(:,1);
    y_analytical = C(:,2);

    % Interpolate on same grid
    if length(unique_L1)>=2
        x_max = max([max(unique_L1),max(x_analytical)]);
        x_grid = linspace(0,x_max,res);
        y_grid_numeric = interp1(unique_L1,C_L1,x_grid,'linear',0);
        y_grid_numeric(y_grid_numeric>1)=1;
        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
        % Calculate difference
        difference(k,1) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
    else
        difference(k,1) = NaN;
    end

    % Interpolate on same grid
    if length(unique_L2)>=2
        x_max = max([max(unique_L2),max(x_analytical)]);
        x_grid = linspace(0,x_max,res);
        y_grid_numeric = interp1(unique_L2,C_L2,x_grid,'linear',0);
        y_grid_numeric(y_grid_numeric>1)=1;
        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
        % Calculate difference
        difference(k,2) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
    else
        difference(k,2) = NaN;
    end

    % Interpolate on same grid
    if length(unique_L3)>=2
        x_max = max([max(unique_L3),max(x_analytical)]);
        x_grid = linspace(0,x_max,res);
        y_grid_numeric = interp1(unique_L3,C_L3,x_grid,'linear',0);
        y_grid_numeric(y_grid_numeric>1)=1;
        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
        % Calculate difference
        difference(k,3) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
    else
        difference(k,3) = NaN;
    end
end

% Deduce diameter
fitted_diameter = zeros(3,1);
idx = zeros(3,1);
for k=1:1:3
    if ~isnan(difference(1,k))
        idx(k) = find( difference(:,k) == min(difference(:,k)) );
        fitted_diameter(k) = 2*Rs(idx(k));
    else
        idx(k) = NaN;
        fitted_diameter(k) = NaN;
    end
end

% Interpolate on same grid
if ~isnan(idx(1))
    [~, C] = function_cumulativefct_chordlength_sphere(Rs(idx(1)),res);
    x_analytical = C(:,1);
    y_analytical = C(:,2);
    x_max = max([max(unique_L1),max(x_analytical)]);
    x_grid = linspace(0,x_max,res);
    y_grid_numeric = interp1(unique_L1,C_L1,x_grid,'linear',0);
    y_grid_numeric(y_grid_numeric>1)=1;
    id0 = find(x_grid<2);
    if ~isempty(id0)
        y_grid_numeric(id0)=1;
    end
    CDF(1).numerical = y_grid_numeric';
    y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
    CDF(1).analytical = y_grid_analytical';
    CDF(1).x = x_grid';
else
    CDF(1).numerical = NaN;
    CDF(1).analytical = NaN;
    CDF(1).x = NaN; 
end

% Interpolate on same grid
if ~isnan(idx(2))
    [~, C] = function_cumulativefct_chordlength_sphere(Rs(idx(2)),res);
    x_analytical = C(:,1);
    y_analytical = C(:,2);
    x_max = max([max(unique_L2),max(x_analytical)]);
    x_grid = linspace(0,x_max,res);
    y_grid_numeric = interp1(unique_L2,C_L2,x_grid,'linear',0);
    y_grid_numeric(y_grid_numeric>1)=1;
    id0 = find(x_grid<2);
    if ~isempty(id0)
        y_grid_numeric(id0)=1;
    end
    CDF(2).numerical = y_grid_numeric';
    y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
    CDF(2).analytical = y_grid_analytical';
    CDF(2).x = x_grid';
else
    CDF(1).numerical = NaN;
    CDF(1).analytical = NaN;
    CDF(1).x = NaN;
end

% Interpolate on same grid
if ~isnan(idx(3))
    [~, C] = function_cumulativefct_chordlength_sphere(Rs(idx(3)),res);
    x_analytical = C(:,1);
    y_analytical = C(:,2);
    x_max = max([max(unique_L3),max(x_analytical)]);
    x_grid = linspace(0,x_max,res);
    y_grid_numeric = interp1(unique_L3,C_L3,x_grid,'linear',0);
    y_grid_numeric(y_grid_numeric>1)=1;
    id0 = find(x_grid<2);
    if ~isempty(id0)
        y_grid_numeric(id0)=1;
    end
    CDF(3).numerical = y_grid_numeric';
    y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
    CDF(3).analytical = y_grid_analytical';
    CDF(3).x = x_grid';
else
    CDF(1).numerical = NaN;
    CDF(1).analytical = NaN;
    CDF(1).x = NaN;
end

%% FUNCTIONS

    function [LL] = get_lengtharray(OneDarray,LL,dir,ka,kb)
        % Identify starting and ending points
        d = diff(OneDarray);
        idstart = find(d==1);
        idend = find(d==-1);
        if ~isempty(idstart) && ~isempty(idend)
            % Remove truncated (edge effect)
            if idend(1) < idstart(1)
                idend(1)=[];
            end
            if length(idstart)>length(idend)
                idstart(end)=[];
            end
            % Length
            current_L = idend-idstart;
            % Location
            if dir==1
                for k=1:1:length(idstart)
                    LL(idstart(k)+1:idend(k),ka,kb) = current_L(k);
                end
            elseif dir==2
                for k=1:1:length(idstart)
                    LL(ka,idstart(k)+1:idend(k),kb) = current_L(k);
                end
            else
                for k=1:1:length(idstart)
                    LL(ka,kb,idstart(k)+1:idend(k)) = current_L(k);
                end
            end
        end
    end

    function [P_L,C_L] = get_cumulativedistfct(L)
        % Calculate probability
        unique_L = unique(L);
        n_L = length(unique_L);
        P_L = zeros(n_L,1);
        for k=1:1:n_L
            P_L(k) = sum(sum(sum(L==unique_L(k))));
        end
        % Normalize
        P_L=P_L/sum(P_L);
        % Deduce cumulative distribution function
        C_L=zeros(n_L,1);
        for k=1:1:n_L
            idx = find(unique_L>=unique_L(k));
            if ~isempty(idx)
                C_L(k) = sum(P_L(idx));
            end
        end
    end

end

