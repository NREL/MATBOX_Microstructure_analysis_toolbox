function [fitted2D] = Function_particle_size_ChordLength_section_Algorithm(BW,direction_todo,size_filter)

% Parameter (analytical mesh)
res = 400;
nR = 1000;

sz = size(BW);
if length(sz)==2
    sz = [sz 1];
end
fitted2D = [];

for kdir = 1:1:3
    if kdir==1 && direction_todo(kdir)
        fitted2D(kdir).x = 1:1:sz(kdir);
        fitted2D(kdir).diameter_ellipse = zeros(sz(kdir),2);
        fitted2D(kdir).area_ellipse = zeros(sz(kdir),1);
        fitted2D(kdir).diameter_eqdisc = zeros(sz(kdir),1);
        for kx=1:1:sz(kdir)
            slice_ = squeeze(BW(kx,:,:));
            n = sum(sum(slice_));
            if n>0
                sz_slice = size(slice_);

                Length_slice_along2 = zeros(sz_slice);
                for kz=1:1:sz_slice(2)
                    OneDarray = reshape(slice_(:,kz),[sz_slice(1),1]);
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
                        for k=1:1:length(idstart)
                            Length_slice_along2(idstart(k)+1:idend(k),kz) = current_L(k);
                        end
                    end
                end

                Length_slice_along3 = zeros(sz_slice);
                for ky=1:1:sz_slice(1)
                    OneDarray = reshape(slice_(ky,:),[sz_slice(2),1]);
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
                        for k=1:1:length(idstart)
                            Length_slice_along3(ky,idstart(k)+1:idend(k)) = current_L(k);
                        end
                    end
                end

                Length_slice_along2(Length_slice_along2==0)=[];
                Length_slice_along3(Length_slice_along3==0)=[];

                % Size filter
                if size_filter>0
                    Length_slice_along2(Length_slice_along2<=size_filter)=[];
                    Length_slice_along3(Length_slice_along3<=size_filter)=[];
                end

                % Probabily and cumulative distribution function
                [~,C_L2] = get_cumulativedistfct(Length_slice_along2);
                [~,C_L3] = get_cumulativedistfct(Length_slice_along3);

                % Comparison with analytical cumulative distribution function for a disc
                unique_L2 = unique(Length_slice_along2);
                unique_L3 = unique(Length_slice_along3);
                Rmax = max( [max(unique_L2),max(unique_L3)] )/2;
                if isempty(Rmax) || Rmax==0 % || isempty(unique_L2) || isempty(unique_L3)
                    fitted2D(kdir).diameter_ellipse(kx,1) = NaN;
                    fitted2D(kdir).diameter_ellipse(kx,2) = NaN;
                    fitted2D(kdir).area_ellipse(kx,1) = NaN;
                    fitted2D(kdir).diameter_eqdisc(kx,1) = NaN;
                else
                    Rs = linspace(Rmax/50,Rmax*5,nR);
                    % Loop over all potential radius
                    difference = zeros(nR,2);
                    for k=1:1:nR
                        [~, C] = function_cumulativefct_chordlength_disc(Rs(k),res);
                        x_analytical = C(:,1);
                        y_analytical = C(:,2);

                        % Interpolate on same grid
                        x_max = max([max(unique_L2),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_L2,C_L2,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,1) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));

                        % Interpolate on same grid
                        x_max = max([max(unique_L3),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_L3,C_L3,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end                        
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,2) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
                    end

                    % Deduce ellipse radius
                    for k=1:1:2
                        idx = find( difference(:,k) == min(difference(:,k)) );
                        fitted2D(kdir).diameter_ellipse(kx,k) = 2*Rs(idx);
                    end
                    fitted2D(kdir).area_ellipse(kx,1) = pi*fitted2D(kdir).diameter_ellipse(kx,1)/2 *fitted2D(kdir).diameter_ellipse(kx,2)/2;
                    fitted2D(kdir).diameter_eqdisc(kx,1) = 2*sqrt(fitted2D(kdir).area_ellipse(kx,1)/pi);

                end
            else
                fitted2D(kdir).diameter_ellipse(kx,1) = NaN;
                fitted2D(kdir).diameter_ellipse(kx,2) = NaN;
                fitted2D(kdir).area_ellipse(kx,1) = NaN;
                fitted2D(kdir).diameter_eqdisc(kx,1) = NaN;
            end
        end

    elseif kdir==2  && direction_todo(kdir)
        fitted2D(kdir).x = 1:1:sz(kdir);
        fitted2D(kdir).diameter_ellipse = zeros(sz(kdir),2);
        fitted2D(kdir).area_ellipse = zeros(sz(kdir),1);
        fitted2D(kdir).diameter_eqdisc = zeros(sz(kdir),1);        
        for ky=1:1:sz(kdir)
            slice_ = squeeze(BW(:,ky,:));

            n = sum(sum(slice_));
            if n>0
                sz_slice = size(slice_);

                Length_slice_along1 = zeros(sz_slice);
                for kz=1:1:sz_slice(2)
                    OneDarray = reshape(slice_(:,kz),[sz_slice(1),1]);
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
                        for k=1:1:length(idstart)
                            Length_slice_along1(idstart(k)+1:idend(k),kz) = current_L(k);
                        end
                    end
                end

                Length_slice_along3 = zeros(sz_slice);
                for kx=1:1:sz_slice(1)
                    OneDarray = reshape(slice_(kx,:),[sz_slice(2),1]);
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
                        for k=1:1:length(idstart)
                            Length_slice_along3(kx,idstart(k)+1:idend(k)) = current_L(k);
                        end
                    end
                end

                Length_slice_along1(Length_slice_along1==0)=[];
                Length_slice_along3(Length_slice_along3==0)=[];

                % Size filter
                if size_filter>0
                    Length_slice_along1(Length_slice_along1<=size_filter)=[];
                    Length_slice_along3(Length_slice_along3<=size_filter)=[];
                end

                % Probabily and cumulative distribution function
                [~,C_L1] = get_cumulativedistfct(Length_slice_along1);
                [~,C_L3] = get_cumulativedistfct(Length_slice_along3);

                % Comparison with analytical cumulative distribution function for a disc
                unique_L1 = unique(Length_slice_along1);
                unique_L3 = unique(Length_slice_along3);
                Rmax = max( [max(unique_L1),max(unique_L3)] )/2;
                if isempty(Rmax) || Rmax==0 % || isempty(unique_L1) || isempty(unique_L3)
                    fitted2D(kdir).diameter_ellipse(kx,1) = NaN;
                    fitted2D(kdir).diameter_ellipse(kx,2) = NaN;
                    fitted2D(kdir).area_ellipse(kx,1) = NaN;
                    fitted2D(kdir).diameter_eqdisc(kx,1) = NaN;
                else
                    Rs = linspace(Rmax/50,Rmax*5,nR);
                    % Loop over all potential radius
                    difference = zeros(nR,2);
                    for k=1:1:nR
                        [~, C] = function_cumulativefct_chordlength_disc(Rs(k),res);
                        x_analytical = C(:,1);
                        y_analytical = C(:,2);

                        % Interpolate on same grid
                        x_max = max([max(unique_L1),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_L1,C_L1,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end                        
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,1) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));

                        % Interpolate on same grid
                        x_max = max([max(unique_L3),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_L3,C_L3,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end                        
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,2) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
                    end

                    % Deduce ellipse radius
                    for k=1:1:2
                        idx = find( difference(:,k) == min(difference(:,k)) );
                        fitted2D(kdir).diameter_ellipse(ky,k) = 2*Rs(idx);
                    end
                    fitted2D(kdir).area_ellipse(ky,1) = pi*fitted2D(kdir).diameter_ellipse(ky,1)/2 * fitted2D(kdir).diameter_ellipse(ky,2)/2;
                    fitted2D(kdir).diameter_eqdisc(ky,1) = 2*sqrt(fitted2D(kdir).area_ellipse(ky,1)/pi);
                end
            else
                fitted2D(kdir).diameter_ellipse(ky,1) = NaN;
                fitted2D(kdir).diameter_ellipse(ky,2) = NaN;
                fitted2D(kdir).area_ellipse(ky,1) = NaN;
                fitted2D(kdir).diameter_eqdisc(ky,1) = NaN;
            end

        end
    elseif kdir==3  && direction_todo(kdir)
        fitted2D(kdir).x = 1:1:sz(kdir);
        fitted2D(kdir).diameter_ellipse = zeros(sz(kdir),2);
        fitted2D(kdir).area_ellipse = zeros(sz(kdir),1);
        fitted2D(kdir).diameter_eqdisc = zeros(sz(kdir),1);             
        for kz=1:1:sz(kdir)
            slice_ = BW(:,:,kz);
            n = sum(sum(slice_));
            if n>0
                sz_slice = size(slice_);

                Length_slice_along1 = zeros(sz_slice);
                for ky=1:1:sz_slice(2)
                    OneDarray = reshape(slice_(:,ky),[sz_slice(1),1]);
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
                        for k=1:1:length(idstart)
                            Length_slice_along1(idstart(k)+1:idend(k),ky) = current_L(k);
                        end
                    end
                end

                Length_slice_along2 = zeros(sz_slice);
                for kx=1:1:sz_slice(1)
                    OneDarray = reshape(slice_(kx,:),[sz_slice(2),1]);
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
                        for k=1:1:length(idstart)
                            Length_slice_along2(kx,idstart(k)+1:idend(k)) = current_L(k);
                        end
                    end
                end

                Length_slice_along1(Length_slice_along1==0)=[];
                Length_slice_along2(Length_slice_along2==0)=[];

                % Size filter
                if size_filter>0
                    Length_slice_along1(Length_slice_along1<=size_filter)=[];
                    Length_slice_along2(Length_slice_along2<=size_filter)=[];
                end

                % Probabily and cumulative distribution function
                [~,C_L1] = get_cumulativedistfct(Length_slice_along1);
                [~,C_L2] = get_cumulativedistfct(Length_slice_along2);

                % Comparison with analytical cumulative distribution function for a disc
                unique_L1 = unique(Length_slice_along1);
                unique_L2 = unique(Length_slice_along2);
                Rmax = max( [max(unique_L1),max(unique_L2)] )/2;
                if isempty(Rmax) || Rmax==0 % || isempty(unique_L1) || isempty(unique_L2)
                    fitted2D(kdir).diameter_ellipse(kx,1) = NaN;
                    fitted2D(kdir).diameter_ellipse(kx,2) = NaN;
                    fitted2D(kdir).area_ellipse(kx,1) = NaN;
                    fitted2D(kdir).diameter_eqdisc(kx,1) = NaN;
                else
                    Rs = linspace(Rmax/50,Rmax*5,nR);
                    % Loop over all potential radius
                    difference = zeros(nR,2);
                    for k=1:1:nR
                        [~, C] = function_cumulativefct_chordlength_disc(Rs(k),res);
                        x_analytical = C(:,1);
                        y_analytical = C(:,2);

                        % Interpolate on same grid
                        x_max = max([max(unique_L1),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_L1,C_L1,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end                        
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,1) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));

                        % Interpolate on same grid
                        x_max = max([max(unique_L2),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_L2,C_L2,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end                        
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,2) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
                    end

                    % Deduce ellipse radius
                    for k=1:1:2
                        idx = find( difference(:,k) == min(difference(:,k)) );
                        fitted2D(kdir).diameter_ellipse(kz,k) = 2*Rs(idx);
                    end
                    fitted2D(kdir).area_ellipse(kz,1) = pi*fitted2D(kdir).diameter_ellipse(kz,1)/2 * fitted2D(kdir).diameter_ellipse(kz,2)/2;
                    fitted2D(kdir).diameter_eqdisc(kz,1) = 2*sqrt(fitted2D(kdir).area_ellipse(kz,1)/pi);
                end
            else
                fitted2D(kdir).diameter_ellipse(kz,1) = NaN;
                fitted2D(kdir).diameter_ellipse(kz,2) = NaN;
                fitted2D(kdir).area_ellipse(kz,1) = NaN;
                fitted2D(kdir).diameter_eqdisc(kz,1) = NaN;
            end
            
        end
    else
        fitted2D(kdir).diameter_eqdisc = [];
    end
end

%% FUNCTIONS
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