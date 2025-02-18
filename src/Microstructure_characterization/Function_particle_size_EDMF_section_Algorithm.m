function [fitted2D] = Function_particle_size_EDMF_section_Algorithm(BW,direction_todo,removeinclusion)

% Parameter (analytical mesh)
res = 400;
nR = 1000;

sz = size(BW);
fitted2D = [];

for kdir = 1:1:3
    if kdir==1 && direction_todo(kdir)
        fitted2D(kdir).x = 1:1:sz(kdir);
        fitted2D(kdir).diameter_eqdisc = zeros(sz(kdir),1);
        for kx=1:1:sz(kdir)
            slice_ = squeeze(BW(kx,:,:));
            n = sum(sum(slice_));
            if n>0

                % Inclusion removal
                slice_true = slice_;
                if removeinclusion
                    L = bwlabel(~slice_,4);
                    [C,~,ic] = unique(L);
                    if length(C)>2
                        counts = accumarray(ic,1);
                        value_counts = [C, counts];
                        idx_binary = find(value_counts(:,1)==0);
                        value_counts(idx_binary,:)=[];
                        value_counts = sortrows(value_counts,-2);
                        for k=2:1:length(C)-1
                            slice_(L==value_counts(k,1)) = 1;
                        end
                    end
                end

                % EDM
                distance_transform=double(bwdist(~slice_,'euclidean')); % Euclidean distance map
                %distance_transform=distance_transform-0.5; % From voxel center to surface
                distance_transform(distance_transform<=0)=0; % Complementary volume
                distance_transform(slice_true~=1)=0; % Complementary volume

                % 1D ARRAY
                distance_transform_vector=reshape(distance_transform,[numel(slice_),1]); % Convert in a vector
                distance_transform_vector(distance_transform_vector<=0)=[]; % % Remove all 0

                % Cumulative function
                [~,C_D] = get_cumulativedistfct(distance_transform_vector);

                % Comparison with analytical cumulative distribution function for a disc
                unique_D = unique(distance_transform_vector);
                Rmax = max(unique_D);
                Rs = linspace(Rmax/50,Rmax*5,nR);
                if Rmax==0
                    fitted2D(kdir).diameter_eqdisc(kx) = NaN;
                else
                    % Loop over all potential radius
                    difference = zeros(nR,1);
                    for k=1:1:nR
                        R=Rs(k);
                        d=linspace(0,R,res);
                        C = (R-d).^2./R.^2;
                        x_analytical = d;
                        y_analytical = C;

                        % Interpolate on same grid
                        x_max = max([max(distance_transform_vector),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_D,C_D,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,1) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
                    end

                    % Deduce disc radius
                    idx = find( difference == min(difference) );
                    fitted2D(kdir).diameter_eqdisc(kx) = 2*Rs(idx);
                end
            else
                fitted2D(kdir).diameter_eqdisc(kx) = NaN;
            end
        end

    elseif kdir==2  && direction_todo(kdir)
        fitted2D(kdir).x = 1:1:sz(kdir);
        fitted2D(kdir).diameter_eqdisc = zeros(sz(kdir),1);
        for ky=1:1:sz(kdir)
            slice_ = squeeze(BW(:,ky,:));
            n = sum(sum(slice_));
            if n>0

                % Inclusion removal
                slice_true = slice_;
                if removeinclusion
                    L = bwlabel(~slice_,4);
                    [C,~,ic] = unique(L);
                    if length(C)>2
                        counts = accumarray(ic,1);
                        value_counts = [C, counts];
                        idx_binary = find(value_counts(:,1)==0);
                        value_counts(idx_binary,:)=[];
                        value_counts = sortrows(value_counts,-2);
                        for k=2:1:length(C)-1
                            slice_(L==value_counts(k,1)) = 1;
                        end
                    end
                end

                % EDM
                distance_transform=double(bwdist(~slice_,'euclidean')); % Euclidean distance map
                %distance_transform=distance_transform-0.5; % From voxel center to surface
                distance_transform(distance_transform<=0)=0; % Complementary volume
                distance_transform(slice_true~=1)=0; % Complementary volume

                % 1D ARRAY
                distance_transform_vector=reshape(distance_transform,[numel(slice_),1]); % Convert in a vector
                distance_transform_vector(distance_transform_vector<=0)=[]; % % Remove all 0

                % Cumulative function
                [~,C_D] = get_cumulativedistfct(distance_transform_vector);

                % Comparison with analytical cumulative distribution function for a disc
                unique_D = unique(distance_transform_vector);
                Rmax = max(unique_D);
                Rs = linspace(Rmax/50,Rmax*5,nR);
                if Rmax==0
                    fitted2D(kdir).diameter_eqdisc(ky) = NaN;
                else
                    % Loop over all potential radius
                    difference = zeros(nR,1);
                    for k=1:1:nR
                        R=Rs(k);
                        d=linspace(0,R,res);
                        C = (R-d).^2./R.^2;
                        x_analytical = d;
                        y_analytical = C;

                        % Interpolate on same grid
                        x_max = max([max(distance_transform_vector),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_D,C_D,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,1) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
                    end

                    % Deduce disc radius
                    idx = find( difference == min(difference) );
                    fitted2D(kdir).diameter_eqdisc(ky) = 2*Rs(idx);
                end
            else
                fitted2D(kdir).diameter_eqdisc(ky) = NaN;
            end
        end

    elseif kdir==3  && direction_todo(kdir)
        fitted2D(kdir).x = 1:1:sz(kdir);
        fitted2D(kdir).diameter_eqdisc = zeros(sz(kdir),1);
        for kz=1:1:sz(kdir)
            slice_ = BW(:,:,kz);
            n = sum(sum(slice_));
            if n>0

                % Inclusion removal
                slice_true = slice_;
                if removeinclusion
                    L = bwlabel(~slice_,4);
                    [C,~,ic] = unique(L);
                    if length(C)>2
                        counts = accumarray(ic,1);
                        value_counts = [C, counts];
                        idx_binary = find(value_counts(:,1)==0);
                        value_counts(idx_binary,:)=[];
                        value_counts = sortrows(value_counts,-2);
                        for k=2:1:length(C)-1
                            slice_(L==value_counts(k,1)) = 1;
                        end
                    end
                end

                % EDM
                distance_transform=double(bwdist(~slice_,'euclidean')); % Euclidean distance map
                %distance_transform=distance_transform-0.5; % From voxel center to surface
                distance_transform(distance_transform<=0)=0; % Complementary volume
                distance_transform(slice_true~=1)=0; % Complementary volume

                % 1D ARRAY
                distance_transform_vector=reshape(distance_transform,[numel(slice_),1]); % Convert in a vector
                distance_transform_vector(distance_transform_vector<=0)=[]; % % Remove all 0

                % Cumulative function
                [~,C_D] = get_cumulativedistfct(distance_transform_vector);

                % Comparison with analytical cumulative distribution function for a disc
                unique_D = unique(distance_transform_vector);
                Rmax = max(unique_D);
                Rs = linspace(Rmax/50,Rmax*5,nR);
                if Rmax==0
                    fitted2D(kdir).diameter_eqdisc(kz) = NaN;
                else
                    % Loop over all potential radius
                    difference = zeros(nR,1);
                    for k=1:1:nR
                        R=Rs(k);
                        d=linspace(0,R,res);
                        C = (R-d).^2./R.^2;
                        x_analytical = d;
                        y_analytical = C;

                        % Interpolate on same grid
                        x_max = max([max(distance_transform_vector),max(x_analytical)]);
                        x_grid = linspace(0,x_max,res);
                        y_grid_numeric = interp1(unique_D,C_D,x_grid,'linear',0);
                        y_grid_numeric(y_grid_numeric>1)=1;
                        id0 = find(x_grid<2);
                        if ~isempty(id0)
                            y_grid_numeric(id0)=1;
                        end
                        y_grid_analytical = interp1(x_analytical,y_analytical,x_grid,'linear',0);
                        % Calculate difference
                        difference(k,1) = trapz(x_grid,abs(y_grid_numeric-y_grid_analytical));
                    end

                    % Deduce disc radius
                    idx = find( difference == min(difference) );
                    fitted2D(kdir).diameter_eqdisc(kz) = 2*Rs(idx);
                end
            else
                fitted2D(kdir).diameter_eqdisc(kz) = NaN;
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