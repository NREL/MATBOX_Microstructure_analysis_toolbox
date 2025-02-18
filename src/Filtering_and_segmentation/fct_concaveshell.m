function [shells,newtype,foo] = fct_concaveshell(Labels,p)
foo=[];
newtype = 'same';

sz = size(Labels);
dimension = length(sz);

shells = zeros(sz);
unique_labels = unique(Labels);
unique_labels(unique_labels==0)=[];
n_labels = length(unique_labels);

for k=1:1:n_labels
    % Select label
    idx = find( Labels==unique_labels(k) );
    [IX, IY, IZ] = ind2sub(sz,idx);
    x_min = min(IX); x_max = max(IX);
    y_min = min(IY); y_max = max(IY);
    z_min = min(IZ); z_max = max(IZ);
    BW = Labels(x_min:x_max, y_min:y_max, z_min:z_max);
    BW(BW~=unique_labels(k))=0;
    BW(BW~=0)=1;

    % Add empty layer
    if length(BW)>1
        sz_tmp = size(BW);
        if length(sz_tmp)==dimension
            tmp = zeros(sz_tmp+2);
            if dimension==2
                tmp(2:end-1,2:end-1) = BW;
            else
                tmp(2:end-1,2:end-1,2:end-1) = BW;
            end
        else
            tmp = zeros(sz_tmp(1)+2,sz_tmp(2)+2,3);
            tmp(2:end-1,2:end-1,2) = BW(:,:,1);
        end
    else
        if dimension==2
            tmp = zeros(3,3);
            tmp(2,2) = BW;
        else
            tmp = zeros(3,3,3);
            tmp(2,2,2) = BW;
        end
    end
    BW = tmp;
    tmp_sz = size(BW);
    tmp_shell = zeros(tmp_sz);


    diffrow = diff(BW,1,1);
    diffcol = diff(BW,1,2);

    if dimension==2
        for y=1:1:tmp_sz(2)
            if length(unique(diffrow(:,y,:)))>1
                ids = find(diffrow(:,y,:)~=0);
                id_start = min(ids)+1;
                id_end = max(ids)+1;
                tmp_shell(id_start:id_end,y,:)=1;
            end
        end
        for x=1:1:tmp_sz(1)
            if length(unique(diffcol(x,:,:)))>1
                ids = find(diffcol(x,:,:) ~=0);
                id_start = min(ids)+1;
                id_end = max(ids)+1;
                tmp_shell(x,id_start:id_end,:)=tmp_shell(x,id_start:id_end,:)+1;
            end
        end
    elseif dimension==3
        diffhei = diff(BW,1,3);

        for y=1:1:tmp_sz(2)
            for z=1:1:tmp_sz(3)
                if length(unique(diffrow(:,y,z)))>1
                    ids = find(diffrow(:,y,z)~=0);
                    id_start = min(ids)+1;
                    id_end = max(ids)+1;
                    tmp_shell(id_start:id_end,y,z)=1;
                end
            end
        end

        for x=1:1:tmp_sz(1)
            for z=1:1:tmp_sz(3)
                if length(unique(diffcol(x,:,z)))>1
                    ids = find(diffcol(x,:,z) ~=0);
                    id_start = min(ids)+1;
                    id_end = max(ids)+1;
                    tmp_shell(x,id_start:id_end,z)=tmp_shell(x,id_start:id_end,z)+1;
                end
            end
        end

        for x=1:1:tmp_sz(1)
            for y=1:1:tmp_sz(2)
                if length(unique(diffhei(x,y,:)))>1
                    ids = find(diffhei(x,y,:)~=0);
                    id_start = min(ids)+1;
                    id_end = max(ids)+1;
                    tmp_shell(x,y,id_start:id_end)=tmp_shell(x,y,id_start:id_end)+1;
                end
            end
        end     

    end

    % Get concave shell
    tmp_shell(tmp_shell~=dimension)=0;
    tmp_shell(tmp_shell~=0)=unique_labels(k);

    % Remove thin layer
    if dimension == 2
        tmp_shell = tmp_shell(2:end-1,2:end-1);
    else
        tmp_shell = tmp_shell(2:end-1,2:end-1,2:end-1);
    end

    % Insert in full array
    old_shells = shells(x_min:x_max, y_min:y_max, z_min:z_max);
    old_shells(tmp_shell==unique_labels(k)) = unique_labels(k); % Overwritte
    shells(x_min:x_max, y_min:y_max, z_min:z_max) = old_shells;
end

[shells] = fct_intconvert(shells);

end