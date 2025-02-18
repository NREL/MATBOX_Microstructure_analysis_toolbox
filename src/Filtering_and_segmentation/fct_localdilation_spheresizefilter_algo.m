function [M,diameter] = fct_localdilation_spheresizefilter_algo(M,p)

sz = size(M);
dimension = length(sz);

if p.extendFOV % Extend FOV to avoid edge effect
    sz = sz + 2*(p.spherediameter+1);
    tmp = zeros(sz)+p.background;
    if dimension==2
        tmp(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1) = M;
    else
        tmp(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1) = M;
    end
    M = tmp;
end

labels = unique(M);
nonbackground_labels = labels;
nonbackground_labels(nonbackground_labels==p.background)=[];

Background = zeros(sz);
Background(M == p.background)=1;
if strcmp(p.choice,'One label') && p.removeotherlabel && length(nonbackground_labels)>1
    cond1 = M ~= p.labeltodilate;
    cond2 = M ~= p.background;
    idx_otherlabels = find(cond1.*cond2==1);
    Background(idx_otherlabels) = 1;
    sav_otherlabel = M(idx_otherlabels);
end

approximate=false;
if p.per_cluster
    if dimension==2
        C = bwlabel(~Background,4);
    else
        C = bwlabeln(~Background,6);
    end
    unis = unique(C);
    unis(unis==0)=[];
    ncluster = length(unis);
    for k=1:1:ncluster
        %k/ncluster

        % if k==4
        %     keyboard
        % end
        %k

        % Extract subdomain centered on initial cluster
        idx = find(C==unis(k));
        if dimension==2
            [IX, IY] = ind2sub(sz,idx);
            z_min = 1; z_max = 1;
        else
            [IX, IY, IZ] = ind2sub(sz,idx);
            z_min = min(min(min(IZ))); z_max = max(max(max(IZ)));
        end
        x_min = min(min(min(IX))); x_max = max(max(max(IX)));
        y_min = min(min(min(IY))); y_max = max(max(max(IY)));
        M_sub = M(x_min:x_max, y_min:y_max, z_min:z_max); % We extract the domain that is modified at each new cluster, not the original one...
        sz_sub = size(M_sub);
        subdimension = length(sz_sub);

        % ... therefore we have to recalculate the background as well
        Background_sub = zeros(sz_sub);
        Background_sub(M_sub == p.background)=1;
        if strcmp(p.choice,'One label') && p.removeotherlabel && length(nonbackground_labels)>1
            cond1 = M_sub ~= p.labeltodilate;
            cond2 = M_sub ~= p.background;
            Background_sub(cond1.*cond2==1) = 1;
        end      

        if subdimension==2
            C_sub = bwlabel(~Background_sub,4);
        else
            C_sub = bwlabeln(~Background_sub,6);
        end
        unis_sub = unique(C_sub);
        unis_sub(unis_sub==0)=[];
        ncluster_sub = length(unis_sub);
        if ncluster_sub>1 % We have to identify the correct cluster, as the domain is modified at each new cluster
            tmp = C(x_min:x_max, y_min:y_max, z_min:z_max);
            tmp(tmp~=unis(k))=0;
            tmp(tmp~=0)=1;
            overlap = tmp.*C_sub; % We will select the cluster calculated in the current subdomain that overlap the most with the original cluster
            [counts, groupnames] = groupcounts(reshape(overlap,[numel(overlap) 1]));
            tmp = [counts groupnames];
            idtmp = find(groupnames==0);
            tmp(idtmp,:)=[];
            tmp = sortrows(tmp,2,'descend');
            Background_sub(C_sub~=tmp(1,2))=1;
        end
        
        sz_sub = sz_sub + 2*(p.spherediameter+1);
        M_sub_ext = zeros(sz_sub)+p.background;    
        Background_sub_ext = ones(sz_sub);  
        if subdimension==2
            M_sub_ext(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1) = M_sub;
            Background_sub_ext(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1) = Background_sub;
        else
            M_sub_ext(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1) = M_sub;
            Background_sub_ext(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1) = Background_sub;
        end
        [diameter_sub,~,IDX_sub] = Function_particle_size_CPSD_Algorithm(Background_sub_ext,p.round_dmap,approximate);
            

        % figure; imagesc(diameter_sub);axis equal; axis tight;
        % unis_ = unique(diameter_sub);
        % col = turbo(length(unis_));
        % col(1,:)=[0.5 0.5 0.5];
        % colormap(col);
        % pause(0.1);



        labels_sub = unique(M_sub);
        nonbackground_labels_sub = labels_sub;
        nonbackground_labels_sub(nonbackground_labels_sub==p.background)=[];

        if strcmp(p.choice,'One label')
            cond1 = diameter_sub<=p.spherediameter;
            cond2 = M_sub_ext(IDX_sub)==p.labeltodilate;
            M_sub_ext(cond1.*cond2==1) = p.labeltodilate;
        elseif strcmp(p.choice,'All non-background labels')
            idx = find(diameter_sub<=p.spherediameter);
            M_sub_ext(idx) = M_sub_ext(IDX_sub(idx));
        end

        % Merge
        if subdimension==2
            Msub_new = M_sub_ext(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1);
        else
            Msub_new = M_sub_ext(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1);
        end
        M(x_min:x_max, y_min:y_max, z_min:z_max) = Msub_new;
    end
    diameter = zeros(sz);

else
    [diameter,~,IDX] = Function_particle_size_CPSD_Algorithm(Background,p.round_dmap,approximate);
    if length(nonbackground_labels)>1
        if strcmp(p.choice,'One label')
            cond1 = diameter<=p.spherediameter;
            cond2 = M(IDX)==p.labeltodilate;
            M(cond1.*cond2==1) = p.labeltodilate;
        elseif strcmp(p.choice,'All non-background labels')
            idx = find(diameter<=p.spherediameter);
            M(idx) = M(IDX(idx));
        end
    else
        M(diameter<=p.spherediameter)=nonbackground_labels;
    end
end

if strcmp(p.choice,'One label') && p.removeotherlabel && length(nonbackground_labels)>1
    M(idx_otherlabels) = sav_otherlabel;
end

if p.extendFOV % Shrink it back
    if dimension==2
        M = M(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1);
        diameter = diameter(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1);
    else
        M = M(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1);
        diameter = diameter(p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1, p.spherediameter+2:end-p.spherediameter-1);
    end
end


end