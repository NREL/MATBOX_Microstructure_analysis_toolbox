function [Labels,newtype] = Assign_voxel_tonearest_cluster(MA,MB)

newtype = 'Segmented (instance)';

sz = size(MB);
dimension = length(sz);

Labels = MB;

% Assign missing voxels to nearest cluster
[~, IDX] = bwdist(MB);
Missingvoxels = ~double(MB).*double(MA);
id_missing = find(Missingvoxels);

Labels(id_missing) = Labels(IDX(id_missing));

% Check contiguity clusters
unis = unique(Labels);
unis(unis==0)=[];
m = max(max(max(unis)));
for k=1:1:length(unis)

    % Sub volume
    idx = find(Labels==unis(k));
    if dimension == 2
        [IX,IY] = ind2sub(sz,idx);
        zmin = 1; zmax=1;
    else
        [IX,IY,IZ] = ind2sub(sz,idx);
        zmin = min(IZ); zmax = max(IZ);
    end
    xmin = min(IX); xmax = max(IX);
    ymin = min(IY); ymax = max(IY);
    sub = Labels(xmin:xmax, ymin:ymax, zmin:zmax);
    sz_sub = size(sub);

    % Connectivity
    Cluster = zeros(sz_sub);
    Cluster(sub==unis(k))=1;
    if dimension==2
        Ctmp = bwlabel(Cluster,4);
    else
        Ctmp = bwlabeln(Cluster,6);
    end
    unistmp = unique(Ctmp);
    unistmp(unistmp==0)=[];
    n = length(unistmp);
    if n>1
        for kk=2:1:n
            m=m+1;
            sub(Ctmp==unistmp(kk))=m;
        end
    end

    % Apply to full volume
    Labels(xmin:xmax, ymin:ymax, zmin:zmax) = sub;




    % Cluster = zeros(sz);
    % Cluster(Labels==unis(k))=1;
    % if dimension==2
    %     Ctmp = bwlabel(Cluster,4);
    % else
    %     Ctmp = bwlabeln(Cluster,6);
    % end
    % unistmp = unique(Ctmp);
    % unistmp(unistmp==0)=[];
    % n = length(unistmp);
    % if n>1
    %     for kk=2:1:n
    %         m=m+1;
    %         Labels(Ctmp==unistmp(kk))=m;
    %     end
    % end
end

[Labels] = fct_intconvert(Labels);

end