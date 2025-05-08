function [Ms,newtype,foo] = fct_calculateconnectivity2(M,p)

newtype = 'Segmented (phase)';
foo = 1;

sz = size(M);
dimension = length(sz);


if strcmp(p.definition,'Isotropic')
    % Id
    MC = 1; % Main cluster
    IC = 2; % Isolated cluster(s)
    UC = 3; % Unknown cluster(s)
    if ~p.backgroundisdomainedges
        MC = MC+1;
        IC = IC+1;
        UC = UC+1;
    end

    if p.backgroundisdomainedges
        Ms = zeros(sz);
        Ms(M~=p.complementarylabel)=IC;
    else
        Ms = M;
        cond1 = double(M==p.complementarylabel);
        cond2 = double(M==p.backgroundlabel);
        Ms(cond1 + cond2 ==0) = IC;
        %Ms(M~=p.complementarylabel)=IC;
        %Ms(M~=p.backgroundlabel)=IC;
        % Ms(M>1)=IC;
    end

    % All clusters and their size
    unis = double(unique(M));
    count = histc(reshape(M,[numel(M) 1]),unis);
    id0 = find(unis==p.complementarylabel); % Remove background
    unis(id0,:)=[];
    count(id0,:)=[];
    if ~p.backgroundisdomainedges
        idc = find(unis==p.backgroundlabel); % Remove background
        unis(idc,:)=[];
        count(idc,:)=[];
    end

    res = [unis count];

    % Identify main cluster
    res = sortrows(res,2,"descend");
    % Ms(M==res(1,1))=MC; % Not yet

    if p.backgroundisdomainedges
        x0 = unique(M(1,:,:));
        x1 = unique(M(end,:,:));
        y0 = unique(M(:,1,:));
        y1 = unique(M(:,end,:));
        x0 = reshape(x0,[length(x0) 1]);
        x1 = reshape(x1,[length(x1) 1]);
        y0 = reshape(y0,[length(y0) 1]);
        y1 = reshape(y1,[length(y1) 1]);        
        if dimension == 2
            attheborders = unique([x0;x1;y0;y1]);
        else
            z0 = unique(M(:,:,1));
            z1 = unique(M(:,:,end));
            z0 = reshape(z0,[length(z0) 1]);
            z1 = reshape(z1,[length(z1) 1]);            
            attheborders = unique([x0;x1;y0;y1;z0;z1]);
        end
        attheborders(attheborders==p.complementarylabel)=[];
    else
        BW = zeros(sz);
        BW(M==p.backgroundlabel)=1;
        tmp = ones(sz+2);
        if dimension == 2
            tmp(2:end-1,2:end-1) = BW;
        else
            tmp(2:end-1,2:end-1,2:end-1) = BW;
        end
        dmap = bwdist(tmp,'chessboard');
        if dimension == 2
            attheborders = unique(M(dmap(2:end-1,2:end-1)==1));
        else
            attheborders = unique(M( dmap(2:end-1,2:end-1,2:end-1)==1 ));
        end
        attheborders(attheborders==p.complementarylabel)=[];
        attheborders(attheborders==p.backgroundlabel)=[];
    end

    if ~isempty(attheborders)
        for k=1:1:length(attheborders)
            Ms(M==attheborders(k))=UC;
        end
    end

    % Overwritte main cluster
    Ms(M==res(1,1))=MC;

elseif strcmp(p.definition,'Anisotropic')
    % Connected
    ax_conn = str2num(p.anisotropy_alongaxe);
    if ax_conn==3 && dimension==2
        Ms = M;
        warning('Connectivity rules are for a 3D volume, not a 2D image')
        return
    end

    % Unknown
    ax_unkn = 1:1:dimension;
    ax_unkn(ax_unkn==ax_conn) = [];

    % Id
    MC = 1; % Connected cluster
    IC = 2; % Isolated cluster(s)
    UC = 3; % Unknown cluster(s)
    if ~p.backgroundisdomainedges
        MC = MC+1;
        IC = IC+1;
        UC = UC+1;
    end

    % % Initialize with isolated clusters
    Ms = zeros(sz);
    if p.backgroundisdomainedges
        Ms(M~=p.complementarylabel)=IC;
    else
        %Ms(M~=p.complementarylabel)=IC;
        %Ms(M~=p.backgroundlabel)=IC;
        cond1 = double(M==p.complementarylabel);
        cond2 = double(M==p.backgroundlabel);
        Ms(cond1 + cond2 ==0) = IC;
    end

    % % Find unknown clusters
    x0 = []; x1 = []; y0 = []; y1 = []; z0 = []; z1 = [];
    if sum(ax_unkn==1)
        x0 = unique(M(1,:,:));
        x1 = unique(M(end,:,:));
    end
    if sum(ax_unkn==2)
        y0 = unique(M(:,1,:));
        y1 = unique(M(:,end,:));
    end
    if sum(ax_unkn==3)
        z0 = unique(M(:,:,1));
        z1 = unique(M(:,:,end));
    end
    attheborders = unique([x0;x1;y0;y1;z0;z1]);
    if ~p.backgroundisdomainedges
        BW = zeros(sz);
        BW(M==p.backgroundlabel)=1;
        dmap = bwdist(BW,'chessboard');
        attheborders = unique([attheborders; unique(M(dmap==1))]);
    end
    attheborders(attheborders==p.complementarylabel)=[];
    if ~p.backgroundisdomainedges
        attheborders(attheborders==p.backgroundlabel)=[];
    end

    if ~isempty(attheborders)
        for k=1:1:length(attheborders)
            Ms(M==attheborders(k))=UC;
        end
    end

    % % Find connected clusters
    if strcmp(p.anisotropy_musttouch,'first slice')
        sl_conn = 1;
    elseif strcmp(p.anisotropy_musttouch,'last slice')
        sl_conn = sz(ax_conn);
    elseif strcmp(p.anisotropy_musttouch,'both opposite slices')
        sl_conn = [1 sz(ax_conn)];
    end

    for k=1:1:length(sl_conn)
        if ax_conn==1
            connected_labels(k).vals = unique(M(sl_conn(k),:,:));
        elseif ax_conn==2
            connected_labels(k).vals = unique(M(:,sl_conn(k),:));
        elseif ax_conn==3
            connected_labels(k).vals = unique(M(:,:,sl_conn(k)));
        end
    end
    if length(sl_conn) > 1
        connected_labels = intersect(connected_labels(1).vals, connected_labels(2).vals);
    else
        connected_labels = connected_labels(1).vals;
    end

    connected_labels(connected_labels==p.complementarylabel)=[];
    if ~p.backgroundisdomainedges
        connected_labels(connected_labels==p.backgroundlabel)=[];
    end

    if ~isempty(connected_labels)
        for k=1:1:length(connected_labels)
            Ms(M==connected_labels(k))=MC;
        end
    end


elseif strcmp(p.definition,'Custom')
    if sum(p.custom_rules(:,3))>0 && dimension==2
        Ms = M;
        warning('Connectivity rules are for a 3D volume, not a 2D image')
        return
    end

    if sum(sum(p.custom_rules)) < 1
        Ms = M;
        warning('No connectivity rules!')
        return
    end

    % Id
    MC = 1; % Connected cluster
    IC = 2; % Isolated cluster(s)
    UC = 3; % Unknown cluster(s)
    if ~p.backgroundisdomainedges
        MC = MC+1;
        IC = IC+1;
        UC = UC+1;
    end

    % % Initialize with isolated clusters
    Ms = zeros(sz);
    if p.backgroundisdomainedges
        Ms(M~=p.complementarylabel)=IC;
    else
        %Ms(M~=p.complementarylabel)=IC;
        %Ms(M~=p.backgroundlabel)=IC;
        cond1 = double(M==p.complementarylabel);
        cond2 = double(M==p.backgroundlabel);
        Ms(cond1 + cond2 ==0) = IC;        
    end

    % % Find unknown clusters: either touch domain's edges or background if any
    x0 = []; x1 = []; y0 = []; y1 = []; z0 = []; z1 = [];
    if ~p.custom_rules(1,1)
        x0 = unique(M(1,:,:));
    end
    if ~p.custom_rules(2,1)
        x1 = unique(M(end,:,:));
    end

    if ~p.custom_rules(1,2)
        y0 = unique(M(:,1,:));
    end
    if ~p.custom_rules(2,2)
        y1 = unique(M(:,end,:));
    end

    if ~p.custom_rules(1,3)
        z0 = unique(M(:,:,1));
    end
    if ~p.custom_rules(2,3)
        z1 = unique(M(:,:,end));
    end      
    attheborders = unique([x0;x1;y0;y1;z0;z1]);
    if ~p.backgroundisdomainedges
        BW = zeros(sz);
        BW(M==p.backgroundlabel)=1;
        dmap = bwdist(BW,'chessboard');
        attheborders = unique([attheborders; unique(M(dmap==1))]);
    end
    attheborders(attheborders==p.complementarylabel)=[];
    if ~p.backgroundisdomainedges
        attheborders(attheborders==p.backgroundlabel)=[];
    end

    if ~isempty(attheborders)
        for k=1:1:length(attheborders)
            Ms(M==attheborders(k))=UC;
        end
    end

    % % Find connected clusters
    k=0;
    if p.custom_rules(1,1)
        k=k+1;
        labels(k).vals = unique(M(1,:,:));
    end
    if p.custom_rules(2,1)
        k=k+1;
        labels(k).vals = unique(M(end,:,:));
    end
    if p.custom_rules(1,2)
        k=k+1;
        labels(k).vals = unique(M(:,1,:));
    end
    if p.custom_rules(2,2)
        k=k+1;
        labels(k).vals = unique(M(:,end,:));
    end
    if p.custom_rules(1,3)
        k=k+1;
        labels(k).vals = unique(M(:,:,1));
    end
    if p.custom_rules(2,3)
        k=k+1;
        labels(k).vals = unique(M(:,:,end));
    end
    if strcmp(p.custom_mustverify,'all connectivity conditions')
        connected_labels = labels(1).vals;
        if k > 1
            for kk=2:1:k
                connected_labels = intersect(labels(kk-1).vals, labels(kk).vals);
            end
        end
        connected_labels = unique(connected_labels);
    elseif strcmp(p.custom_mustverify,'at least one connectivity condition')
        connected_labels = labels(1).vals;
        for kk=2:1:k
            connected_labels = union(connected_labels, labels(kk).vals);
        end
        connected_labels = unique(connected_labels);
    end

    connected_labels(connected_labels==p.complementarylabel)=[];
    if ~p.backgroundisdomainedges
        connected_labels(connected_labels==p.backgroundlabel)=[];
    end

    if ~isempty(connected_labels)
        for k=1:1:length(connected_labels)
            Ms(M==connected_labels(k))=MC;
        end
    end

end

[Ms] = fct_intconvert(Ms);

end