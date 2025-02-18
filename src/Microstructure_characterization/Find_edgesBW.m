function [Edges] = Find_edgesBW(BW_original)

sz_original = size(BW_original);
dimension = length(size(BW_original));

BW = zeros(sz_original+2);
sz = size(BW);

if dimension == 2
    BW(2:end-1,2:end-1) = BW_original;
    idx = find(BW==1);

    % Edges = zeros(sz);
    % Edges(:,1:end-1) = BW(:,1:end-1)+BW(:,2:end);
    % Edges(:,2:end) = Edges(:,1:end-1);
    % Edges(1:end-1,:) = BW(1:end-1,:)+BW(2:end,:);
    % Edges(2:end,:) = Edges(1:end-1,:);
    % Edges(BW==0)=0;

    [IX, IY] = ind2sub(sz,idx);
    IXp = IX+1;
    IXm = IX-1;
    IYp = IY+1;
    IYm = IY-1;

    idx_IYp = sub2ind(sz,IX,IYp);
    idx_IYm = sub2ind(sz,IX,IYm);
    idx_IXp = sub2ind(sz,IXp,IY);
    idx_IXm = sub2ind(sz,IXm,IY);

    tmp = zeros(sz);
    tmp(idx_IYp)=1;
    id_edge_yp = BW+tmp==1;
    tmp = zeros(sz);
    tmp(idx_IYm)=1;
    id_edge_ym = BW+tmp==1;
    tmp = zeros(sz);
    tmp(idx_IXp)=1;
    id_edge_xp = BW+tmp==1;
    tmp = zeros(sz);
    tmp(idx_IXm)=1;
    id_edge_xm = BW+tmp==1;    

    Edges = BW;
    Edges(id_edge_yp)=2;
    Edges(id_edge_ym)=2;
    Edges(id_edge_xp)=2;
    Edges(id_edge_xm)=2;
    Edges = Edges.*BW;

    Edges=Edges(2:end-1,2:end-1);

else
    BW(2:end-1,2:end-1,2:end-1) = BW_original;
    idx = find(BW==1);

    [IX, IY, IZ] = ind2sub(sz,idx);
    IXp = IX+1;
    IXm = IX-1;
    IYp = IY+1;
    IYm = IY-1;
    IZp = IZ+1;
    IZm = IZ-1;    

    idx_IYp = sub2ind(sz,IX,IYp,IZ);
    idx_IYm = sub2ind(sz,IX,IYm,IZ);
    idx_IXp = sub2ind(sz,IXp,IY,IZ);
    idx_IXm = sub2ind(sz,IXm,IY,IZ);
    idx_IZp = sub2ind(sz,IX,IY,IZp);
    idx_IZm = sub2ind(sz,IX,IY,IZm);    

    tmp = zeros(sz);
    tmp(idx_IYp)=1;
    id_edge_yp = BW+tmp==1;

    tmp = zeros(sz);
    tmp(idx_IYm)=1;
    id_edge_ym = BW+tmp==1;

    tmp = zeros(sz);
    tmp(idx_IXp)=1;
    id_edge_xp = BW+tmp==1;

    tmp = zeros(sz);
    tmp(idx_IXm)=1;
    id_edge_xm = BW+tmp==1;    

    tmp = zeros(sz);
    tmp(idx_IZp)=1;
    id_edge_zp = BW+tmp==1;

    tmp = zeros(sz);
    tmp(idx_IZm)=1;
    id_edge_zm = BW+tmp==1;   

    Edges = BW;
    Edges(id_edge_yp)=2;
    Edges(id_edge_ym)=2;
    Edges(id_edge_xp)=2;
    Edges(id_edge_xm)=2;
    Edges(id_edge_zp)=2;
    Edges(id_edge_zm)=2;    
    Edges = Edges.*BW;

    Edges=Edges(2:end-1,2:end-1,2:end-1);

end

end