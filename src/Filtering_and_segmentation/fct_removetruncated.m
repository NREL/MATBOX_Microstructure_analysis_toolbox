function [M,newtype,foo] = fct_removetruncated(M,p)

newtype = 'same';
foo = 1;

sz = size(M);
dimension = length(sz);

if p.background_is_edges
    x0 = unique(M(1,:,:));
    x1 = unique(M(end,:,:));
    y0 = unique(M(:,1,:));
    y1 = unique(M(:,end,:));
    z0 = unique(M(:,:,1));
    z1 = unique(M(:,:,end));
    attheborders = unique([x0;x1;y0;y1;z0;z1]);
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
end

if ~isempty(attheborders)
    for k=1:1:length(attheborders)
        M(M==attheborders(k))=p.complimentarylabel;
    end
end




[M] = fct_intconvert(M);

end