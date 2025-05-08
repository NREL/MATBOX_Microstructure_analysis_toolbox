function [Ms,newtype,foo] = fct_checksurroundings_percluster(M,p)

p.choice = 'Assign to a phase';
p.label = 2; % SSE
p.assigntolabel = 1; % Defect
p.threshold = 1.0;
p.ignore = [];

foo=[];
newtype = 'same';

sz = size(M);
dimension = length(sz);

Ms = M;

BW = zeros(sz);
BW(M==p.label)=1;
if dimension == 2
    C = bwlabel(BW,4);
else
    C = bwlabeln(BW,6);
end

unis = unique(C);
unis(unis==0)=[];
n = length(unis);
for k=1:n
    idx = find(C==unis(k));
    if dimension == 2
        [IX,IY]=ind2sub(sz,idx);
        z0=1; z1=1;
    else
        [IX,IY,IZ]=ind2sub(sz,idx);
        z0 = max([1, min(IZ)-1]);
        z1 = min([sz(3), max(IZ)+1]);
    end
    x0 = max([1, min(IX)-1]);
    y0 = max([1, min(IY)-1]);
    x1 = min([sz(1), max(IX)+1]);
    y1 = min([sz(2), max(IY)+1]);

    subBW = C(x0:x1,y0:y1,z0:z1);
    subBW(subBW~=unis(k))=0;
    sub = M(x0:x1,y0:y1,z0:z1);

    dmap = bwdist(subBW,'chessboard');
    adjacent_labels = sub(dmap==1);

    [counts, groupnames] = groupcounts(adjacent_labels);
    tmp = [double(groupnames), double(counts)];

    if ~isempty(p.ignore)
        for kk=1:length(p.ignore)
            id0 = find(tmp(:,1)==p.ignore(kk));
            if ~isempty(id0)
                tmp(id0,:)=[];
            end
        end
    end

    tmp = sortrows(tmp,2,'descend');

    if ~isempty(tmp)

        tot = sum(tmp(:,2));
        if strcmp(p.choice,'Assign to most surrounded phase')
            if tmp(1,2) >= p.threshold * tot
                Ms(idx) = tmp(1,1);
            end
        elseif strcmp(p.choice,'Assign to a phase')
            idlabel = find(tmp(:,1)==p.assigntolabel);
            if ~isempty(idlabel)
                if tmp(idlabel,2) >= p.threshold * tot
                    Ms(idx) = p.assigntolabel;
                end
            end
        end

    end


end

[Ms] = fct_intconvert(Ms);

end