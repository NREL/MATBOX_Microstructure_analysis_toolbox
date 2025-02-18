function [Ms,newtype,foo] = fct_calculateconnectivity_simplified(M,p)

newtype = 'Segmented (phase)';
foo = 1;

sz = size(M);
dimension = length(sz);

if p.backgroundisdomainedges
    Ms = zeros(sz);
else
    Ms = M;
end
   

% Id
MC = 1; % Main cluster
IC = 2; % Isolated cluster(s)
UC = 3; % Unknown cluster(s)
if ~p.backgroundisdomainedges
    MC = MC+1;
    IC = IC+1;
    UC = UC+1;
end

% All clusters and their size
unis = double(unique(M));
count = histc(reshape(M,[numel(M) 1]),unis);
id0 = find(unis==0); % Remove background
unis(id0,:)=[];
count(id0,:)=[];
if ~p.backgroundisdomainedges
    idc = find(unis==1); % Remove background
    unis(idc,:)=[];
    count(idc,:)=[];
end

res = [unis count];

% Identify main cluster
res = sortrows(res,2,"descend");
Ms(M==res(1,1))=MC;

% Discriminate between isolated and unknown clusters
[n,~] = size(res);
if n>=2
    for k=2:1:n
        idx = find(M==res(k,1));
        if dimension == 3
            [IX,IY,IZ] = ind2sub(sz,idx);
            zmin = min(IZ); zmax = max(IZ);
        else
            [IX,IY] = ind2sub(sz,idx);
        end
        xmin = min(IX); xmax = max(IX);
        ymin = min(IY); ymax = max(IY);

        ICcase = false;
        if dimension == 3
            if xmin==1 || ymin==1 || zmin==1 || xmax==sz(1) || ymax==sz(2) || zmax==sz(3)
                Ms(M==res(k,1))=UC;
            else
                Ms(M==res(k,1))=IC; ICcase=true;
            end
        else
            if xmin==1 || ymin==1 || xmax==sz(1) || ymax==sz(2)
                Ms(M==res(k,1))=UC;
            else
                Ms(M==res(k,1))=IC; ICcase=true;
            end
        end

        if ~p.backgroundisdomainedges && ICcase  % Some isolated clusters may be actually unknown clusters
            xxmin = max([xmin-1 1]); xxmax = min([xmax+1 sz(1)]);
            yymin = max([ymin-1 1]); yymax = min([ymax+1 sz(2)]);
            if dimension == 3
                zzmin = max([zmin-1 1]); zzmax = min([zmax+1 sz(3)]);
            else
                zzmin = 1; zzmax = 1;
            end
            sub = M(xxmin:xxmax,yymin:yymax,zzmin:zzmax);
            unis_sub = unique(sub);
            if max(unis_sub==0) % Possible contact
                BW = zeros(size(sub));
                BW(sub==res(k,1))=1;
                dmap=bwdist(BW,'chessboard');
                idxsub = find(dmap<=1.01);
                allids = unique(sub(idxsub));
                if max(allids==0)
                    Ms(idx)=UC; % Actually unknown
                end
            end

        end

    end
end

[Ms] = fct_intconvert(Ms);

end