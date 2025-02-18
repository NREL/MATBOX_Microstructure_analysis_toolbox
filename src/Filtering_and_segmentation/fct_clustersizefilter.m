function [M,newtype,foo] = fct_clustersizefilter(M,p)

sz = size(M);
dimension = length(sz);
foo=[];
newtype = 'same';

labels = unique(M);
BW = zeros(sz);
BW(M == p.labeltofilter)=1;

if dimension==2
    C=bwlabel(BW,4);
else
    C=bwlabeln(BW,6);
end

if length(labels)==2
    otherlabel = labels;
    otherlabel(otherlabel==p.labeltofilter)=[];
end

filter_size = p.filter_size_length^dimension;
uni = unique(C);
uni(uni==0)=[]; % Remove background
n = length(uni); % Number of clusters
% Consider only small and isolated clusters
for k=1:1:n
    idx = find(C==uni(k));
    [IX,IY,IZ]=ind2sub(sz,idx);
    if dimension==2
        if (min(IX)>1 && min(IY)>1 && max(IX)<sz(1) && max(IY)<sz(2)) || p.include_edgeclusters
            if (length(IX)<=filter_size && strcmp(p.suporinf,'<=')) || (length(IX)>filter_size && strcmp(p.suporinf,'>'))
                if length(labels)==2
                    M(idx)=otherlabel;
                else
                    if strcmp(p.reassign_rule,'Assign to another label')
                        M(idx)=p.reassignlabel;
                    else
                        warning('TO DO: fct_clustersizefilter.m: rule to be added');
                        % Rule of assign to add
                    end
                end
            end
        end
    else
        if (min(IX)>1 && min(IY)>1 && min(IZ)>1 && max(IX)<sz(1) && max(IY)<sz(2) && max(IZ)<sz(3)) || p.include_edgeclusters
            if (length(IX)<=filter_size && strcmp(p.suporinf,'<=')) || (length(IX)>filter_size && strcmp(p.suporinf,'>'))
                if length(labels)==2
                    M(idx)=otherlabel;
                else
                    if strcmp(p.reassign_rule,'Assign to another label')
                        M(idx)=p.reassignlabel;
                    else
                        warning('TO DO: fct_clustersizefilter.m: rule to be added');
                        % Rule of assign to add
                    end
                end
            end
        end
    end
end


end