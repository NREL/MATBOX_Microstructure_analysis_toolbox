function [M,newtype,foo] = fct_sortpersize(M,p)

keyboard

uni = unique(M);
if p.keepunchangedsomelabels
    cst  = str2num(p.cst);
    n_cst = length(cst);
    for k=1:1:n_cst
        uni(uni==cst(k))=[]; % remove constant labels
    end
end
n = length(uni);

newtype = 'same';
foo = 1;

M_1D = reshape(M,[numel(M) 1]);
[counts, groupnames] = groupcounts(M_1D);
res = [counts double(groupnames)];
idx = find(groupnames==0);
res(idx,:)=[]; % Remove background

% Sort
if strcmp(p.order,'decreasing')
    res = sortrows(res,1,"descend");
elseif strcmp(p.order,'increasing')
    res = sortrows(res,1,"ascend");
end
[n,~] = size(res);

Msav = M;
for k=1:1:n
    M(Msav == res(k,2))=k;
end

[M] = fct_intconvert(M);

end