function [Ms,newtype,foo] = fct_calculatesize(M,p)

newtype = 'Channel';
foo = 1;

sz = size(M);
dimension = length(sz);

% Same
% M_1D = reshape(M,[numel(M) 1]);
% [counts, groupnames] = groupcounts(M_1D);
% res = [counts double(groupnames)];
% idx = find(groupnames==0);
% res(idx,:)=[]; % Remove background
%
% Ms = zeros(sz);
% [n,~] = size(res);
% for k = 1:1:n
%     Ms(M==res(k,2))=res(k,1);
% end


% % Slightly faster
% Ms = zeros(sz);
% labels = unique(M);
% labels(labels==0)=[];
% for k = 1:1:length(labels)
%     idx = M==labels(k);
%     Ms(idx)=sum(sum(sum(idx)));
%     %Ms(idx)=length(idx);
% end

M = M+1; % no zero for regionprops

if dimension == 2
    idx_tmp = regionprops(M,"PixelIdxList");
    res = regionprops("table",M,"Area");
    res = res.Area;
else
    idx_tmp = regionprops3(M,"VoxelIdxList");
    res = regionprops3(M,"Volume");
    res = res.Volume;
end

Ms = zeros(sz,'single');

if strcmp(p.choice,'One label')
    labels = p.label + 1;
else
    labels = unique(M);
    if strcmp(p.choice,'All labels except some labels')
        excluded_labels = str2num(p.excluded_labels);
        if ~isempty(excluded_labels)
            excluded_labels = excluded_labels + 1;
            for k=1:1:length(excluded_labels)
                labels(labels==excluded_labels(k))=[];
            end
        end
    end
end

n = length(labels);
if dimension == 2
    for k = 1:1:n
        idx = idx_tmp(labels(k)).PixelIdxList;
        Ms(idx) = res(labels(k));
    end
else
    for k = 1:1:n
        idx = idx_tmp.VoxelIdxList{labels(k),1};
        Ms(idx) = res(labels(k));
    end
end

Ms = round(Ms);
[Ms] = fct_intconvert(Ms);

end


