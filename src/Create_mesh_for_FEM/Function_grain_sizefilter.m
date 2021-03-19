function [Ma] = Function_grain_sizefilter(M,background_code,grain_critical_size)

backgdound_location = find(M==background_code);
grains = unique(M);
grains(grains==background_code)=[];
n_grain = length(grains);
sz = size(M);
for k=1:1:n_grain
    grain = grains(k);
    idx = find(M==grain);
    size_grain = length(idx);
    if size_grain<grain_critical_size
        M(idx) = background_code;
    end
end
% Assign background voxel to the nearest grain
BW = zeros(sz);
BW(M~=background_code)=1;
[~,idx] = bwdist(BW,'Euclidean');
Ma=M;
for k=1:1:numel(idx)
    Ma(k)=Ma(idx(k));
end
Ma(backgdound_location)=background_code;

end