function [M,achieved_volume_fraction,achieved_depth] = function_generate_coating(M,additive_id,deposit_sens,deposit_phase,target_volume_fraction,target_depth,target_choice)

BW = zeros(size(M));
BW(M==deposit_phase)=1;
depositphase_volumefraction = sum(sum(sum(BW)))/numel(BW);
complementary_volumefraction = 1-depositphase_volumefraction;

if strcmp(deposit_sens,'From particle surface to background')
    dmap = bwdist(BW); r=depositphase_volumefraction;
elseif strcmp(deposit_sens,'From particle surface to particle interior')
    dmap = bwdist(~BW); r=complementary_volumefraction;
end

dall = unique(dmap);
dall(dall==0)=[];
for k=1:1:length(dall)
    coating_depth = dall(k);
    if strcmp(target_choice,'Volume fraction')
        coating_vol = (sum(sum(sum(dmap <= dall(k) ))) / numel(BW))-r;
        if coating_vol>=target_volume_fraction
            BW(dmap <= dall(k))=2;
            break
        end
    elseif strcmp(target_choice,'Coating depth in voxel')
        if coating_depth>=target_depth
            BW(dmap <= dall(k))=2;
            break
        end
    end
end

% Overwritte
if strcmp(deposit_sens,'From particle surface to background')
    BW(M == deposit_phase)=1;
elseif strcmp(deposit_sens,'From particle surface to particle interior')
    BW(M ~= deposit_phase)=0;
end
M(BW==2)=additive_id;

achieved_volume_fraction = (sum(sum(sum(M == additive_id )))) / numel(M);
achieved_depth = double(coating_depth);

end