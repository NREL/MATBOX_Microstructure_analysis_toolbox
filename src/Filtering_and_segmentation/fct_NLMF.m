function [M,newtype,estDos] = fct_NLMF(M,p)
sz = size(M);
dimension = length(sz);
newtype = 'same';

if dimension==2
    sz = [sz 1];
end

for k=1:1:sz(3)
    slice_ = M(:,:,k);
    if p.estimatedegreeofsmoothing
        [M(:,:,k), estDos] = imnlmfilt(slice_,'SearchWindowSize',p.searchwindowssize,'ComparisonWindowSize',p.comparizonwindowssize);
        %estDos % Estimated degree of smoothing
    else
        estDos = [];
        M(:,:,k) = imnlmfilt(slice_,'SearchWindowSize',p.searchwindowssize,'ComparisonWindowSize',p.comparizonwindowssize,'DegreeOfSmoothing',p.degreeofsmoothing);
    end
end

[M] = fct_intconvert(M);

end