function [M,newtype,foo] = fct_spheresizefilter(M,p)

sz = size(M);
foo=[];
newtype = 'same';

labels = unique(M);
BW = zeros(sz);
BW(M == p.labeltofilter)=1;

[diameter] = Function_particle_size_CPSD_Algorithm(BW,p.smooth_surface);

% cond1 = diameter <= 6;
% cond2 = BW;
% idx = find(cond1.*cond2 == 1);
% diameter(idx)=1e9;
% tmp = zeros(sz);
% tmp(idx)=1;
% Fig = figure; imagesc(tmp); axis equal; axis tight; colormap gray;

cmap = turbo(length(unique(diameter)));
cmap(1,:) = [0.5 0.5 0.5];
Fig = figure; imagesc(diameter); axis equal; axis tight; colormap(cmap);

M(diameter<=p.filter_size_diameter)=p.reassign_to;

end