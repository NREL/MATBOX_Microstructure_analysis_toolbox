function [BWscaled] = function_stochastic_surfaceroughness(BW,ratio_outside_sr,ratio_inside_sr,niter,forceconnectivity)

for iter=1:1:niter
    if ratio_outside_sr>0
        dmap_0to1=bwdist(BW);
        cond1=dmap_0to1< (2-0.001);
        cond2=BW==0;
        id_0t01 = find(cond1.*cond2);
        n=length(id_0t01);
        p = randperm(n,round(ratio_outside_sr*n));
        BW(id_0t01(p))=1;
    end

    if ratio_inside_sr>0
        dmap_1to0=bwdist(~BW);
        cond1=dmap_0to1< (2-0.001);
        cond2=BW==0;
        id_1t00 = find(cond1.*cond2);
        n=length(id_1t00);
        p = randperm(n,round(ratio_outside_sr*n));
        BW(id_1t00(p))=0;
    end

    % Set parameters
    parameters_scaling.scaling_factor = 1/2;
    parameters_scaling.label_or_greylevel = 'Label';
    parameters_scaling.background = 0;
    % Scale
    BWscaled = function_scaling(BW,parameters_scaling);

    % Prepare next iteration
    BW=BWscaled;
end

if forceconnectivity
    C=bwlabeln(~BWscaled,6);
    [GC,GR] = groupcounts(reshape(C,[numel(C), 1]));
    tmp = [GC,GR];
    id0 = find(GR==0);
    tmp(id0,:)=[];
    tmp=sortrows(tmp,1,"descend");
    n_cluster = length(tmp(:,1));
    if n_cluster>1
        BWscaled = ones(size(BWscaled));
        BWscaled(C==tmp(1,2))=0;
    end

    C=bwlabeln(BWscaled,6);
    [GC,GR] = groupcounts(reshape(C,[numel(C), 1]));
    tmp = [GC,GR];
    id0 = find(GR==0);
    tmp(id0,:)=[];
    tmp=sortrows(tmp,1,"descend");
    n_cluster = length(tmp(:,1));
    if n_cluster>1
        BWscaled = zeros(size(BWscaled));
        BWscaled(C==tmp(1,2))=1;
    end
end

end