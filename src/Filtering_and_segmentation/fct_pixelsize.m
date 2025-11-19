function [M2,newtype,foo] = fct_pixelsize(M,p)

sz = size(M);
foo = [];
newtype = 'Channel';
dimension = length(sz);

if strcmp(p.choice,'One label')
    labels = p.label;
else
    labels = unique(M);
    if strcmp(p.choice,'All labels except some labels')
        excluded_labels = str2num(p.excluded_labels);
        for k=1:1:length(excluded_labels)
            labels(labels==excluded_labels(k))=[];
        end
    end
end

n = length(labels);
M2 = zeros(sz);
for k=1:1:n
    idx = M==labels(k);
    BW = zeros(sz);
    BW(idx)=1;
    if strcmp(p.metric,'c-PSD')
        [Mtmp,~,~] = Function_particle_size_CPSD_Algorithm(BW,p.round,p.round_dmap_digit);
    elseif strcmp(p.metric,'EDM')
        [Mtmp, ~] = bwdist(~BW,p.distance);
        if p.round
            Mtmp = round(Mtmp,p.round_dmap_digit);
        end
    end
    if strcmp(p.choice,'One label')
        M2 = Mtmp;
    else
        M2(idx) = Mtmp(idx);
    end
end

if p.round
    M2 = round(M2);
    [M2] = fct_intconvert(M2);
end

end
