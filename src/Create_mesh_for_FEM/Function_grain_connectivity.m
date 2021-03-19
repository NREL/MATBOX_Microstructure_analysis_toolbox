function [new_M] = Function_grain_connectivity(M,background_code,connectivity)

grains = unique(M);
grains(grains==background_code)=[];
n_grain = length(grains);
sz = size(M);
new_M = zeros(sz) + background_code;
nn_grain = 0;
for k=1:1:n_grain
    grain = grains(k);
    BW = zeros(sz);
    idx = find(M==grain);
    BW(idx)=1;
    L=bwlabeln(BW,connectivity);
    uniqueL = unique(L);
    uniqueL(uniqueL==0)=[]; % Remove background
    n = length(uniqueL);    
    if n>1 % Grain is discontinous
        for kk=1:1:n
            idx = find(L==uniqueL(kk));
            nn_grain=nn_grain+1;
            new_M(idx) = nn_grain;
        end
    else
        nn_grain = nn_grain+1;
        new_M(idx) = nn_grain;
    end
end

end