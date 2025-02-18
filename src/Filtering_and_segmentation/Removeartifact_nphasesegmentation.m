function [BW] = Removeartifact_nphasesegmentation(BW,width_threshold,illustrate)
% Remove thin layer of witdh <=width_threshold for phase i, between phase i-1 and phase i+1

sz = size(BW);
labels = unique(BW);
numberphase = length(labels);
idx_artifact = [];
% Loop over intermediate phase
if numberphase>2
    for k=1:1:numberphase
        tmp = zeros(sz);
        tmp(BW==labels(k))=1;
        dmaps(k).array = bwdist(tmp);
    end
    for current=2:1:numberphase-1
        darker = current-1;
        brighter = current+1;
        dmapdarker_tmp = dmaps(darker).array;
        dmapbrighter_tmp = dmaps(brighter).array;
        idcurrent = find(BW==labels(current));

        if illustrate
            if current==2 % Only for illustration
                dmapdarker = NaN(sz);
                dmapdarker(idcurrent) = dmapdarker_tmp(idcurrent);
                dmapbrighter = NaN(sz);
                dmapbrighter(idcurrent) = dmapbrighter_tmp(idcurrent);
                cmap = turbo; cmap(1,:) = [1 1 1];
                figure;
                imagesc(dmapdarker(:,:,1)); axis equal; colormap(cmap); title ('First intermediate phase. Distance map: to darker phase');
                figure;
                imagesc(dmapbrighter(:,:,1)); axis equal; colormap(cmap); title ('First intermediate phase. Distance map: to brigther phase');
            end
        end

        dmapdarker = zeros(sz)+1e9; % Very high value
        dmapdarker(idcurrent) = dmapdarker_tmp(idcurrent);
        dmapbrighter = zeros(sz)+1e9; % Very high value
        dmapbrighter(idcurrent) = dmapbrighter_tmp(idcurrent); 
        %cond1 = dmapdarker<=width_threshold;
        %cond2 = dmapbrighter<=width_threshold;
        %idtmp = find( (cond1.*cond2) == 1);
        idtmp = find( (dmapdarker+dmapbrighter) <= width_threshold);
        if ~isempty(idtmp)
            idx_artifact = [idx_artifact; idtmp];
        end
    end
end
artifacts = zeros(sz);
artifacts(idx_artifact) = 1;
if illustrate
    figure;
    imagesc(artifacts(:,:,1)); axis equal; colormap gray; title ('Detected artifacts');
end

% Remove artifacts
[~, idx_nearest] = bwdist(~artifacts);
BW(idx_artifact) = BW(idx_nearest(idx_artifact));

if illustrate
    figure;
    imagesc(BW(:,:,1)); axis equal; colormap gray; title ('Artifacts removed');
end

end