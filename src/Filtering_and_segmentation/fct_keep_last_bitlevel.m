function [S5to8,newtype,unis] = fct_keep_last_bitlevel(M,p)
sz = size(M);
dimension = length(sz);
newtype = 'same';

%
M = double(M);
M = round(M*255 / max(max(max(M))));
M = uint8(M);
%

unis(1) = length(unique(M));

% Get bit level
if dimension == 2
    M8 = zeros([sz(1) sz(2) 8],"uint8");
    for k=1:1:8
        M8(:,:,k) = bitget(M,k,"uint8");
    end
else
    M8 = zeros([sz(1) sz(2) sz(3) 8],"uint8");
    for k=1:1:8
        M8(:,:,:,k) = bitget(M,k,'uint8');
    end
end
if p.show_bitlevel
    figure;
    for i=1:8
        subplot(2,4,i)
        if dimension == 2
            strt_title = 'Bit levels';
            imagesc(M8(:,:,i)); axis equal; axis tight; colormap gray;
        else
            imagesc(M8(:,:,round(sz(3)/2),i)); axis equal; axis tight; colormap gray;
            strt_title = 'Bit levels (middle slice)';
        end
        title(['Bit level: ' num2str(i)])
    end
    sgtitle(strt_title);
end

% Keep only last bits (contains most of the image information, removing noise)
S5to8 = zeros(sz,"uint8");
for k=8:-1:8-p.lastbit_to_keep+1
    if dimension==2
        S5to8 = S5to8 + M8(:,:,k)*2^(k-1);
    else
        S5to8 = S5to8 + M8(:,:,:,k)*2^(k-1);
    end
end
S5to8 = uint8(S5to8);

if p.show_bitlevel
    if dimension==2
        figure; imagesc(S5to8); axis equal; axis tight; colormap gray;
    else
        figure; imagesc(S5to8(:,:,round(sz(3)/2))); axis equal; axis tight; colormap gray;
    end
    title(['Grey level with only the ' num2str(p.lastbit_to_keep,'%i') ' last bits'])
end
unis(2) = length(unique(S5to8));

end