function [M,newtype,foo] = fct_removeinclusion(M,p)

newtype = 'same';
foo = 1;

sz = size(M);
dimension = length(sz);

unis = unique(M);
unis(unis==p.complementarylabel) = [];
if ~p.background_is_edges
    unis(unis==p.backgroundlabel) = [];
end

% I could use imfill, but that would require to loop on clusters, and get
% id of each cluster which is very expensive

BW = ones(sz);
BW(M==p.complementarylabel) = 0;
if ~p.background_is_edges
    BW(M==p.backgroundlabel) = 0;
end

labels = [unis ones(length(unis),1)]; % By default, all to be reassigned

conv_uni = zeros(max(unis),1);
conv_uni(1:unis(1)) = 1;
for k=2:1:length(unis)
    x0 = unis(k-1)+1;
    x1 = unis(k);
    conv_uni(x0:x1) = ones(unis(k)-unis(k-1),1)*k;
end


% BW2 = BW;
% BW2(1,:,:) = 0;
% BW2(end,:,:) = 0;
% BW2(:,1,:) = 0;
% BW2(:,end,:) = 0;
% if dimension == 3
%     BW2(:,:,1) = 0;
%     BW2(:,:,end) = 0;
% end
% idx = find(BW2);
% if dimension == 2
%     [IX,IY] = ind2sub(sz,idx);
% else
%     [IX,IY, IZ] = ind2sub(sz,idx);
% end
% 
% for k=1:1:length(idx)
%     current_label = M(idx(k));
%     if dimension == 2
%         label_x0 = M(IX(k)-1,IY(k));
%         label_x1 = M(IX(k)+1,IY(k));
%         label_y0 = M(IX(k),IY(k)-1);
%         label_y1 = M(IX(k),IY(k)+1);
%         surrounding_labels = [label_x0 label_x1 label_y0 label_y1];
%     else
%         label_x0 = M(IX(k)-1,IY(k),IZ(k));
%         label_x1 = M(IX(k)+1,IY(k),IZ(k));
%         label_y0 = M(IX(k),IY(k)-1,IZ(k));
%         label_y1 = M(IX(k),IY(k)+1,IZ(k));
%         label_z0 = M(IX(k),IY(k),IZ(k)-1);
%         label_z1 = M(IX(k),IY(k),IZ(k)+1);
%         surrounding_labels = [label_x0 label_x1 label_y0 label_y1 label_z0 label_z1];
%     end
%     if sum(surrounding_labels==p.complementarylabel)
%         labels(conv_uni(current_label),2) =0; % No reassignment for this one
%     end
% end


idx = find(BW);
if dimension == 2
    [IX,IY] = ind2sub(sz,idx);
else
    [IX,IY, IZ] = ind2sub(sz,idx);
end

n = length(idx);
x0s = max([IX-1 ones(n,1)],[],2);
x1s = min([IX+1 ones(n,1)*sz(1)],[],2);
y0s = max([IY-1 ones(n,1)],[],2);
y1s = min([IY+1 ones(n,1)*sz(2)],[],2);
if dimension == 3
    z0s = max([IZ-1 ones(n,1)],[],2);
    z1s = min([IZ+1 ones(n,1)*sz(3)],[],2);
end

for k=1:1:length(idx)
    current_label = M(idx(k));
    if dimension == 2
        label_x0 = M(x0s(k),IY(k));
        label_x1 = M(x1s(k),IY(k));
        label_y0 = M(IX(k),y0s(k));
        label_y1 = M(IX(k),y1s(k));
        surrounding_labels = [label_x0 label_x1 label_y0 label_y1];
    else
        label_x0 = M(x0s(k),IY(k),IZ(k));
        label_x1 = M(x1s(k),IY(k),IZ(k));
        label_y0 = M(IX(k),y0s(k),IZ(k));
        label_y1 = M(IX(k),y1s(k),IZ(k));
        label_z0 = M(IX(k),IY(k),z0s(k));
        label_z1 = M(IX(k),IY(k),z1s(k));
        surrounding_labels = [label_x0 label_x1 label_y0 label_y1 label_z0 label_z1];
    end
    if sum(surrounding_labels==p.complementarylabel)
        labels(conv_uni(current_label),2) =0; % No reassignment for this one
    end
end

id_tobereassinged = find(labels(:,2));
labels_reassign = labels(id_tobereassinged,1);

if p.Assigntosurrundinglabel
    BW3 = zeros(sz);
    for k=1:1:length(idx)
        if sum(M(idx(k))==labels_reassign)
            BW3(idx(k))=1;
        end
    end

    % Assign missing voxels to nearest cluster
    id_missing = find(BW3);
    BW(id_missing)=0;
    [~, IDX] = bwdist(BW);
    M(id_missing) = M(IDX(id_missing));
else
    for k=1:1:length(idx)
        if sum(M(idx(k))==labels_reassign)
            M(idx(k))=p.Assigntolabel;
        end
    end
end

[M] = fct_intconvert(M);

end