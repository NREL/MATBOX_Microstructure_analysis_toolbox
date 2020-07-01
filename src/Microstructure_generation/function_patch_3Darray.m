function [f,v] = function_patch_3Darray(BW)
% Build faces f and vertices v array only for the voxels at the interface
% of the binary array BW

% Add empty layer at domain's extremity
domain_size = size(BW);
tmp=zeros(domain_size+2);
tmp(2:end-1,2:end-1,2:end-1)=BW;
BW=tmp;
clear tmp;
domain_size = size(BW);

% Detect interface along x
BW_tmp = zeros(domain_size);
BW_tmp(2:end,:,:) = BW(1:end-1,:,:);
sub = BW - BW_tmp; clear BW_tmp;
idx = find(sub==-1);
[x,y,z]=ind2sub(domain_size,idx); clear idx;
x=x-2; % Correct added empty layer: -1, correct -1 are located in the next voxel: -1
y=y-1; % Correct added empty layer: -1
z=z-1; % Correct added empty layer: -1
% Create vertices and facets
n_facets = length(x);
f_xplus = zeros(n_facets,4); % Square
v_xplus = zeros(n_facets*4,3); % 4 vertices per square
for n=1:1:n_facets
    f_xplus(n,:)= (n-1)*4+1 + [0,1,2,3];
    v_xplus( f_xplus(n,1),: ) = [x(n) y(n)-1 z(n)-1];
    v_xplus( f_xplus(n,2),: ) = [x(n) y(n)-1 z(n)];
    v_xplus( f_xplus(n,3),: ) = [x(n) y(n)   z(n)];
    v_xplus( f_xplus(n,4),: ) = [x(n) y(n)   z(n)-1];
end
idx = find(sub==1); clear sub;
[x,y,z]=ind2sub(domain_size,idx); clear idx;
x=x-2; % Correct added empty layer: -1
y=y-1; % Correct added empty layer: -1
z=z-1; % Correct added empty layer: -1
% Create vertices and facets
f_xminus = zeros(n_facets,4); % Square
v_xminus = zeros(n_facets*4,3); % 4 vertices per square
for n=1:1:n_facets
    f_xminus(n,:)= (n-1)*4+1 + [0,1,2,3];
    v_xminus( f_xminus(n,1),: ) = [x(n) y(n)-1 z(n)-1];
    v_xminus( f_xminus(n,2),: ) = [x(n) y(n)-1 z(n)];
    v_xminus( f_xminus(n,3),: ) = [x(n) y(n)   z(n)];
    v_xminus( f_xminus(n,4),: ) = [x(n) y(n)   z(n)-1];
end
n_ = length(v_xplus);
f_xminus = f_xminus + n_;

% Detect interface along y
BW_tmp = zeros(domain_size);
BW_tmp(:,2:end,:) = BW(:,1:end-1,:);
sub = BW - BW_tmp; clear BW_tmp;
idx = find(sub==-1);
[x,y,z]=ind2sub(domain_size,idx); clear idx;
x=x-1; % Correct added empty layer: -1
y=y-2; % Correct added empty layer: -1, correct -1 are located in the next voxel: -1
z=z-1; % Correct added empty layer: -1
% Create vertices and facets
n_facets = length(y);
f_yplus = zeros(n_facets,4); % Square
v_yplus = zeros(n_facets*4,3); % 4 vertices per square
for n=1:1:n_facets
    f_yplus(n,:)= (n-1)*4+1 + [0,1,2,3];
    v_yplus( f_yplus(n,1),: ) = [x(n)-1 y(n) z(n)-1];
    v_yplus( f_yplus(n,2),: ) = [x(n)-1 y(n) z(n)];
    v_yplus( f_yplus(n,3),: ) = [x(n)   y(n) z(n)];
    v_yplus( f_yplus(n,4),: ) = [x(n)   y(n) z(n)-1];
end
n_ = n_ + length(v_xminus);
f_yplus = f_yplus + n_;
idx = find(sub==1); clear sub;
[x,y,z]=ind2sub(domain_size,idx); clear idx;
x=x-1; % Correct added empty layer: -1
y=y-2; % Correct added empty layer: -1
z=z-1; % Correct added empty layer: -1
% Create vertices and facets
f_yminus = zeros(n_facets,4); % Square
v_yminus = zeros(n_facets*4,3); % 4 vertices per square
for n=1:1:n_facets
    f_yminus(n,:)= (n-1)*4+1 + [0,1,2,3];
    v_yminus( f_yminus(n,1),: ) = [x(n)-1 y(n) z(n)-1];
    v_yminus( f_yminus(n,2),: ) = [x(n)-1 y(n) z(n)];
    v_yminus( f_yminus(n,3),: ) = [x(n)   y(n) z(n)];
    v_yminus( f_yminus(n,4),: ) = [x(n)   y(n) z(n)-1];
end
n_ = n_ + length(v_yplus);
f_yminus = f_yminus + n_;

% Detect interface along z
BW_tmp = zeros(domain_size);
BW_tmp(:,:,2:end) = BW(:,:,1:end-1);
sub = BW - BW_tmp; clear BW_tmp;
idx = find(sub==-1);
[x,y,z]=ind2sub(domain_size,idx); clear idx;
x=x-1; % Correct added empty layer: -1
y=y-1; % Correct added empty layer: -1
z=z-2; % Correct added empty layer: -1, correct -1 are located in the next voxel: -1
% Create vertices and facets
n_facets = length(z);
f_zplus = zeros(n_facets,4); % Square
v_zplus = zeros(n_facets*4,3); % 4 vertices per square
for n=1:1:n_facets
    f_zplus(n,:)= (n-1)*4+1 + [0,1,2,3];
    v_zplus( f_zplus(n,1),: ) = [x(n)-1 y(n)-1 z(n)];
    v_zplus( f_zplus(n,2),: ) = [x(n)-1 y(n)   z(n)];
    v_zplus( f_zplus(n,3),: ) = [x(n)   y(n)   z(n)];
    v_zplus( f_zplus(n,4),: ) = [x(n)   y(n)-1 z(n)];
end
n_ = n_ + length(v_yminus);
f_zplus = f_zplus + n_;
idx = find(sub==1); clear sub;
[x,y,z]=ind2sub(domain_size,idx); clear idx;
x=x-1; % Correct added empty layer: -1
y=y-1; % Correct added empty layer: -1
z=z-2; % Correct added empty layer: -1, correct -1 are located in the next voxel: -1
% Create vertices and facets
f_zminus = zeros(n_facets,4); % Square
v_zminus = zeros(n_facets*4,3); % 4 vertices per square
for n=1:1:n_facets
    f_zminus(n,:)= (n-1)*4+1 + [0,1,2,3];
    v_zminus( f_zminus(n,1),: ) = [x(n)-1 y(n)-1 z(n)];
    v_zminus( f_zminus(n,2),: ) = [x(n)-1 y(n)   z(n)];
    v_zminus( f_zminus(n,3),: ) = [x(n)   y(n)   z(n)];
    v_zminus( f_zminus(n,4),: ) = [x(n)   y(n)-1 z(n)];
end
n_ = n_ + length(v_zplus);
f_zminus = f_zminus + n_;

% Concatenate all facets and vertices
v = [v_xplus; v_xminus; v_yplus; v_yminus; v_zplus; v_zminus];
f = [f_xplus; f_xminus; f_yplus; f_yminus; f_zplus; f_zminus];

% Remove doublons verteces
v_unique = unique(v,'rows');
f_tmp = zeros(size(f));
for k=1:1:length(v_unique)
    RowIdx = find(ismember(v, v_unique(k,:),'rows'));
    f_tmp(ismember(f,RowIdx))=k;
end

% Rename
v=v_unique; clear v_unique
f=f_tmp; clear f_tmp
end