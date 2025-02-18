function [angles] = Charact_CrackParticle_orientation(Semantic, idcrack, Instance, centroids, connectivity)

sz = size(Semantic);
dimension = length(sz);

idx = find(Semantic == idcrack);
n = length(idx);
if dimension == 2
    [IX,IY]=ind2sub(sz,idx);
    xyzC = [IX, IY];
else
    [IX,IY,IZ]=ind2sub(sz,idx);
    xyzC = [IX, IY, IZ];
end

angles = zeros(n*10*dimension,3);

i=0;
for k=1:1:n

    idA = Instance(idx(k));
    xyzA = centroids(idA,:);

    conn = connectivity(idA).adjacentinstances;
    if ~isempty(conn)
        for kb = 1:1:length(conn)
            idB = conn(kb);
            xyzB = centroids(idB,:);

            u = xyzA - xyzB;
            v = xyzA - xyzC(k,:);

            if sum(abs(v)) > 0
                % https://www.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab

                % Method 1
                % CosTheta = max(min(real(dot(u,v))/(norm(u)*norm(v)),1),-1);
                % theta = rad2deg(acos(CosTheta));

                %Method 2
                nu = norm(u);
                nv = norm(v);
                uu = u*nv;
                vv = v*nu;
                a = uu-vv;
                b = uu+vv;
                theta = rad2deg(2*atan(sqrt(sum(a.^2)/sum(b.^2))));
                theta = min([abs(theta-180),theta]);

                i = i + 1;
                angles(i,1) = idA; % Particle A id
                angles(i,2) = idB; % Partilce B id
                angles(i,3) = k;   % Crack id
                angles(i,4) = theta; % Angle
            end
         
        end
    end
end

angles(i+1:end,:)=[];

% For each crack id that belong to A, find B that minimize the angle
angles = sortrows(angles,[3 4]);
unis = unique(angles(:,3));
[~,idx] = ismember(unis,angles(:,3));
angles = angles(idx,:);
angles = sortrows(angles,[1 2 3 4]);



%%
%angles = sortrows(angles,[1 2 3 4]);
%d = double(logical(diff(angles(:,1)))) + double(logical(diff(angles(:,2))));
%t = find(d>=1);
%t = [0; t; i];

% Plot distribtion per particle
% Plot distribution of all angles
% for k=1:1:length(t)-1
%     n0 = t(k)+1;
%     n1 = t(k+1);
%     a = angles(n0:n1,3);
% 
%     % Calculate probability density function
%     histogram(a)
% 
%     % Calculate distribution metric 
% end

end