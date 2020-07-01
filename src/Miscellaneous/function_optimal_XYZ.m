function [optimal_threshold_1,optimal_threshold_2,optimal_threshold] = function_optimal_XYZ(x,y,z)
%UNTITLED7 Summary of this function goes here

unique_x = unique(x);
n_1=length(unique_x);
optimal_threshold_1=zeros(n_1,3);
optimal_threshold_1(:,1)=unique_x;
for k=1:1:n_1
    idx = find(x==unique_x(k));
    iddxx = find(z(idx)==max(z(idx)));
    if ~isempty(iddxx)
        optimal_threshold_1(k,2) = y(idx(iddxx(1)));
        optimal_threshold_1(k,3) = z(idx(iddxx(1)));
    else
        optimal_threshold_1(k,2)=-1;
    end
end
idx = find(optimal_threshold_1(:,2)==-1);
if ~isempty(idx)
    optimal_threshold_1(idx,:)=[];
end
unique_y = unique(y);
n_2=length(unique_y);
optimal_threshold_2=zeros(n_2,3);
optimal_threshold_2(:,1)=unique_y;
for k=1:1:n_2
    idx = find(y==unique_y(k));
    iddxx = find(z(idx)==max(z(idx)));
    if ~isempty(iddxx)
        optimal_threshold_2(k,2) = x(idx(iddxx(1)));
        optimal_threshold_2(k,3) = z(idx(iddxx(1)));
    else
        optimal_threshold_2(k,2)=-1;
    end
end
idx = find(optimal_threshold_2(:,2)==-1);
if ~isempty(idx)
    optimal_threshold_2(idx,:)=[];
end
iter=0;
[n_1,~]=size(optimal_threshold_1);
for k=1:1:n_1
    x=optimal_threshold_1(k,1);
    y1 = optimal_threshold_1(k,2);
    z1 = optimal_threshold_1(k,3);
    idx = find(optimal_threshold_2(:,2)==x);
    if ~isempty(idx)
        iter=iter+1;
        y2 = optimal_threshold_2(idx(1),1);
        z2 = optimal_threshold_2(idx(1),3);
        if z1>=z2
            optimal_threshold(iter,1)=x;
            optimal_threshold(iter,2)=y1;
            optimal_threshold(iter,3)=z1;
        else
            optimal_threshold(iter,1)=x;
            optimal_threshold(iter,2)=y2;
            optimal_threshold(iter,3)=z2;
        end
    end
end
[n_2,~]=size(optimal_threshold_2);
for k=1:1:n_2
    y=optimal_threshold_2(k,1);
    x2 = optimal_threshold_2(k,2);
    z2 = optimal_threshold_2(k,3);
    idx = find(optimal_threshold_1(:,2)==y);
    if ~isempty(idx)
        iter=iter+1;
        x1 = optimal_threshold_1(idx(1),1);
        z1 = optimal_threshold_1(idx(1),3);
        if z1>=z2
            optimal_threshold(iter,1)=x1;
            optimal_threshold(iter,2)=y;
            optimal_threshold(iter,3)=z1;
        else
            optimal_threshold(iter,1)=x2;
            optimal_threshold(iter,2)=y;
            optimal_threshold(iter,3)=z2;
        end
    end
end
optimal_threshold=unique(optimal_threshold,'rows');
[~,idx] = sort(optimal_threshold(:,1)); % sort just the first colum
optimal_threshold = optimal_threshold(idx,:);   % sort the whole matrix using the sort indices

end

