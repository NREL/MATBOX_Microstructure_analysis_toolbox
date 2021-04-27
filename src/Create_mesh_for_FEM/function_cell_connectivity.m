function cellcluster = function_cell_connectivity(elem)

% Initialize
not_assigned = elem;

cluster_id = 0;
while ~isempty(not_assigned)
    current_cluster = not_assigned(1,:); % Initialize
    cluster_id = cluster_id+1;
    numel_old = 0; numel_new = numel(current_cluster);
    while numel_new>numel_old
        numel_old = numel(current_cluster);
        Lia = ismember(not_assigned,current_cluster);
        idx = find(sum(Lia,2)>=1);
        if ~isempty(idx)
            current_cluster = [current_cluster; not_assigned(idx,:)];
            not_assigned(idx,:) = [];
        end
        numel_new = numel(current_cluster);
        cellcluster(cluster_id).cells=current_cluster;
    end
end

for k=1:1:cluster_id
    Lia = ismember(elem,cellcluster(k).cells);
    idx = find(sum(Lia,2)==4);        
    cellcluster(k).index=idx;
end

end

