function [x0,x1] = find_extremities(idx,x,n)

if ~isempty(idx)
    idx_before = idx;
    idx_before(idx > x) = [];
    if ~isempty(idx_before)
        x0 = max(idx_before)+1;
    else
        x0 = 1;
    end
    idx_after = idx;
    idx_after(idx < x) = [];
    if ~isempty(idx_after)
        x1 = min(idx_after)-1;
    else
        x1 = n;
    end
else
    x0 = 1;
    x1 = n;
end

end