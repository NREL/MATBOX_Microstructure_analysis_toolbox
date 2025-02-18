function [x0, x1, y0, y1, maxarea] = fct_findlargestrectangle2D(BW)

sz = size(BW);
maxarea = 0;
hist = zeros(1,sz(2));
for r=1:1:sz(1)
    row = BW(r,:);
    hist = (row + hist);
    hist(row==0)=0;
    % Calculate largest rectangle in histogram
    % At least on bar is fully included in the largest rectangle in the histogram
    area = 0;
    for c=1:1:sz(2)
        width = 1;
        id0 = find(hist<hist(c));
        id0_before = id0;
        id0_before(id0_before>c)=[];
        id0_after = id0;
        id0_after(id0_after<c)=[];
        if isempty(id0_before)
            %width = width + c-1;
            ymin = 1;
        else
            %width = width + c-max(id0_before)-1;
            ymin = max(id0_before)+1;
        end
        if isempty(id0_after)
            %width = width + sz(2)-c;
            ymax = sz(2);
        else
            %width = width + min(id0_after)-c-1;
            ymax = min(id0_after)-1;
        end
        height = hist(c);
        width = ymax-ymin+1;
        current_area = width*height;
        if current_area>maxarea
            x0 = r-height+1;
            x1 = r;
            y0 = ymin;
            y1 = ymax;
        end
        area = max([area current_area]);
    end
    maxarea = max([maxarea area]);
end

end