function [x,y] = Find_allindex_2Dline(x0,x1,y0,y1)

% x is always increasing
% y can be increasing or decreasing

nx = x1-x0+1;
if y1>y0
    ny = y1-y0+1;
    delta = 1;
else
    ny = y0-y1+1;
    delta = -1;
end
if nx==ny
    x = x0:1:x1;
    y = y0:delta:y1;
elseif nx<ny
    y = y0:delta:y1;
    a = (y1-y0)/(x1-x0);
    b = y1-a*x1;
    x = round((y-b)./a);
else
    x = x0:1:x1;
    a = (y1-y0)/(x1-x0);
    b = y1-a*x1;
    y = round(a.*x+b);
end

if max(abs(diff(x)))>1 || max(abs(diff(y)))>1
    keyboard
end

end