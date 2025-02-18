function [P, C] = function_cumulativefct_chordlength_sphere(R,res)
% Calculate chord length probability and cumulative function for a sphere of radius R
%   R: sphere radius
%   res: mesh along sphere diameter 

% Calculate probability
L=linspace(0,2*R,res);
dL=L(2)-L(1);
Lm=L+dL;
Lm(end)=[];
P=zeros(res-1,2);
P(:,1) = Lm;
for k=1:1:res-1
    La = L(k);
    Lb = L(k+1);
    Ra = sqrt(R^2 - (La/2)^2);
    Rb = sqrt(R^2 - (Lb/2)^2);
    P(k,2) = abs((pi*Rb^2 - pi*Ra^2)*(La+Lb)/2 / (4/3*pi*R^3));    
end

% Deduce cumulative distribution function
C=zeros(res,2);
C(:,1) = L;
for k=1:1:res
    idx = find(Lm>L(k));
    if ~isempty(idx)
        C(k,2) = sum(P(idx,2));
    end
end

end