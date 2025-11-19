function [Deff,Mc,Tau,Bruggeman,eps] = Call_TauFactor2_binary(M,direction,label)
% [Deff,Mc,Tau,Bruggeman,eps] = Call_TauFactor2_binary(M,direction,label)
% Inputs:
%   Mandatory inputs:
%   - M: array with integer >=0
%   - direction: 1, 2 or 3
%   Optional input:
%   - label: integer >=0
%     if not provided, homogenization is performed on domain M==1
%     if provided, homogenization is performed on domain M==label
% Outputs:
%  Deff is calculated with an homogenization calculation performed by TauFactor2
%  Kench et al. (2023). TauFactor 2: A GPU accelerated python tool for microstructural analysis. Journal of Open Source Software, 8(88), 5358.
%  https://doi.org/10.21105/joss.05358. 
%  eps is the volume fraction of the domain homogenized
%  Tau and Bruggeman are deduced from: Deff/1 = eps/tau = eps^Bruggeman
%  Mc is MacMullin number

nargin; % Number of input variable when the function is call
if nargin==2
    label = 1;
end

% TauFactor 2 solves along first axis
if direction == 2
    p.action = 'Swap axis 1 with axis 2';
    [M,~,~] = fct_flipswap(M,p);
elseif direction == 3
    p.action = 'Swap axis 1 with axis 3';
    [M,~,~] = fct_flipswap(M,p);
end

% Binary volume
sz = size(M);
BW = zeros(sz);
BW(M==label) = 1;
clear M

% Volume fraction
eps = sum(sum(sum(BW))) / numel(BW);

if eps>0
    try
        BWpython = py.numpy.array(BW);
        res = pyrunfile("TauFactor2_binary.py", "res", img=BWpython, verbose=false); % Call TauFactor2
        Deff = double(res(1));
    catch ME
        if (strcmp(ME.identifier,'MATLAB:Python:PyException'))
            warning('Python was not loaded correctly, please restart MATLAB.')
        end
        Deff = NaN;
    end

    if isnan(Deff)
        Tau = NaN;
        Bruggeman = NaN;
        Mc = NaN;
    else
        if Deff<=0
            Deff = 0;
            Tau = NaN;
            Bruggeman = NaN;
            Mc = NaN;
        else
            Tau = double(res(2));
            Bruggeman = 1 - log(Tau)/log(eps);
            Mc = 1/Deff;
        end
    end

else
    Deff = 0;
    Tau = NaN;
    Bruggeman = NaN;
    Mc = NaN;
end

end