function [Deff,Mc,Tau,Bruggeman] = Call_TauFactor2_multiphase(Dbulk,direction,eps)
% [Deff,Tau,Bruggeman] = Call_TauFactor2_multiphase(Dbulk,direction,eps)
% Inputs:
%   Mandatory inputs:
%   - Dbulk: array with values ranging from 0 to 1
%   - direction: 1, 2 or 3
%   Optional input:
%   - eps: total porosity
% Outputs:
%  Deff is calculated with an homogenization calculation performed by TauFactor2
%  Kench et al. (2023). TauFactor 2: A GPU accelerated python tool for microstructural analysis. Journal of Open Source Software, 8(88), 5358.
%  https://doi.org/10.21105/joss.05358. 
%  - if eps is provided, Tau and Bruggeman are deduced from: Deff/1 = eps/tau = eps^Bruggeman
%  - if eps is not provided, Tau is deduced using Dmean = mean(Dbulk), Tau = Dmean/Deff
%    that is, how much we deviate from a rule of mixture (upper bound for Deff)
%    Note that this defintion of Tau is not the tortuosity factor.
%    It is then recommended to provide eps if known.

% There is an infinity of Dbulk value between 0 and 1
% However, TauFactor2 has a limit for the number of different labels.
% So we round Dbulk: e.g.: 0.374678... becomes 0.3747, assigned to label 3747
nargin; % Number of input variable when the function is call
if nargin==2
    eps = []; % We do not know the true porosity. We will use the tortuosity definition provided by TauFactor 2 in this case
end

prec = 4;
scale = 10^prec;

labels_taufactor2 = 1:1:scale;
Dbulk_taufactor2 = labels_taufactor2/scale; % [0,1] 

BW_Taufactor2 = round(Dbulk*scale);

labels_nonzero = labels_taufactor2*0;
labels_nonzero(BW_Taufactor2(BW_Taufactor2~=0))=1;

labels_taufactor2=labels_taufactor2(labels_nonzero==1);
Dbulk_taufactor2=Dbulk_taufactor2(labels_nonzero==1);

% TauFactor 2 solves along first axis
if direction == 2
    p.action = 'Swap axis 1 with axis 2';
    [BW_Taufactor2,~,~] = fct_flipswap(BW_Taufactor2,p);
elseif direction == 3
    p.action = 'Swap axis 1 with axis 3';
    [BW_Taufactor2,~,~] = fct_flipswap(BW_Taufactor2,p);
end

try
    BWpython = py.numpy.array(BW_Taufactor2);
    Deff = pyrunfile("TauFactor2_multiphase.py", "Deff", img=BWpython, labels=labels_taufactor2, dbulk=Dbulk_taufactor2, verbose=false);  % Call TauFactor2
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
        if ~isempty(eps)
            Bruggeman = log(Deff)/log(eps);
            Tau = eps^(1-Bruggeman);
        else
            % We do not know the true porosity. We will use the tortuosity definition provided by TauFactor 2 in this case
            Dmean = mean(mean(mean( Dbulk  ))); % Weight average
            Tau = Dmean/Deff;
            Bruggeman = NaN;
        end
        Mc = 1/Deff;
    end
end

end