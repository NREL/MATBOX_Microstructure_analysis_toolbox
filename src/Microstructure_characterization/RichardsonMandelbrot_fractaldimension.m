function [fractal_dimension, fractal_propensity, logs] = RichardsonMandelbrot_fractaldimension(topology_dimension,resolution,values)
% Fractal dimension according to Bertei et al., https://doi.org/10.1016/j.nanoen.2017.06.028 (Richardson/Mandelbrot formula)
% Log(property) = m + (topology_dimension-fractal dimension)*log(voxel size)

% Check for 0
id0 = find(resolution==0);
if ~isempty(id0)
    resolution(id0)=[];
    values(id0)=[];
end
id0 = find(values==0);
if ~isempty(id0)
    resolution(id0)=[];
    values(id0)=[];
end

logV = log(resolution);
logP = log(values);

logs = [reshape(logV,[numel(logV),1]) reshape(logP,[numel(logP),1])];

pf = polyfit(logV,logP,1);
fractal_dimension = topology_dimension-pf(1);
fractal_propensity = abs( topology_dimension - fractal_dimension );

end
