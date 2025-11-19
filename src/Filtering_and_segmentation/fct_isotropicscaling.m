function [M,newtype,foo] = fct_isotropicscaling(M,p)

sz = size(M);
foo=[];
newtype = 'same';

if strcmp(p.datatype,'Grey level') || strcmp(p.datatype,'Channel')
    p.label_or_greylevel = 'Grey level';
else
    p.label_or_greylevel = 'Label';
end

p.scaling_factor=1/p.scaling_factor; % I have to harmonize this will all the different scaling in this toolbox...

[M] = function_scaling(M,p);

[M] = fct_intconvert(M);

end