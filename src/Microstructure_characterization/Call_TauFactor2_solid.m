function [solid_Deff,solid_Mc,solid_Tau,solid_p,solid_eps,n_voxel] = Call_TauFactor2_solid(Nanoporosity,Bulkconductivity,direction)

idsolid = find(Nanoporosity<1);  % Check binary or dual scale case
n_voxel = length(idsolid);
solid_eps = sum( 1-Nanoporosity(idsolid) )/numel(Nanoporosity);
Ks = unique(Bulkconductivity(idsolid));
cond0 = Ks==0;
cond1 = Ks==1;
condK01 = cond0+cond1; % OR

if ~isempty(idsolid)
    if min(condK01)==0 % Multi scale D=[0,1]
        [solid_Deff,solid_Mc,solid_Tau,solid_p] = Call_TauFactor2_multiphase(Bulkconductivity,direction,solid_eps);
    else % Single scale D=0 or D=1
        [solid_Deff,solid_Mc,solid_Tau,solid_p,~] = Call_TauFactor2_binary(Bulkconductivity,direction);
    end
else
    solid_Deff = 0;
    solid_Tau = NaN;
    solid_p = NaN;
    solid_Mc = NaN;
end

end