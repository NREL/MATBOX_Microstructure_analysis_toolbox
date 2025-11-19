function [pore_Deff,pore_Mc,pore_Tau,pore_p,pore_eps,n_voxel] = Call_TauFactor2_pore(Nanoporosity,Wetting,Bulkdiffusivity,direction)

idpore = find(Nanoporosity.*Wetting>0); % Check binary or dual scale case
n_voxel = length(idpore);
pore_eps = sum(Nanoporosity(idpore).*Wetting(idpore))/numel(Nanoporosity);
Ds = unique(Bulkdiffusivity(idpore));
cond0 = Ds==0;
cond1 = Ds==1;
condD01 = cond0+cond1; % OR

if ~isempty(idpore)
    if min(condD01)==0 % Multi scale D=[0,1]
        [pore_Deff,pore_Mc,pore_Tau,pore_p] = Call_TauFactor2_multiphase(Bulkdiffusivity,direction,pore_eps);
    else % Single scale D=0 or D=1
        [pore_Deff,pore_Mc,pore_Tau,pore_p,~] = Call_TauFactor2_binary(Bulkdiffusivity,direction);
    end
else
    pore_Deff = 0;
    pore_Tau = NaN;
    pore_p = NaN;
    pore_Mc = NaN;
end

end