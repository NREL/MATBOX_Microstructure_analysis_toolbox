import taufactor as tau
from taufactor.utils import flux_direction
import tifffile
import torch

# Call TauFactor2
# Kench et al. (2023). TauFactor 2: A GPU accelerated python tool for microstructural analysis. Journal of Open Source Software, 8(88), 5358.
# https://doi.org/10.21105/joss.05358.


img = tifffile.imread('C:/Users/fussegli/Documents/GitHub/Test_TauFactor2/Test_TauFactor2.tif')
#flux_direction(img)
#img = torch.permute(torch.tensor(img), (1,2,0))
#img = torch.permute(torch.tensor(img), (2,0,1))


#s = tau.Solver(img)
#s = tau.Solver(img, device='cuda')

s = tau.PeriodicSolver(img)


#s.solve(verbose=verbose,conv_crit=conv_crit_std)

#s.solve(verbose=verbose,conv_crit=1e-6)


s.solve(verbose=True)




if s.D_eff==0:
    Deff = s.D_eff
    tortuosity = s.tau
else:
    Deff = s.D_eff.item()
    tortuosity = s.tau.item()

res = [Deff, tortuosity]