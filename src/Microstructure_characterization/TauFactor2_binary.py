import taufactor as tau

# Call TauFactor2
# Kench et al. (2023). TauFactor 2: A GPU accelerated python tool for microstructural analysis. Journal of Open Source Software, 8(88), 5358.
# https://doi.org/10.21105/joss.05358.

#s = tau.Solver(img)
s = tau.Solver(img, device='cuda')
s.solve(verbose=verbose,conv_crit=conv_crit_std)

if s.D_eff==0:
    Deff = s.D_eff
    tortuosity = s.tau
else:
    Deff = s.D_eff.item()
    tortuosity = s.tau.item()

res = [Deff, tortuosity]