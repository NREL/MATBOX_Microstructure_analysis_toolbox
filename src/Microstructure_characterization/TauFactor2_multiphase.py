import taufactor as tau

# Call TauFactor2
# Kench et al. (2023). TauFactor 2: A GPU accelerated python tool for microstructural analysis. Journal of Open Source Software, 8(88), 5358.
# https://doi.org/10.21105/joss.05358.

cond = {}
for i in range(len(labels)):
    cond[labels[i]] = dbulk[i]

s = tau.MultiPhaseSolver(img, cond=cond, device='cuda')
s.solve(verbose=verbose)

Deff = s.D_eff.item()
tortuosity = s.tau.item()

res = [Deff, tortuosity]



