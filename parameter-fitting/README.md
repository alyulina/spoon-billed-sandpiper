This folder contains the code used to generate figures in `../figures`, as well as the `fitted_s_parameters.txt` in `../data`.

To determine which parameter combinations explain the levels of heterozygosity observed in the spoon-billed sandpiper population, we first simulate heterozygosity $H(s,h)$ across a range of fitness costs $s$ and dominance coefficients $h$ (see `../simulations`). In our simulations, we consider a freely recombining panmictic diploid population with mutation rate $\mu = 1.2 \cdot 10^-8$ per site per generation, whose size dynamics are described in the manuscript. We then assume that a fraction $p$ of fitness costs of new mutations is drawn from a gamma distribution with a shape parameter $\alpha$ and a rate parameter $\beta$, such that $p(s) \propto s^{\alpha-1}e^{-\beta s}$, and that these mutations are codominant (i.e. with the dominance coefficient $h=1/2$). We further assume that the rest of new mutations are strongly deleterious ($s=10^{-2}$) and recessive ($h=0$). On a grid of parameters, we then compute the heterozygosity numerically as  
```math
H=(1-p)\cdot \int_{0}^{1}dsH(s,1/2)p(s)+p\cdot H(10^{-2},0),
```
by approximating the integral as a Riemann sum with step size $ds$ on a logarithmic scale; we use linear interpolation to find the parameter values that would match the observed spoon-billed-sandpiper heterozygosity. 
