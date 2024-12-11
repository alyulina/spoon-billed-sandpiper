This folder contains the code used to generate _Figure 3_ in the main text, as well as _Supplemental Figure 7_.

To determine which parameter combinations explain the levels of heterozygosity observed in the spoon-billed sandpiper population, we first simulate heterozygosity $H(s,h)$ across a range of fitness costs $s$ and dominance coefficients $h$ (see `../simulations`). In our simulations, we consider a freely recombining panmictic diploid population with mutation rate $\mu = 1.2 \cdot 10^-8$ per site per generation, whose size dynamics are described in the manuscript. We then assume that a fraction $p$ of fitness costs of new mutations is drawn from a gamma distribution with a shape parameter $\alpha$ and a rate parameter $\beta$, such that $p(s) \propto s^{\alpha-1}e^{-\beta s}$, and that these mutations are codominant (i.e. with the dominance coefficient $h=1/2$). We further assume that the rest of new mutations are strongly deleterious ($s=10^{-2}$) and recessive ($h=0$). On a grid of parameters, we then compute the heterozygosity numerically as  
```math
H=(1-p)\cdot \int_{0}^{1}H(s,1/2)p(s)ds+p\cdot H(10^{-2},0),
```
by approximating the integral as a Riemann sum with step size $ds$ on a logarithmic scale; we use linear interpolation to find the parameter values that would match the observed spoon-billed-sandpiper heterozygosity. 

Running the code requires the web-based interactive computational environment `jupyter notebook`, `python 3.8.3` or above, as well as the following libraries: `numpy 1.23.1`, `scipy 1.10.1`, `statsmodels 0.13.5`, `matplotlib 3.7.1`, `seaborn 0.12.2`, which can be installed on any operation system. For installation instructions, see [python.org/downloads](https://www.python.org/downloads/) and [jupyter.org/install](https://jupyter.org/install). To launch the notebook, simply type `jupyter notebook fitting_selection_coefficients.ipynb` from the command line. Running the code in this notebook should take about half an hour. 

**Required input data** (which can be obtained by running the code in `../simulations`):  
`../data/H_s_recessive.txt` and `../data/H_s_codominant.txt`, which contain the simulated heterozygocities under the inferred demography for $h=0$ and $h=1/2$, respectively;  
`../data/H_s_recessive_bottleneck.txt` and `../data/H_s_codominant_bottleneck.txt`, which contain the simulated heterozygocities under the simpler demography with just a bottleneck for $h=0$ and $h=1/2$, respectively;  
`../data/fitness_diff.txt` and `fitness_diff_bottleneck.txt`, which contain the simulated genotype fitness values under the two demographies that we considered;  
`../data/inbred_fitness_diff.txt` and `inbred_fitness_diff_bottleneck.txt`, same as above but for inbred individuals.  

**Output data**:  
`data/fitted_s_parameters.txt` contains the model parameters that fit the observed levels of spoon-billed sandpiper heterozygosity;  
_Figure 3_ in the main text, `../figures/H_p=0.4.pdf`;  
_Supplemental Figure 7_, `../figures/H.pdf`.  
