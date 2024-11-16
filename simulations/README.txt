In this directory, there are several models written in Rust that run individual based models. SBS folder contains the model that has a gamma distribution of selection coefficients, with lethal alleles are fully recessive and all other alleles in the distribution having a dominance coefficient of 0.5. The results of these simulations, summarized by the get_averages.pl Perl script, were used to create Figures 3e,f. All other folders contain individual models where for each run all of the alleles have the same selection and dominance coefficients - these simulations were used to create Figures 3a-d.

The folders *_bottleneck contain models where the population was not subjected to complex demographic changes except for a decline of the population at the end of the simulation, mimicking a genetic bottleneck.

.sh bash files were used to run the models on a cluster.

In the SBS folder the parameters.txt file contains the gamma distributions (reported in Figure 3b) that were used for the individual simulations to make Figures 3e,f. The results in files averages_demography.txt and averages_bottleneck.txt report the differences in the fitness of the invididuals in the population before population size changes and after, with averages_demography.txt reporting the difference in fitness after complex demogrpahic changes and the file averages_bottleneck.txt reporting the difference in fitness after a simple bottleneck.

To run the models use the Rust cargo environment with "cargo build --release" to make the runables and "cargo run" to run the model, or use the .sh scripts.
