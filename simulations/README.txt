In this directory, there are several models written in Rust that run individual based models. SBS folder contains the model that has a gamma distribution of selection coefficients, with lethal alleles are fully recessive and all other alleles in the distribution having a dominance coefficient of 0.5. The results of these simulations, summarized by the get_averages.pl Perl script, were used to create Figures 3e,f. All other folders contain individual models where for each run all of the alleles have the same selection and dominance coefficients - these simulations were used to create Figures 3a-d.

The SBS_by* folders contain Rust simulations that model SBS-like demography on a population where all sites have the same selection coefficent. The results of these models were used to crate Figures 3a-d, and Supplementary Figure 7a-d. We ran these simulations 30 times, and the script getAVEfromfiles.pl is a Perl script that was used to average the output of multiple runs that takes a text file as input of the list of files to be averaged.

The folders *_bottleneck contain models where the population was not subjected to complex demographic changes except for a decline of the population at the end of the simulation, mimicking a genetic bottleneck.

.sh bash files were used to run the models on a cluster.

In the SBS folder the parameters.txt file contains the gamma distributions (reported in Figure 3b and Supplementary Figure 7b) that were used for the individual simulations to make Figures 3e-h (Supplementary Figure 7e-h). The output of this simulation is averaged across 100 runs of the simulations by the get_average.pl Perl script. The results in files averages_demography.txt and averages_bottleneck.txt report the differences in the fitness of the invididuals in the population before population size changes and after, with averages_demography.txt reporting the difference in fitness after complex demogrpahic changes and the file averages_bottleneck.txt reporting the difference in fitness after a simple bottleneck. The files averages_bottleneck_inbreeding.txt and averages_demography_inbreeding.txt report the differences in fitness between individuals created as an inbreeding of full siblings (averages_demography_inbreeding.txt file will report the difference in fitness of inbred individuals after SBS-like demography and at equilibrium and the averages_bottleneck_inbreeding.txt reports the difference in fitness of inbred individuals after a bottleneck and at equilibrium). To run the scimulation download the SBS folder and do:

>cargo build --release

if you run the the Slurm environment you can run:
>sbatch --array=1-190 SBS_run.sh

Otherwise, run individual parameters from paramters.txt file, for example:

>/target/release/SBS 0.0     300.0000000000001       2.293806075758451e-07

To average the fitness values across multiple runs you can make a list of the files:

>ls SBS_fitness_before_0.* > list.txt

you can then make output files (averages_demography.txt, averages_bottleneck.txt, averages_bottleneck_inbreeding.txt and averages_demography_inbreeding.txt) by:

>perl get_average.pl list.txt

The Rust environment should be compatible with all systems, our environment was a UNIX system and UNIX emulator on Windows 11 (Cygwin). To run the simulation copy a folder that includes the .toml file and the src folder. Compite and download the dependebles with cargo build --release command. Run a specific simulation using the cargo run --release command. The parameters of the gamma distribution used to run the individual based files are found in the parameters.txt file in the SBS folder.
