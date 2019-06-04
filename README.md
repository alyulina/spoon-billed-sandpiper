#### Final versions of the scripts used to compute spoon-billed sandpiper's and red-necked stint's density of polymorphisms and substitutions:

`pols.m` and `pols-h.m`\
heterozygosity is defined as the expected value of 2x(1-x) given the Wright's distribution of allele frequency x and the truncated skew normal distribution of selection coefficients with relative fitness 2s, s, and 0 or s, hs, and 0 assigned to AA, Aa, and aa, respectively;

`subs.m` and `subs-h.m`\
the number of substitutions per site is defined using Kimura's fixation probability with relative fitness 2s, s, and 0 or s, hs, and 0 assigned to AA, Aa, and aa, respectively, and assuming N = Ne;

`distr.m`\
calculates the proportion of deleterious mutations for the skew normal distributed selection coefficients on a grid of distribution paremeters;

`parser.ipynb`\
parses the output of `pols.m`, `subs.m`, and `distr.m` to pass it to `plots.nb`;

`plots.nb`\
plots plots

