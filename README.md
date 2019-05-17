#### Final versions of the scripts used to compute spoon-billed sandpiper's and red-necked stint's density of polymorphisms and substitutions:

`pols.m`\
heterozygosity is defined as the expected value of 2x(1-x) given the Wright's distrion of allele frequency x and the truncated skew normal distribution of selection coefficients with relative fitness 2s, s, and 0 assigned to AA, Aa, and aa, respectively;

`subs.m`\
the number of substitutions per site is defined using Kimura's fixation probability with relative fitness 2s, s, and 0 assigned to AA, Aa, and aa, respectively;

`parser.ipynb`\
just a notebook to parse the output of `pols.m` and `subs.m` to pass it to `plots.m`;

`plots.m`
