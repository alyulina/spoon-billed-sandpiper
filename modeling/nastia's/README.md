#### Final versions of the scripts used to compute spoon-billed sandpiper's and red-necked stint's density of polymorphisms and substitutions:

`pols.m` and `pols-h.m`\
heterozygosity is defined as the expected value of _2x(1-x)_ given the Wright's distribution of allele frequency x and the truncated skew normal distribution of selection coefficients with relative fitness _2s_, _s_, and _0_ __or__ _s_, _hs_, and _0_ assigned to _AA_, _Aa_, and _aa_, respectively; uses Duffy transformation and global adaptive algorithm for numerical integration or quasi Monte Carlo method;

`subs.m` and `subs-h.m`\
the number of substitutions per site is defined using Kimura's fixation probability with relative fitness _2s_, _s_, and _0_ __or__ _s_, _hs_, and _0_ assigned to _AA_, _Aa_, and _aa_, respectively, and assuming _N_ = _Ne_; uses trapezoid rule for numerical integration;

`distr.m`\
calculates the proportion of deleterious mutations for the skew normal distributed selection coefficients on a grid of distribution paremeters;

`parser.ipynb`\
parses the output of `pols.m`, `subs.m`, and `distr.m` to pass it to `plots.nb`;

`plots.nb`\
plots plots

Note: see [papers](../../papers) for references.
