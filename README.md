To do:

(right now)
- [ ] once the jobs are finished: make plots of mu-sigma grids to compare the results of Rai's and Fedya's codes; also compute errors from fluctuations in Fedya's code; the path to files on the cluster is `/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver20/rai` and locally `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver20/Rai's` see `parsing again.ipynb`
- [ ] do what Fedya asked to do analytically
- [ ] wait for Rai's reply to Fedya's email

(eventually)
- [ ] have a proper description here and in other READMEs
- [ ] upload final code / plots
- [ ] write down paths to files on the cluster

For the paper:
- [ ] make a repo for all the code (with Rai and possibly Fedya)
- [ ] start writing methods 
- [ ] rewrite Mathematica code?




**Current state of things:**

ver18, ver19 not scaling sigma in Fedya's code, but now it does not work for large abs mu and small sigma because we have a cut-off on the distribution of s; change it to something else sinse s is not limited by -1 and 1 (maybe we should start taking exp(s))?

ver20: taking exp(s) instead of s for fitness + no limits on the distr of s; still looks bad (Rai's results do not match Fedya's) – what's wrong?

trying to understand why Fedya's results do not match Rai's? either something is wrong with interpolation / me parsing things (running stuff to see if this is the case) or with Rai's code / math (currently focusing on this) or with Fedya's code 
