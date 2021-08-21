A local version of this repo is in `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/github-repo`.

**To do:**

(right now)
- [ ] take points near the border of 1% of s larger than 1/4Ne (where there are some positive s but not more than 1% is above the 1/4Ne treshold); discard everything with abs(mu) or sigma > 0.01; take Rai's point + add Fedya's if needed, redo interpolation in that triangle (we do not believe in other distributions anyway), and then run Fedya's code on the line where the interpolated values match Mateusz' results (Fedya talked to Mateusz and he said that there are better methods for calculating dens/het and that he would use them and send us new numbers); update: looked at points where there are positive s, but no more than 1% larger than the cut-off – these are the same points; results of running Fedya's code on points near the border are in `/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver22/pos` (pos stands for positive); have not finished the rest yet
- [x] maybe also look at the sh<1/4Ne treshold, which is equal to 1/2Ne for the typical values of h of 1/2? and compare it with the 1/4Ne case; did that - the points where the results do not match are the same, so it does not really matter for the calculations (but might generally make more sense since most h are 1/2?); eventually started looking at points with more than 5% of the probability density for the distribution of s was positive
- [x] run Fedya's code for the same points as Rai where it does not match Rai's results and use this to do interpolation + infer points where distributions match Mateusz' data; for that, need to check the % of s>1/4Ne for each point and took those above a treshold that Fedya set (what to do with points that are very slow? can I discard them and run interpolation not on a square grid in Mathematica/python?)
- [x] run the latest version of Fedya's code (fixed: sigmoid instead of exp + beta is scaled inside now) on Rai's points and compare it with the 0 to 1 Rai's frequency version – this is in `/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver22/rai`; write to Fedya once done – hopefully it all works everywhere except for where there are beneficial allelel because the theory can't describe them
- [x] once the jobs are finished: make plots of mu-sigma grids to compare the results of Rai's and Fedya's codes; also compute errors from fluctuations in Fedya's code; the path to files on the cluster is `/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver21/rai` and locally `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver21/Rai's` see `parsing again.ipynb`
- [x] do what Fedya asked to do analytically
- [x] wait for Rai's reply to Fedya's email

(eventually)
- [ ] have a proper description here and in other READMEs
- [ ] upload final code / plots
- [ ] write down paths to files on the cluster

For the paper:
- [ ] make a repo for all the code (with Rai and possibly Fedya)
- [ ] rewrite Mathematica code? maybe even have it in python?
- [ ] start writing methods 
- [ ] have all lit references (what we know about the shape of the dfe in theory + h(s) dependancy + any studies that inferred that from data)


**Current state of things:**

ver18, ver19: not scaling sigma in Fedya's code, but now it does not work for large abs mu and small sigma because we have a cut-off on the distribution of s; change it to something else sinse s is not limited by -1 and 1 (maybe we should start taking exp(s))?

ver20: taking exp(s) instead of s for fitness + no limits on the distr of s; still looks bad (Rai's results do not match Fedya's) – what's wrong?

trying to understand why Fedya's results do not match Rai's? either something is wrong with interpolation / me parsing things (running stuff to see if this is the case) or with Rai's code / math (currently focusing on this) or with Fedya's code 

Fedya: maybe we should be scaling sigma after all

ver21: yes!!

ver22: running the latest version of Fedya's code (fixed: sigmoid instead of exp + beta is scaled inside now) on Rai's points to compare it with the 0 to 1 Rai's frequency version (data in ver14?) – hopefully it all works everywhere except for where there are beneficial allelel because the theory can't describe them (is this because there is no equilibrium freq?)

talked to Rai: sigma is time in the heat equation so this dependence makes sense (same het with diff mu but same sigma)

issues: Rai's neutral het is very different from Hexp when only including frequencies up to 1-1/2N (which makes sense but is potentially a problem);
theory does not match experiments when there are many positive s in the distribution (theory only works for s < 1/2N?) – this is ok, but as of now they also do not match when beta > 0 (should be fixed in ver22) and when het is very low (which also makes sense – need to try to estimate errors)

talked to Fedya: decided to discard everything with abs(mu) or sigma > 0.01 and where there are too many positive s; take Rai's points + add Fedya's if needed, redo interpolation in that triangle (we do not believe in other distributions anyway), and then run Fedya's code on the line where the interpolated values match Mateusz' results (Fedya talked to Mateusz and he said that there are better methods for calculating dens/het and that he would use them and send us new numbers)

used Mathematica for every point on Rai's grid to find how much of the probability density of s is above some treshhold (tried 0, 1/2Ne, and 1/4Ne, where Ne = 500,000 is the red-necked stint's Ne) - this is in `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver22/area.nb` - and looked at how this looks like in python - `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver22/Rai's/notebook-ver22.ipynb` - chose points where there are some positive but not more than 0.05 of the density is above 0 + three points below that (the treshhold turned out to not be very important and maybe it might be a good idea to just switch to 1/2Ne because this is the drift barrier for h = 1/2) - for viz see `Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver22/Rai's/to_run.png` and for points see `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver22/pos/RNS_input.txt`

ran Fedya's code on the points above (the cluster folder for that is `/nfs/scistore08/kondrgrp/alyulina/sandpiper/dynamics/h=sigm/ver22/pos`); am looking at the results now (in `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver22/pos/out/`) - need to add a section to `/Users/alyulina/Projects/Kondrashov/Sandpiper/h=sigm/ver22/Rai's/notebook-ver22.ipynb`

talked to Fedya and he said that the results are okay and that we should merge them and interpolate and run the simulation for the sandpiper

decided to run some more of Fedya's code to include abs(mu) and sigma of 0.01; once it's done: need to concatenate Fedya's and Rai's (saving as `to_use_rai.txt`) results 


**How to do things:**
- generate points (`make points.ipynb`) + an input file (`RNS_input.txt`) with mu, sigma, alpha, beta, and the output file name and location (`make shs.ipynb`)
- make a batch script (`run_RNS.sh`) to run Fedya's code (`g++ -o RNSv8 RNSv8.cpp`) on the cluster + make a folder for output (`./rns`) and slurm output (`./rns/outs`)
- use a python notebook on the cluster (`notebook-ver22.ipynb`) to calculate averages across the ten runs (`rns-hets-points_rns_avg-het.txt`, `rns-hets-points_rns_pols.txt`, and split into separate files with different alpha and beta in `./out`); alternatively, if jupyterhub is down, convert it to a pythin script and run that 
- download data on the local machine
- run a Mathematica notebook to make plots (`plots-ver22-rns.nb`)
