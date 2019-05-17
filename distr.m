Mathematica 11.1.0
Off[General::stop]

n1=50000;
n2=500000;
\[Nu]=1.2*10^-8;

distr[\[Mu]_,\[Sigma]_,\[Alpha]_]:=NIntegrate[PDF[TruncatedDistribution[{\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},SkewNormalDistribution[\[Mu],\[Sigma],\[Alpha]]],s],{s,\[Mu]-5*\[Sigma],0}]

d0=Table[distr[Chop[\[Mu]],\[Sigma],\[Alpha]],{\[Alpha],-10,1,0.56},{\[Sigma],10^-8,2.5*10^-4,1.3*10^-5},{\[Mu],-1*10^-4,5*10^-5,0.76*10^-5}];

Export["distr-20.csv",d0]