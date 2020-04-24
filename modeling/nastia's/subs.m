Mathematica 11.1.0
Off[General::stop]

n1=50000;
n2=500000;
\[Nu]=1.2*10^-8;

subs[\[Mu]_,\[Sigma]_,\[Alpha]_,n_]:=200000*2*n*\[Nu]*NIntegrate[(1-Exp[-2*s])*PDF[TruncatedDistribution[{\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},SkewNormalDistribution[\[Mu],\[Sigma],\[Alpha]]],s]/(1-Exp[-4*n*s]),{s,\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},Method->"Trapezoidal"]

subs1=Table[subs[Chop[\[Mu]],\[Sigma],\[Alpha],n1],{\[Alpha],-10,1,0.56},{\[Sigma],10^-8,2.5*10^-4,1.3*10^-5},{\[Mu],-1*10^-4,5*10^-5,0.76*10^-5}];
Export["subs_SkewNormal_50,000-20.csv",subs1]

subs2=Table[subs[Chop[\[Mu]],\[Sigma],\[Alpha],n2],{\[Alpha],-10,1,0.56},{\[Sigma],10^-8,2.5*10^-4,1.3*10^-5},{\[Mu],-1*10^-4,5*10^-5,0.76*10^-5}];
Export["subs_SkewNormal_500,000-20.csv",subs2]
