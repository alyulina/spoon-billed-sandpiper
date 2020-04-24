Mathematica 11.1.0
Off[General::stop]

n1=50000;
n2=500000;
\[Nu]=1.2*10^-8;

Mx[s_,x_,h_]:=s*x*(1-x)*(h+(1-2*h)*x)
Vx[x_,n_]:=x*(1-x)/(2*n)
G[s_,n_,h_,x_]:=Exp[2*n*s*x*((2*h-1)*x-2*h)]
fix[s_,n_,h_]:=Integrate[G[s,n,h,x],{x,0,1/(2*n)}]/Integrate[G[s,n,h,x],{x,0,1}]
subs[\[Mu]_,\[Sigma]_,\[Alpha]_,n_,h_]:=200000*2*n*\[Nu]*NIntegrate[fix[s,n,h]*PDF[TruncatedDistribution[{\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},SkewNormalDistribution[\[Mu],\[Sigma],\[Alpha]]],s],{s,\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},Method->"Trapezoidal"]

subs00=Table[subs[\[Mu],\[Sigma],\[Alpha],n2,0],{\[Alpha],-10,1,1.2},{\[Sigma],10^-8,2.5*10^-4,2.75*10^-5},{\[Mu],-1*10^-4,5*10^-5,1.65*10^-5}];
Export["subs_SkewNormal_500,000_0-Trapezoidal-10.csv",subs00]

subs05=Table[subs[\[Mu],\[Sigma],\[Alpha],n2,0.5],{\[Alpha],-10,1,1.2},{\[Sigma],10^-8,2.5*10^-4,2.75*10^-5},{\[Mu],-1*10^-4,5*10^-5,1.65*10^-5}];
Export["subs_SkewNormal_500,000_0.5-Trapezoidal-10.csv",subs05]
