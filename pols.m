Mathematica 11.1.0
Off[General::stop]

n1=50000;
n2=500000;
\[Nu]=1.2*10^-8;

const[n_,s_,v_,u_]:=Integrate[Exp[4*n*s*y]*y^(4*n*v-1)*(1-y)^(4*n*u-1),{y,0,1}]
\[Phi][x_,n_,s_,v_,u_]:=Exp[4*n*s*x]*x^(4*n*v-1)*(1-x)^(4*n*u-1)
pols[\[Mu]_,\[Sigma]_,\[Alpha]_,n_]:=NIntegrate[2*x*(1-x)*\[Phi][x,n,s,\[Nu],\[Nu]]/const[n,s,\[Nu],\[Nu]]*PDF[TruncatedDistribution[{\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},SkewNormalDistribution[\[Mu],\[Sigma],\[Alpha]]],s],{x,0,1},{s,\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},Method->"QuasiMonteCarlo"]//Re

pols1=Table[pols[\[Mu],\[Sigma],\[Alpha],n1],{\[Alpha],-10,1,0.56},{\[Sigma],10^-8,2.5*10^-4,1.3*10^-5},{\[Mu],-1*10^-4,5*10^-5,0.76*10^-5}];
Export["pols_SkewNormal_50,000-QuasiMonteCarlo-20.csv",pols1]

pols2=Table[pols[\[Mu],\[Sigma],\[Alpha],n2],{\[Alpha],-10,1,0.56},{\[Sigma],10^-8,2.5*10^-4,1.3*10^-5},{\[Mu],-1*10^-4,5*10^-5,0.76*10^-5}];
Export["pols_SkewNormal_500,000-QuasiMonteCarlo-20.csv",pols2]
