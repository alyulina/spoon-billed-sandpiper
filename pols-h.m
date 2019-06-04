Mathematica 11.1.0
Off[General::stop]

n1=50000;
n2=500000;
\[Nu]=1.2*10^-8;

a[s_,x_,h_]:=s*x^2+2*s*h*x*(1-x)
const[n_,s_,v_,u_,h_]:=Integrate[Exp[2*n*a[s,y,h]]*y^(4*n*v-1)*(1-y)^(4*n*u-1),{y,0,1}]
\[Phi][x_,n_,s_,v_,u_,h_]:=Exp[2*n*a[s,x,h]]*x^(4*n*v-1)*(1-x)^(4*n*u-1)
pols[\[Mu]_,\[Sigma]_,\[Alpha]_,n_,h_]:=NIntegrate[2*x*(1-x)*\[Phi][x,n,s,\[Nu],\[Nu],h]/const[n,s,\[Nu],\[Nu],h]*PDF[TruncatedDistribution[{\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},SkewNormalDistribution[\[Mu],\[Sigma],\[Alpha]]],s],{x,0,1},{s,\[Mu]-5*\[Sigma],\[Mu]+5*\[Sigma]},Method->"DuffyCoordinates"]

pols00=Table[pols[\[Mu],\[Sigma],\[Alpha],n2,0],{\[Alpha],-10,1,1.2},{\[Sigma],10^-8,2.5*10^-4,2.75*10^-5},{\[Mu],-1*10^-4,5*10^-5,1.65*10^-5}];
Export["pols_SkewNormal_500,000_0-DuffyCoordinates-10.csv",pols00]

pols05=Table[pols[\[Mu],\[Sigma],\[Alpha],n2,0],{\[Alpha],-10,1,1.2},{\[Sigma],10^-8,2.5*10^-4,2.75*10^-5},{\[Mu],-1*10^-4,5*10^-5,1.65*10^-5}];
Export["pols_SkewNormal_500,000_0.5-DuffyCoordinates-10.csv",pols05]
