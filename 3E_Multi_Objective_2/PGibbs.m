function Pg_sw=PGibbs(T_C ,w)

% Calculating derivative of seaawater Gibbs function ( ∂g_sw/ ∂w )
% Calculating derivatives of seawater enthalpy ( ∂h_sw/∂w) and entropy
% ( ∂s_sw/ ∂w)
% ∂g_sw/ ∂w = ∂h_sw/∂w - (T+273.15).∂s_sw/ ∂w

% Inputs:
%               T_C : Temperature                [degree Celsius]
%               w   : Salinity                   [ kg salt / kg solution]

a1=-2.348e4;
a2=3.152e5;
a3=2.803e6;
a4=-1.446e7;
a5=7.826e3;
a6=-4.417e1;
a7=2.139e-1;
a8=-1.991e4;
a9=2.778e4;
a10=9.728e1;

 Ph_sw=  -( a1 + 2*a2*w + 3*a3*w^2 + 4*a4*w^3 + a5*T_C + a6*T_C^2 ...
         + a7*T_C^3 + 2*a8*w*T_C + 3*a9*w^2*T_C + 2*a10*w*T_C^2)/1000;
  
  
b1=-4.231e2;
b2=1.463e4;
b3=-9.88e4;
b4=3.095e5;
b5=2.562e1;
b6=-1.443e-1;
b7=5.879e-4;
b8=-6.111e1;
b9=8.041e1;
b10=3.035e-1;


Ps_sw=  -( b1 + 2*b2*w + 3*b3*w^2 + 4*b4*w^3 + b5*T_C + b6*T_C^2 ...
        + b7*T_C^3 + 2*b8*w*T_C + 3*b9*w^2*T_C + 2*b10*w*T_C^2)/1000;


  
  Pg_sw=Ph_sw-(273.15+T_C)*Ps_sw;
 
end