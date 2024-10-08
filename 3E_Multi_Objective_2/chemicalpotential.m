function [mu_w, mu_S]=chemicalpotential(T0_C, w)

% Calculating   Chemical Potential (μ) for water and salt
%  μ_w = g_sw - w* ∂g_sw/ ∂w        : water chemical potential
%  μ_s = g_sw + (1-w)*∂g_sw/ ∂w     : Salt Chemical Potential

% Inputs:
%               T0 : Temperature         [degree Celsius]
%               w  : Salinity            [ kg salt / kg solution]


[h_sw,s_sw]=seawaterprop(T0_C,w);                    %[kj/kg , kj/kg.K]
g_sw=Gibbs(h_sw,s_sw,T0_C);                          %[kj/kg]
Pg_sw=PGibbs(T0_C,w);                                %[kj/kg]

mu_s= g_sw*1000 + (1000-w*1000)*Pg_sw;             %[j/kg]
mu_S=1000*w*mu_s;                                  %[j/kg] 
mu_w=g_sw*1000 - w*1000*Pg_sw;                     %[j/kg]

end