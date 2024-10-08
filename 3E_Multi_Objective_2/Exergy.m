function e=Exergy(T_C,w,T0_C,w0)



%Calculating total specific exergy of state
% InPut:
%           T_C  : state temperature               [degree Celsius] 
%           w    : Salinity                        [ kg salt / kg solution]
%           T0_C : Reference State Temperature     [degree Celsius]
%           w0   : Reference State Salinity        [ kg salt / kg solution]

% Enthalpy and Entropy calculation of desired state
[h_sw, s_sw]=seawaterprop(T_C,w);                 %[kj/kg , kj/kg.K]
h_sw=h_sw*1000;                                   %[j/kg]
s_sw=s_sw*1000;                                   %[j/kg.K]
% Enthalpy and Entropy calculation of Reference state
% μ*_w and  μ
[h_sw_star, s_sw_star] = seawaterprop(T0_C,w);   
h_sw_star=h_sw_star*1000;
s_sw_star=s_sw_star*1000;


% calculating Restricted Dead State Chemical Potential (μ*_w and  μ*_s)
[mu_w_star, Smu_s_star ] = chemicalpotential(T0_C,w);


% calculating Total Dead State Chemical Potential (μ*_w and  μ*_s)
[mu_w_0,S0mu_s_0 ]    = chemicalpotential(T0_C,w0);
Smu_s_0 = (S0mu_s_0/w0).*w;
w=w*1000;
     e = (h_sw - h_sw_star)-(T0_C+273.15).*(s_sw-s_sw_star)...
        + (1-0.001*w).*(mu_w_star-mu_w_0)+0.001*(Smu_s_star-Smu_s_0);
end