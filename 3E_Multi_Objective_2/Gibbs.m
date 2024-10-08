function [g_sw ]=Gibbs(h_sw, s_sw, T_C)


% Calculating Gibbs Function
% Inputs:
%               T_C  : Temperature                  %[degree Celsius]
%               h_sw : seawater enthalpy            %[kj/kg]
%               s_sw : seawater entropy             %[kj/kg.K]


g_sw=h_sw-(T_C+273.15)*s_sw;                        %[kj/kg]
end
