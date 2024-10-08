clc
clear 
close all


    


    fluid = 'water';
    T0=25+273;                                           %[k]
    T0_C=25;                                             %[C]
    P0=101;                                              %[kPa]
    w0=0.032;
    T=zeros(1,15);
    T_C=zeros(1,15);
    P=zeros(1,15);
    Q=zeros(1,15);
    m=zeros(1,15);
    w=zeros(1,15);
%% Data Input From Reference Article



T_Geo=200:5:300;
P2_Fig=linspace(1,1.5,15)*1000;
P5_Fig=[0.04,0.064,0.079,0.097,0.119,0.147,0.183,0.226,0.277,0.344,...	
0.427,0.529,0.655,0.809,1]*1000;
i_Fig=linspace(0.02,0.12,15);
n_Fig=linspace(20,35,15);


Fig={T_Geo, P2_Fig,P5_Fig,i_Fig,n_Fig};
X_axis_labels={'Geo Fluid Temperatue (C)'; 'Seperator Pressure (kPa)';...
    'Turbine Outlet Pressure (kPa)' ; 'Interest Rate %';...
    'System Lifetime (year)'};

T_C(1)=200;
P(2)=1300;                           %[kPa] expansion valve outlet
P(5)=100;                             %[kPa] Turbine outlet  pressure
int_rat=0.03;                                                  % Intrest rate 
n=25;                                                    %[year]

for q=1:5
    ExD_Fig=zeros(length(cell2mat(Fig(q))),6);
    Overall_Eff_Fig=zeros(length(cell2mat(Fig(q))),6);
    C_tot_Fig=zeros(length(cell2mat(Fig(q))),6);
 
for j=1:length(cell2mat(Fig(q)))
    
    
%% Inputs 
    % GEothermal Steam Turbine Input
    if q==1
        clear T_C(1)
        T_C(1)=T_Geo(j);                          %[C]   Geo-fluid temperature
    elseif q==2
        clear P(2)
        P(2)=P2_Fig(j);
        T_C(1)=200;
    elseif q==3
        clear P(5)
        P(5)=P5_Fig(j);
        T_C(1)=200;
        P(2)=1300;
    elseif q==4
        clear int_rat
        int_rat=i_Fig(j);
        T_C(1)=200;
        P(2)=1300;
        P(5)=100;
        
    elseif q==5
        clear n
        n=n_Fig(j);
        T_C(1)=200;
        P(2)=1300;
        P(5)=100;
        int_rat=0.03;
    end


    T(1)=T_C(1)+273;                     %[k] 
    m(1)=150;                            %[kg/s]Geo-fluid mass flow rate 
    Q(1)=0;                              %      Quality of production well
    eta_turbine=0.75;                    %      Turbine efficiency
    T_C(7)=45;                           %[C]   Heat exchanger outlet Temp.
    T(7)=T_C(7)+273;                     %[k]
    
    % MED Input
    T_C(8)=25;                           %[C]   Intake feed seawater Temp.
    T(8)=T_C(8)+273;                     %[k]
    P(8)=101;                            %[kPa] Intake feed seawater pres.
    w(8)=32/1000;                        %[g salt/g water] Intake Salinity
    w(9)=70/1000;                        %[g salt/g water] Brine Sal.
    w(10)=0;                             %[g salt/g water] Distillate Sal.
    T_effI_C=60;                         %[C]   First effect Temp.
    T_effI=T_effI_C+273;                 %[k]             
    T_diff=3.2;                          %[k]   Temp. difference between effects
    Effect_no=6;                                 %      Effects No.
    T_effn_C=(T_effI_C-(Effect_no-1)*T_diff);    %[C]   Last effect Temp. 
    T_ave=(T_effI_C+T_effn_C)/2;         %[C]   Average effect Temp.
    T_C(9)=T_ave;                        %[C]   Brine outlet Temp.
    T(9)=T_C(9)+273;                     %[k]
    T_C(10)=T_ave;                       %[C]   Distillate outlet Temp.
    T(10)=T_C(10)+273;                   %[k]


    % RO Input
    T_C(11)=25;                          %[C]   Intake Ro Feedwater Temp
    T(11)=T_C(11)+273;                   %[k]
    P(11)=101;                           %[kPa] Intake Pres.
    w(11)=32/1000;                       %[g salt/g water] Intake Salinity
    w(12)=70/1000;                       %[g salt/g water] Brine Sal.
    w(13)=0;                             %[g salt/g water] Distillate Sal.
    
         % Incompressible fluid assumption
    T(12)=T(11);                         %[k]   Brine outlet Temp.
    T(13)=T(11);                         %[k]   Distillate outlet Temp.
    T_C(12)=T_C(11);                     %[C]
    T_C(13)=T_C(11);                     %[C]
    RO_HP=5500;                          %[kPa] Pressure of RO module high pressure side 
    RO_LP=101;                           %[kPa] Pressure of RO module low pressure side 
    
    T_C(14)=25;                          %[C]   Cooling water Inlet Temp.
    T_C(15)=52;                          %[C]   Cooling water outlet Temp
    T(14)=T_C(14)+273;                   %[k]
    T(15)=T_C(15)+273;                   %[k]
    
%% Fisrt Law Analysis
h=zeros(1,15);
s=zeros(1,15);
 
  %Geothermal System 
  % 'state 1'
  h(1)= refpropm('H','T',T(1),'Q',Q(1),fluid)/1000;         %[kj/kg]
  P(1)= refpropm('P','T',T(1),'Q',Q(1),fluid);              %[kpa]
  s(1)= refpropm('S','T',T(1),'Q',Q(1),fluid)/1000;         %[kj/kg.k]
  
  % 'state 2'
  h(2)=h(1);
  T(2)= refpropm('T','P',P(2),'H',h(2)*1000,fluid);         %[k]
  T_C(2)=T(2)-273;                                          %[C]
  Q(2)= refpropm('Q','P',P(2),'H',h(2)*1000,fluid);
  s(2)= refpropm('S','P',P(2),'H',h(2)*1000,fluid)/1000;    %[kj/kg.k]
  m(2)=m(1); 
  
  % 'state 3'
  P(3)=P(2);
  T(3)=T(2);
  T_C(3)=T_C(2);
  Q(3)=0;
  h(3)= refpropm('H','P',P(3),'Q',Q(3),fluid)/1000;         %[kj/kg]
  s(3)= refpropm('S','P',P(3),'Q',Q(3),fluid)/1000;         %[kj/kg.k]
  m(3)=(1-Q(2))*m(2);                                       %[kg/s]
  
  % 'state 4'
  P(4)=P(2);
  T(4)=T(2);
  T_C(4)=T_C(2);
  Q(4)=1;
  h(4)= refpropm('H','P',P(4),'Q',Q(4),fluid)/1000;         %[kj/kg]
  s(4)= refpropm('S','P',P(4),'Q',Q(4),fluid)/1000;         %[kj/kg.k]
  m(4)=m(2)-m(3);                                           %[kg/s]
  
  % 'state 5'
  s5s=s(4)*1000;
  h5s= refpropm('H','P',P(5),'S',s5s,fluid)/1000;          %[kj/kg]
  h(5)=h(4)-eta_turbine*(h(4)-h5s);                        %[kj/kg]
  s(5)= refpropm('S','P',P(5),'H',h(5)*1000,fluid)/1000;   %[kj/kg]
  Q(5)= refpropm('Q','P',P(5),'H',h(5)*1000,fluid);
  T(5)= refpropm('T','P',P(5),'H',h(5)*1000,fluid);        %[k]
  T_C(5)=T(5)-273;                                         %[C]
  m(5)=m(4);                                               %[kg/s]
  
  % 'state 6'
  P(6)=P(5);
  Q(6)=0;
  h(6)= refpropm('H','P',P(6),'Q',Q(6),fluid)/1000;        %[kj/kg]
  s(6)= refpropm('S','P',P(6),'Q',Q(6),fluid)/1000;        %[kj/kg.k]
  T(6)= refpropm('T','P',P(6),'H',h(6)*1000,fluid);        %[k]
  T_C(6)=T(6)-273;                                         %[C]
  m(6)=m(5);                                               %{kg/s]
  
  % 'state 7'
  P(7)=P(6);
  h(7)= refpropm('H','T',T(7),'P',P(7),fluid)/1000;        %[kj/kg]
  s(7)= refpropm('S','T',T(7),'P',P(7),fluid)/1000;        %[kj/kg.k]
  m(7)=m(6);                                               %[kg/s]


  % 'state 8 to 13'
  mus=zeros(1,15);                       % [kj/kg] Salt Chemical Potential 
  muw=zeros(1,15);                       % [kj/kg] Water Chemical Potential
  rho=zeros(1,15);                       % [kg/m^3] Density
  v=zeros(1,15);                         % [m^3/kg] Specific volume
clear i
  for i=8:13
      P(i)=P(8);
      [h(i), s(i)]=seawaterprop(T_C(i),w(i)); 
      [mus(i), muw(i)]=chemicalpotential(T_C(i),w(i)); 
      rho(i)=Density(T_C(i),w(i));                         %
      v(i)=1/rho(i);
  end

 
  % Mass & Energy Balances for MED
  X1=[1 -1 -1;w(8) -w(9) -w(10);h(8) -h(9) -h(10)];
  Y1=[0;0;m(6)*h(6)-m(5)*h(5)];
  Z1=X1\Y1;
  m(8)=Z1(1);                         %[kg/s] MED Feedwater mass flow rate
  m(9)=Z1(2);                         %[kg/s] MED Brine mass flow rate
  m(10)=Z1(3);                        %[kg/s] MED Distillate mass flow rate
  
  % Mass & Energy Balance for Ro
  Pp=0.33;                            %       Power Portion used in pump
  w_P=v(11)*(RO_HP-RO_LP);            %[kj/kg]
  W_ST=m(4)*h(4)-m(5)*h(5);           %[kW]   Turbine Power output
  W_P=Pp*W_ST;                        %[kW]   Pump Power used
  m(11)=2*W_P/w_P;                    %[kg/s] RO Feedwater mass flow rate
  m(12)=(m(11)*w(11))/w(12);          %[kg/s] RO Brine mass flow rate
  m(13)=m(11)-m(12);                  %[kg/s] RO Distillate mass flow rate
  
  % 'Cooling Water inlet'  ('State 14')
  h(14)= refpropm('H','T',T(14),'P',P0,fluid)/1000;     %[kj/kg]
  s(14)= refpropm('S','T',T(14),'P',P0,fluid)/1000;     %[kj/kg.k]
 
  % 'Cooling Water outlet' ('State 15')
  h(15)= refpropm('H','T',T(15),'P',P0,fluid)/1000;     %[kj/kg]
  s(15)= refpropm('S','T',T(15),'P',P0,fluid)/1000;     %[kj/kg.k]

  Q_HEX=m(6)*h(6)-m(7)*h(7);           %[kW] HeatExchanger Q in
  m(14)=Q_HEX/(h(15)-h(14));           %[kg/s] cooling water mass flow rate
  m(15)=m(14);
  
%% Second Law Analysis

  % 'Dead State'
  h0= refpropm('H','T',T0,'P',P0,fluid)/1000;     %[kj/kg]
  s0= refpropm('S','T',T0,'P',P0,fluid)/1000;     %[kj/kg.k]

  ex=zeros(1,17);                                 % [kj/kg] Specific exergy
  Ex=zeros(1,17);                                 % [kW]    Exergy Rate
clear i
  for i=1:17
      if i<=7
          ex(i)=(h(i)-h0)-T0*(s(i)-s0);
      elseif i>=8 && i<=13
          ex(i)=Exergy(T_C(i),w(i),T0_C,w0)/1000;                                                         
      elseif i==16
          ex(i)=W_ST;
      elseif i==17
          ex(i)=W_P;
      else
          ex(i)=(h(i)-h0)-T0*(s(i)-s0);                    
      end
    
  end
  for i=1:17
      if i<=15
        Ex(i)=m(i)*ex(i);
      elseif i==16
        Ex(i)=m(5)*ex(i);
      elseif i==17 
        Ex(i)=m(11)*ex(i);
      end
  end
  
                        %<<Exergy Destruction>>

  % 'Expansion Valve'
  Ex_dest_expvalve=Ex(1)-Ex(2);                           %[kW]
 
  % 'Separator'
  Ex_dest_Sep=Ex(2)-Ex(3)-Ex(4);                          %[kW]
 
  % 'Steam Turbine'
  Ex_dest_ST=Ex(4)-Ex(5)-W_ST;                            %[kW]
 
  % 'Heat Exchanger'
  Ex_dest_HEX=Ex(6)+Ex(14)-Ex(7)-Ex(15);                  %[kW]

  % 'MED'
  Ex_dest_MED=Ex(5)-Ex(6)+Ex(8)-Ex(9)-Ex(10);             %[kW]
 
  % 'RO'
  Ex_dest_RO=W_P+Ex(11)-Ex(13)-Ex(12);                    %[kW]
  
%% Exergo-Economic Analysis


 CRF=(int_rat*(int_rat+1)^n/((int_rat+1)^n-1));
 phi=1.06; 
 t=7000;                                                  %[h]
 
 % 'Expansion Valve'
 Z_C_expvalve=200;                                        %[$]
 Z_expvalve=(CRF*phi*Z_C_expvalve)/t;                     %[$/h]
 
 % 'Steam Turbine'
 Z_C_ST=6000*(W_ST)^0.7;                                  %[$]
 Z_ST=(CRF*phi*Z_C_ST)/t;                                 %[$/h]
 
 % 'Separator'
 Z_C_Sep=0.07*Z_C_ST;                                     %[$]
 Z_Sep=(CRF*phi*Z_C_Sep)/t;                               %[$/h]
 
 % 'Heat Exchanger'
 DeltaT1=T(6)-T(7);                                       %[C]
 DeltaT2=T(15)-T(14);                                     %[C]
 T_LMTD=(DeltaT1-DeltaT2)/(log(DeltaT1/DeltaT2));         %[C]
 U=2;                                                     %[kW/m^2.C]
 A_HEX=Q_HEX/(T_LMTD*U);                                  %[m^2]
 Z_C_HEX=588*(A_HEX)^0.8;                                 %[$]
 Z_HEX=(CRF*phi*Z_C_HEX)/t;                               %[$/h]
 
 % 'MED'
 Z_C_MED=1500*m(10)*3600*24/997;                          %[kW]
 Z_MED=(CRF*phi*Z_C_MED)/t;                               %[$/h]
 
 % 'RO'
 Z_C_RO=900*m(13)*3600*24/rho(13) ;                       %[kW]
 Z_RO=(CRF*phi*Z_C_RO)/t;                                 %[$/h]
 
 % 'Expansion Valve'
 c_geo=1.3;                                               %[$/Gj]


                    %<< Calculatin cost per unit exergy (c) >>

 %    1    2      3      4      5      6      7      8     9      10        
 A2=[Ex(1),-Ex(2),0     ,0     ,0     ,0     ,0     ,0    ,0     ,0        %1   Number of Equation    
     1    ,0     ,0     ,0     ,0     ,0     ,0     ,0    ,0     ,0        %2
     0    ,Ex(2) ,-Ex(3),-Ex(4),0     ,0     ,0     ,0    ,0     ,0        %3
     0    ,0     ,1     ,-1    ,0     ,0     ,0     ,0    ,0     ,0        %4
     0    ,0     ,0     ,Ex(4) ,-Ex(5),0     ,0     ,0    ,0     ,0        %5 
     0    ,0     ,0     ,1     ,-1    ,0     ,0     ,0    ,0     ,0        %6
     0    ,0     ,0     ,0     ,Ex(5) ,-Ex(6),0     ,Ex(8),-Ex(9),-Ex(10)  %7
     0    ,0     ,0     ,0     ,0     ,0     ,0     ,1    ,0     ,0        %8
     0    ,0     ,0     ,0     ,0     ,0     ,0     ,0    ,1     ,0        %9
     0    ,0     ,0     ,0     ,1     ,-1    ,0     ,0    ,0     ,0        %10
     0    ,0     ,0     ,0     ,0     ,Ex(6) ,-Ex(7),0    ,0     ,0        %11
     0    ,0     ,0     ,0     ,0     ,0     ,0     ,0    ,0     ,0        %12
     0    ,0     ,0     ,0     ,0     ,1     ,-1    ,0    ,0     ,0        %13
     0    ,0     ,0     ,0     ,0     ,0     ,0     ,0    ,0     ,0        %14
     0    ,0     ,0     ,0     ,0     ,0     ,0     ,0    ,0     ,0        %15
     0    ,0     ,0     ,0     ,0     ,0     ,0     ,0    ,0     ,0        %16
     0    ,0     ,0     ,0     ,0     ,0     ,0     ,0    ,0     ,0];      %17


%    11    12      13    Wat_in  Wat_out  Turb.  Ro 
%                          14       15
 A3=[0     ,0      ,0      ,0      ,0      ,0    ,0                        %1
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %2
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %3
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %4
     0     ,0      ,0      ,0      ,0      ,-W_ST,0                        %5
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %6
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %7
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %8
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %9
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %10
     0     ,0      ,0      ,Ex(14) ,-Ex(15),0    ,0                        %11
     0     ,0      ,0      ,1      ,0      ,0    ,0                        %12
     0     ,0      ,0      ,0      ,0      ,0    ,0                        %13
     Ex(11),-Ex(12),-Ex(13),0      ,0      ,0    ,W_P                      %14
     1     ,0      ,0      ,0      ,0      ,0    ,0                        %15
     0     ,1      ,0      ,0      ,0      ,0    ,0                        %16
     0     ,0      ,0      ,0      ,0      ,1    ,-1];                     %17
 A1=[A2,A3];

 B1=[-Z_expvalve                        %1                        
     1.3                                %2
     -Z_Sep*(1e6/3600)                  %3
     0                                  %4
     -Z_ST*(1e6/3600)                   %5
     0                                  %6
     -Z_MED*(1e6/3600)                  %7
     0                                  %8
     0                                  %9
     0                                  %10
     -Z_HEX*(1e6/3600)                  %11
     0                                  %12 
     0                                  %13
     -Z_RO*(1e6/3600)                   %14
     0                                  %15
     0                                  %16
     0];                                %17

  D1=A1\B1;

  c=zeros(1,17);                        % [$/Gj] cost per unit exergy
  C=zeros(1,17);                        % [$/h]  cost rate
  clear i
  for i=1:17
      c(i)=D1(i);
  end
  for i=1:17
      C(i)=Ex(i)*c(i)*3600/1e6;
  end
  


  ExD=[Ex_dest_expvalve, Ex_dest_Sep, Ex_dest_ST, Ex_dest_MED,...
           Ex_dest_HEX,      Ex_dest_RO];
    ExD_Fig(j,:)=ExD;


  Z_dot=[Z_expvalve, Z_Sep, Z_ST, Z_MED,  Z_HEX,  Z_RO]; 



                    %<< Calculating cost Rates >>

% 1: Expansion Valve      
Cf(1)=C(1);
Exf(1)=m(1)*ex(1);
Exp(1)=m(2)*ex(2);

% 2: Seperator
Cf(2)=C(2);
Exf(2)=m(2)*ex(2);
Exp(2)=m(3)*ex(3)+m(4)*ex(4);

% 3: Steam Turbine
Cf(3)=C(4)-C(5);
Exf(3)=m(4)*ex(4)-m(5)*ex(5);
Exp(3)=W_ST;

% 4: MED Unit
Cf(4)=C(5)-C(6)+C(8)-C(9);
Exf(4)=(m(5)*ex(5)-m(6)*ex(6)+m(8)*ex(8)-m(9)*ex(9));
Exp(4)=m(10)*ex(10);


% 5: Heat Exchanger
Cf(5)=C(15)-C(14);
Exf(5)=m(15)*ex(15)-m(14)*ex(14);
Exp(5)=m(6)*ex(6)-m(7)*ex(7);

% 6: RO Unit
Cf(6)=C(13)-Z_RO;
Exf(6)=(W_P+ m(11)*ex(11));
Exp(6)=m(13)*ex(13)+m(12)*ex(12);

cf=zeros(1,6);
for i=1:6
    
        cf(i)=Cf(i)/(Exf(i)*3600/1e6);
    
end


% Exergy losses
Exl(1)=0;               % Expansion Valve 
Exl(2)=m(3)*ex(3);      % Seperator
Exl(3)=0;               % Turbine
Exl(4)=m(9)*ex(9);      % MED Unit
Exl(5)=m(7)*ex(7);      % Heat Exchanger    
Exl(6)=m(12)*ex(12);    % RO Unit


C_DES=zeros(1,6);
C_loss=zeros(1,6);
C_Tot=zeros(1,6);
f=zeros(1,6);
clear i
for i=1:6
C_DES(i)=cf(i)*ExD(i)*3600/1e6;
C_loss(i)=cf(i)*Exl(i)*3600/1e6;
C_Tot(i)=Z_dot(i)+C_DES(i)+C_loss(i);
f(i)=Z_dot(i)/C_Tot(i);
end
C_tot_Fig(j,:)=C_Tot;
Overall_Eff=   ( (W_ST-W_P) + Ex(10) + Ex(13) -Ex(12)-Ex(9)) /...
           ( (Ex(1)-Ex(7)-Ex(3)) + (Ex(8)+ Ex(11)) );
Overall_Eff_Fig(j,1)=Overall_Eff;

end
Ex_Des_Figs{q}=ExD_Fig;
Overall_Eff_Figs{q}=Overall_Eff_Fig;
C_Tot_Figs{q}=C_tot_Fig;
end


%% Figures
color=[1 0 0 ; 0 1 0 ; 0 0 1 ; 0 1 1; 1 0 1 ; 1 1 0];
shape=["square","+","o","diamond","x","^"];
Compon={'Expansion Valve'; 'Steam Turbine'; 'MED'; 'Heat Exchanger';...
    'RO';'Overall Exergy'};
for q=1:3

clear i

Y=Ex_Des_Figs{q};
X=cell2mat(Fig(q));
Overall_Eff_Fig=Overall_Eff_Figs{q};

figure(q)

    for i=1:6
        if i==2
            continue

        else
        yyaxis left
        plot(X,Y(:,i),'Color',color(i,:),'Marker',shape(i)...
            ,'LineWidth',1,'LineStyle','-','MarkerEdgeColor','k')
        ylabel('Exergy Destruction (kW)')
        hold on

        yyaxis right

        plot(X,Overall_Eff_Fig(:,1),'k'...
            ,'LineWidth',1,'MarkerEdgeColor','k')
        hold on
        ylabel('Overall Efficiency')
        end
    end
    legend(Compon);
    
    xlabel(X_axis_labels(q))
    
end
hold off


for q=4:5
    Y=C_Tot_Figs{q};
    X=cell2mat(Fig(q));

    figure(q)
    for i=1:6
        if i==2
            continue
        else
        yyaxis left
     plot(X,Y(:,i),'Color',color(i,:),'Marker',shape(i)...
            ,'LineWidth',1,'LineStyle','-','MarkerEdgeColor','k')
        hold on
        
        end
    end
    legend(Compon(1:end-1));
    xlabel(X_axis_labels(q))
    ylabel('C_T $/h')
    hold off
end

%% Results



Fig_3=readmatrix('Fig3.xlsx');
figure (6)
j=1;
for i=1:3:14
    yyaxis left
    plot(Fig_3(:,i)*1000,Fig_3(:,i+1),'Color',color(j,:),'Marker',shape(j)...
            ,'LineWidth',1,'LineStyle','--','MarkerEdgeColor','k')
        hold on
    hold on
    j=j+1;
end
%       state     T       P      h      s     ex    Q
Art_1=  [1       200    1.55   852.3  2.33   162.3  0
         2       191.6  1.3    852.3  2.33   162    0.02
         3       191.6  1.3    814.6  2.25   148.4  0
         4       191.6  1.3    2786   6.5    856    1
         5       99.62  0.1    2461   6.8    443.6  0.9
         6       99.62  0.1    417.5  1.3    33.8   0
         7       45     0.1    188.5  0.64   2.72   0
         8       25     0.101  100.2  0.35   0      0
         9       53     0.101  202.2  0.66   6.1    0
         10      53     0.101  222    0.74   7.5    0
         11      25     0.101  100.2  0.35   0      0
         12      25     0.101  94.31  0.32   1.35   0
         13      25     0.101  98.53  0.32   2.36   0]';


%     state  Ex(kW)   c($/Gj)    C_dot($/h)
Art_2= [1    24341    1.3        112
        2    24281    1.303      112
        3    21828    1.31       101.1
        4    2453     1.31       11.35
        5    1271     1.31       5.88
        6    96.83    1.31       0.45
        7    7.8      1.31       0.036
        8    0        0          0
        9    144.6    0          0
        10   212.2    49.68      37.32
        11   0        0          0
        12   72.58    0          0
        13   150.7    86.8       46.3]';


%      Z_dot($/h)  c($/Gj)   C_des($/h   C_loss($/h)   C_tot($/h)     f  
Art_3=[0.002       1.3       0.27        0             0.28         0.62   %Expansion valve
       0.34        1.3       0           100.76        101.1        0.33   % Separator
       6.25        1.3       1.15        0             7.41         84.4   %Steam turbine
       0.032       4.34      0.9         0.5           0.93         3.5    %Heat exchanger
       31.9        1.47      4.3         0.76          37           86     %MED unit
       48.13       5.65      0.2         1.32          44.86        96]';  %RO unit
   


   %{
% Table 1 validation data
 Val1=[1:13; T_C(1:13); P(1:13)/1000; h(1:13); s(1:13); ex(1:13); Q(1:13)];
 Tab1=table((1:13)', T_C(1:13)', P(1:13)'/1000, h(1:13)',...
     s(1:13)', ex(1:13)', Q(1:13)','VariableNames',...
     {'State','Temp.','Press.','Enthalpy','Entropy', 'Exergy','Quality'})
  
% Table 2 validation data
 Val2=[1:13; Ex(1:13); c(1:13); C(1:13)];
 Tab2=table((1:13)', Ex(1:13)', c(1:13)', C(1:13)','VariableNames',...
     {'State','Exergy Rate','c','C'})
% Table 3 validation data
Compon={'Expansion Valve'; 'Seperator'; 'Steam Turbine'; 'MED'; 'Heat Exchanger';...
    'RO'};
 Val3=[Z_dot ; cf; C_DES; C_loss; C_Tot; f];
 Tab3=table(Z_dot' , cf', C_DES', C_loss', C_Tot', f','VariableNames',...
     {'Z_dot','cf','C Destruction','C loss',' C Total', 'f'},...
     'RowNames',Compon)


 error1=(Art_1-Val1)./Art_1.*100;error1(1,:)=1:13;
 error2=(Art_2-Val2)./Art_2.*100;error2(1,:)=1:13;
%}