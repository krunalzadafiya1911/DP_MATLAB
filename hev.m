function [X, C, I, out] = hev(inp,~)
%function [X C I out] = hev(inp,par)
%HEV Computes the resulting state-of-charge based on current state-
%   of-charge, inputs and drive cycle demand.
%   
%   [X C I out] = HEV(INP,PAR)
%
%   INP   = input structure
%   PAR   = user defined parameters
%
%   X     = resulting state-of-charge
%   C     = cost matrix
%   I     = infeasible matrix
%   out   = user defined output signals

wg = inp.W{1};           %angular velocity
dwg = inp.W{2};          %angualar acceleration
Ttot = inp.W{3};         %torque 

% Dynamic parameters
u1 = (Ttot>0).*(inp.U{1}<1);     %engine start-stop
	                 %split ratio

T_EM_max = [60;60;60;44.5714285714286;32.5714285714286;24;16.2857142857143;11.1428571428571;7.71428571428571];
w_EM_max = [0,100,200,300,400,500,600,700,800];

theta_EM = 0.1;              % define motor inertia
epsilon = 0.01;              % define epsilon (cf. Slide 3-8 and 3-10)

Tm = Ttot.*(inp.U{1});          % totque of motor
Te = Ttot.*(1 - inp.U{1});

Ttot = Ttot*ones(size(Ttot));
inc = ~((~((Ttot>0).*(u1==0))+(Tm==Ttot)).*(~(Ttot<=0)+(Tm>=(Ttot))));

% TORQUE SPLIT
engine_inertia = 0.2;
motor_inertia  = 0.1;
wg_idle = 105;                              % Engine minimum speed

Te = Te + engine_inertia*dwg;
Tm = Tm + motor_inertia*dwg;

%engine block
P_CE_idle = 8000;
T_CE_idle = P_CE_idle/wg_idle;              % Engine minimum torque
Te = Te.*(Te>T_CE_idle) + T_CE_idle.*(Te<=T_CE_idle);
wg = wg.*(wg>wg_idle) + wg_idle.*(wg<=wg_idle);

wg = ones(size(Te)).*wg;

load OM_622;
%use with interp2
V_CE = H_u* interpn(w_CE_row,T_CE_col,V_CE_map, wg, Te,'linear',0.000142);

P_CE = P_CE_idle.*(wg<=wg_idle).*(Te<=T_CE_idle) + V_CE - V_CE.*((Te<=T_CE_idle).*(wg<=wg_idle));

P_CE = P_CE.*(Te >= (T_CE_idle + 0.01)) + 0.*(Te < (T_CE_idle + 0.01));

w_CE_upper = max(w_CE_max);                 % Upper limit engine speed [rad/s]

ine = ((interp1(w_CE_max,T_CE_max, wg) - abs(Te))<0) | ((w_CE_upper - wg)<0);

%motor block
load EM;

%efficiency of motor

efficiency = interpn(w_EM_row,T_EM_col,eta_EM_map, wg, Tm,'linear',0.000142);

%power used from battery
P_EM = Tm.*wg.*efficiency;

% motor feasibility
% T_EM_max = [60.0 60.0 60.0 44.57142857142858 32.57142857142858 24.0 16.28571428571429 11.142857142857144 7.714285714285714];
% w_EM_upper = 800;               %max speed of motor
w_EM_upper = max(w_EM_max);         % Upper limit motor speed
inm = (interp1(w_EM_row,T_EM_max, wg,'linear') - abs(Tm))<0 | (w_EM_upper - wg)<0 ;

% battery block
load BT;

I_BT_max = 300;              % Max battery current
Q_BT_0 = 36000;              % Battery charge
Q_BT = inp.X{1}*3600*I_0;    % Initial battery charge
P_BT = P_EM;                 % Power of battery

e_C = 1.*((-1.*(P_BT.*(P_BT<0) + 0.*(P_BT>=0)))>0) - 1.*((-1.*(P_BT.*(P_BT<0) + 0.*(P_BT>=0)))<0) + 0.*((-1.*(P_BT.*(P_BT<0) + 0.*(P_BT>=0)))==0);
e_D = 1.*((P_BT.*(P_BT > 0) + 0.*(P_BT<=0))>0) - 1.*((P_BT.*(P_BT > 0) + 0.*(P_BT<=0))<0) + 0.*((P_BT.*(P_BT > 0) + 0.*(P_BT<=0))==0);
e_I = 1 - abs(1.*(P_BT > 0) - 1.*(P_BT < 0) + 0.*(P_BT == 0));

U_BT_L = (e_C > 0).* (c_BT_L3.*Q_BT./Q_BT_0 + c_BT_L1 + sqrt((c_BT_L3.*Q_BT./Q_BT_0 + c_BT_L1).^2 + 4*(c_BT_L4.*Q_BT./Q_BT_0 + c_BT_L2).*(-P_BT./I_0)))/2;
U_BT_E = (e_D > 0).* (c_BT_E3.*Q_BT./Q_BT_0 + c_BT_E1 + sqrt((c_BT_E3.*Q_BT./Q_BT_0 + c_BT_E1).^2 + 4*(c_BT_E4.*Q_BT./Q_BT_0 + c_BT_E2).*(-P_BT./I_0)))/2;
U_BT_Ei = (e_I > 0).* (Q_BT.*c_BT_L3./Q_BT_0 + c_BT_L1);

U1 = (e_D >= 0.5).* U_BT_E +  (e_D < 0.5).* U_BT_Ei;

U_BT = (e_C >= 0.5).* U_BT_L + (e_C < 0.5).* U1;

I_BT = P_BT./U_BT;

% battery infeasibility
inb1 = (2*U_BT-(Q_BT.*c_BT_E3/Q_BT_0 + c_BT_E1))<0;
inb2 = ((I_BT_max - abs(I_BT))<0);
inb = inb1 | inb2;

% Update State variable
X{1} = (Q_BT - I_BT)/(3600*I_0);
X{1} = (conj(X{1})+X{1})/2;

% COST

% Summarize infeasible matrix
I = (inc+ine+inm+inb~=0);

% Calculate cost matrix (fuel mass flow)
C{1}  = P_CE;

if numel(find(I==0))==0
    %condition check for infeasible system.
    keyboard
end

% SIGNALS
out.Te = Te;
out.Tm = Tm;

end