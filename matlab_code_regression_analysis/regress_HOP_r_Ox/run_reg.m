%driver code to run regression analysis and fit oxalic acid yields for
%sonochemical benzyl alcohol oxidation experiment.


clc
clear
global k_HOP r_BA_frag r_BA_cons k_OOH_O2m k_2OOH K_OOH_AB t_exp OA_exp tspan out_prev BA_0


%kinetic data
t_exp = [3600	7200	10800];%	time in seconds
OA_exp = [0.024502	0.0472725	0.065469];%oxalic acid concentration in mM

OA_exp = OA_exp /1000;%convert mM to M

tspan = linspace(0,22000,5000);

%specify the factor to increase K_DFT
K_inc = 1;
%DFT value for O2- addition to HOP intermediate
K_EQ = 4.0E-04;
%diffusion rate constant
k_d = 1.159E10
%calculate pseudo-first-order rate constant for HOP oxidation
k_HOP = K_EQ*K_inc*k_d;
%benzyl alcohol consumption rate
r_BA_cons = 1.897E-8;
%rate constant for OOH + O2m reaction
k_OOH_O2m = 9.7E7;
%rate constant for OOH + OOH reaction
k_2OOH = 8.3E5;
%OOH acid dissociation constant at 315 K
K_OOH_AB = 10^(-4.62);

%initial benzyl alcohol concentration
BA_0 = 5/1000;

%lower and upper bounds for regression analysis
lb = 0;
ub = .91;
%fit the rate of BA -> HOP reaction, in units of x10^(-8) M s-1 
beta = .5;

%run regression with run_ode function
options = optimoptions('lsqnonlin', 'DiffMinChange',.0001,'Display','iter');
[beta_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@model_HOP,beta,lb,ub, options);
ci = nlparci(beta_fit,residual,'jacobian',jacobian);% 95% confidence intervales

%initialize concentrations of HOP, oxalic acid, and formic acid
%initialize with small concentration
C0 = [1,1,1].*1E-10;

%simulate time-dependent yields for HOP, OA, and FA
options=odeset('RelTol',10^-6,'AbsTol',10^-6);
[t,y]=ode23s(@(t,y) rxn_network(t,y,beta_fit)...
            ,tspan,C0);

MAE = mean(abs(residual)./OA_exp')

%save output
save("reg_0")