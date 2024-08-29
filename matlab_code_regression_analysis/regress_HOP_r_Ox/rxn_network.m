function rates = rxn_network(t,y,r_BA_frag)
% input time, concentrations, BA_fragmentation rate (to fit)

    global k_HOP r_BA_cons k_OOH_O2m k_2OOH K_OOH_AB

    HOP = y(1); % hydroxy oxo pentenal
    % y2 and 3 are oxalic acid and formic acid
    OA = y(2);
    FA = y(3);

    
    %initialize pKa values and concentrations of intermediates for input to
    %acid_base.m function. Ox dissociation 1, 2, formic acid dissociation, formates conc, oxalates conc
    p = [10^(-1.23),10^(-4.35),10^(-3.77), OA,	FA]; %last data point
    C0 = [p(end-1)/3,p(end-1)/3,p(end-1)/3,p(end),p(end)/2,p(end)/2];
    %Ox, Ox-, Ox2-, H+, formic acid, formate
    out = fsolve(@(x) acid_base(x,p),C0);
    
    %proton concentration
    Hp = out(4);

    %Pseudo-steady state O2m concentration
    O2m = sqrt(r_BA_cons/(Hp*k_OOH_O2m/K_OOH_AB + Hp*Hp*k_2OOH/K_OOH_AB ));
    
    % rate of BA -> HOP
    r1 = r_BA_frag*1E-8;

    % forward rate of HOP consumption
    r2 = k_HOP*O2m*HOP*Hp;
    
    %changes in concentrations of each product
    r_HOP = r1-r2;
    r_OA = r2;
    r_FA = 2*r_OA;
    rates = [r_HOP;r_OA;r_FA];
end