function out = model_HOP(beta)
    %solve the kinetic model for benzyl alcohol oxidation for a given beta
    %value (rate of BA -> HOP reaction) and calculate the residuals

    global k_HOP r_BA_cons k_OOH_O2m k_2OOH K_OOH_AB t_exp OA_exp tspan out_prev BA_0

    % C1 = HOP, C2 = OA, C3 = FA
    %initialize with small concentration
    C0 = [1,1,1].*1E-10;
   
    %simulate time-dependent yields for HOP, OA, and FA
    options=odeset('RelTol',10^-6,'AbsTol',10^-6);
    [t,y]=ode23s(@(t,y) rxn_network(t,y,beta)...
                ,tspan,C0);

    C_OA_outs = [];

    % collect the oxalic acid yields at times of experimental measurements
    for i=1:length(t_exp)
        t_sel = t_exp(i);
        sel_ind = sum(t<t_sel); %get the row for that data point
        C_OA_out = y(sel_ind,2);

        C_OA_outs = [C_OA_outs;C_OA_out];
    end

    %calculate residuals renormalizing units
    out = (C_OA_outs-OA_exp')./BA_0;
end