function out = acid_base(x,p)
% run acid-base equilibrium calculation for oxalic acid and formic acid in
% aqueous solution

%parameters (p) Ox dissociation 1, 2, formic acid dissociation, formic acid conc, ox conc
%x = 2, Ox-; 3, Ox--; 4, H+; 5, formic acid
    eq0 = x(2)*x(4)-p(1)*x(1); % Ox dissociation 1
    eq1 = x(3)*x(4)-p(2)*x(2); % Ox dissociation 2
    eq2 = x(5)*x(4)-p(3)*x(6); % formic acid dissociation 1
    eq3 = x(2)+2*x(3)+x(5)-x(4); % charnge balance (neglect OH- concentration)
    eq4 = p(4)-x(5)-x(6); %balance glyox
    eq5 = p(5)-x(1)-x(2)-x(3); % balance Ox
    out = [eq0,eq1,eq2,eq3,eq4,eq5]; % concentrations of Ox, Ox-, Ox2-, H+, formic acid, formate-