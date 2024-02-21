function [S_perturbed,idx_perturbed,stab, Try, L_pert] = wsc(gp, Wn, zeta, Wi, var_percent,variations, C_optm, Cff)

s=tf('s');
stability=zeros(8,1);

for i = 1:8
    % calcolo i parametri con la variazione
    gp_perturbed = gp * (1 + variations(i, 1) * var_percent);
    Wn_perturbed = Wn * (1 + variations(i, 2) * var_percent);
    zeta_perturbed = zeta * (1 + variations(i, 3) * var_percent);

    % calcolo il processo i-esimo con i parametri variati
    P_perturbed=(gp_perturbed*Wn_perturbed^2)/(s*(s^2 + 2*zeta_perturbed*Wn_perturbed*s + Wn_perturbed^2));

    % calcolo il loopgain con la perturbazione i-esima
    L_perturbed(i)=C_optm * P_perturbed;

    T_perturbed=feedback(L_perturbed(i),1);
    stability(i)=double(isstable(T_perturbed));

    % calcolo Try i-esimo
    Try_perturbed(i) = Cff * T_perturbed;

    sensibility(i)=min(freqresp(L_perturbed(i),Wi));
end

% calcolo la variazione massima
[S_perturbed, idx_perturbed]=max(sensibility);
stab=stability;
Try = Try_perturbed(idx_perturbed);
L_pert = L_perturbed(idx_perturbed);