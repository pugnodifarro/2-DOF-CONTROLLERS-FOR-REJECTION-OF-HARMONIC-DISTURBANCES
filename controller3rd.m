s = tf('s');    % variabile complessa per tf

% -------------------------- DATI DI PROGETTO -----------------------------
gp = 0.10;      % ms^-1*V^-1
Wn = 97;        % rad/s
zeta = 0.25;    % dumping

Wi = (40:5:80);                     % frequenze insieme delta
Weta = [(10:5:40), (80:5:1000)];    % frequenze insieme eta
M_min = 0.165;                      % modulus margin
alpha_value = 0.08;
P = (gp*Wn^2)/(s*(s^2 + 2*zeta*Wn*s + Wn^2));

lb  = [100, 0.00001];
ub  = [2000,  0.1];

% ----------------------------- 3RD ORDER ---------------------------------
n_k1 = 10; % numero di elementi di k1
n_k2 = 60; % numero di elementi di k2

% intervalli di k1 e k2
k1 = linspace(lb(1), ub(1), n_k1);
k2 = linspace(lb(2), ub(2), n_k2);

optimized_k1 = 0;
optimized_k2 = 0;
optimized_k3 = alpha_value;
obj=-inf;

% metodo brute force per calcolare k1 e k2 optm
for i = n_k1:-1:1
    for j = 1:n_k2
        L = loopgain(k1(i), k2(j), alpha_value, P);
        stab = all(isstable(feedback(L,1), Wi));
        modulus = min(abs(squeeze(1+freqresp(L, Weta))));
        qq = min(abs(1+ squeeze(freqresp(L, Wi))));
        if stab == 1 && modulus > M_min && qq > obj
            obj = qq;
            optimized_k1 = k1(i)
            optimized_k2 = k2(j)
        end
    end
end

%%

% stampa le combinazioni valide
disp("Combinazioni valide:");
disp("k1            k2        alpha_value");
disp([optimized_k1 optimized_k2 optimized_k3]);

% calcolo funzione ottimizzata
C_optm=optimized_k1*(1 + optimized_k2*s)^3/(1 + optimized_k2*optimized_k3*s)^3;
L_optm=C_optm*P;

% calcolo l'effetto minimo di reiezione al disturbo
dummy=abs(squeeze(1+freqresp(L_optm,Wi)));
S=min(dummy);
S_1=1/S

T=feedback(L_optm,1);
T_1=T^-1;

%Controllore FEEDFORWARD 
Cff=T_1/((1+1/5*s)*(1+1/10*s)*(1+1/15*s)); 
Try=minreal(Cff*T);

%%

%WORST CASE SCENARIO (VERIFICA DI ROBUSTEZZA)

stability = zeros(8,1);
sensibility = zeros(8,1);

% variazione percentuale
var_percent = 0.02;

% matrice per calcolare i vertici dell'insieme wcs
variations =  [ 1,  1,  1;     %1° variazione
                1,  1, -1;     %2°
                1, -1,  1;     %3°
                1, -1, -1;     %4°
               -1,  1,  1;     %5°
               -1,  1, -1;     %6°
               -1, -1,  1;     %7°
               -1, -1, -1];    %8° variazione

[S_perturbed,idx_perturbed,stab, Try_perturbed, L_perturbed] = ...
    wsc(gp, Wn, zeta, Wi, var_percent,variations, C_optm, Cff);

S_1_perturbed=1/abs(S_perturbed)

% controllo la differenza sulle prestazioni temporali con gradino in ingresso 
step(Try_perturbed);
stepinfo(Try_perturbed)

%%

%STAMPA DELLE IMMAGINI

% controllo il margine di ampiezza e di fase per L_optm
figure
margin(L_optm)

% stampa luogo delle radici e nyquist ottimizzati
figure 
nyquist(L_optm);
hold on
cerchio(-1,0,M_min);
axis equal
title('Diagramma di Nyquist ottimizzato')
hold off

% controllo il margine di ampiezza e di fase per L_perturbed
figure
margin(L_perturbed)

%%

%FUNZIONE C*P=L
function L = loopgain(k1, k2, alpha_value, P)
s = tf('s');
C = k1*(1 + k2*s)^3/(1 + alpha_value*k2*s)^3;
L = C*P; 
end