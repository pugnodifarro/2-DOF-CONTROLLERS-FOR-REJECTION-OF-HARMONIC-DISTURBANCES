s = tf('s');    % variabile complessa per tf

% -------------------------- DATI DI PROGETTO -----------------------------
gp = 0.10;      % ms^-1*V^-1
Wn = 97;        % rad/s
zeta = 0.25;    % dumping

Wi = (40:5:80);                    % frequenze insieme delta
Weta = [(10:5:35), (85:5:1000)];   % frequenze insieme eta
M_min = 0.2;                       % modulus margin
alpha_value = 0.01;
P = (gp*Wn^2)/(s*(s^2 + 2*zeta*Wn*s + Wn^2));

% Definisco i valori che mi serviranno per fminimax
A   = [];
b   = [];
Aeq = [];
beq = [];
lb  = [100, 0.001];
ub  = [500,  50];
x0  = [300, 0.02];
options = optimset('Display', 'iter');

% --------------------------- PD CONTROLLER -------------------------------
n_k1 = 6000; % numero di elementi di k1 [VARIABILE]
n_k2 = 6000; % numero di elementi di k2 [VARIABILE]

% intervalli di k1 e k2
k1 = linspace(lb(1), ub(1), n_k1);
k2 = linspace(lb(2), ub(2), n_k2);

% funzione obiettivo
objFunc = @(params) -min(abs(squeeze(freqresp(1 + loopgain(params(1), ...
    params(2), alpha_value, P), Wi))));

% vincoli non lineari
nonlincon = @(params) nonlinconstr(params(1), params(2), alpha_value, ...
    P, Weta, M_min);

% ottimizzazione con fminimax
optimizedParams = fminimax(objFunc, x0, A, b, Aeq, beq, lb, ub, ...
    nonlincon, options);

% valori ottimi di k1 e k2
optimized_k1 = optimizedParams(1);
optimized_k2 = optimizedParams(2);
optimized_k3 = alpha_value;

% stampa le combinazioni valide
disp("Combinazioni valide:");
disp("k1       k2        k3");
disp([optimized_k1 optimized_k2 optimized_k3]);

%%

% calcolo funzione ottimizzata
C_optm=optimized_k1*(1 + optimized_k2*s)/(1 + optimized_k2*optimized_k3*s);
L_optm=C_optm*P;

figure
margin(L_optm);

%%

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
% ( non uso più l'indice del massimo valore perché passo già il valore
% massimo usando la funzione wsc.m
step(Try_perturbed);
stepinfo(Try_perturbed)

%%

%STAMPA DELLE IMMAGINI

% controllo margine di ampiezza e di fase su L_optm
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

% controllo il margine di ampiezza e di fase su L_perturbed
figure
margin(L_perturbed)

%%

%FUNZIONI

% funzione per calcolare L = C*P
function L = loopgain(k1, k2, alpha_value, P)
s = tf('s');
C = k1*(1 + k2*s)/(1 + alpha_value*k2*s);
L = C*P; 
end

% funzione per calcolare i vincoli non lineari e stabilità
function [ineqcons, eqcons] = nonlinconstr(k1, k2, alpha_value, P, Weta, M_min)
   
    s = tf('s');
    L = (loopgain(k1, k2, alpha_value, P));
    
    % Estrazione del denominatore del sistema (polinomio caratteristico)
    [~, den] = tfdata(feedback(L, 1), 'v');

    % Verifica dei vincoli non lineari
    ineqcons = M_min - min(abs(squeeze(freqresp(1 + L, Weta))));
    eqcons = 1 - routh(den,0.01);
end