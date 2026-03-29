% TAIPS_02A - Simulator pentru identificarea modelului stocastic (AR/ARMA)

% Autori: Andrei BLEJAN
% Creat: 7 Ianuarie 2026

clear; clc; close all;

% Alegere fisier
while true
    idx = input('Introduceti indexul seriei (ex: 20 pentru V20.mat): ');
    nume_fisier = sprintf('V%02d.mat', idx);

    if exist(nume_fisier, 'file'), break; end

    fprintf('Fisierul nu exista. Incercati alt fisier\n');
end

% Incarcare date
load(nume_fisier, 'V');

% Alegere canal (daca este cazul)

if size(V.y,2) > 1
    fprintf('Fisierul contine %d canale de date.\nAlegeti canalul dorit: ', size(V.y,2));
    while true
        ch = input('');

        if ch > 0 && ch <= size(V.y,2), break; end

        fprintf('Alegeti un canal valid!(1-%d): ', size(V.y,2));
    end
    v_total = V.y(:, ch);
else
    v_total = V.y;
end

N = length(v_total);

% Configurare model
limita_ordin = floor(N / 3);

fprintf('\n--- Selectie Model ---\n');
fprintf('1. Model AR\n');
fprintf('2. Model ARMA\n');
while true
    tip_model = input('Alegeti tipul (1/2): ');
    if tip_model==1 || tip_model==2, break; end
end

nc = 0; nalpha = [];

if tip_model == 1
    fprintf('Ati ales Model AR.\n');
    while true
        na = input(sprintf('Ordin na (1...%d): ', limita_ordin));
        if na > 0 && na <= limita_ordin, break; end
    end
else
    fprintf('Ati ales Model ARMA.\n');
    while true
        na = input(sprintf('Ordin na (1...%d): ', limita_ordin));
        if na > 0 && na <= limita_ordin, break; end
    end
    while true
        nc = input(sprintf('Ordin nc (1...%d): ', limita_ordin));
        if nc > 0 && nc <= limita_ordin, break; end
    end

    nalpha = 3 * (na + nc);
    nalpha_aux = input(sprintf('Ordin nalpha (Enter pt %d): ', nalpha), 's');

    if ~isempty(nalpha_aux), nalpha = str2double(nalpha_aux); end
end

% 3. Identificare si predictie
K = 5;
N_masura = N - K;

v_masura = v_total(1:N_masura);       
v_validare = v_total(N_masura+1:end); 

% --- Identificare ---
[yAR, e, lambda2, theta_AR, theta_ARMA] = stochastic(v_masura, na, nc, nalpha);

% --- Predictie ---
y_pred = zeros(K, 1);       
v_masura_aux = v_masura; % buffer pentru istoricul datelor                      
ordin = length(theta_AR) - 1; % ordinul polinomului de predictie

for k = 1:K
    esantioane = v_masura_aux(end:-1:end-ordin+1); % extragem esantioanele necesare pentru calcul ( y[n-1] -> y[n-na] )
    valoare_predictata = - (theta_AR(2:end).' * esantioane); % valoarea predictata la momentul k
    y_pred(k) = valoare_predictata;
    v_masura_aux = [v_masura_aux; valoare_predictata]; % adaugam valoarea predictata in istoric
end

% --- Calcul erori si SNR ---

% Orizontul de masura
var_v_mas = var(v_masura, 1);
var_e_mas = var(e, 1);    
SNR_masura = sqrt(var_v_mas) / sqrt(var_e_mas);

% Orizontul de predictie
epsilon_pred = v_validare - y_pred; % eroarea de predictie
var_v_val = var(v_validare, 1);
var_eps_pred = var(epsilon_pred, 1);
SNR_pred = sqrt(var_v_val) / sqrt(var_eps_pred);

% --- Afisare coeficienti in consola ---
fprintf('\n========================================\n');
fprintf(' Rezultate Identificare\n');
fprintf('========================================\n');
if nc == 0
    % Cazul AR
    fprintf('Model AR(%d):.\n', na);
    fprintf('Coeficienti polinom A(q^-1):\n');
    disp(theta_AR');
else
    % Cazul ARMA
    fprintf('Model ARMA(%d, %d):\n', na, nc);

    A_coeffs = theta_ARMA(1:(na+1), 1);
    C_coeffs = theta_ARMA(1:(nc+1), 2);
    
    fprintf('Coeficienti polinom A(q^-1) (Auto-Regresiv):\n');
    disp(A_coeffs');
    
    fprintf('Coeficienti polinom C(q^-1) (Medie Mobila):\n');
    disp(C_coeffs');
end

% --- Grafice ---

figure('Name', ['Analiza ' nume_fisier], 'NumberTitle', 'off');

% Axele de timp
t_total = 0:N-1;
t_masura = 0:N_masura-1;
t_pred = N_masura:N-1;

% Grafic 1: Zgomot colorat complet + suprapunere model
subplot(3, 2, [1 2]);

plot(t_masura, v_masura, 'b', 'LineWidth', 1); hold on;
plot(t_pred, v_validare, 'c', 'LineWidth', 1); 
plot(t_masura, yAR, 'g', 'LineWidth', 2); 
plot(t_pred, y_pred, 'r', 'LineWidth', 2); 

% Separator vertical
line([N_masura-0.5 N_masura-0.5], ylim, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);

% Titlu in functie de model
if nc == 0
    titlu = sprintf('Model AR(%d): Zgomot Colorat vs Model', na);
else
    titlu = sprintf('Model ARMA(%d, %d): Zgomot Colorat vs Model', na, nc);
end
title(titlu);

% Legenda
legend('Real (Masura)', 'Real (Validare)', 'Model (Simulare)', 'Model (Predictie)', 'Separator', 'Location', 'best');
xlabel('Timp (esantioane)'); ylabel('Amplitudine');
grid on; axis tight;

% Grafic 2: Eroarea de model
subplot(3, 2, 3);

plot(t_masura, e, 'y');
title(['Eroare masura: \lambda^2=' num2str(var_e_mas, '%.4f') ' SNR=' num2str(SNR_masura, '%.2f')]);
xlabel('Timp'); grid on; axis tight;

% Grafic 3: Eroarea de predictie
subplot(3, 2, 4);

plot(t_pred, epsilon_pred, 'm', 'LineWidth', 1.5);
title(['Eroare predictie: \lambda^2=' num2str(var_eps_pred, '%.4f') ' SNR=' num2str(SNR_pred, '%.2f')]);
xlabel('Timp'); grid on; axis tight;

% Grafic 4: Zoom Predictie + Tub de incredere
subplot(3, 2, [5 6]);

plot(t_pred, v_validare, 'c--o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
plot(t_pred, y_pred, 'r--o', 'LineWidth', 2, 'MarkerSize', 6);

% Calculare Tub de incredere (de tip 3sigma):
limita_superioara = y_pred + 3 * sqrt(lambda2);
limita_inferioara = y_pred - 3 * sqrt(lambda2);

plot(t_pred, limita_superioara, 'g', 'LineWidth', 1);
plot(t_pred, limita_inferioara, 'g', 'LineWidth', 1);

title('Zoom Predictie si Tub de Incredere');
legend('Real (Validare)', 'Predictie', 'Tub de incredere', 'Location', 'best');
xlabel('Timp'); ylabel('Amplitudine');
grid on; axis tight;