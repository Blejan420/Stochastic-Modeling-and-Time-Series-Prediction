% TAIPS_02C - Simulator comparativ AR vs ARMA

% Autori: Andrei BLEJAN
% Creat: 13 Ianuarie 2026

clear; close all; clc;

% --- Incarcare si pregatire date ---
n = input('Introduceti numarul seriei de timp (ex: 3 pentru Y03) (max 20): ');
nume_fisier = sprintf('Y%02d', n);

if exist([nume_fisier '.mat'], 'file')
    load([nume_fisier '.mat']);
else
    error('Fisierul %s.mat nu a fost gasit.', nume_fisier);
end

% Selectie canal (daca exista mai multe)
nr_canale = size(Y.y, 2);
if nr_canale > 1
    fprintf('Aceasta serie are %d canale. ', nr_canale);
    ch = input(['Alegeti canalul dorit (1-' num2str(nr_canale) '): ']);
    if isempty(ch) || ch < 1 || ch > nr_canale, error('Canal invalid.'); end
else
    ch = 1;
end

y_total = Y.y(:, ch);
N = length(y_total);

% Selectie orizont de predictie
K = input(sprintf('Introduceti orizontul de predictie K (total %d date): ', N));

% K nu poate depasi 2/3*N deoarece ar rezulta 0 date de antrenament
% Ne asiguram ca avem minim 10 date de antrenament
if isempty(K) || K < 1 || K > (floor(2/3 * N) - 10)
    fprintf('================================================================================================\n');
    fprintf('K invalid. Se alege valoarea implicita K = 5\nK trebuie ales intre 1 si %d (minim 10 date rezervate pentru antrenament)\n',floor(2/3 * N)-10); 
    fprintf('================================================================================================\n');
    pause(3);
    K = 5;
end

% --- Separare Masura de Predictie ---
N_masura = N - K;
y_masura = y_total(1:N_masura); 
y_extrapolare = y_total(N_masura+1:end);

% --- Separare Antrenare de Validare Interna ---
K_uitat = ceil(K/2); 
N_antrenare = N_masura - K_uitat;

y_antrenare = y_total(1:N_antrenare);
y_validare = y_total(N_antrenare+1 : N_masura);

% Vectori timp pentru validare
t_trend_valid = (N_antrenare : N_antrenare + K_uitat - 1)';
t_seasonal_valid  = (N_antrenare + 1 : N_antrenare + K_uitat)'; 

% Parametri cautare
na_max = floor(N/9);

% --- AR optimal (STRATEGIA 2) ---

PQ_opt_2 = -inf;
p_opt_2 = 0;
na_opt_2 = 1;

% Initializare matrice cu valorile PQ pentru Grafic 2
P_vector = 0:5;
% Ordonam crescator valorile aleatoare pentru claritatea graficului
if na_max <= 24, na_valori_2 = 1:na_max; else, na_valori_2 = sort(randperm(na_max, 24)); end
PQ_matrice = -inf(length(P_vector), length(na_valori_2));

for i_p = 1:length(P_vector)
    p = P_vector(i_p);

    % Calculare model determinist pe antrenare
    [ysta_2, YT_antrenare_2, theta_trend_2] = trend(y_antrenare, p);
    [YS_antrenare_2, v_2, theta_seas_2, P_vec_2] = seasonal(ysta_2);
    
    % Calculare model determinist pe validare
    YT_opt_2 = zeros(K_uitat, 1);
    for i=1:length(theta_trend_2)
        YT_opt_2 = YT_opt_2 + theta_trend_2(i)*(t_trend_valid.^(i-1)); 
    end
    
    YS_opt_2 = zeros(K_uitat, 1);
    if ~isempty(P_vec_2)
        for k = 1 : length(P_vec_2)
            omega = 2*pi/P_vec_2(k);
            YS_opt_2 = YS_opt_2 + (theta_seas_2(2*k-1)*sin(omega*t_seasonal_valid) + theta_seas_2(2*k)*cos(omega*t_seasonal_valid));
        end
    end
    
    Y_det_opt_2 = YT_opt_2 + YS_opt_2;
    
    for i_na = 1:length(na_valori_2)
        na = na_valori_2(i_na);

        [y_AR_2, ~, ~, th_AR, ~] = stochastic(v_2, na);
        
        y_AR_predictie_2 = zeros(K_uitat, 1);
        v_aux = v_2;
        ordin = length(th_AR) - 1;
        
        for k = 1 : K_uitat
            esantioane = v_aux(end:-1:end-ordin+1); 
            valoare_predictata = - (th_AR(2:end).' * esantioane);
            y_AR_predictie_2(k) = valoare_predictata;
            v_aux = [v_aux; valoare_predictata];
        end
        
        % Asamblare model final pe validare
        y_model_validare_2 = Y_det_opt_2 + y_AR_predictie_2;
        
        % Asamblare model final pe antrenare
        y_model_antrenare_2 = (YT_antrenare_2 + YS_antrenare_2) + y_AR_2;
        
        % Protectie instabilitate
        if any(isnan(y_model_antrenare_2)) || any(isinf(y_model_antrenare_2)) || ...
           any(isnan(y_model_validare_2)) || any(isinf(y_model_validare_2))
            pq = -1;
        else
            id_m = iddata(y_antrenare, [], 1);
            id_v = iddata(y_validare, [], 1);
            id_s = iddata(y_model_antrenare_2, [], 1);
            id_p = iddata(y_model_validare_2, [], 1);
            
            [pq, ~, ~] = pred_qual(id_m, id_v, id_s, id_p, th_AR);
        end

        % Salvare in matrice pentru Grafic 2
        PQ_matrice(i_p, i_na) = pq;

        if pq > PQ_opt_2
            PQ_opt_2 = pq;
            p_opt_2 = p;
            na_opt_2 = na;
        end
    end
end

% --- ARMA optimal (STRATEGIA 2) ---
fprintf('<TAIPS_02C>: Se cauta ARMA optimal. Procesul poate dura un timp! \n');
PARMA_opt_val = -inf;

% Intervalul de cautare
P_vector_parma = 0:3;
nr_perechi = 36;

% Structura cell pentru a salva datele celor 4 suprafete (pentru grafice)
rezultate_PARMA = cell(length(P_vector_parma), 1);

for i_p = 1:length(P_vector_parma)
    p = P_vector_parma(i_p);
    
    % Antrenare model determinist
    [ysta_parma, YT_antr_parma, th_trend_parma] = trend(y_antrenare, p);
    [YS_antr_parma, v_parma, th_seas_parma, P_vec_parma] = seasonal(ysta_parma);
    
    % Validare model determinist
    YT_valid_parma = zeros(K_uitat, 1);
    for i=1:length(th_trend_parma)
        YT_valid_parma = YT_valid_parma + th_trend_parma(i)*(t_trend_valid.^(i-1)); 
    end
    
    YS_valid_parma = zeros(K_uitat, 1);
    if ~isempty(P_vec_parma)
        for k = 1:length(P_vec_parma)
            omega = 2*pi/P_vec_parma(k);
            YS_valid_parma = YS_valid_parma + (th_seas_parma(2*k-1)*sin(omega*t_seasonal_valid) + th_seas_parma(2*k)*cos(omega*t_seasonal_valid));
        end
    end

    % Model determinist complet
    Y_det_antr_parma = YT_antr_parma + YS_antr_parma;
    Y_det_valid_parma = YT_valid_parma + YS_valid_parma;
    
    % --- Generare perechi ---

    % Se foloseste randi pentru generarea perechilor [na,nc] deoarece randperm 
    % nu permite refolosirea aceleasi valori deja alese pentru o pereche noua 
    % (util la generarea unui singur vector de valori, nu la generarea a perechi de vectori)
    % Pentru nc, se dau valori incepand de la 1, pentru a evita cazul in 
    % care nc = 0, unde stochastic ar calcula un model AR
    valori_na = randi([0, na_max], nr_perechi, 1);
    valori_nc = randi([1, na_max], nr_perechi, 1);
    perechi = [valori_na, valori_nc];

    % Matrice pt salvare: [na, nc, PQ]
    date_PARMA = zeros(nr_perechi, 3);
    
    % Construire model nedeterminist ARMA
    for k = 1:nr_perechi
        na = perechi(k, 1);
        nc = perechi(k, 2);
        
        % Apelare stochastic
        try
            [y_arma_antr, e, ~, ~, theta_ARMA] = stochastic(v_parma, na, nc); % nalpha e calculat implicit in stochastic
            
            % Predictie ARMA
            y_arma_valid = zeros(K_uitat, 1);
            
            % Extragere polinoame din theta_ARMA
            A_arma = theta_ARMA(:, 1); 
            C_arma = theta_ARMA(:, 2);
    
            % Vectorii de valori pe tot orizontul (antrenare + validare)
            v_total = [v_parma; zeros(K_uitat, 1)];
            e_total = [e; zeros(K_uitat, 1)];
            
            % Predictie pas cu pas
            for idx = 1:K_uitat
                 n = length(v_parma) + idx;
                 
                 % Termen AR: - sumi(a_i * v(t-i))
                 termen_AR = 0;
                 % A_arma(1) = 1. Coeficientii incep de la index 2.
                 % Pentru a evita erorile de indexare care pot fi cauzate de
                 % lungimea polinomului A returnat de stochastic (deoarece la
                 % final se umple vectorul cu zerouri), indexarea se parcurge
                 % pana la minimul dintre na si lungimea lui A
                 for ia = 1:min(na, length(A_arma)-1)
                     if n-ia > 0
                        termen_AR = termen_AR - A_arma(ia+1) * v_total(n-ia);
                     end
                 end
                 
                 % Termen MA: + sumj(c_j * e(t-j))
                 termen_MA = 0;
                 for ic = 1:min(nc, length(C_arma)-1)
                     if n-ic > 0
                        termen_MA = termen_MA + C_arma(ic+1) * e_total(n-ic);
                     end
                 end
    
                 y_arma_valid(idx) = termen_AR + termen_MA; % se populeaza valorile predictate
                 v_total(n) = termen_AR + termen_MA; % se populeaza feedback-ul pentru predictii viitoare
            end
            
            % Evaluare PQ
            y_mod_valid = Y_det_valid_parma + y_arma_valid;
            y_mod_antr = Y_det_antr_parma + y_arma_antr;
            
            % Protectie la instabilitate (cand modelul diverge)
            if any(isnan(y_mod_antr)) || any(isinf(y_mod_antr)) || ...
               any(isnan(y_mod_valid)) || any(isinf(y_mod_valid))
                 pq = -1;
            else
                 id_m = iddata(y_antrenare,[],1);
                 id_v = iddata(y_validare,[],1);
                 id_sa = iddata(y_mod_antr,[],1);
                 id_sv = iddata(y_mod_valid,[],1);
                 
                 % Apelare pred_qual
                 [pq, ~, ~] = pred_qual(id_m, id_v, id_sa, id_sv, theta_ARMA);
            end
        catch % cazul in care na+nc (alese aleator) > length(v_parma) -> eroare in functia armax
            pq = -1;
        end
        
        % Salvare date pentru grafic
        date_PARMA(k, :) = [na, nc, pq];
        
        % Retinere optim global
        if pq > PARMA_opt_val
            PARMA_opt_val = pq;
            parma_p = p; parma_na = na; parma_nc = nc;
        end
    end
    
    % Salvare rezultate pentru gradul p curent
    rezultate_PARMA{i_p} = date_PARMA;
end

% --- Afisare grafica ---
axa_x = Y.TimeUnit;
axa_y = Y.OutputUnit{ch};
if size(Y.y, 2) == 1
    fig1 = figure('Name', sprintf('Analiza Model %s | Pagina 1', nume_fisier), 'NumberTitle', 'off'); % pagina 1
    fig2 = figure('Name', sprintf('Analiza Model %s | Pagina 2', nume_fisier), 'NumberTitle', 'off'); % pagina 2
    fig3 = figure('Name', sprintf('Analiza Model %s | Pagina 3', nume_fisier), 'NumberTitle', 'off'); % pagina 3
    fig4 = figure('Name', sprintf('Analiza Model %s | Pagina 4', nume_fisier), 'NumberTitle', 'off'); % pagina 4
    nume_serie = Y.Notes;
else 
    fig1 = figure('Name', sprintf('Analiza Model %s, Canal %d | Pagina 1', nume_fisier, ch), 'NumberTitle', 'off');
    fig2 = figure('Name', sprintf('Analiza Model %s, Canal %d | Pagina 2', nume_fisier, ch), 'NumberTitle', 'off');
    fig3 = figure('Name', sprintf('Analiza Model %s, Canal %d | Pagina 3', nume_fisier, ch), 'NumberTitle', 'off');
    fig4 = figure('Name', sprintf('Analiza Model %s, Canal %d | Pagina 4', nume_fisier, ch), 'NumberTitle', 'off');
    nume_serie = Y.OutputName{ch};
end

t_masura = 0:N_masura-1;
t_extrapolare = N_masura:N-1;

% --- Reconstructie PAR ---

% Componenta determinista
[ysta_par, YT_par_m, th_tr_par] = trend(y_masura, p_opt_2);
[YS_par_mas, v_par, th_ss_par, P_vec_par] = seasonal(ysta_par);

% Extrapolare
t_trend_extra = (N_masura : N - 1)';
t_seas_extra  = (N_masura + 1 : N)';

YT_par_extra = zeros(K, 1);
for i=1:length(th_tr_par), YT_par_extra = YT_par_extra + th_tr_par(i)*(t_trend_extra.^(i-1)); end

YS_par_e = zeros(K, 1);
if ~isempty(P_vec_par)
    for k = 1:length(P_vec_par)
        om = 2*pi/P_vec_par(k);
        YS_par_e = YS_par_e + (th_ss_par(2*k-1)*sin(om*t_seas_extra) + th_ss_par(2*k)*cos(om*t_seas_extra));
    end
end
Y_det_par_m = YT_par_m + YS_par_mas;
Y_det_par_e = YT_par_extra + YS_par_e;

% Componenta nedeterminista
[y_nedet_par_mas, e_par, lambda2_par, th_nedet_par] = stochastic(v_par, na_opt_2, 0);

% Predictie
y_nedet_par_extra = zeros(K, 1);
v_aux = v_par;
ordin = length(th_nedet_par) - 1;
for k = 1:K
    esantioane = v_aux(end:-1:end-ordin+1);
    valoarea_predicata = -(th_nedet_par(2:end).' * esantioane);
    y_nedet_par_extra(k) = valoarea_predicata;
    v_aux = [v_aux; valoarea_predicata];
end

% Model complet
Y_tot_par_mas = Y_det_par_m + y_nedet_par_mas;
Y_tot_par_extra = Y_det_par_e + y_nedet_par_extra;

% Calcul PQ final
id_m = iddata(y_masura,[],1); id_e = iddata(y_extrapolare,[],1);
id_sim = iddata(Y_tot_par_mas,[],1); id_pred = iddata(Y_tot_par_extra,[],1);
[PQ_final_par, SNR_N_par, SNR_K_par] = pred_qual(id_m, id_e, id_sim, id_pred, th_nedet_par);

% --- Afisare consola PAR ---
fprintf('\n');
fprintf('======================================================\n');
fprintf('REZULTATE PENTRU: PAR (AR Optimal)\n');
fprintf('Parametrii optimi: p = %d, na = %d\n', p_opt_2, na_opt_2);
fprintf('PQ Final: %.2f%% | SNR Masura: %.2f | SNR Pred: %.2f\n', PQ_final_par, SNR_N_par, SNR_K_par);
fprintf('======================================================\n');
fprintf('Coeficientii polinomului tendinta (p=%d):\n', p_opt_2);
disp(th_tr_par');

if ~isempty(P_vec_par)
    fprintf('Componenta sezoniera:\n');
    fprintf('   | %-10s | %-15s | %-15s |\n', 'Perioada', 'Ak (Sin)', 'Bk (Cos)');
    fprintf('   ------------------------------------------------------\n');
    for k = 1:length(P_vec_par)
        ak_val = th_ss_par(2*k-1);
        bk_val = th_ss_par(2*k);
        fprintf('   | %-10d | %-15.4f | %-15.4f |\n', P_vec_par(k), ak_val, bk_val);
    end
else
    fprintf('Fara componenta sezoniera detectata.\n');
end

% --- Grafice PAR ---
figure(fig1);

% Grafic 1: Seria de Timp + Media    
subplot(2,2,1);
plot(0:N-1, y_total, 'b', 'LineWidth', 1); hold on;
line([0 N-1], [mean(y_total) mean(y_total)], 'Color', 'y', 'LineWidth', 1.5);
title(['Seria ' nume_fisier ': ' nume_serie ' | Media = ' num2str(mean(y_total), '%.2f')]);
xlabel(axa_x); ylabel(axa_y);
legend('Seria de timp', 'Media', 'Location', 'best');
grid on; axis tight;

% Grafic 2: Suprafata valorilor PQ
subplot(2,2,2);
[x_grid, y_grid] = meshgrid(na_valori_2, P_vector);
surf(x_grid, y_grid, PQ_matrice);
shading interp; hold on;
plot3(na_opt_2, p_opt_2, PQ_opt_2, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
title(['Suprafata PQ (Maxim: p=' num2str(p_opt_2) ', na=' num2str(na_opt_2) ', PQ=' num2str(PQ_opt_2, '%.2f') '%)']);
xlabel('Ordin AR (na)'); ylabel('Grad trend (p)'); zlabel('PQ [%]');
view(-30, 30); grid on; axis tight; colorbar;

% Grafic 3: Model determinist optim suprapus
subplot(2,2,3);
plot(t_masura, y_masura, 'b', 'LineWidth', 1); hold on;
plot(t_extrapolare, y_extrapolare, 'c', 'LineWidth', 1); 
plot(t_masura, Y_det_par_m, 'g', 'LineWidth', 1.2); 
plot(t_extrapolare, Y_det_par_e, 'r', 'LineWidth', 1.2);
line([N_masura-0.5 N_masura-0.5], ylim, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1);
title(['PAR: Model determinist optim (p=' num2str(p_opt_2) ')']);
xlabel(axa_x); ylabel(axa_y);
legend('Real (Masura)', 'Real (Extrapolare)', 'Model (Simulare)', 'Model (Predictie)', 'Location', 'best');
grid on; axis tight;

% Grafic 4: Model complet suprapus
subplot(2,2,4);
plot(t_masura, y_masura, 'b', 'LineWidth', 1); hold on;
plot(t_extrapolare, y_extrapolare, 'c', 'LineWidth', 1); 
plot(t_masura, Y_tot_par_mas, 'g', 'LineWidth', 1); 
plot(t_extrapolare, Y_tot_par_extra, 'r', 'LineWidth', 1.5);
line([N_masura-0.5 N_masura-0.5], ylim, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1);
title(['PAR: Model complet (p=' num2str(p_opt_2) ', na=' num2str(na_opt_2) ')']);
xlabel(axa_x); ylabel(axa_y);
legend('Real (Masura)', 'Real (Extrapolare)', 'Model (Simulare)', 'Model (Predictie)', 'Location', 'best');
grid on; axis tight;

figure(fig2);

% Calcul eroare predictie
epsilon_par = y_extrapolare - Y_tot_par_extra;
var_eps = var(epsilon_par, 1);

% Grafic 5: Eroarea de masura
subplot(3,1,1); 
plot(t_masura, e_par, 'y');
title(['PAR: Eroare masura | \lambda^2=' num2str(lambda2_par, '%.4f') ' | SNR_N=' num2str(SNR_N_par, '%.2f')]);
xlabel(axa_x); ylabel(axa_y); 
grid on; axis tight;

% Grafic 6: Eroarea de predictie
subplot(3,1,2);
plot(t_extrapolare, epsilon_par, 'm', 'LineWidth', 1.5);
title(['PAR: Eroare predictie | \lambda^2=' num2str(var_eps, '%.4f') ' | SNR_K=' num2str(SNR_K_par, '%.2f')]);
xlabel(axa_x); ylabel(axa_y);  
grid on; axis tight;

% Grafic 7: Zoom pe predictie + Tub de incredere
subplot(3,1,3);
% Calcul Tub de incredere (tip 3sigma)
tub_sus = Y_tot_par_extra + 3*sqrt(lambda2_par);
tub_jos = Y_tot_par_extra - 3*sqrt(lambda2_par);

plot(t_extrapolare, y_extrapolare, 'c--o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
plot(t_extrapolare, Y_tot_par_extra, 'r--x', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(t_extrapolare, tub_sus, 'g', 'LineWidth', 1);
plot(t_extrapolare, tub_jos, 'g', 'LineWidth', 1);

title(sprintf('PAR: Zoom Predictie cu Tub de Incredere (\\pm3\\sigma) | PQ_{optim} = %.2f%%', PQ_final_par));
xlabel(axa_x); ylabel(axa_y);
legend('Real', 'Predictie', 'Tub Incredere', 'Location', 'best');
grid on; axis tight;

% --- Reconstructie PARMA ---

% Componenta determinista
[ysta_parma, YT_parma_mas, th_tr_parma] = trend(y_masura, parma_p);
[YS_parma_mas, v_parma, th_ss_parma, P_vec_parma] = seasonal(ysta_parma);

% Extrapolare
YT_parma_extra = zeros(K, 1);
for i=1:length(th_tr_parma), YT_parma_extra = YT_parma_extra + th_tr_parma(i)*(t_trend_extra.^(i-1)); end

YS_parma_e = zeros(K, 1);
if ~isempty(P_vec_parma)
    for k = 1:length(P_vec_parma)
        om = 2*pi/P_vec_parma(k);
        YS_parma_e = YS_parma_e + (th_ss_parma(2*k-1)*sin(om*t_seas_extra) + th_ss_parma(2*k)*cos(om*t_seas_extra));
    end
end
Y_det_parma_mas = YT_parma_mas + YS_parma_mas;
Y_det_parma_extra = YT_parma_extra + YS_parma_e;

% Componenta nedeterminista
[y_nedet_parma_mas, e_parma, lambda2_parma, ~, theta_ARMA_opt] = stochastic(v_parma, parma_na, parma_nc);

% Predictie 
y_nedet_parma_extra = zeros(K, 1);
A_opt = theta_ARMA_opt(:, 1);
C_opt = theta_ARMA_opt(:, 2);

v_total = [v_parma; zeros(K, 1)];
e_total = [e_parma; zeros(K, 1)];

for idx = 1:K
     n = length(v_parma) + idx;
     
     % Termen AR
     termen_AR = 0;
     for ia = 1:min(parma_na, length(A_opt)-1)
         if n-ia > 0, termen_AR = termen_AR - A_opt(ia+1) * v_total(n-ia); end
     end
     
     % Termen MA
     termen_MA = 0;
     for ic = 1:min(parma_nc, length(C_opt)-1)
         if n-ic > 0, termen_MA = termen_MA + C_opt(ic+1) * e_total(n-ic); end
     end

     y_nedet_parma_extra(idx) = termen_AR + termen_MA;
     v_total(n) = termen_AR + termen_MA;
end

% Model complet PARMA
Y_parma_mas = Y_det_parma_mas + y_nedet_parma_mas;
Y_parma_extra = Y_det_parma_extra + y_nedet_parma_extra;

% Calcul PQ
id_sim_parma = iddata(Y_parma_mas,[],1); 
id_pred_parma = iddata(Y_parma_extra,[],1);

[PQ_final_parma, SNR_N_parma, SNR_K_parma] = pred_qual(id_m, id_e, id_sim_parma, id_pred_parma, theta_ARMA_opt);

% --- Afisare consola PARMA ---
fprintf('\n');
fprintf('======================================================\n');
fprintf('REZULTATE PENTRU: PARMA (ARMA Optimal)\n');
fprintf('Parametrii optimi: p = %d, na = %d, nc = %d\n', parma_p, parma_na, parma_nc);
fprintf('PQ Final: %.2f%% | SNR Masura: %.2f | SNR Pred: %.2f\n', PQ_final_parma, SNR_N_parma, SNR_K_parma);
fprintf('======================================================\n');
fprintf('Coeficientii polinomului tendinta (p=%d):\n', parma_p);
disp(th_tr_parma');

if ~isempty(P_vec_parma)
    fprintf('Componenta sezoniera:\n');
    fprintf('   | %-10s | %-15s | %-15s |\n', 'Perioada', 'Ak (Sin)', 'Bk (Cos)');
    fprintf('   ------------------------------------------------------\n');
    for k = 1:length(P_vec_parma)
        ak_val = th_ss_parma(2*k-1);
        bk_val = th_ss_parma(2*k);
        fprintf('   | %-10d | %-15.4f | %-15.4f |\n', P_vec_parma(k), ak_val, bk_val);
    end
else
    fprintf('Fara componenta sezoniera detectata.\n');
end

% --- Afisare grafica PARMA ---

figure(fig3);

for i_p = 1:4
    subplot(2,2,i_p);
    data = rezultate_PARMA{i_p}; % Matricea [na, nc, pq]
    
    % Functia surf genereaza erori cand o apelam cu grila incompleta in care nu
    % sunt definite toate combinatiile [na,nc,pq] (din cauza randi pe perechile 
    % [na,nc], deci apelam la interpolare pentru generarea suprafetelor

    % Filtrare puncte invalide (pq = -1)
    idx_valid = data(:,3) ~= -1;
    pct_valid = data(idx_valid, :);
    
    x = pct_valid(:,1); % na
    y = pct_valid(:,2); % nc
    z = pct_valid(:,3); % pq
    
    % Creare grila pentru interpolare (100x100 puncte)
    [X_coord, Y_coord] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
    
    % --- Interpolare ---

    % Construim o functie F care "invata" forma suprafetei trecand prin punctele valabile
    % ('natural' = suprafata neteda; 'none' = fara valori interpolate in afara ariei)
    F = scatteredInterpolant(x, y, z, 'natural', 'none');
    
    Z_coord = F(X_coord, Y_coord);
    
    % Apelare surf
    surf(X_coord, Y_coord, Z_coord); 
    shading interp;
    hold on;
    
    % Marcare maxim
    [max_loc, id_max] = max(z);
    plot3(x(id_max), y(id_max), max_loc, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    
    title(['Suprafata PQ, p=' num2str(P_vector_parma(i_p)) ' | Maxim: ' num2str(max_loc, '%.1f') '% (na=' num2str(x(id_max)) ', nc=' num2str(y(id_max)) ')']);
    xlabel('na'); ylabel('nc'); zlabel('PQ [%]'); 
    view(-30, 30); grid on; axis tight; colorbar;
end

figure(fig4);

% Grafic 1: Seria de timp + Media
subplot(3,2,1);
plot(0:N-1, y_total, 'b', 'LineWidth', 1); hold on;
line([0 N-1], [mean(y_total) mean(y_total)], 'Color', 'y', 'LineWidth', 1.5);
title(['Seria ' nume_fisier ': ' nume_serie ' | Media = ' num2str(mean(y_total), '%.2f')]);
xlabel(axa_x); ylabel(axa_y);
legend('Seria', 'Media', 'Location', 'best');
grid on; axis tight;

% Grafic 2: Model determinist suprapus
subplot(3,2,2);
plot(t_masura, y_masura, 'b', 'LineWidth', 1); hold on;
plot(t_extrapolare, y_extrapolare, 'c', 'LineWidth', 1); 
plot(t_masura, Y_det_parma_mas, 'g', 'LineWidth', 1.2); 
plot(t_extrapolare, Y_det_parma_extra, 'r', 'LineWidth', 1.2);
line([N_masura-0.5 N_masura-0.5], ylim, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1);
title(['PARMA: Model determinist optim (p=' num2str(parma_p) ')']);
xlabel(axa_x); ylabel(axa_y);
legend('Masura', 'Extra', 'Model M', 'Model E', 'Location', 'best');
grid on; axis tight;

% Grafic 3: Model complet suprapus
subplot(3,2,3);
plot(t_masura, y_masura, 'b', 'LineWidth', 1); hold on;
plot(t_extrapolare, y_extrapolare, 'c', 'LineWidth', 1); 
plot(t_masura, Y_parma_mas, 'g', 'LineWidth', 1); 
plot(t_extrapolare, Y_parma_extra, 'r', 'LineWidth', 1.5);
line([N_masura-0.5 N_masura-0.5], ylim, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1);
title(['PARMA: Model complet (p=' num2str(parma_p) ', na=' num2str(parma_na) ', nc=' num2str(parma_nc) ')']);
xlabel(axa_x); ylabel(axa_y);
legend('Masura', 'Extra', 'Model M', 'Model E', 'Location', 'best');
grid on; axis tight;

% Calcul erori
epsilon_parma = y_extrapolare - Y_parma_extra;
var_eps_parma = var(epsilon_parma, 1);

% Grafic 4: Eroarea de masura
subplot(3,2,4);
plot(t_masura, e_parma, 'y');
title(['PARMA: Eroare masura | \lambda^2=' num2str(lambda2_parma, '%.4f') ' | SNR_N=' num2str(SNR_N_parma,'%.2f')]);
xlabel(axa_x); ylabel('Eroare'); 
grid on; axis tight;

% Grafic 5: Eroarea de predictie
subplot(3,2,5);
plot(t_extrapolare, epsilon_parma, 'm', 'LineWidth', 1.5);
title(['PARMA: Eroare predictie |\lambda^2=' num2str(var_eps_parma, '%.4f') ' | SNR_K=' num2str(SNR_K_parma,'%.2f')]);
xlabel(axa_x); ylabel('Eroare'); 
grid on; axis tight;

% Grafic 6: Zoom predictie + Tub de incredere
subplot(3,2,6);
tub_sus_parma = Y_parma_extra + 3*sqrt(lambda2_parma);
tub_jos_parma = Y_parma_extra - 3*sqrt(lambda2_parma);

plot(t_extrapolare, y_extrapolare, 'c--o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
plot(t_extrapolare, Y_parma_extra, 'r--x', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(t_extrapolare, tub_sus_parma, 'g', 'LineWidth', 1);
plot(t_extrapolare, tub_jos_parma, 'g', 'LineWidth', 1);

title(sprintf('PARMA: Zoom Predictie + Tub (\\pm3\\sigma) | PQ_{optim} = %.2f%%', PQ_final_parma));
xlabel(axa_x); ylabel(axa_y);
legend('Real', 'Predictie', 'Tub Incredere', 'Location', 'best');
grid on; axis tight;