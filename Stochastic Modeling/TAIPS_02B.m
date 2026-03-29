% TAIPS_02B - Simulator pentru alegerea predictorului optimal (AR) prin
%             doua strategii distincte

% Autori: Andrei BLEJAN
% Creat: 10 Ianuarie 2026

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

% Parametri cautare
P_max = 11; 
na_max = floor(N/3);

% --- STRATEGIA 1 ---

% Maximizare EQ
EQ_maxim = -inf;
p_opt = 0;

% Vectori de timp pentru extrapolare in validare interna
t_trend_valid = (N_antrenare : N_antrenare + K_uitat - 1)';
t_seasonal_valid  = (N_antrenare + 1 : N_antrenare + K_uitat)'; 

for p = 0:P_max
    % Identificare trend pe antrenare
    [ysta, YT_antrenare, theta_trend] = trend(y_antrenare, p);
    
    % Identificare sezonalitate
    [ys_antrenare, ~, theta_seas, P_vec] = seasonal(ysta);
    
    % Extrapolare trend pe orizontul uitat
    YT_validare = zeros(K_uitat, 1);
    for i = 1:length(theta_trend)
        YT_validare = YT_validare + theta_trend(i) * (t_trend_valid.^(i-1));
    end
    
    % Extrapolare sezonalitate
    YS_validare = zeros(K_uitat, 1);
    if ~isempty(P_vec)
        for k = 1:length(P_vec)
            omega = 2 * pi / P_vec(k);
            ak = theta_seas(2*k-1);
            bk = theta_seas(2*k);
            YS_validare = YS_validare + (ak * sin(omega * t_seasonal_valid) + bk * cos(omega * t_seasonal_valid));
        end
    end
    
    % Model determinist total
    Y_det_antrenare = YT_antrenare + ys_antrenare;
    Y_det_validare = YT_validare + YS_validare;
    
    % Calcul EQ
    id_y_a = iddata(y_antrenare, [], 1);
    id_y_v = iddata(y_validare, [], 1);
    id_ym_a = iddata(Y_det_antrenare, [], 1);
    id_ym_v = iddata(Y_det_validare, [], 1);
    
    [eq, ~, ~] = extra_qual(id_y_a, id_y_v, id_ym_a, id_ym_v);
    
    if eq > EQ_maxim
        EQ_maxim = eq;
        p_opt = p;
    end
end

% Recalculare reziduu pentru p optim pe setul de antrenare
[ysta_opt, YT_antrenare_final, theta_opt] = trend(y_antrenare, p_opt);
[YS_antrenare_final, v, theta_seas_opt, P_vec_opt] = seasonal(ysta_opt);

% Reconstruire model determinist optim
YT_opt = zeros(K_uitat, 1);
YS_opt = zeros(K_uitat, 1);

for i=1:length(theta_opt) 
    YT_opt = YT_opt + theta_opt(i)*(t_trend_valid.^(i-1)); 
end

if ~isempty(P_vec_opt)
    for k=1:length(P_vec_opt)
        omega = 2*pi/P_vec_opt(k);
        YS_opt = YS_opt + (theta_seas_opt(2*k-1)*sin(omega*t_seasonal_valid) + theta_seas_opt(2*k)*cos(omega*t_seasonal_valid));
    end
end
Y_det_opt = YT_opt + YS_opt;

% Folosim functia randperm pentru a extrage 12 valori pseudo-aleatoare cu
% distributie uniforma in scopul maximizarii PQ

if na_max <= 12
    % Daca intervalul posibil (N/3) e mai mic de 12, luam toate valorile
    na_valori = 1:na_max; 
else
    na_valori = randperm(na_max, 12);
end

PQ_opt = -inf;
na_opt = 1;

for na = na_valori
    
    [y_AR, ~, ~, theta_AR, ~] = stochastic(v, na);

    y_AR_predictie = zeros(K_uitat, 1);
    v_aux = v;
    ordin = length(theta_AR) - 1;

    for k = 1:K_uitat
        esantioane = v_aux(end:-1:end-ordin+1); 
        valoare_predictata = - (theta_AR(2:end).' * esantioane); 
        y_AR_predictie(k) = valoare_predictata;
        v_aux = [v_aux; valoare_predictata]; 
    end
    
    % Asamblare model complet pe validare
    y_model_validare = Y_det_opt + y_AR_predictie;
    
    % Asamblare model pe antrenare
    y_model_antrenare = (YT_antrenare_final + YS_antrenare_final) + y_AR;
    
    % Protectie pentru instabilitate (cand AR calculat diverge, rezultand valori NaN in vectori)
    if any(isnan(y_model_antrenare)) || any(isinf(y_model_antrenare)) || ...
       any(isnan(y_model_validare)) || any(isinf(y_model_validare))
        pq = -1; % penalizare maxima
    else
        % Modelul nu diverge, deci se calculeaza PQ
        id_meas_real = iddata(y_antrenare, [], 1);
        id_valid_real = iddata(y_validare, [], 1);
        id_sim_train = iddata(y_model_antrenare, [], 1);
        id_pred_valid = iddata(y_model_validare, [], 1);
        
        [pq, ~, ~] = pred_qual(id_meas_real, id_valid_real, id_sim_train, id_pred_valid, theta_AR);
    end
    
    if pq > PQ_opt
        PQ_opt = pq;
        na_opt = na;
    end
end

% --- STRATEGIA 2 ---

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

% --- Afisare grafica ---
axa_x = Y.TimeUnit;
axa_y = Y.OutputUnit{ch};
if size(Y.y, 2) == 1
    fig1 = figure('Name', sprintf('Analiza Model %s | Pagina 1', nume_fisier), 'NumberTitle', 'off');
    nume_serie = Y.Notes;
else 
    fig1 = figure('Name', sprintf('Analiza Model %s, Canal %d | Pagina 1', nume_fisier, ch), 'NumberTitle', 'off');
    nume_serie = Y.OutputName{ch};
end

% Grafic 1: Seria de Timp + Media    
subplot(3,2,1);
plot(0:N-1, y_total, 'b', 'LineWidth', 1); hold on;
line([0 N-1], [mean(y_total) mean(y_total)], 'Color', 'y', 'LineWidth', 1.5);
title(['Seria ' nume_fisier ': ' nume_serie ' | Media = ' num2str(mean(y_total), '%.2f')]);
xlabel(axa_x); ylabel(axa_y);
legend('Seria de timp', 'Media', 'Location', 'best');
grid on; axis tight;

% Grafic 2: Suprafata valorilor PQ (pentru Strategia 2)
subplot(3,2,2);
% Construire grila de coordonate pentru surf
[x, y] = meshgrid(na_valori_2, P_vector);

% Desenare suprafata cu surf
surf(x, y, PQ_matrice);
shading interp;
hold on;

% Marcare maxim
plot3(na_opt_2, p_opt_2, PQ_opt_2, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');

title(['Suprafata PQ (Maxim: p=' num2str(p_opt_2) ', na=' num2str(na_opt_2) ', PQ=' num2str(PQ_opt_2, '%.2f') '%)']);
xlabel('Ordin AR (na)'); ylabel('Grad trend (p)'); zlabel('PQ [%]');
view(-30, 30); grid on; axis tight;
colorbar;

% Pagina 2 pentru afisarea graficelor
if size(Y.y, 2) == 1
    fig2 = figure('Name', sprintf('Analiza Model %s | Pagina 2', nume_fisier), 'NumberTitle', 'off');
    nume_serie = Y.Notes;
else 
    fig2 = figure('Name', sprintf('Analiza Model %s, Canal %d | Pagina 2', nume_fisier, ch), 'NumberTitle', 'off');
    nume_serie = Y.OutputName{ch};
end

% --- Generare grafice unice fiecarei strategii ---

% Folosim o matrice cu valorile optime gasite pentru a folosi o bucla
% pentru generarea graficelor 3-7 si afisarea la consola
valori_optime = [p_opt, na_opt; p_opt_2, na_opt_2];
titlu = {'Strategia 1', 'Strategia 2'};

for i = 1:2
    p_final = valori_optime(i, 1);
    na_final = valori_optime(i, 2);    

    % --- Indentificare model pe tot orizontul de masura ---
    % Se foloseste functia din MATLAB evalc pentru a stoca textul afisat in
    % consola de rutinele trend si seasonal intr-o variabila text la care
    % renuntam. Scopul este ca mesajele rutinelor sa nu interfereze cu
    % afisarile de interes din consola
    [~] = evalc('[ysta_final, YT_final, th_tr_f] = trend(y_masura, p_final)');
    [~] = evalc('[YS_final, v_final, th_ss_f, P_v_f] = seasonal(ysta_final)');
    [y_AR_final, e_final, ~, theta_AR_final, ~] = stochastic(v_final, na_final);
    
    % --- Predictie model pe orizontul de extrapolare ---
    t_trend_final = (N_masura : N_masura + K - 1)';
    YT_extra_final = zeros(K, 1);
    for j = 1:length(th_tr_f)
        YT_extra_final = YT_extra_final + th_tr_f(j)*(t_trend_final.^(j-1)); 
    end

    t_seas_final = (N_masura + 1 : N_masura + K)';
    YS_extra_final = zeros(K, 1);
    if ~isempty(P_v_f)
        for k=1:length(P_v_f)
            omega = 2*pi/P_v_f(k);
            YS_extra_final = YS_extra_final + (th_ss_f(2*k-1)*sin(omega*t_seas_final) + th_ss_f(2*k)*cos(omega*t_seas_final));
        end
    end
    
    % --- Predictie AR ---
    y_AR_extra_final = zeros(K, 1);
    v_aux = v_final;
    ordin = length(theta_AR_final) - 1;
    for k = 1:K
        esantioane = v_aux(end:-1:end-ordin+1);
        valoare_predictata = - (theta_AR_final(2:end).' * esantioane);        
        y_AR_extra_final(k) = valoare_predictata;
        v_aux = [v_aux; valoare_predictata];
    end
    
    % Asamblare
    Y_total_masura = YT_final + YS_final + y_AR_final;
    Y_total_extrapolare = YT_extra_final + YS_extra_final + y_AR_extra_final;
    
    % --- Calcul PQ_final ---
    id_m_f = iddata(y_masura, [], 1);
    id_e_f = iddata(y_extrapolare, [], 1);
    id_sim_f = iddata(Y_total_masura, [], 1);
    id_pred_f = iddata(Y_total_extrapolare, [], 1);
    
    [PQ_final, SNR_N, SNR_K] = pred_qual(id_m_f, id_e_f, id_sim_f, id_pred_f, theta_AR_final);
    
    % --- Afisare consola ---
    fprintf('\n');
    fprintf('======================================================\n');
    fprintf('REZULTATE PENTRU: %s\n', titlu{i});
    fprintf('Parametrii optimi: p = %d, na = %d\n', p_final, na_final);
    fprintf('PQ Final: %.2f%% | SNR Masura: %.2f | SNR Pred: %.2f\n', PQ_final, SNR_N, SNR_K);
    fprintf('======================================================\n');

    fprintf('Coeficientii polinomului tendinta (p=%d):\n', p_final);
    disp(th_tr_f');
    
    if ~isempty(P_v_f)
        fprintf('Componenta sezoniera:\n');
        fprintf('   | %-10s | %-15s | %-15s |\n', 'Perioada', 'Ak (Sin)', 'Bk (Cos)');
        fprintf('   ------------------------------------------------------\n');
        for k = 1:length(P_v_f)
            ak_val = th_ss_f(2*k-1);
            bk_val = th_ss_f(2*k);
            fprintf('   | %-10d | %-15.4f | %-15.4f |\n', P_v_f(k), ak_val, bk_val);
        end
    else
        fprintf('Fara componenta sezoniera detectata.\n');
    end

    figure(fig1);

    % Grafic 3: Model determinist optim suprapus (pentru fiecare strategie)
    subplot(3,2,1+2*i);

    t_masura = 0:N_masura-1;
    t_extrapolare = N_masura:N-1;
   
    % Seria de timp
    plot(t_masura, y_masura, 'b', 'LineWidth', 1); hold on;
    plot(t_extrapolare, y_extrapolare, 'c', 'LineWidth', 1); 
    
    % Modelul determinist optim
    plot(t_masura, YT_final + YS_final, 'g', 'LineWidth', 1); 
    plot(t_extrapolare, YT_extra_final + YS_extra_final, 'r', 'LineWidth', 1)

    line([N_masura-0.5 N_masura-0.5], ylim, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);

    title(sprintf('Model Determinist Optim (S%d): p=%d', i, p_final));
    xlabel(axa_x); ylabel(axa_y);
    legend('Real (Masura)', 'Real (Extrapolare)', 'Model (Simulare)', 'Model (Predictie)', 'Location', 'best');
    grid on; axis tight;

    % Grafic 4: Model complet suprapus (pentru fiecare strategie)
    subplot(3,2,2+2*i);
   
    % Seria de timp
    plot(t_masura, y_masura, 'b', 'LineWidth', 1); hold on;
    plot(t_extrapolare, y_extrapolare, 'c', 'LineWidth', 1); 
    
    % Modelul complet (determinist + nedeterminist)
    plot(t_masura, Y_total_masura, 'g', 'LineWidth', 1); 
    plot(t_extrapolare, Y_total_extrapolare, 'r', 'LineWidth', 1);

    line([N_masura-0.5 N_masura-0.5], ylim, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);

    title(sprintf('Model Complet (S%d): p=%d, na=%d', i, p_final, na_final));
    xlabel(axa_x); ylabel(axa_y);
    legend('Real (Masura)', 'Real (Extrapolare)', 'Model (Simulare)', 'Model (Predictie)', 'Location', 'best');
    grid on; axis tight;
    
    figure(fig2);
    
    % Calcul Dispersii 
    lambda2_e = var(e_final, 1); % zgomotul alb
    
    epsilon = y_extrapolare - Y_total_extrapolare; % eroarea de predictie
    lambda2_eps = var(epsilon, 1);
    
    % Grafic 5: Eroarea de masura
    subplot(3,2,i*2-1); 
    
    plot(t_masura, e_final, 'y');
    title([titlu{i} ' | Eroare masura: \lambda^2=' num2str(lambda2_e, '%.4f') ' SNR=' num2str(SNR_N, '%.2f')]);
    xlabel(axa_x); ylabel(axa_y); 
    grid on; axis tight;
    
    % Grafic 6: Eroarea de predictie
    subplot(3,2,i*2);
    
    plot(t_extrapolare, epsilon, 'm', 'LineWidth', 1.5);
    title([titlu{i} ' | Eroare predictie: \lambda^2=' num2str(lambda2_eps, '%.4f') ' SNR=' num2str(SNR_K, '%.2f')]);
    xlabel(axa_x); ylabel(axa_y); 
    grid on; axis tight;

    % Grafic 7: Zoom pe predictie + Tub de incredere
    subplot(3,2,4+i)

    % Calcul Tub de incredere (tip 3sigma)
    tub_sus = Y_total_extrapolare + 3*sqrt(lambda2_e);
    tub_jos = Y_total_extrapolare - 3*sqrt(lambda2_e);
    
    % Seria de timp pe extrapolare
    plot(t_extrapolare, y_extrapolare, 'c--o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
    
    % Predictie model
    plot(t_extrapolare, Y_total_extrapolare, 'r--o', 'LineWidth', 1.5, 'MarkerSize', 6);
    
    % Tub de Incredere
    plot(t_extrapolare, tub_sus, 'g', 'LineWidth', 1);
    plot(t_extrapolare, tub_jos, 'g', 'LineWidth', 1);

    title(sprintf('(S%d) Zoom Predictie cu Tub de Incredere (\\pm3\\sigma) | PQ_{optim} = %.2f%%', i, PQ_final));
    xlabel(axa_x); ylabel(axa_y);
    legend('Real', 'Predictie', 'Tub Incredere', 'Location', 'best');
    grid on; axis tight;

end