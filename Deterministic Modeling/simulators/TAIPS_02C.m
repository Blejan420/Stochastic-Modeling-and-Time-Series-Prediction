% TAIPS_02C - Simulator pentru determinarea automata a gradului optim p
clear; clc; close all;

% --- Incarcare si pregatire date ---
n = input('Introduceti numarul seriei de timp (ex: 3 pentru Y3) (max 20): ');
nume_fisier = sprintf('Y%02d', n);
cale = fullfile('TAIPS/TAIPSMAT', [nume_fisier '.mat']);

if ~exist(cale, 'file')
    if exist([nume_fisier '.mat'], 'file')
        load([nume_fisier '.mat']);
    else
        error('Fisierul %s nu a fost gasit.', nume_fisier);
    end
else
    load(cale); 
end

% Selectie canal (pentru obiectele Y16-Y20)
nr_canale = size(Y.y, 2);
if nr_canale > 1
    fprintf('Aceasta serie are %d canale. ', nr_canale);
    ch = input(['Alegeti canalul dorit (1-' num2str(nr_canale) '): ']);
    if isempty(ch) || ch < 1 || ch > nr_canale, error('Canal invalid.'); end
else
    ch = 1;
end

% Extragere date pentru canalul ales
y_total = Y.y(:, ch);
N = length(y_total);

% Selectie orizont de extrapolare
K = input(sprintf('Introduceti orizontul de extrapolare K (total %d date in serie): ', N));

% K nu poate depasi 2/3*N deoarece ar rezulta 0 date de antrenament
% Ne asiguram ca avem minim 10 date de antrenament
if isempty(K) || K < 1 || K > (floor(2/3 * N) - 10)
    fprintf('================================================================================================\n');
    fprintf('K invalid. Se alege valoarea implicita K = 5\nK trebuie ales intre 1 si %d (minim 10 date rezervate pentru antrenament)\n',floor(2/3 * N)-10); 
    fprintf('================================================================================================\n');
    pause(3);
    K = 5;
end

% Separare masura si extrapolare
N_masura = N - K;
y_masura = y_total(1:N_masura);
y_extra  = y_total(N_masura+1:end);

% --- Determinarea gradului optim ---

% Definim "datele uitate" pentru validare interna
K_uitat = ceil(K / 2); 
N_antrenare = N_masura - K_uitat;

% Datele pentru antrenament
y_antrenare = y_total(1:N_antrenare);
% Datele pentru validare
y_validare = y_total(N_antrenare+1 : N_masura); 

P_max = 10;
% Vector pentru stocarea EQ pentru fiecare grad p
EQ_val = zeros(P_max + 1, 1); % index 1 pt p=0

for p = 0:P_max
    % Identificare model pe setul de antrenament
    [ysta_antrenare, YT_antrenare, theta_trend] = trend(y_antrenare, p);
    [ys_antrenare, ~, theta_sez, P_vec] = seasonal(ysta_antrenare);
    
    % Model complet pe antrenament
    Y_model_antrenare = YT_antrenare + ys_antrenare;
    
    % Extrapolare model pe setul de validare (datele "uitate")
    
    % Extrapolare trend
    t_trend_validare = (N_antrenare : N_antrenare + K_uitat - 1)';
    YT_validare = zeros(K_uitat, 1);
    for i = 1:p+1
        YT_validare = YT_validare + theta_trend(i) * (t_trend_validare.^(i-1));
    end
    
    % Extrapolare seasonal
    t_seasonal_validare = (N_antrenare + 1 : N_antrenare + K_uitat)'; % indexarea difera intre rutinele trend si seasonal
    ys_validare = zeros(K_uitat, 1);
    if ~isempty(P_vec)
        for k = 1:length(P_vec)
            omega = 2 * pi / P_vec(k);
            ak = theta_sez(2*k-1);
            bk = theta_sez(2*k);
            ys_validare = ys_validare + (ak * sin(omega * t_seasonal_validare) + bk * cos(omega * t_seasonal_validare));
        end
    end
    
    Y_model_validare = YT_validare + ys_validare;
    
    % Calculare EQ prin rutina extra_qual
    % Construire obiecte iddata (doar OutputData pentru simplitate)
    id_y_a = iddata(y_antrenare, [], 1);
    id_y_v = iddata(y_validare, [], 1);
    id_ym_a = iddata(Y_model_antrenare, [], 1);
    id_ym_v = iddata(Y_model_validare, [], 1);
    
    % Apelare rutina
    [eq, ~, ~] = extra_qual(id_y_a, id_y_v, id_ym_a, id_ym_v);
    
    EQ_val(p+1) = eq;
end

% Gasirea maximului
[max_EQ, idx_max] = max(EQ_val);
p_0 = idx_max - 1; % convertim index (1..11) la grad (0..10)

% --- Construirea modelului determinist final ---

% Refacem calculele folosind y_masura (N-K puncte) si p_0

% Determinare model pe orizontul intreg masura
[ysta_m, YT_m, theta_trend_final] = trend(y_masura, p_0);
[ys_m, v_m, theta_seasonal_final, P_vec_final] = seasonal(ysta_m);

Y_m_total = YT_m + ys_m;

% Extrapolare pe orizontul real de extrapolare
t_trend_extra = (N_masura : N-1)';
YT_extra = zeros(K, 1);

for i = 1:p_0+1
    YT_extra = YT_extra + theta_trend_final(i) * (t_trend_extra.^(i-1));
end

t_seasonal_extra = (N_masura+1 : N)';
ys_extra = zeros(K, 1);

if ~isempty(P_vec_final)
    for k = 1:length(P_vec_final)
        omega = 2 * pi / P_vec_final(k);
        ak = theta_seasonal_final(2*k-1);
        bk = theta_seasonal_final(2*k);
        ys_extra = ys_extra + (ak * sin(omega * t_seasonal_extra) + bk * cos(omega * t_seasonal_extra));
    end
end

Y_extra_total = YT_extra + ys_extra;

% --- Evaluare EQ real pe ultimele K date ---

v_p = y_masura - Y_m_total;
epsilon_p = y_extra - Y_extra_total;

id_y_m_final = iddata(y_masura, [], 1);
id_y_e_final = iddata(y_extra, [], 1);
id_ym_m_final = iddata(Y_m_total, [], 1);
id_ym_e_final = iddata(Y_extra_total, [], 1);

[EQ_final, snr_n, snr_k] = extra_qual(id_y_m_final, id_y_e_final, id_ym_m_final, id_ym_e_final);

% --- Afisare Consola ---

fprintf('\n--- Rezultate %s (Canal %d) ---\n', nume_fisier, ch);
fprintf('Coeficientii polinomului tendinta (p=%d):\n', p_0);
disp(theta_trend_final');

if ~isempty(P_vec_final)
    fprintf('Componenta sezoniera:\n');
    fprintf('   | %-10s | %-15s | %-15s |\n', 'Perioada', 'Ak (Sin)', 'Bk (Cos)');
    fprintf('   ------------------------------------------------------\n');
    for k = 1:length(P_vec_final)
        fprintf('   | %-10d | %-15.4f | %-15.4f |\n', P_vec_final(k), theta_seasonal_final(2*k-1), theta_seasonal_final(2*k));
    end
else
    fprintf('Fara componenta sezoniera.\n');
end

% --- Grafice ---

if size(Y.y, 2) == 1
    figure('Name', ['Analiza Model ' nume_fisier], 'NumberTitle', 'off');
    nume_serie = Y.Notes;
else 
    figure('Name', ['Analiza Model ' nume_fisier ', Canal ' num2str(ch)], 'NumberTitle', 'off');
    nume_serie = Y.OutputName{ch};
end

t_axa = 0:N-1;
axa_x = Y.TimeUnit;
axa_y = Y.OutputUnit{ch};

i = 0;
if (p_0 > 0) || (~isempty(P_vec_final)), j = 5; else, j = 4; end

% Grafic 1: Seria completa + Media
subplot(j,2,[i+1 i+2]);
plot(t_axa, y_total, 'b', 'LineWidth', 1); hold on;
line([0 N-1], [mean(y_total) mean(y_total)], 'Color', 'y', 'LineWidth', 2);
title(['Seria ' nume_fisier ': ' nume_serie ' | Media = ' num2str(mean(y_total))]);
xlabel(axa_x); ylabel(axa_y);
legend('Seria de timp', 'Media', 'Location', 'best');
grid on; axis tight;

% Grafic 2: Valorile criteriului EQ
subplot(j,2,[i+3 i+4]);
stem(0:P_max+1, [EQ_val; EQ_final], 'filled', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;

plot(p_0, max_EQ, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
title(['Valorile EQ, unde p_{0}=' num2str(p_0) ' cu EQ_{0}=' num2str(max_EQ, '%.2f') '%']);
xlabel('grad polinom (p)');
ylabel('EQ [%]');
ylim([0 100]);
grid on;
% Etichetare
xticks(0:P_max+1); 
etichete = string(0:P_max);
etichete(end+1) = "real";
xticklabels(etichete);
legend('Valori EQ', 'Optim', 'Location', 'best');

if j == 5

    % Grafic 3: Model determinist optim suprapus
    subplot(j,2,[i+5 i+6]);

    plot(t_axa(1:N_masura), y_masura, 'b', 'LineWidth', 1); hold on;
    plot(t_axa(N_masura+1:end), y_extra, 'c', 'LineWidth', 1); 

    plot(t_axa(1:N_masura), Y_m_total, 'g', 'LineWidth', 2); 
    plot(t_axa(N_masura+1:end), Y_extra_total, 'r', 'LineWidth', 2); 

    yl = ylim;
    line([N_masura-0.5 N_masura-0.5], yl, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);

    if ~isempty(P_vec_final)
        title(['Model determinist optim (trend p=' num2str(p_0) ' + componenta sezoniera). ']);
    else
        title(['Model determinist optim (trend p=' num2str(p_0) ' fara componenta sezoniera). ']);
    end

    xlabel(axa_x);
    ylabel(axa_y);
    legend('Masura Real', 'Extrapolare Real', 'Model Determinist', 'Predictie', 'Separator', 'Location', 'best');
    grid on; axis tight;
    
    i = 2;
end

% Grafic 4: Eroarea de masura
subplot(j,2,i+5);
plot(0:N_masura-1, v_p, 'm');
title(['Eroare masura: \lambda^2=' num2str(var(v_p, 1), 3) ' SNR=' num2str(snr_n, 3)]);
xlabel(axa_x); grid on; axis tight;

% Grafic 5: Eroarea de extrapolare
subplot(j,2,i+6);
plot(N_masura:N-1, epsilon_p, 'r');
title(['Eroare extrapolare: \lambda^2=' num2str(var(epsilon_p, 1), 3) ' SNR=' num2str(snr_k, 3)]);
xlabel(axa_x); grid on; axis tight;

% Grafic 6: Zoom pe Extrapolare
subplot(j,2,[i+7 i+8]);
t_zoom = N_masura:N-1;
plot(t_zoom, y_extra, 'c', 'LineWidth', 1.5); hold on;
plot(t_zoom, Y_extra_total, 'r', 'LineWidth', 2);
title(['Zoom Extrapolare (EQ = ' num2str(EQ_final, '%.2f') '%)']);
legend('Real', 'Predictie', 'Location', 'best');
xlabel(axa_x); ylabel(axa_y);
grid on; axis tight;
