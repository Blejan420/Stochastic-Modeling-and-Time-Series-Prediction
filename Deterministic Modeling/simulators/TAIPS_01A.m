% TAIPS_01A - Simulator pentru testarea rutinei 'trend'
clear; clc; close all;

% 1. Incarcarea datelor
n = input('Introduceti numarul seriei de timp (ex: 3 pentru Y3) (max 20): ');
nume_fisier = sprintf('Y%02d', n);
cale = fullfile('TAIPS/TAIPSMAT', [nume_fisier '.mat']);

if ~exist(cale, 'file')
    error('Fisierul %s nu a fost gasit pe calea TAIPS/TAIPSMAT.', cale);
end
load(cale); % Incarcam obiectul IDDATA

% 2. Parametri de simulare
p = input('Introduceti gradul polinomului tendinta (0-10): ');
K = 5; % Orizontul de extrapolare implicit
 
% Luam canalele de iesire la rand
ch = 1;
while ch <= size(Y.y, 2) 
    y = Y.y(:,ch); 
    N = length(y);
    N_masura = N - K; % Orizontul de masura
    
    y_masura = y(1:N_masura);
    y_extrapolare = y(N_masura+1:end);
    
    % Apelarea functiei TREND pe orizontul de masura
    [ysta_m, YT_m, theta] = trend(y_masura, p);
    
    % Extrapolarea tendintei pe ultimele K puncte
    t_extra = (N_masura:N-1)';
    YT_extra = zeros(K, 1);
    for i = 1:p+1
        YT_extra = YT_extra + theta(i) * (t_extra.^(i-1));
    end
    
    % Calcule statistice (Erori si SNR)
    % Eroarea pe masura (v_p)
    v_p = ysta_m; 
    lambda_v = var(v_p, 1);
    snr_m = std(y_masura, 1) / sqrt(lambda_v);
    
    % Eroarea pe extrapolare (epsilon_p)
    epsilon_p = y_extrapolare - YT_extra;
    lambda_e = var(epsilon_p, 1);
    snr_e = std(y_extrapolare, 1) / sqrt(lambda_e);
    
    % Reprezentare Grafica
    if size(Y.y, 2) == 1
        figure('Name', ['Analiza Trend ' nume_fisier], 'NumberTitle', 'off');
        nume_serie = Y.Notes;
    else 
        figure('Name', ['Analiza Trend ' nume_fisier ', Canal ' num2str(ch)], 'NumberTitle', 'off');
        nume_serie = Y.OutputName{ch};
    end
    
    t_axa = 0:N-1;
    axa_x = Y.TimeUnit;
    axa_y = Y.OutputUnit{ch};
    
    % Indici pentru generarea subplot-urilor
    i = 0;
    if p, j = 4; else, j = 3; end

    % Grafic 1: Seria completa si Media
    subplot(j,2,[i+1 i+2]);
    plot(t_axa, y, 'b', 'LineWidth', 1); hold on;
    line([0 N-1], [mean(y) mean(y)], 'Color', 'y', 'LineWidth', 1.5); % Media
    title(['Seria ' nume_fisier ': ' nume_serie ' | Media = ' num2str(mean(y))]);
    xlabel(axa_x);
    ylabel(axa_y);
    legend('Seria de timp', 'Media seriei', 'Location', 'northwest');
    grid on;
    
    if j == 4
        % Grafic 2: Seria completa si tendinta (p nenul)
        subplot(j,2,[i+3 i+4]);
        plot(t_axa(1:N_masura), y_masura, 'b', 'LineWidth', 1); hold on;
        plot(t_axa(N_masura+1:end), y_extrapolare, 'c', 'LineWidth', 1); % Ultimele K date
        plot(t_axa(1:N_masura), YT_m, 'g', 'LineWidth', 2); % Valori simulate
        plot(t_axa(N_masura+1:end), YT_extra, 'r', 'LineWidth', 2); % Valori extrapolate
    
        % Linie verticala de separare
        yl = ylim;
        line([N_masura-0.5 N_masura-0.5], yl, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);
    
        title(['Seria ' nume_fisier ': ' nume_serie ' | Grad p = ' num2str(p)]);
        xlabel(axa_x);
        ylabel(axa_y);
        legend('Masura', 'Real extrapolare', 'Valori simulate', 'Valori extrapolate', 'Separator', 'Location', 'northwest');
        grid on;
        i = 2;
    end

    % Grafic 3 si 4: Eroarea de model (Masura) si Eroarea de extrapolare
    subplot(j,2,i+3);
    plot(0:N_masura-1, v_p, 'm');
    title(['Eroare model: \lambda^2=' num2str(lambda_v, 4) ' SNR=' num2str(snr_m, 3)]);
    xlabel(axa_x); grid on;
    
    subplot(j,2,i+4);
    plot(N_masura:N-1, epsilon_p, 'r');
    title(['Eroare extrapolare: \lambda^2=' num2str(lambda_e, 4) ' SNR=' num2str(snr_e, 3)]);
    xlabel(axa_x); grid on;
    
    % Grafic 5: Extrapolarea (Zoom pe ultimele K puncte)
    subplot(j,2,[i+5 i+6]);
    t_zoom = N_masura:N-1;
    plot(t_zoom, y_extrapolare, 'co-', 'LineWidth', 1); hold on; % Date reale
    plot(t_zoom, YT_extra, 'ro-', 'LineWidth', 2); % Predictia
    title('Extrapolare (Zoom)');
    legend('Real', 'Trend', 'Location', 'best');
    xlabel(axa_x); ylabel(axa_y); grid on;

    % Afisare coeficienti in consola
    fprintf('\n--- Rezultate %s - Canal %d ---\n', nume_fisier, ch);
    fprintf('Coeficientii polinomului tendinta:\n');
    disp(theta');
    ch = ch + 1;
end