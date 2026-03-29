% TAIPS_01B - Simulator pentru testarea trendului si sezonalitatii
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
    [ysta_m, YT_m, theta_trend] = trend(y_masura, p);
    
    % Extragerea componentei sezoniere
    [ys_m, v_m, theta_sez, P_vec] = seasonal(ysta_m);
    
    % Extrapolare
    t_seasonal_extrapolare = (N_masura+1:N)'; % Indicii conform rutinei seasonal
    t_trend_extrapolare   = (N_masura:N-1)'; % Indicii conform rutinei trend
    
    % Extrapolare trend
    YT_extra = zeros(K, 1);
    for i = 1:p+1
        YT_extra = YT_extra + theta_trend(i) * (t_trend_extrapolare.^(i-1));
    end
    
    % Extrapolare seasonal
    ys_extra = zeros(K, 1);
    if ~isempty(P_vec)
        for k = 1:length(P_vec)
            per = P_vec(k);
            omega = 2 * pi / per;
            ak = theta_sez(2*k-1);
            bk = theta_sez(2*k);

            ys_extra = ys_extra + (ak * sin(omega * t_seasonal_extrapolare) + bk * cos(omega * t_seasonal_extrapolare));
        end
    end
    
    % Calculare model complet
    Y_m_total = YT_m + ys_m; 
    Y_extra_total  = YT_extra + ys_extra;

    % Eroarea de masura
    v_p = y_masura - Y_m_total; 

    lambda_vp = var(v_p, 1);
    snr_vp = std(y_masura, 1) / sqrt(lambda_vp);
    
    % Eroarea de extrapolare
    epsilon_p = y_extrapolare - Y_extra_total;

    lambda_ep = var(epsilon_p, 1);
    snr_ep = std(y_extrapolare, 1) / sqrt(lambda_ep);
    
    % Reprezentare Grafica
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
    if p
        % Grafic 1: Seria completa si Media
        subplot(4,2,[1 2]);
        plot(t_axa, y, 'b', 'LineWidth', 1); hold on;
        line([0 N-1], [mean(y) mean(y)], 'Color', 'y', 'LineWidth', 1.5); % Media
        title(['Seria ' nume_fisier ': ' nume_serie ' | Media = ' num2str(mean(y))]);
        xlabel(axa_x);
        ylabel(axa_y);
        legend('Seria de timp', 'Media seriei', 'Location', 'northwest');
        grid on;
        
        if ~isempty(P_vec)
            % Grafic 2: Model determinist complet
            subplot(4,2,[3 4]);
            % Datele reale (masura + extrapolare)
            plot(t_axa(1:N_masura), y_masura, 'b', 'LineWidth', 1); hold on;
            plot(t_axa(N_masura+1:end), y_extrapolare, 'c', 'LineWidth', 1); 
            
            % Modelul determinist (trend + sezonalitate)
            plot(t_axa(1:N_masura), Y_m_total, 'g', 'LineWidth', 2); 
            plot(t_axa(N_masura+1:end), Y_extra_total, 'r', 'LineWidth', 2); 
        
            % Linie verticala de separare
            yl = ylim;
            line([N_masura-0.5 N_masura-0.5], yl, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);
    
            title(['Model determinist (trend p=' num2str(p) ' + componenta sezoniera). ']);
            xlabel(axa_x);
            ylabel(axa_y);
            legend('Masura Real', 'Extrapolare Real', 'Model Determinist', 'Predictie', 'Separator', 'Location', 'best');
            grid on; axis tight;
        else
            % Cazul fara sezonalitate
            subplot(4,2,[3 4]);
            % Datele reale (masura + extrapolare)
            plot(t_axa(1:N_masura), y_masura, 'b', 'LineWidth', 1); hold on;
            plot(t_axa(N_masura+1:end), y_extrapolare, 'c', 'LineWidth', 1); 
            
            % Grafic 2: Doar Tendinta (YT_m si YT_extra)
            plot(t_axa(1:N_masura), YT_m, 'g', 'LineWidth', 2); 
            plot(t_axa(N_masura+1:end), YT_extra, 'r', 'LineWidth', 2); 
        
            % Linie verticala de separare
            yl = ylim;
            line([N_masura-0.5 N_masura-0.5], yl, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);
    
            title(['Model tendinta (p=' num2str(p) '). Nu exista componenta sezoniera.']);
            xlabel(axa_x);
            ylabel(axa_y);
            legend('Masura Real', 'Extrapolare Real', 'Tendinta Simulata', 'Tendinta Extrapolata', 'Separator', 'Location', 'best');
            grid on; axis tight;
        end
    
        % Grafic 3: Eroarea de model v_p
        subplot(4,2,5);
        plot(0:N_masura-1, v_p, 'm'); 
        title(['Eroare model: \lambda^2=' num2str(lambda_vp, 4) ' SNR=' num2str(snr_vp, 3)]);
        xlabel(axa_x); grid on; axis tight;
        
        % Grafic 4: Eroarea de extrapolare epsilon_p
        subplot(4,2,6);
        plot(N_masura:N-1, epsilon_p, 'r');
        title(['Eroare extrapolare: \lambda^2=' num2str(lambda_ep, 4) ' SNR=' num2str(snr_ep, 3)]);
        xlabel(axa_x); grid on; axis tight;
        
        % Grafic 5: Extrapolarea (Zoom)
        subplot(4,2,[7 8]);
        t_zoom = N_masura:N-1;
        plot(t_zoom, y_extrapolare, 'co-', 'LineWidth', 1); hold on; 
        plot(t_zoom, Y_extra_total, 'ro-', 'LineWidth', 2); 
        title('Extrapolare (Zoom)');
        legend('Real', 'Predictie', 'Location', 'best');
        xlabel(axa_x); ylabel(axa_y); grid on;
    
    else
        i = 0; 
        if ~isempty(P_vec), j = 4; else, j = 3; end
        
        % Grafic 1: Seria completa si Media
        subplot(j,2,[i+1 i+2]);
        plot(t_axa, y, 'b', 'LineWidth', 1); hold on;
        line([0 N-1], [mean(y) mean(y)], 'Color', 'y', 'LineWidth', 1.5); % Media
        title(['Seria ' nume_fisier ': ' nume_serie ' | Media = ' num2str(mean(y))]);
        xlabel(axa_x);
        ylabel(axa_y);
        legend('Seria de timp', 'Media seriei', 'Location', 'northwest');
        grid on;
        
        if j == 4 % daca avem componenta sezoniera, modelul difera de media seriei
            % Grafic 2: Model determinist complet
            subplot(j,2,[i+3 i+4]);
            % Datele reale (masura + extrapolare)
            plot(t_axa(1:N_masura), y_masura, 'b', 'LineWidth', 1); hold on;
            plot(t_axa(N_masura+1:end), y_extrapolare, 'c', 'LineWidth', 1); 
            
            % Modelul determinist (trend + sezonalitate)
            plot(t_axa(1:N_masura), Y_m_total, 'g', 'LineWidth', 2); 
            plot(t_axa(N_masura+1:end), Y_extra_total, 'r', 'LineWidth', 2); 
        
            % Linie verticala de separare
            yl = ylim;
            line([N_masura-0.5 N_masura-0.5], yl, 'Color', 'w', 'LineStyle', '-.', 'LineWidth', 1.5);
    
            title(['Model determinist (trend p=' num2str(p) ' + componenta sezoniera). ']);
            xlabel(axa_x);
            ylabel(axa_y);
            legend('Masura Real', 'Extrapolare Real', 'Model Determinist', 'Predictie', 'Separator', 'Location', 'best');
            grid on; axis tight; 
            i = 2;
        end

        % Grafic 3: Eroarea de model v_p
        subplot(j,2,i+3);
        plot(0:N_masura-1, v_p, 'm'); 
        title(['Eroare model: \lambda^2=' num2str(lambda_vp, 4) ' SNR=' num2str(snr_vp, 3)]);
        xlabel(axa_x); grid on; axis tight;
        
        % Grafic 4: Eroarea de extrapolare epsilon_p
        subplot(j,2,i+4);
        plot(N_masura:N-1, epsilon_p, 'r');
        title(['Eroare extrapolare: \lambda^2=' num2str(lambda_ep, 4) ' SNR=' num2str(snr_ep, 3)]);
        xlabel(axa_x); grid on; axis tight;
        
        % Grafic 5: Extrapolarea (Zoom)
        subplot(j,2,[i+5 i+6]);
        t_zoom = N_masura:N-1;
        plot(t_zoom, y_extrapolare, 'co-', 'LineWidth', 1); hold on; 
        plot(t_zoom, Y_extra_total, 'ro-', 'LineWidth', 2); 
        title('Extrapolare (Zoom)');
        legend('Real', 'Predictie', 'Location', 'best');
        xlabel(axa_x); ylabel(axa_y); grid on;
        
    end

    % Afisare coeficienti in consola
    fprintf('\n--- Rezultate %s - Canal %d ---\n', nume_fisier, ch);
    fprintf('Coeficientii polinomului tendinta:\n');
    disp(theta_trend');

    % Afisare perioade ale modelului si perechile de coef Fourier estimati
    if ~isempty(P_vec)
        fprintf('   Sirul de perioade si coeficientii Fourier estimati:\n');
        fprintf('   ------------------------------------------------------\n');
        fprintf('   | %-10s | %-15s | %-15s |\n', 'Perioada', 'Ak (Sin)', 'Bk (Cos)');
        fprintf('   ------------------------------------------------------\n');
        
        for k = 1:length(P_vec)
            Ak = theta_sez(2*k-1);
            Bk = theta_sez(2*k);

            fprintf('   | %-10d | %-15.4f | %-15.4f |\n', P_vec(k), Ak, Bk);
        end
        fprintf('   ------------------------------------------------------\n');
    else
        fprintf('Nu s-a detectat nicio componenta sezoniera.\n');
    end

    ch = ch + 1;
end