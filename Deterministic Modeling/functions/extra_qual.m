function [EQ, snr_n, snr_k] = extra_qual(id_y_masura, id_y_extra, id_ymodel_masura, id_ymodel_extra)
% extra_qual - Calculeaza criteriul de calitate a extrapolarii
% Date de intrare:
%   id_y_masura      - Obiect IDDATA al seriei de timp pe orizontul de masura
%   id_y_extra      - Obiect IDDATA al seriei de timp pe orizontul de extrapolare
%   id_ymodel_masura   - Obiect IDDATA al modelului determinist pe orizontul de masura
%   id_ymodel_extra   - Obiect IDDATA al modelului determinist pe orizontul de extrapolare
% Date de iesire:
%   EQ          - Valoarea criteriului de calitate a extrapolarii [%]
%   snr_n       - Raportul Semnal-Zgomot pe orizontul de masura
%   snr_k       - Raportul Semnal-Zgomot pe orizontul de extrapolare

% Folosiri: TAIPS_02C
% Autor: Andrei BLEJAN
% Creat: Decembrie 27, 2025

    if nargin ~= 4
        error('<extra_qual> necesita patru argumente de intrare: [EQ, snr_n, snr_k] = extra_qual(id_y_masura, id_y_extra, id_yhat_masura, id_yhat_extra');
    end
    
    % Validare argumente intrare de tip IDDATA
    if ~isa(id_y_masura, 'iddata') || ~isa(id_y_extra, 'iddata') || ...
       ~isa(id_ymodel_masura, 'iddata') || ~isa(id_ymodel_extra, 'iddata')
        error('Toate cele 4 argumente de intrare trebuie sa fie obiecte IDDATA.');
    end

    % Constante de ponderare
    w_alpha = 0.75; 
    w_beta  = 0.25; 

    % Extragerea datelor
    y_m = id_y_masura.y;       
    y_e = id_y_extra.y;        
    ymodel_m = id_ymodel_masura.y; 
    ymodel_e = id_ymodel_extra.y;  

    N = length(y_m); % orizont de masura
    K = length(y_e); % orizont de extrapolare

    % --- Orizontul de masura ---
    % Deviatia standard a semnalului
    sigma_y_N = std(y_m, 1); 
    var_y_N   = sigma_y_N^2;

    % Eroarea de masura
    v_p = y_m - ymodel_m;
    lambda_v_sq = var(v_p, 1); % disperia erorii
    lambda_v = std(v_p, 1); % deviatia standard a erorii

    % SNR N
    snr_n = sigma_y_N / lambda_v;


    % --- Orizontul de extrapolare ---
    % Deviatia standard a semnalului
    sigma_y_K = std(y_e, 1);
    var_y_K   = sigma_y_K^2;

    % Eroarea de extrapolare
    epsilon_p = y_e - ymodel_e;
    lambda_eps_sq = var(epsilon_p, 1); 
    lambda_eps    = std(epsilon_p, 1);

    % SNR K
    snr_k = sigma_y_K / lambda_eps;

    % Calculare SNR Global (SNR_NK)
    snr_global = sqrt(N * var_y_N + K * var_y_K / N * lambda_v_sq + K * lambda_eps_sq);

    % Calculare EQ
    
    % Metoda 1 (EQ_alpha) [%]
    
    eq_alpha = 100 / (1 + 1 / (snr_n * snr_k));

    % Metoda 2 (EQ_beta) [%]
    
    eq_beta = 100 / (1 + 1 / snr_global);

    % --- EQ ponderat [%] ---

    EQ = w_alpha * eq_alpha + w_beta * eq_beta;

end