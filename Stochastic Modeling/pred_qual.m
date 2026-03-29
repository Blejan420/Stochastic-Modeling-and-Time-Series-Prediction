function [PQ, SNR_N, SNR_K] = pred_qual(id_masura, id_validare, id_simulare, id_predictie, theta)
% PRED_QUAL - Evalueaza criteriul de calitate a predictiei (PQ)

% Date de intrare:
%   id_masura    - Datele reale pe orizontul de masura (N)
%   id_validare  - Datele reale pe orizontul de predictie (K)
%   id_simulare  - Iesirea modelului pe orizontul de masura
%   id_predictie - Iesirea modelului pe orizontul de predictie
%   theta - Coeficientii modelului nedeterminist, necesari pentru calculul 
%           lui alpha in PQ_gamma.
% Date de iesire:
%   PQ           - Calitatea de predictie [%]
%   SNR_N        - SNR pe masura
%   SNR_K        - SNR pe predictie

% Folosiri: TAIPS_02B, TAIPS_02C
% Autori: Andrei BLEJAN
% Creat: 8 Ianuarie 2026

    if nargin ~= 5
        error('<pred_qual> necesita cinci argumente de intrare: [PQ_global, SNR_mas, SNR_pred] = pred_qual(id_masura, id_validare, id_simulare, id_predictie, theta)');
    end
    
    % Validare argumente intrare de tip IDDATA
    if ~isa(id_masura, 'iddata') || ~isa(id_validare, 'iddata') || ...
       ~isa(id_simulare, 'iddata') || ~isa(id_predictie, 'iddata')
        error('Primele 4 argumente de intrare trebuie sa fie obiecte IDDATA.');
    end

    % Definire ponderi
    w_alpha = 0.45; 
    w_beta  = 0.15;
    w_gamma = 0.4; 

    % Extragere date din IDDATA
    yN_real = id_masura.y;
    yN_sim  = id_simulare.y;
    
    yK_real = id_validare.y;
    yK_pred = id_predictie.y;

    % Lungimile orizonturilor
    N = length(yN_real);
    K = length(yK_real);

    % --- Calcul erori si SNR ---
    
    % Orizont de masura
    eN = yN_real - yN_sim;              
    sigma_y_N = std(yN_real, 1);        
    lambda_N = sqrt(sum(eN.^2) / N); 
    SNR_N = sigma_y_N / lambda_N; 

    % Orizont de predictie
    eK = yK_real - yK_pred;             
    sigma_y_K = std(yK_real, 1);        
    lambda_K = sqrt(sum(eK.^2) / K);
    SNR_K = sigma_y_K / lambda_K;
    
    % SNR Global     
    SNR_NK = sqrt((N * (sigma_y_N^2) + K * (sigma_y_K^2)) / (N * (lambda_N^2)  + K * (lambda_K^2)));

    % --- Calcul PQ_alpha ---
    PQ_alpha = 100 / (1 + (1 / (SNR_N * SNR_K)));

    % --- Calcul PQ_beta ---
    PQ_beta = 100 / (1 + (1 / SNR_NK));

    % --- Calcul PQ_gamma ---
    % Calculul recursiv al vectorului sigma
    sigma2 = zeros(K, 1);
    
    % Parametrii de initializare
    lambda_p_na = sqrt(mean(eK.^2));

    % Calculul alpha pe baza coeficientiilor modelului
    if size(theta,2) == 2
        % Model ARMA
        A = theta(:, 1);
        C = theta(:, 2);
    else
        % Model AR
        A = theta;
        C = 1;
    end
    
    % Impartirea polinoamelor este echivalenta cu aplicarea unui filtru C/A
    % pe semnalul Dirac
    delta = zeros(K, 1);
    delta(1) = 1; % impuls Dirac
    alpha = filter(C, A, delta);
    
    % Initializare pas 1
    sigma2(1) = 0;
    
    % Bucla recursiva pentru calculul lui sigma
    for k = 2:K
        sigma2(k) = sigma2(k-1) + (lambda_p_na^2) * (alpha(k-1)^2);
    end

    sigma_k = sqrt(sigma2); 
    abs_eK = abs(eK);
    
    % Calculare PQ_gamma
    
    % Identificare multimi (in interiorul sau exteriorul tubului)
    idx_in  = find(abs_eK <= 3 * sigma_k);
    idx_out = find(abs_eK > 3 * sigma_k);
    
    % Termen 1 al PQ_gamma
    if isempty(idx_in)
        termen_in = 0;
    else                
        termen_in = sqrt((sum((sigma_k(idx_in).^2) .* (eK(idx_in).^2))) ...
                    / (sum(sigma_k(idx_in).^2) * sum(eK(idx_in).^2)));
    end
    
    % Termen 2 al PQ_gamma
    if isempty(idx_out)
        termen_out = 0;
    else
        termen_out = sum(abs_eK(idx_out)./(3*sigma_k(idx_out)));
    end
    
    PQ_gamma = 100 / (1 + termen_in + termen_out);

    % --- Calcul PQ Global ---
    PQ = w_alpha * PQ_alpha + w_beta * PQ_beta + w_gamma * PQ_gamma;

end