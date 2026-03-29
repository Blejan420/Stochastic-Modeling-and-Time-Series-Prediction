function [P, ys1P, V] = WiRo(ysta)
% WIRO - Algoritmul Whittaker-Robinson pentru estimarea unei perioade
% Date de intrare:
%   ysta - vectorul de date care contine seria de timp stationarizata
% Date de iesire:
%   P    - perioada componentei sezoniere
%   ys1P - componenta sezoniera unimodala detectata
%   V    - vectorul valorilor criteriului de adecvanta

% Folosiri: TAPIS_01B, TAIPS_02C
% Autor: Andrei BLEJAN
% Creat: Decembrie 21, 2025

    if nargin ~= 1
        error('<WiRo> necesita un argument de intrare: [P, ys1P, V] = WiRo(ysta)');
    end

    N = length(ysta);
    mu = 0.4; % Alegem valoarea mediana de 40%
    M = N; % Esantionare uniforma
    V = inf(1, floor(M/2)); % initializam valorile criteriului de adecvanta

    for P = 2:floor(M/2)
        % Determinarea numarului de perioade intregi
        nr_perioade = floor(N / P);
        
        % Construirea matricii de semnal
        % Se iau primele nr_perioade*P puncte si se aseaza matricea
        % de dimensiune P x nr.perioade, dupa care le transpunem, deoarece
        % MATLAB prin reshape pune valorile pe coloane, nu pe linii
        Y_sta = reshape(ysta(1:nr_perioade*P), P, nr_perioade)';
        
        % Estimarea coeficientilor sezonieri
        y_Sp = mean(Y_sta, 1); % media pe coloane
        
        % Construirea componentei sezoniere prin prelungire
        % Se repeta coeficientii pentru a acoperi toata lungimea N
        y_S_P_aux = repmat(y_Sp, 1, nr_perioade + 1);
        y_S_P = y_S_P_aux(1:N)';
        
        % Calculul criteriului de adecvanta V[P]
        err = ysta - y_S_P;
        V(P) = sum(err.^2);
    end
    
    % Determinarea minimului criteriului de adecvanta
    [min_val, P_0] = min(V);
    
    % Calculul energiei totale a seriei stationarizate
    E_ysta = sum(ysta.^2);
    
    % Verificarea pragului energetic pentru selectia perioadei optime
    if min_val <= mu * E_ysta
        P = P_0;
        % Reconstructia componentei sezoniere unimodale pentru optim
        M_opt = floor(N / P);
        Y_sta_opt = reshape(ysta(1:M_opt*P), P, M_opt)';
        ys1P = mean(Y_sta_opt, 1);
    else
        P = 0;
        ys1P = [];
    end

end