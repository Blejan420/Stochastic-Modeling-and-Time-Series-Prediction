function [ysta, YT, theta] = trend(y, p)
% trend - Estimeaza tendinta polinomiala de grad p a seriei de timp y
% Date de intrare:
%   y     - vectorul care contine seria de timp
%   p     - ordinul polinomului de aproximare a tendintei
% Date de iesire:
%   ysta  - seria de timp stationarizata (y - YT)
%   YT    - valorile polinomului tendinta la fiecare moment de timp
%   theta - vectorul coeficientilor polinomului tendinta

% Folosiri: TAIPS_01A, TAPIS_01B, TAIPS_02C
% Autor: Andrei BLEJAN
% Creat: Decembrie 20, 2025

    % Prevenirea erorilor
    if (nargin < 2), p = 0; end
    if isempty(p), p = 0; end
    p = abs(round(p));
    if (p > 10)
        p = 10;
        disp('<trend>: Gradul polinomului a fost limitat la 10.');
    end
    if (p > 3)
        disp('<trend>: Atentie! Gradul polinomului depaseste valoarea 3.');
    end

    y = y(:); % Vector coloana
    N = length(y); % Lungimea seriei
    t = (0:N-1)'; % Momente de esantionare uniforme
    
    % Constructia recursiva a matricii RN si vectorului rN
    % RN(i,j) contine media puterilor timpului t^(i+j-2)
    % rN(i) contine media produselor t^(i-1) * y
    RN = zeros(p+1, p+1);
    rN = zeros(p+1, 1);
    
    % Calculam toate puterile necesare ale timpului t (pana la 2*p)
    t_puteri = zeros(2*p + 1, 1);
    for k = 0:2*p
        t_puteri(k+1) = mean(t.^k);
    end
    
    for i = 1:p+1
        for j = 1:p+1
            RN(i,j) = t_puteri(i+j-1);
        end
        rN(i) = mean((t.^(i-1)) .* y);
    end

    % Balansare numerica
    % Folosim o matrice de balansare BM ca in curs
    M = N; % Numarul de balansare pentru serii uniforme (Ts = 1) (Tmax = N)
    BM = diag(M .^ (-0.5:-1:-p-0.5)); % BM = diag(1/M^(1/2), 1/M^(3/2),...)
    
    % Estimarea parametrilor: theta = BM * (BM*RN*BM)^-1 * BM * rN
    theta = BM * ((BM * RN * BM) \ (BM * rN));

    % Evaluarea tendintei pe orizontul de masura
    YT = zeros(N, 1);
    for i = 1:p+1
        YT = YT + theta(i) * (t.^(i-1));
    end

    % Stationarizarea si corectia
    ysta = y - YT;
    theta(1) = theta(1) + mean(ysta);
    ysta = ysta - mean(ysta);

end