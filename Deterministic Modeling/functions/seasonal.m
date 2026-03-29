function [ys, v, theta_F, P_vec] = seasonal(ysta)
% SEASONAL - Identifica componenta sezoniera folosind Analiza Fourier
%            si Testul Fisher cu prag de netezire.
% Date de intrare:
%   ysta - Seria de timp stationarizata (fara trend)
% Date de iesire:
%   ys      - Componenta sezoniera reconstruita (suma armonicelor valide)
%   v       - Zgomotul colorat (ysta - ys)
%   theta_F - Vectorul coeficientilor Fourier [A1, B1, A2, B2...] pentru frecventele detectate
%   P_vec   - Vectorul perioadelor detectate corespunzatoare coeficientilor

% Folosiri: TAIPS_01B, TAIPS_02C
% Autor: Andrei BLEJAN
% Creat: Decembrie 22, 2025
% Actualizat: Decembrie 26, 2025

    if nargin ~= 1
        error('<seasonal> necesita un argument de intrare: [ys, v, theta_F, P_vec] = seasonal(ysta)');
    end

    N = length(ysta);
    t = (1:N)'; % axa timpului pentru calculul armonicelor  
    % Pragul energetic se verifica la apelarea rutinei WiRo, deci stabilim
    % doar pragul de netezire
    niu = 12.5;
    
    ysta_wiro = ysta; % Initial lucram cu seria completa
    P_wiro = zeros(1, floor(N/2)); % Lista perioadelor detectate
    i = 0;
    
    % Determinarea setului de perioade
    % Se apeleaza WiRo in bucla pana cand nu se mai detecteaza nicio perioada
    while true
        % Apelam WiRo pana cand perioada componentei sezoniere este 0
        [P_crt, ys1P_crt, ~] = WiRo(ysta_wiro);

        if P_crt == 0
            break;
        end
        
        % Adaugam perioada in lista
        i = i + 1;
        P_wiro(i) = P_crt;
        
        % Reconstruim componenta pe perioada data ca sa o scadem din reziduu
        j = mod(0:N-1, P_crt) + 1;
        ysta_aux = ys1P_crt(j)';
        ysta_wiro = ysta_wiro - ysta_aux;
    end

    % Eliminam zerourile ramase de la initializare
    P_wiro = P_wiro(1:i);
    
    % Procesare si filtrare
    
    if isempty(P_wiro)
        % Cazul in care nu exista sezonalitate
        ys = zeros(N, 1);
        v = ysta;
        theta_F = [];
        P_vec = [];
        disp('<seasonal>: Nu s-au detectat perioade.');
    else
        % Se ordoneaza descrescator setul de perioade
        P_wiro = sort(P_wiro, 'descend');
        
        % Determinarea benzii de trecere a filtrului
        % Filtrul trebuie sa lase sa treaca toate frecventele sezoniere detectate.
        % Cea mai inalta frecventa corespunde celei mai mici perioade detectate.
        P_min = P_wiro(end); 
        
        % Selectia parametrilor conform tabelului
        if P_min <= 3     
            na = 35; fc = 0.5;
        elseif P_min == 4  
            na = 34; fc = 0.234;
        elseif P_min == 5
            na = 32; fc = 0.21;
        elseif P_min == 6
            na = 27; fc = 0.17;
        elseif P_min == 7
            na = 23; fc = 0.148;
        elseif P_min == 8
            na = 21; fc = 0.13;
        elseif P_min == 9
            na = 20; fc = 0.12;
        elseif P_min == 10
            na = 19; fc = 0.108;
        elseif P_min == 11
            na = 17; fc = 0.096;
        elseif P_min == 12
            na = 16; fc = 0.088;
        elseif P_min >= 13 && P_min <= 15
            na = 15; fc = 0.08;
        elseif P_min >= 16 && P_min <= 18
            na = 14; fc = 0.066;
        elseif P_min >= 19 && P_min <= 21
            na = 13; fc = 0.055;
        elseif P_min >= 22 && P_min <= 25
            na = 12; fc = 0.046;
        elseif P_min >= 26 && P_min <= 32
            na = 11; fc = 0.039;
        elseif P_min >= 33 && P_min <= 43
            na = 10; fc = 0.031;
        elseif P_min >= 44 && P_min <= 62
            na = 9; fc = 0.023;
        elseif P_min >= 63 && P_min <= 100
            na = 8; fc = 0.016;
        elseif P_min >= 101 && P_min <= 250
            na = 7; fc = 0.01;
        elseif P_min >= 251 && P_min <= 285
            na = 6; fc = 0.004;
        elseif P_min >= 286 && P_min <= 714
            na = 5; fc = 0.0035;
        elseif P_min >= 715 && P_min <= 2000
            na = 4; fc = 0.0014;
        else % P_min > 2000
            na = 3; fc = 0.0005;
        end
  
        % Apelarea functiei Chebysev si filtrarea seriei

        [b, a] = cheby2(na, 50, fc);
        
        ysta_f = filter(b,a,ysta);
        
        % Construirea modelului Fourier

        L = length(P_wiro);
        % Se construieste matricea M cu 2*L coloane pentru sin si cos
        % corespunzatoare perioadelor din P_wiro
        M = zeros(N, 2*L);
        for k = 1:L
            omega = 2 * pi / P_wiro(k);
            M(:, 2*k-1) = sin(omega * t);
            M(:, 2*k)   = cos(omega * t);
        end

        R_NL = M' * M; 
        r_NL = M' * ysta_f;
        % Coeficientii Fourier conform MCMMP
        theta_F = R_NL \ r_NL;

        % Echilibrare Energetica
        % Calcularea modelului temporar folosind coeficientii actuali
        ys_tilda = M * theta_F;
        
        % Determinarea amplitudinii de echilibrare energetica
        A_tilda = (ys_tilda' * ysta) / (ys_tilda' * ys_tilda);
        
        % Corectarea coeficientilor Fourier
        theta_F = A_tilda * theta_F;
        
        % Periodograma Schuster si Pragul de netezire       
        % Calcularea puterii pentru fiecare perioada si puterea maxima
        Putere = zeros(L, 1);
        for k = 1:L
            Putere(k) = sqrt(theta_F(2*k-1)^2 + theta_F(2*k)^2);
        end
        P_max = max(Putere);
        
        % Selectam indicii care trec de pragul de netezire
        % Conditia: P_i >= (v/100) * P_max
        i_valid = find(Putere >= (niu/100) * P_max);
        
        % Refiltrare si recalculare date
    
        % Verificam daca armonica de perioada minima (ultima din P_wiro) a fost eliminata
        % P_wiro este sortat descrescator, deci P_min este la ultimul index (L)

        if ~isempty(i_valid) && ~ismember(L, i_valid)
            % 1. Determinam noul P_min din perioadele ramase valide
            P_wiro2 = P_wiro(i_valid);
            P_min2 = min(P_wiro2);
            
            % 2. Reproiectam filtrul pentru noul P_min
            if P_min2 <= 3,      na = 35; fc = 0.5;
            elseif P_min2 == 4,  na = 34; fc = 0.234;
            elseif P_min2 == 5,  na = 32; fc = 0.21;
            elseif P_min2 == 6,  na = 29; fc = 0.189;
            elseif P_min2 == 7,  na = 23; fc = 0.148;
            elseif P_min2 == 8,  na = 22; fc = 0.14;
            elseif P_min2 == 9,  na = 20; fc = 0.12;
            elseif P_min2 == 10, na = 19; fc = 0.108;
            elseif P_min2 == 11, na = 17; fc = 0.096;
            elseif P_min2 == 12, na = 16; fc = 0.088;
            elseif P_min2 >= 13 && P_min2 <= 15, na = 15; fc = 0.08;
            elseif P_min2 >= 16 && P_min2 <= 18, na = 14; fc = 0.066;
            elseif P_min2 >= 19 && P_min2 <= 21, na = 13; fc = 0.055;
            elseif P_min2 >= 22 && P_min2 <= 25, na = 12; fc = 0.046;
            elseif P_min2 >= 26 && P_min2 <= 32, na = 11; fc = 0.039;
            elseif P_min2 >= 33 && P_min2 <= 43, na = 10; fc = 0.031;
            elseif P_min2 >= 44 && P_min2 <= 62, na = 9;  fc = 0.023;
            elseif P_min2 >= 63 && P_min2 <= 100, na = 8; fc = 0.016;
            elseif P_min2 >= 101 && P_min2 <= 250, na = 7; fc = 0.01;
            elseif P_min2 >= 251 && P_min2 <= 285, na = 6; fc = 0.004;
            elseif P_min2 >= 286 && P_min2 <= 714, na = 5; fc = 0.0035;
            elseif P_min2 >= 715 && P_min2 <= 2000, na = 4; fc = 0.0014;
            else, na = 3; fc = 0.0005;
            end

            [b, a] = cheby2(na, 50, fc);
        
            ysta_f2 = filter(b,a,ysta);

            L2 = length(P_wiro2);
            M2 = zeros(N, 2*L2);
            for k = 1:L2
                omega = 2 * pi / P_wiro2(k);
                M2(:, 2*k-1) = sin(omega * t);
                M2(:, 2*k)   = cos(omega * t);
            end
            
            % MCMMP
            theta2 = (M2' * M2) \ (M2' * ysta_f2);
            
            % Re-echilibrare Energetica
            ys_tilda2 = M2 * theta2;
            A_tilda = (ys_tilda2' * ysta) / (ys_tilda2' * ys_tilda2);
            theta2 = A_tilda * theta2;
            
            P_final = P_wiro2;     % lista perioadelor finale
            theta_final = theta2;  % lista coeficientilor finali
        else
            
            P_final = P_wiro(i_valid);
            theta_final = zeros(2 * length(i_valid), 1);
            for k = 1:length(i_valid)
                idx = i_valid(k);
                theta_final(2*k-1) = theta_F(2*idx-1);
                theta_final(2*k)   = theta_F(2*idx);
            end
        end

        % Reconstructia finala a componentei sezoniere
        ys = zeros(N, 1);
        theta_F = zeros(2 * length(P_final), 1); % cate 2 coeficienti per perioada
        P_vec = zeros(1, length(P_final));

        for k = 1:length(P_final)

            P_val = P_final(k);
            
            % Reconstruim componenta sezoniera
            ak = theta_final(2*k-1);
            bk = theta_final(2*k);
            omega = 2 * pi / P_val;

            ys = ys + (ak * sin(omega * t) + bk * cos(omega * t));
            
            % Salvam in vectori doar componentele pastrate
            theta_F(2*k-1) = ak;
            theta_F(2*k)   = bk;
            P_vec(k) = P_val;
        end
        
        % Calculam zgomotul colorat
        v = ysta - ys;

    end
end
    