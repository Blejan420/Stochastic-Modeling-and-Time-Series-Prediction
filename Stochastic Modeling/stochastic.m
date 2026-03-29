function [yAR, e, lambda2, theta_AR, theta_ARMA] = stochastic(v, na, nc, nalpha)
% STOCHASTIC - Estimeaza componenta nedeterminista (AR sau ARMA)
% Date de intrare:
%   v      - Zgomotul colorat (reziduul determinist)
%   na     - Ordinul partii Auto-Regresive (AR)
%   nc     - (optional) Ordinul partii de medie mobila (MA). 
%            Daca lipseste, se identifica un model AR simplu.
%   nalpha - (optional) Ordinul modelului AR aproximant (folosit doar la ARMA).
%            Daca lipseste, se calculeaza implicit.
% Date de iesire:
%   yAR         - Componenta stocastica reconstruita 
%   e          - Zgomotul alb estimat
%   lambda2    - Dispersia zgomotului alb
%   theta_AR   - Coeficientii modelului AR (principal sau aproximant)
%   theta_ARMA - Coeficientii modelului ARMA (2 coloane: A si C)

% Folosiri: TAIPS_02A, TAIPS_02B, TAIPS_02C
% Autori: Andrei BLEJAN
% Creat: 6 Ianuarie 2026

    % Validare argumente de intrare
    if nargin < 2
        error('<stochastic> necesita cel putin argumentele v si na.');
    end
    
    % Initializare implicita a argumentelor optionale
    if nargin < 3 || isempty(nc)
        nc = 0;
        nalpha = 0;
    end
    
    % Determinare mod de lucru: AR sau ARMA
    arma = (nc > 0);
    
    % --- Identificare AR (Levinson-Durbin) ---

    % Mod AR -> ordinul este na
    % Mod ARMA -> model AR aproximant de ordin nalpha

    if arma
        % Daca nalpha lipseste, se foloseste regula din cursul de TAIPS
        if nargin < 4
            nalpha = 3 * max(na, nc); 
        end
        
        ordin_levinson = nalpha;
    else
        ordin_levinson = na;
    end
    
    % --- Indentificare AR ---

    % xcorr returneaza valori de la -ordin_levinson la +ordin_levinson. 
    % Se alege doar partea pozitiva (inclusiv 0)
    R_aux = xcorr(v, ordin_levinson, 'biased');
    R = R_aux((ordin_levinson+1):end);
    
    % Apelare Levinson pentru coeficientii theta_AR si dispersia lambda2
    [theta_aux, lambda2] = levinson(R, ordin_levinson);
    theta_AR = theta_aux(:); % vector coloana

    % Calculare reziduu e prin filter: e(t) = A(q) * v(t)
    e = filter(theta_AR, 1, v);
    
    % --- Identificare ARMA (daca e cazul) ---

    theta_ARMA = [];
    
    if arma
        data = iddata(v, [], 1);

        % Apelare armax
        model = armax(data, [na nc]);
        
        % Extragere polinoame
        A_polinom = model.A;
        C_polinom = model.C;
        
        % Construire coeficienti theta_ARMA 
        % Vectorii pentru coeficientii A si C, completati cu zero daca lungimile difera
        l = max(length(A_polinom), length(C_polinom));
        theta_ARMA = zeros(l, 2);
        theta_ARMA(1:length(A_polinom), 1) = A_polinom(:);
        theta_ARMA(1:length(C_polinom), 2) = C_polinom(:);
        
        % Calculare si suprascriere reziduu si dispersie zgomot alb rezidual
        e = resid(model, data).y;
        lambda2 = model.NoiseVariance;
    end
    
    % --- Reconstructia componentei stocastice ---

    % Indiferent de model: v(t) = yAR(t) + e(t)
    yAR = v - e;

end