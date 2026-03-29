function [y_uniform, N_uniform] = splineB3(t, y)
% SPLINEB3 - Interpoleaza si re-esantioneaza seria de timp
% Date de intrare:
%   t - vectorul momentelor de timp originale
%   y - vectorul valorilor seriei
% Date de iesire:
%   y_uniform - seria re-esantionata uniform
%   N - noua lungime a seriei

% Folosiri: -
% Autor: Andrei BLEJAN
% Creat: Decembrie 23, 2025

    if nargin ~= 2
        error('<splineB3> necesita doua argumente de intrare: [y_uniform, N_uniform] = splineB3(t, y)');
    end

    % Determinarea pasului de esantionare uniform Ts = min(delta_t)
    dt = diff(t);
    Ts = min(dt);

    if Ts <= 0
        Ts = 1; 
    end

    % Construirea uniforma a momentelor de esantionare
    t_start = t(1);
    t_end = t(end);
    t_uniform = t_start : Ts : t_end;
    
    % Interpolarea
    y_uniform = interp1(t, y, t_uniform, 'spline');
    N_uniform = length(y_uniform);

end