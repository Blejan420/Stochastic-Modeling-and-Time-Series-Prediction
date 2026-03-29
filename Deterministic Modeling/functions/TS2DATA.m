function DATA = TS2DATA(nume_fisier)
% TS2DATA - Converteste un fisier .m de tip serie de timp in obiect IDDATA
% Date de intrare:
%   nume_fisier - numele fisierului .m (ex: 'Y16' sau 'Y16.mat')
% Date de iesire:
%   DATA - obiectul IDDATA

% Folosiri: in afara simulatoarelor
% Autor: Andrei BLEJAN
% Creat: Decembrie 19, 2025

    % 1. Executare .M pentru a incarca variabilele in workspace
    try
        run(nume_fisier);
    catch
        error('Eroare la deschiderea fisierului %s.m', nume_fisier);
    end

    % 2. Verificarea datelor obligatorii
    if ~exist('y', 'var')
        error('Fisierul nu contine variabila de date ''y''.');
    end

    % 3. Crearea obiectului IDDATA
    % Verificam daca Ts este valid (uniform) sau lipseste (neuniform)
    if exist('Ts', 'var') && ~isnan(Ts)
        DATA = iddata(y(:), [], Ts);
    else
        DATA = iddata(y(:), []);
    end

    % 4. Includerea informatiilor auxiliare
    DATA.Name = nume_fisier;
    
    if exist('ntime', 'var')
        DATA.SamplingInstants = ntime(:);
    end
    
    if exist('label', 'var')
        DATA.Notes = label;
    end
    
    if exist('yunit', 'var')
        DATA.OutputUnit = {yunit};
    end

    if exist('unit', 'var')
        try
            % Extragem doar cuvantul din paranteze drepte, daca exista
            % Ex: din 'Time [months]' extragem 'months'
            extras = regexp(unit, '\[([^\]]+)\]', 'tokens');
            if ~isempty(extras)
                DATA.TimeUnit = lower(extras{1}{1});
            else
                DATA.TimeUnit = lower(unit);
            end
        catch
            % Daca unitatea nu este standard MATLAB (nanoseconds, seconds, months, years etc.), o stocam in UserData
            DATA.UserData = unit;
        end
    end

    % 5. Completarea informatiilor lipsa de la consola
    if ~exist('label', 'var') || isempty(label)
        DATA.Notes = input(['Introduceti note pentru ' nume_fisier ': '], 's');
    end

    % 6. Salvarea pe disc in format .MAT
    % Se salveaza intr-un subdirector "TAIPSMAT" care se afla in directorul "TAIPS" pentru organizare

    cale_fisier = fullfile('TAIPS/TAIPSMAT/', [nume_fisier '.mat']);
    Y = DATA;
    save(cale_fisier, 'Y');

    
    fprintf('<TS2DATA>: Fisierul %s.m a fost convertit in %s.mat.\n', nume_fisier, nume_fisier);
end