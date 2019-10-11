% Created by Artemio - May 2019

function addToolbarExplorationButtons(varargin)
    % This function prevents an error on older MATLAB versions (previous to
    % R2018b)
    % If you are using R2018 or older, remove this file.

    v = ver('matlab');
    % Check if version is older than 2018a
    if str2double(v.Version) >= 9.4
        run('C:\Program Files\MATLAB\R2019a\toolbox\matlab\plottools\addToolbarExplorationButtons');
    else
        warning('Using an older version of MATLAB. The function ''addToolbarExplorationButtons'' has been hardcoded here.');
    end

    return;
end