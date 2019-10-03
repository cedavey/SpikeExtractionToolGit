% Created by Artemio - May 2019

function addToolbarExplorationButtons(varargin)
    % This function prevents an error on older MATLAB versions (previous to
    % R2018b)
    % If you are using R2018 or older, remove this file.

    v = ver('matlab');
    % Check if version is older than 2018b
    if str2double(v.Version) >= 9.5
        error('Remove the file ''addToolbarExplorationButtons.m from'' your current folder. You are using a recent MATLAB version and don''t need to have that file.');
    else
        warning('Using an older version of MATLAB. The function ''addToolbarExplorationButtons'' has been hardcoded here. Find the file in the current folder.');
    end

    return;
end