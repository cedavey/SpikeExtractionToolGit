% Created by Artemio - May 2019

function addToolbarExplorationButtons(varargin)
    % This function prevents an error on older MATLAB versions (previous to
    % R2018a)
    % If you are using R2018a or older, you can remove this file.

    v = ver('matlab');
    % Check if version is equal or older than 2018a
    if str2double(v.Version) >= 9.4
        try
           path = fullfile( matlabroot, 'toolbox', 'matlab', 'plottools' );
           file = fullfile( path, 'addToolbarExplorationButtons.p' );
           if exist( path, 'file' ) && exist( file, 'file' )
              run( file );
           end
        catch
           error('The function ''addToolbarExplorationButtons'' has been overwritten for compatibility purposes. If you are using a recent version of MATLAB (2018a or higher), you can remove the file ''addToolbarExplorationButtons.m'' from the SpikeExtractionTool directory.');
        end
    else
        warning('Using an older version of MATLAB. The function ''addToolbarExplorationButtons'' has been hardcoded here.');
    end

    return;
end