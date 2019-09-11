% BATCHPROCESSING Rescales and/or extracts spikes from several files
%
%  Syntax:  batchProcessing(); % Opens a dialog to chose the folder or
%                                file(s) to process. Only does rescale.
%           batchProcessing(handles); % Called from GUI
%           batchProcessing(handles, opts); % Includes options
%
%  Inputs: 
%     handles  -  (Optional) GUI's handles
%     opts     -  (Optional) Struct with the following fields:
%                    + tool - list of tools to run. If not given, the
%                       option choses from a dialog window.
%                    + files - list of files to process. If not given, the
%                       user choses from a dialog window.
%                    + saveFolder - location to save the results. If not
%                       given, it is current active folder.
%                    + saveName - Name of the files to save. Default:
%                       'SET_result_i.mat'
%
% Artemio Soto-Breceda | 5 September 2019
function opts = batchProcessing(varargin)
   % Check inputs
   switch nargin
      case 0
         [opts.files, opts.path] = userChoice();
      case 1
         handles = varargin{1};
         [opts.files, opts.path] = userChoice();
      case 2
         handles = varargin{1};
         opts = varargin{2};
         
         if ~isfield (opts, 'files') || isempty(opts.files)
            [opts.files, opts.path] = userChoice();
         end
      otherwise
         error('Wrong number of inputs. Seek help.');
   end
   
   if ~isfield (opts, 'tool') || isempty(opts.tool)
      opts.tool = selectTools();
   end
   
   if ~isfield (opts, 'saveFolder') || isempty(opts.tool)
      opts.saveFolder = [opts.path filesep 'batchProcessing' filesep];
      opts.saveName = 'result';
   end
end

function [files, path] = userChoice()

   ftypes = {'*.mat;*.smr;*.txt', 'Voltage files (*.mat, *.smr, *.txt)'; ...
        '*.mat', 'M-files (*.mat)'; ...
        '*.smr','Analyze image files (*.img)'; ...
        '*.txt','Text files (*.txt)'; ...
        '*.*',  'All Files (*.*)'};

   files = [];
   [file, path] = uigetfile(ftypes,'Chose a folder to process', 'MultiSelect', 'on');
   if iscell(file) || ischar(file)
      if ischar(file)
         files = {file};
      else
         for i = 1:numel(file)
            files{i} = file{i};
         end
      end
   else
      fprintf('\tThe user didn''t chose a file location. No file was loaded.\n');
      error('No_choice');
   end
end

function tool = selectTools()
   toolList = getBatchToolList();
   idx = listdlg('ListString', toolList);
   tool = toolList(idx);
end