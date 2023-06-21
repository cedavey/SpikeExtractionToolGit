% EXPORTTOEXCEL From the spike view, it exports the number of spikes per
% family, and average spike rates.
%
%  Syntax:  exportToExcel(tseries);
%
%  Inputs: 
%     tseries  -  Should be a spike type tseries
%
% GitHub example ()delete this line
%
% Artemio Soto-Breceda | 10 September 2019
function exportToExcel(tseries) 
   Ntemplates = numel(tseries.APstimes);
   Nfamilies = ones(Ntemplates,1);
   % Create a cell variable to save into an Excel file
   Nspikes = cell(1,4);
   % Title row
   Nspikes{1} = 'Template';
   Nspikes{2} = 'Familiy';
   Nspikes{3} = 'Spikes [No.]';
   Nspikes{4} = 'Mean firing rate [Hz]';
   
   famCount = 1; % Total family counter
   for i = 1:Ntemplates
      Nfamilies(i) = numel(tseries.APstimes{i}); % Families in current template
      for ii = 1:Nfamilies(i)
         famCount = famCount + 1; % Starting on row 2, because row 1 is the titles
         Nspikes{famCount,1} = i; % Template
         Nspikes{famCount,2} = ii; % Family
         Nspikes{famCount,3} = numel(tseries.APstimes{i}{ii}); % Spike count
         Nspikes{famCount,4} = Nspikes{famCount,3}/tseries.time(end); % Spike rate
      end
   end
   
   % User choses file location
   [file,path] = uiputfile(['spikes_info' filesep 'spikes_info.xlsx'], 'Save file name');
   if file
      file_name = [path file];
      try
         % Save Excel file
         writecell(Nspikes, file_name);
      catch E
         str = sprintf('\tExcel file couldn''t be saved. Please try again.\n\tError: %s\n', E.message);
         runtimeErrorHandler(E, 'message', str);
      end
      
   else
      fprintf('\tUser didn''t chose a file location. The Excel file wasn''t saved.\n');
   end
   
end