% updateGUIParams(handles, tseries)
% Currently a max of 4 parameter boxes - may need to expand on this!
% updateGUIParams(handles, tseries)
% Updates the Generation Parameters panel by putting parameter names and 
% values in the parameter/value display boxes
function updateGUIParams(handles, tseries)
   % To Do: 
   % - start adjusting font size based on parameter/value string lengths
   
   % set parameter & value text box font sizes - allows me to experiment
   ph = findobj(handles.figure1, '-regexp','tag','param[0-9]');
   set(ph, 'fontsize', 12); 
   ph = findobj(handles.figure1, '-regexp','tag','val[0-9]');
   set(ph, 'fontsize', 12); 


   maxParams = 8;
   if isfield( tseries, 'params' ) && ~isempty(tseries.params)
      params = tseries.params;
      % get tool & method name & then remove from params so only
      % param/value pairs remain in the struct
      toolstr   = lower( params.tool );   % tool name
      methodstr = lower( params.method ); % method used to implement tool 
      params    = rmfield( params, {'tool', 'method'} ); 
      genstr    = sprintf('Generation parameters for %s using %s', toolstr, methodstr);
      fontsize  = ternaryOp( length(genstr)>40, 14, 16 );
      set( handles.gen_params, 'String', genstr );
      
      names  = fieldnames(params);
      % don't populate with nested parameters 
      if any( strcmpi(names, 'nested_params') )
         names( find( strcmpi(names, 'nested_params') ) ) = [];
      end
      Np = min( length(names), maxParams ); % currently max of 8 param boxes
      
   else
      Np     = 0;
      i      = 0;
      set( handles.gen_params, 'String', 'Generation parameters' );
   end

   for i=1:Np
      value       = params.(names{i}).value;
      descript    = params.(names{i}).descript;
      name        = params.(names{i}).name;
      units       = params.(names{i}).units;
      type        = params.(names{i}).type;
      
      % don't include unit type when it's something not useful, like
      % units=wavelet, or units=integer
      switch units
         case {'integer', 'wavelet', 'float'}
            units = '';
         otherwise
            % some units don't need to be printed
            if any( strcmpi(units, {'distribution', 'true or false'}) )
               units = '';
            else
               units = sprintf('(%s)', units);
            end
      end
      switch isType(value)
         case 'string'
            value = sprintf('%s %s', value, units);
            % change all occurences of ' to '' for printing
            value = insertAfter(value, '''', '''');
         case 'float'
            value = sprintf('%g %s', value, units);
         case 'boolean'
            value = sprintf('%s %s', ternaryOp(value,'true','false'), units);
         otherwise
            value = sprintf('%g %s', value, units);
      end
      
      eval(['set(handles.param' num2str(i) ', ''String'', ''' name ''');']);
      eval(['set(handles.param' num2str(i) ', ''Visible'',  ''on'');']);
      eval(['set(handles.val'   num2str(i) ', ''String'', ''' value ''');']);
      eval(['set(handles.val'   num2str(i) ', ''Visible'',  ''on'');']);

   end
   % if i is 0 the above for loop makes it empty - gotta fix this shit!
   i = ternaryOp(isempty(i), 0, i);
   for j=(i+1):maxParams
      eval(['set(handles.param' num2str(j) ', ''String'', ''off'');']);
      eval(['set(handles.param' num2str(j) ', ''Visible'', ''off'');']);
      eval(['set(handles.val'   num2str(j) ', ''String'', ''off'');']);
      eval(['set(handles.val'   num2str(j) ', ''Visible'', ''off'');']);
   end

%          voltage.denoise.threshold.minposduration.value      = 5;
%          voltage.denoise.threshold.minposduration.name       = 'min pos duration';
%          voltage.denoise.threshold.minposduration.descript   = 'remove positive voltages that become negative too quickly (units: milliseconds)';
%          voltage.denoise.threshold.minposduration.type       = 'positive float';
%          voltage.denoise.threshold.minposduration.units      = 'ms'; % was \mu s for us coz dialogue has tex interpreter
% 
%          voltage.denoise.threshold.minnegduration.value      = 3;
%          voltage.denoise.threshold.minnegduration.name       = 'min neg duration';
%          voltage.denoise.threshold.minnegduration.descript   = 'remove negative voltages that become positive too quickly (units: milliseconds)';
%          voltage.denoise.threshold.minnegduration.type       = 'positive float';
%          voltage.denoise.threshold.minnegduration.units      = 'ms'; % was \mu s for us coz dialogue has tex interpreter
% 
%          voltage.denoise.threshold.positivethreshold.value   = 5;
%          voltage.denoise.threshold.positivethreshold.name    = 'positive threshold';
%          voltage.denoise.threshold.positivethreshold.descript= 'positive duration must have peak at least this big (units: std dev)';
%          voltage.denoise.threshold.positivethreshold.type    = 'positive float';
%          voltage.denoise.threshold.positivethreshold.units   = 'std dev';
% 
%          voltage.denoise.threshold.negativethreshold.value   = 1;
%          voltage.denoise.threshold.negativethreshold.name    = 'negative threshold';
%          voltage.denoise.threshold.negativethreshold.descript= 'negative duration must have peak magnitude at least this big (units: std dev, 0 to ignore)';
%          voltage.denoise.threshold.negativethreshold.type    = 'positive float';
%          voltage.denoise.threshold.negativethreshold.units   = 'std dev';

end