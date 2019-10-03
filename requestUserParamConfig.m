% params = requestUserParamConfig(params, method_name)
% 
% Get user to set configuration parameters for the tool and implementation
% method. 
% Inputs:
%  params - current parameter values, defaults from getDefaultToolParams
%  title  - string with tool/method name for dialogue title
function [params, cancel] = requestUserParamConfig(params, dlg_title)
   % each parameter in the method_params struct has the following fields:
   %  name     - name of parameter, but it can have spaces & other punctuation
   %  value    - current value of parameter
   %  descript - longer sentence explaining the parameter & its units
   %  type     - data type
   %  list     - only when type is list, this gives cell array of list options

   % possible data type options are:
   %  float, positive float, negative float
   %  integer, positive integer, negative integer
   %  boolean
   %  list - where the type is list there's an additional field with list items
   % generate user dialogue to set parameters

   % E.g. for method_params.voltage.rescale.variance - voltage data type,
   % rescale tool, and variance as the method to implement the tool, which
   % has the following parameter (plus others, but only 1 shown here)
   % avg_window.value     = 10; % MA window in seconds
   % avg_window.name      = 'moving average window';
   % avg_window.descript  = 'duration of moving average window (units: seconds)';
   % avg_window.type      = 'positive float';

   names  = fieldnames(params);
   Np     = length(names);
   def    = cell(Np,1);
   prompt = cell(Np,2);

   % different types of input dialogues
   sz = 0;
   pos_float = struct('type','edit','format','float',  'style', 'edit', 'size', sz, 'limits', [   0 Inf], 'items', []);
   neg_float = struct('type','edit','format','float',  'style', 'edit', 'size', sz, 'limits', [-Inf   0], 'items', []);
   norm_float= struct('type','edit','format','float',  'style', 'edit', 'size', sz, 'limits', [   0   1], 'items', []);
   float     = struct('type','edit','format','float',  'style', 'edit', 'size', sz, 'limits', [-Inf Inf], 'items', []);
   percentage= struct('type','edit','format','float',  'style', 'edit', 'size', sz, 'limits', [   0 100], 'items', []);
   pos_int   = struct('type','edit','format','integer','style', 'edit', 'size', sz, 'limits', [   0 Inf], 'items', []);
   neg_int   = struct('type','edit','format','integer','style', 'edit', 'size', sz, 'limits', [-Inf   0], 'items', []);
   int       = struct('type','edit','format','integer','style', 'edit', 'size', sz, 'limits', [-Inf Inf], 'items', []);
   list      = struct('type','list','format','integer','style', 'popupmenu','size',sz);

   for i=1:Np
      defvalue    = params.(names{i}).value;
      def{i,1}    = defvalue;                   % default val is current val
      descript    = params.(names{i}).descript;
      name        = params.(names{i}).name;
      units       = params.(names{i}).units;
      prompt{i,1} = descript;                   % user prompt (use name or descript?!)
      prompt{i,2} = names{i};                   % name of struct field for result
      prompt{i,3} = units;                      % units of parameter
      type        = params.(names{i}).type;

      switch lower(strip(type))
         % float options
         case {'float'}
            formats(i,1) = float;
         case {'positive float'}
            formats(i,1) = pos_float;
         case {'negative float'}
            formats(i,1) = neg_float;
         case {'normalised float'}
            formats(i,1) = norm_float;
         % integer options
         case {'integer'}
            formats(i,1) = int;
         case {'positive integer'}
            formats(i,1) = pos_int;
         case {'negative integer'}
            formats(i,1) = neg_int;
            
         case {'boolean'}
            tmp          = list;
            tmp.limits   = [1 2];
            tmp.items    = {'false','true'};
            formats(i,1) = tmp;
            
            % update default value so instead of being 0/1 it's 1/2
%             defvalue     = find(strcmpi(defvalue, options));
            def{i,1}     = defvalue + 1;

         % list options
         case {'list'}
            options      = params.(names{i}).list;
            tmp          = list;
            tmp.limits   = [1 length(options)];
            tmp.items    = options;
            formats(i,1) = tmp;
            
            % need to update default value for lists from being string to
            % being index into list, else inputsdlg has a hissy
             def{i,1}    = find(strcmpi(defvalue, options));
            
         case {'percentage'}
            formats(i,1) = percentage;
            
         otherwise
            str = sprintf('Error made when setting data type for %s (type: %s)\n', names{i}, lower(strip(type)));
            displayErrorMsg( str );
            cprintf('Errors*', str);
            return;
      end
         
   end

   other.Resize      = 'on';
   other.WindowStyle = 'normal';
   other.Interpreter = 'tex';
   
   % denoising by wavelets - plot mother wavelets
   if strcmpi(names{1}, 'mother_wavelet')
      % if user is selecting mother wavelet, plot options to choose from
      % mother_wavelet.list = {'symlets','daubechies','biorthogonal','coiflets','biorsplines','meyer'};
      wavelets = params.(names{1}).list; % list of possible wavelets
      nr = length(wavelets);
      nc = 8;
      nfilt = cell(nr,1); % number of valid wavelets in family
      [nfilt{:}] = deal(0);
      figure; 
      for j=1:nr
         for i=1:nc
            filter = getMotherWavelet(wavelets{j}, i);
            if ~isempty(filter)
               nfilt{j} = nfilt{j} + 1;
               title_str  = sprintf('%s: %d', wavelets{j}, i);
               subplot(nr,nc,(j-1)*nc+nfilt{j}),
                  plot(filter); 
                  title(title_str,'fontsize',10);
            end
         end
      end
   end
 
   dlg_title = ['Set parameters for ' dlg_title];
   try
      % Fixing bug about wrong parameter format on APs view
      if strcmp('timeseries name',prompt{end}) && ~isinteger(def{end})
         def{end} = [];
      end
      
      [userparams, cancel] = inputsdlg(prompt, dlg_title, formats, def, other);
   catch ME
      str = getCatchMEstring( ME, 'Error setting parameters', true );
      if contains(str, 'Cell array can contain only non-empty character vectors')
         for ii = 1:length(formats(end).items)
            if iscell(formats(end).items{ii})
               formats(end).items{ii} = formats(end).items{ii}{1};
            end
         end
         [userparams, cancel] = inputsdlg(prompt, dlg_title, formats, def, other);
      else
         runtimeErrorHandler(ME,'message',str);
         displayErrorMsg( 'Error setting parameters, reverting to old values' );
         return;
      end
   end
   if cancel
      return;
   end
   
   %% Integrity test user parameter config
   % if using wavelets, check that mother wavelet choice is valid
   str = [];
   if strcmpi( names{1}, 'mother_wavelet' )
      filter = getMotherWavelet(wavelets{userparams.mother_wavelet}, userparams.wavelet_number);
      if isempty(filter)
         str = 'Invalid mother wavelet choice - check wavelet figure for options';
      end
      
   % denoising by filtering
   elseif any( strcmpi( 'filter_order', names ) )
      if userparams.high_cutoff_frequency < userparams.low_cutoff_frequency
         str = 'High cut-off frequency must be greater than low cut-off frequency (duh!), try again';
      end
   
   % rescaling by variance
   elseif any( strcmpi( 'avg_window', names ) ) && any( strcmpi( 'skip_window', names ) ) && userparams.avg_window>0
      if userparams.avg_window < userparams.skip_window
         str = 'You''re skipping ahead by more than the averaging window so you''re ignoring samples, try again';
      end
      
   % avg firing rate
   elseif any( strcmpi( 'moving_avg_window', names ) ) && any( strcmpi( 'skip_window', names ) ) && userparams.moving_avg_window>0
      if userparams.moving_avg_window < userparams.skip_window
         str = 'You''re skipping ahead by more than the averaging window so you''re ignoring samples, try again';
      end
      
   % rescaling (not by variance)
   elseif any( strcmpi( 'glitch_magnitude', names ) )
      if userparams.glitch_magnitude < userparams.voltage_magnitude
         str = 'Glitch threshold is lower than spike threshold so you''ll never identify spike peaks, try again';
      end
   end
   if ~isempty(str)
      displayErrorMsg(str);
      return;
   end
   
   % now update param values with user selections
   for i=1:Np
      type = params.(names{i}).type;

      switch type
         % float options
         case {'percentage', 'float', 'normalised float', 'positive float', 'negative float', 'integer', 'positive integer', 'negative integer'}
            params.(names{i}).value = userparams.(names{i});
            
         case {'boolean'}
            % for the list types extract value from index into list
            params.(names{i}).value = ternaryOp(userparams.(names{i})==1, false, true); % convert from 0/1 to 1/2

         case {'list'}
            % for the list types extract value from index into list
            params.(names{i}).value = formats(i).items{userparams.(names{i})};

      end

   end
end
