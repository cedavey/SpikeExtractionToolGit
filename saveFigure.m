% saveFigure(figs,name,fillscreen,varargin)
% 
% Expands figure to take up entire screen, & then saves in multiple
% formats. List desired formats in varargin - one string input per format.
% If no filename provided, use filename in figure properties (gcf.filename)
%
% Inputs:
%  figs     - array of figure handles to save
%  names    - cell array of figure names, unless only 1 figure then its a string
%             Can be empty to use previously saved filename in fig properties
%  fillScreen - true if you want the axes to fill the screen before saving
%  varargin - string for every format you want the figure saved in
function saveFigure(figs,names,fillScreen,varargin)
   if nargin==0
      help saveFigure;
      return;
   end
   % If only 1 figure, names is character array, else in cell array
   if nargin<3, fillScreen = true; end
   if length(figs)==1
      if ~isempty(names) && ischar(names), names = {names}; end
   elseif ~isempty(figs) && length(figs) ~= length(names)
      str = sprintf('saveFigure:inputError: need one name for every figure \n');
      cprintf('*Errors', str);
      return;
   end
   if ~all(isvalid(figs))
      str = sprintf('Deleted figure handle (%d deleted), exiting without saving...\n', sum(~isvalid(figs)));
      cprintf('*Errors', str);
      return;
   end
   nargs = length(varargin);
   if nargs==0
      % if fillscreen is cell instead of logical, user has forgotten to set it
      % & has started passing in the variable arguments from here instead
      if iscell(fillScreen)
         varargin   = fillScreen; 
         fillScreen = true; 
      % if no optional arguments provided & fillScreen is appropriately
      % boolean, set varargin to default figure formats
      else
         varargin{1} = 'pdf'; varargin{2} = 'fig'; 
      end
   elseif nargs==1 && iscell(varargin{1})
      tmp = cell(length(varargin{1}),1);
      [tmp{:}] = varargin{1}{:};
      varargin = tmp;
   end
   units    = get(figs,'units');
   curr_pos = get(figs,'position');
   if fillScreen
      set(figs,'units','normalized','position',[0 0 1 1]);
   end
   for ff=1:length(figs)
      if isempty(names)
         figname = get(figs(ff),'filename');
         % if fig name has no fig type, append .fig
         if ~any(figname=='.')
            figname = [figname '.fig'];
         end
         s = regexp(figname,'[.]','split'); % split in to name & file type
         if length(s)==1 && (ischar(s) && isempty(s) || iscell(s) && isempty(s{:}))
            str = sprintf('\nsaveFigureError: no filename input, and no filename found, exiting...\n');
            cprintf('systemcommands',str);
            return;
         end
         
         [fname,ftype]=s{:};
      end
%       set(figs(ff),'units','normalized','outerposition',[0 0 1 1]);
      for vi=1:length(varargin)
         set(figs(ff), 'units', 'points');
         size = get(figs(ff),'Position');
         size = size(3:4); % the last two elements are width and height of the figure
%          size = size(3:4) + size(1:2); % the last two elements are width and height of the figure
         set(figs(ff),'PaperUnit','points'); % unit for the property PaperSize
         set(figs(ff),'PaperSize',size); 
%          try
         switch varargin{vi}
            case 'eps'
               format  = '-depsc'; % convert eps to colour with 'c' suffix
               fn      = @print;
               options = {};
            case 'fig'
               format  = 'fig';
               fn      = @saveas;
               options = {};
            case {'jpg'}
               format  = '-djpeg';
               fn      = @print;
               options = {};
            case {'png'}
               format  = '-dpng';
               fn      = @print;
               options = {};
            case {'pdf'}
               format  = ['-d' varargin{vi}];
               fn      = @print;
               options = {};
%                options = {'-bestfit'};
            otherwise
               format = varargin{vi};
               fn = @print;
         end
         try
            if ~isempty(names)
               fn(figs(ff), names{ff}, format, options{:});
            else
               fn(figs(ff), fname, format, options{:});
            end

         catch ME
            str = sprintf('Error saving in format %s (%s), ignoring \n', varargin{vi}, ME.message);
            cprintf('keywords', str);
         end
      end
      if length(figs)==1
         set(figs,'units',units,'position',curr_pos);
      else
         arrayfun(@(i) set(figs(i),'units',units{i},'position',curr_pos{i}), 1:length(figs));
      end
   end
end