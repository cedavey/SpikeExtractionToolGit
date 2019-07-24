% SETTOOLTIPS Assigns each item in the application its tooltip (description
% on mouse hover).
%
% Syntax: 
%     h = setTooltips(h); % Sets all the default tooltips to all the
%                               default elements
%     h = setTooltips(h, {'object1' 'object2'}, {'String 1', 'String 2'});
%
% Inputs:
%     h        - GUI handles
%     object   - (Optional) updates object's tooltip with input string
%                 'str'. Can be a single object's name or a vector of names
%     str      - (Optional) Input string to update object's tooltip. It can
%                 be a single string or a cell vector.
%     
% Outputs:
%     h        - GUI handles
%
% Artemio - 11/July/2019
function h = setTooltips(h, varargin)
   if nargin == 1
      object = {'load_voltage', 'add_voltage', 'clear_voltage',...
         'toggleZoomButton', 'time_slider', 'voltage_slider',...
         'new_figure', 'access_voltage', 'save_voltage',...
         'set_tool_params', 'run_tool', 'tool_list', 'method_list',...
         'scroll_axes', 'automatic_params', 'curr_signal'};
      
      str = getTooltips(object);

   elseif nargin ~= 3
      str = sprintf('\tInvalid number of parameters ''getTooltips()''\n');
      printMessage('off', 'Errors', str);
      return
      
   else
      object = varargin{1};
      str = varargin{2};
   end
   
   for i = 1:min(length(object), length(str))
      set(h.(object{i}), 'TooltipString', str{i});
   end
end