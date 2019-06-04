% [unit_scale, unit_label] = getUnitScale(value, label)
%
% Get the scale of an input unit. E.g. if scaling time by dt, then return 
% the timescale in s, ms, us, ns, etc, and then the label 's', 'ms', ...
%
% Inputs:
%  dt    - unit to determine the scale of
%  label - unit type, e.g. 's' for seconds
% Outputs:
%  scale - scale from as small as pica to as large as tera
%  label - unit label including the scale, e.g. KB, ms, ...
function [scale, label] = getUnitScale(value, label)
   if nargin==0
      help getUnitScale; return;
   end
   if nargin<2, label = ''; end
   
   if any(strcmpi(label, {'s','sec','secs','second','seconds'}))
      if value/1e12 > 1
         scale = 1e12;
         label = [' T' label];
      elseif value/(60*60*24) >= 1
         scale = 60*60*24;
         label = 'days';
      elseif value/(60*60) >= 1
         scale = 60*60;
         label = 'hours';
      elseif value/60 >= 1
         scale = 60;
         label = 'minutes';
      elseif value/1 > 1
         scale = 1;
         label = label;
      elseif value/1e-3 >= 1
         scale = 1e-3;
         label = [' m' label];
      elseif value/1e-6 >= 1
         scale = 1e-6;
         label = [' \mu' label];
      elseif value/1e-9 >= 1
         scale = 1e-9;
         label = [' \eta' label];
      else % dt/1e-12 >= 1
         scale = 1e-12;
         label = [' p' label];  
      end
      return;
   end
   
   if value/1e12 > 1
      scale = 1e12;
      label = [' T' label];
   elseif value/1e9 >= 1
      scale = 1e9;
      label = [' G' label];
   elseif value/1e6 >= 1
      scale = 1e6;
      label = [' M' label];
   elseif value/1e3 >= 1
      scale = 1e3;
      label = [' K' label];
   elseif value/1 > 1
      scale = 1;
      label = label;
   elseif value/1e-3 >= 1
      scale = 1e-3;
      label = [' m' label];
   elseif value/1e-6 >= 1
      scale = 1e-6;
      label = [' \mu' label];
   elseif value/1e-9 >= 1
      scale = 1e-9;
      label = [' \eta' label];
   else % dt/1e-12 >= 1
      scale = 1e-12;
      label = [' p' label];  
   end
end