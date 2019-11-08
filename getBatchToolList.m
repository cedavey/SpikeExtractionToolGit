% GETBATCHTOOLLIST Returns the list of tools available for batch processing
%
% Created: 11-September-2019
function [toolList, methodList] = getBatchToolList()
   % TODO: Return the list according to the loaded data type.
   toolList = {'rescale' 'denoise' 'identify ap templates'...
      'extract spikes' 'firing rate' 'export to excel'};
   
   methodList = {'particle filter' 'wavelets' 'threshold'...
      'matched filter' 'moving average' 'spike rate and count'};
end