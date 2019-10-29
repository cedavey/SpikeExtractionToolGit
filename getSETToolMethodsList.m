% methods = getSETToolMethodsList(tool, data_type)
% 
% Get list of available methods for given tool and data type. Returns a
% cell array of strings.
function methods = getSETToolMethodsList(tool, data_type)
   switch lower(data_type)
      case 'voltage'
         
         switch lower(tool)
            case 'rescale'
               methods = {'Recursive least squares', 'Particle filter', 'Recursive mean', 'Variance'};

            case 'denoise'
               methods = {'Wavelets', 'Threshold', 'Filter'};
               
            case 'identify ap templates'
               methods = {'Threshold', 'Wavelets', 'K means'};
               
            case 'extract spikes'
               methods = {'Matched Filter', 'K means'};

            case 'utilities'
               methods = {'Downsample', 'Truncate'};

         end
         
      case 'ap'
         switch lower(tool)
            case 'extract spikes'
%                methods = {'Matched filter', 'Wavelets'};
               methods = {'Matched filter'};
               
            case 'merge templates'
               methods = {'User selection'};
               
            case 'delete templates'
               methods = {'User selection'};
         end
         
      case 'spike'
         switch lower(tool)
            case 'firing rate'
               methods = {'Moving average'};
      
            case 'statistics'
               methods = {'Interspike interval', 'Amplitude change', 'PSTH', 'Raster'};
               
            case 'spike operations'
               methods = {'Merge spikes'};
               
            case 'export to excel'
               methods = {'Spike rate and count'};
         end
         
      case 'rate'
         switch lower(tool)
            case 'statistics'
               methods = {'Autocorrelation'};

         end

   end
end
