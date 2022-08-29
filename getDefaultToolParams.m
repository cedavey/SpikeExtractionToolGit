% params = getDefaultToolParams
%
% Get all default parameters for all tools and possible implementations
function params = getDefaultToolParams
   % Possible parameters that are not included at the moment
   % voltage.denoise.wavelets - method (currently set to 'modwtsqtwolog')
   % format of param struct is: 
      % data_type.tool.method.param.value;    % default value of param
      % data_type.tool.method.param.name;     % name of param for input dlg
      % data_type.tool.method.param.descript; % help on param for user
      % data_type.tool.method.param.type;     % data type for when we request user values
      % possible types: 'positive float', 'boolean', 'percentage',
      %                 'normalised float', 'positive integer'
      
   %% Voltage data 
   
      % Rescale - based on variance (implicit Gaussian assumption)
         % In the interests of creating as few configurable parameters together I
         % think only allow a MA estimation of sigma, as opposed to allowing
         % different filters, such as lowpass, butterworth etc
         voltage.preprocess.downsample.keep_every_n.value                       = 5; % MA window in seconds
         voltage.preprocess.downsample.keep_every_n.name                        = 'keep every N';
         voltage.preprocess.downsample.keep_every_n.descript                    = 'keep every Nth sample to reduce voltage size';
         voltage.preprocess.downsample.keep_every_n.type                        = 'positive integer';
         voltage.preprocess.downsample.keep_every_n.units                       = 'sample';
         
      % Rescale - using variance
         voltage.rescale.variance.moving_avg_window.value                       = 10; % MA window in seconds
         voltage.rescale.variance.moving_avg_window.name                        = 'moving avg window';
         voltage.rescale.variance.moving_avg_window.descript                    = 'duration of moving average window';
         voltage.rescale.variance.moving_avg_window.type                        = 'positive float';
         voltage.rescale.variance.moving_avg_window.units                       = 'seconds';
         
         voltage.rescale.variance.skip_window.value                             = 5; % skip ahead window in seconds
         voltage.rescale.variance.skip_window.name                              = 'skip window';
         voltage.rescale.variance.skip_window.descript                          = 'skip ahead time to next averaging window';
         voltage.rescale.variance.skip_window.type                              = 'positive float';
         voltage.rescale.variance.skip_window.units                             = 'seconds';
                    
         voltage.rescale.variance.separate_pos_neg.value                        = false; % consider +/- voltages separately
         voltage.rescale.variance.separate_pos_neg.name                         = 'separate pos neg';
         voltage.rescale.variance.separate_pos_neg.descript                     = 'consider positive and negative voltages separately';
         voltage.rescale.variance.separate_pos_neg.type                         = 'boolean';
         voltage.rescale.variance.separate_pos_neg.units                        = 'true or false';
        
         voltage.rescale.variance.outlier.value                                 = 10; % remove large voltages
         voltage.rescale.variance.outlier.name                                  = 'outlier';
         voltage.rescale.variance.outlier.descript                              = 'remove voltages larger than this many std devs from 0';
         voltage.rescale.variance.outlier.type                                  = 'positive float';
         voltage.rescale.variance.outlier.units                                 = 'std dev';
         
         voltage.rescale.variance.distribution.value                            = 'gaussian'; % or student's-t
         voltage.rescale.variance.distribution.name                             = 'data distribution';
         voltage.rescale.variance.distribution.descript                         = 'assumed distribution of data';
         voltage.rescale.variance.distribution.type                             = 'list';
         voltage.rescale.variance.distribution.list                             = {'gaussian','student''s-t'};
         voltage.rescale.variance.distribution.units                            = 'distribution';
         
         voltage.rescale.variance.initial.value                                 = 0;
         voltage.rescale.variance.initial.name                                  = 'Start rescaling at';
         voltage.rescale.variance.initial.descript                              = 'Starting rescale point';
         voltage.rescale.variance.initial.type                                  = 'positive float';
         voltage.rescale.variance.initial.units                                 = 'seconds or samples';
         
         voltage.rescale.variance.initial_type.value                            = 'seconds';
         voltage.rescale.variance.initial_type.name                             = 'type';
         voltage.rescale.variance.initial_type.descript                         = 'Starting rescale point in seconds or samples';
         voltage.rescale.variance.initial_type.type                             = 'list';
         voltage.rescale.variance.initial_type.list                             = {'seconds', 'samples'};
         voltage.rescale.variance.initial_type.units                            = '';
         
         voltage.rescale.variance.final.value                                   = 0;
         voltage.rescale.variance.final.name                                    = 'Stop rescaling at';
         voltage.rescale.variance.final.descript                                = 'Final rescale point';
         voltage.rescale.variance.final.type                                    = 'positive float';
         voltage.rescale.variance.final.units                                   = 'seconds or samples';
         
         voltage.rescale.variance.final_type.value                              = 'seconds';
         voltage.rescale.variance.final_type.name                               = 'type';
         voltage.rescale.variance.final_type.descript                           = 'Final rescale point in seconds or samples';
         voltage.rescale.variance.final_type.type                               = 'list';
         voltage.rescale.variance.final_type.list                               = {'seconds', 'samples'};
         voltage.rescale.variance.final_type.units                              = '';

      % Rescale - based on a particle filter 
         voltage.rescale.particle_filter.voltage_magnitude.value                = 5; 
         voltage.rescale.particle_filter.voltage_magnitude.name                 = 'voltage magnitude';
         voltage.rescale.particle_filter.voltage_magnitude.descript             = 'Threshold at which voltage is considered a spike rather than noise';
         voltage.rescale.particle_filter.voltage_magnitude.type                 = 'positive float';
         voltage.rescale.particle_filter.voltage_magnitude.units                = 'std dev';
         
         voltage.rescale.particle_filter.glitch_magnitude.value                 = 15; 
         voltage.rescale.particle_filter.glitch_magnitude.name                  = 'glitch magnitude';
         voltage.rescale.particle_filter.glitch_magnitude.descript              = 'Threshold at which voltage is considered a glitch of some sort';
         voltage.rescale.particle_filter.glitch_magnitude.type                  = 'positive float';
         voltage.rescale.particle_filter.glitch_magnitude.units                 = 'std dev';
         
         voltage.rescale.particle_filter.forgetting_factor.value                = 0.99;
         voltage.rescale.particle_filter.forgetting_factor.name                 = 'forgetting factor';
         voltage.rescale.particle_filter.forgetting_factor.descript             = 'Controls how quickly we forget past samples (0 - forget fastest, 1 - forget slowest)';
         voltage.rescale.particle_filter.forgetting_factor.type                 = 'normalised float';
         voltage.rescale.particle_filter.forgetting_factor.units                = '\in [0, 1]';
         
         voltage.rescale.particle_filter.select_peaks.value                     = 'separate';
         voltage.rescale.particle_filter.select_peaks.name                      = 'select peaks';
         voltage.rescale.particle_filter.select_peaks.descript                  = 'Select positive, negative, or both positive and negative peaks ';
         voltage.rescale.particle_filter.select_peaks.type                      = 'list';
         voltage.rescale.particle_filter.select_peaks.list                      = {'all', 'positive', 'negative','separate'};
         voltage.rescale.particle_filter.select_peaks.units                     = 'voltage';
         
         voltage.rescale.particle_filter.jump_ahead.value                       = 10;
         voltage.rescale.particle_filter.jump_ahead.name                        = 'jump ahead';
         voltage.rescale.particle_filter.jump_ahead.descript                    = 'Get largest spike in this amount of time';
         voltage.rescale.particle_filter.jump_ahead.type                        = 'positive float';
         voltage.rescale.particle_filter.jump_ahead.units                       = 'seconds';
         
         voltage.rescale.particle_filter.initial.value                          = 0;
         voltage.rescale.particle_filter.initial.name                           = 'Start rescaling at';
         voltage.rescale.particle_filter.initial.descript                       = 'Starting rescale point';
         voltage.rescale.particle_filter.initial.type                           = 'positive float';
         voltage.rescale.particle_filter.initial.units                          = 'seconds or samples';
         
         voltage.rescale.particle_filter.initial_type.value                     = 'seconds';
         voltage.rescale.particle_filter.initial_type.name                      = 'type';
         voltage.rescale.particle_filter.initial_type.descript                  = 'Starting rescale point in seconds or samples';
         voltage.rescale.particle_filter.initial_type.type                      = 'list';
         voltage.rescale.particle_filter.initial_type.list                      = {'seconds', 'samples'};
         voltage.rescale.particle_filter.initial_type.units                     = '';
         
         voltage.rescale.particle_filter.final.value                            = 0;
         voltage.rescale.particle_filter.final.name                             = 'Stop rescaling at';
         voltage.rescale.particle_filter.final.descript                         = 'Final rescale point';
         voltage.rescale.particle_filter.final.type                             = 'positive float';
         voltage.rescale.particle_filter.final.units                            = 'seconds or samples';
         
         voltage.rescale.particle_filter.final_type.value                       = 'seconds';
         voltage.rescale.particle_filter.final_type.name                        = 'type';
         voltage.rescale.particle_filter.final_type.descript                    = 'Final rescale point in seconds or samples';
         voltage.rescale.particle_filter.final_type.type                        = 'list';
         voltage.rescale.particle_filter.final_type.list                        = {'seconds', 'samples'};
         voltage.rescale.particle_filter.final_type.units                       = '';
         
      % Rescale voltage - least squares
         voltage.rescale.recursive_least_squares.voltage_magnitude.value        = 5; 
         voltage.rescale.recursive_least_squares.voltage_magnitude.name         = 'voltage magnitude';
         voltage.rescale.recursive_least_squares.voltage_magnitude.descript     = 'Size of voltage when it is considered a spike rather than noise';
         voltage.rescale.recursive_least_squares.voltage_magnitude.type         = 'positive float';
         voltage.rescale.recursive_least_squares.voltage_magnitude.units        = 'std dev';

         voltage.rescale.recursive_least_squares.glitch_magnitude.value         = 15; 
         voltage.rescale.recursive_least_squares.glitch_magnitude.name          = 'glitch magnitude';
         voltage.rescale.recursive_least_squares.glitch_magnitude.descript      = 'Threshold at which voltage is considered a glitch of some sort';
         voltage.rescale.recursive_least_squares.glitch_magnitude.type          = 'positive float';
         voltage.rescale.recursive_least_squares.glitch_magnitude.units         = 'std dev';
         
         voltage.rescale.recursive_least_squares.forgetting_factor.value        = 0.99;
         voltage.rescale.recursive_least_squares.forgetting_factor.name         = 'forgetting factor';
         voltage.rescale.recursive_least_squares.forgetting_factor.descript     = 'Controls how quickly we forget past samples (0 - forget fastest, 1 - forget slowest)';
         voltage.rescale.recursive_least_squares.forgetting_factor.type         = 'normalised float';
         voltage.rescale.recursive_least_squares.forgetting_factor.units        = '\in [0, 1]';
         
         voltage.rescale.recursive_least_squares.select_peaks.value             = 'separate';
         voltage.rescale.recursive_least_squares.select_peaks.name              = 'select peaks';
         voltage.rescale.recursive_least_squares.select_peaks.descript          = 'Select positive, negative, or both positive and negative peaks ';
         voltage.rescale.recursive_least_squares.select_peaks.type              = 'list';
         voltage.rescale.recursive_least_squares.select_peaks.list              = {'positive', 'negative', 'all','separate'};
         voltage.rescale.recursive_least_squares.select_peaks.units             = 'voltage';
         
         voltage.rescale.recursive_least_squares.jump_ahead.value               = 10;
         voltage.rescale.recursive_least_squares.jump_ahead.name                = 'jump ahead';
         voltage.rescale.recursive_least_squares.jump_ahead.descript            = 'Jump ahead to get largest spike in this amount of time';
         voltage.rescale.recursive_least_squares.jump_ahead.type                = 'positive float';
         voltage.rescale.recursive_least_squares.jump_ahead.units               = 'seconds';
         
         voltage.rescale.recursive_least_squares.initial.value                  = 0;
         voltage.rescale.recursive_least_squares.initial.name                   = 'Start rescaling at';
         voltage.rescale.recursive_least_squares.initial.descript               = 'Starting rescale point';
         voltage.rescale.recursive_least_squares.initial.type                   = 'positive float';
         voltage.rescale.recursive_least_squares.initial.units                  = 'seconds or samples';
         
         voltage.rescale.recursive_least_squares.initial_type.value             = 'seconds';
         voltage.rescale.recursive_least_squares.initial_type.name              = 'type';
         voltage.rescale.recursive_least_squares.initial_type.descript          = 'Starting rescale point in seconds or samples';
         voltage.rescale.recursive_least_squares.initial_type.type              = 'list';
         voltage.rescale.recursive_least_squares.initial_type.list              = {'seconds', 'samples'};
         voltage.rescale.recursive_least_squares.initial_type.units             = '';
         
         voltage.rescale.recursive_least_squares.final.value                    = 0;
         voltage.rescale.recursive_least_squares.final.name                     = 'Stop rescaling at';
         voltage.rescale.recursive_least_squares.final.descript                 = 'Final rescale point';
         voltage.rescale.recursive_least_squares.final.type                     = 'positive float';
         voltage.rescale.recursive_least_squares.final.units                    = 'seconds or samples';
         
         voltage.rescale.recursive_least_squares.final_type.value               = 'seconds';
         voltage.rescale.recursive_least_squares.final_type.name                = 'type';
         voltage.rescale.recursive_least_squares.final_type.descript            = 'Final rescale point in seconds or samples';
         voltage.rescale.recursive_least_squares.final_type.type                = 'list';
         voltage.rescale.recursive_least_squares.final_type.list                = {'seconds', 'samples'};
         voltage.rescale.recursive_least_squares.final_type.units               = '';

      % Rescale voltage - recursive mean
         voltage.rescale.recursive_mean.voltage_magnitude.value                 = 5; 
         voltage.rescale.recursive_mean.voltage_magnitude.name                  = 'voltage magnitude';
         voltage.rescale.recursive_mean.voltage_magnitude.descript              = 'Size of voltage when it is considered a spike rather than noise';
         voltage.rescale.recursive_mean.voltage_magnitude.type                  = 'positive float';
         voltage.rescale.recursive_mean.voltage_magnitude.units                 = 'std dev';

         voltage.rescale.recursive_mean.glitch_magnitude.value                  = 15; 
         voltage.rescale.recursive_mean.glitch_magnitude.name                   = 'glitch magnitude';
         voltage.rescale.recursive_mean.glitch_magnitude.descript               = 'Threshold at which voltage is considered a glitch of some sort';
         voltage.rescale.recursive_mean.glitch_magnitude.type                   = 'positive float';
         voltage.rescale.recursive_mean.glitch_magnitude.units                  = 'std dev';
         
         voltage.rescale.recursive_mean.forgetting_factor.value                 = 0.99;
         voltage.rescale.recursive_mean.forgetting_factor.name                  = 'forgetting factor';
         voltage.rescale.recursive_mean.forgetting_factor.descript              = 'Controls how quickly we forget past samples (0 - forget fastest, 1 - forget slowest)';
         voltage.rescale.recursive_mean.forgetting_factor.type                  = 'normalised float';
         voltage.rescale.recursive_mean.forgetting_factor.units                 = '\in [0, 1]';
         
         voltage.rescale.recursive_mean.select_peaks.value                      = 'positive';
         voltage.rescale.recursive_mean.select_peaks.name                       = 'select peaks';
         voltage.rescale.recursive_mean.select_peaks.descript                   = 'Select positive, negative, or both positive and negative peaks ';
         voltage.rescale.recursive_mean.select_peaks.type                       = 'list';
         voltage.rescale.recursive_mean.select_peaks.list                       = {'positive', 'negative', 'all','separate'};
         voltage.rescale.recursive_mean.select_peaks.units                      = 'voltage';
         
         voltage.rescale.recursive_mean.initial.value                           = 0;
         voltage.rescale.recursive_mean.initial.name                            = 'Start rescaling at';
         voltage.rescale.recursive_mean.initial.descript                        = 'Starting rescale point';
         voltage.rescale.recursive_mean.initial.type                            = 'positive float';
         voltage.rescale.recursive_mean.initial.units                           = 'seconds or samples';
         
         voltage.rescale.recursive_mean.initial_type.value                      = 'seconds';
         voltage.rescale.recursive_mean.initial_type.name                       = 'type';
         voltage.rescale.recursive_mean.initial_type.descript                   = 'Starting rescale point in seconds or samples';
         voltage.rescale.recursive_mean.initial_type.type                       = 'list';
         voltage.rescale.recursive_mean.initial_type.list                       = {'seconds', 'samples'};
         voltage.rescale.recursive_mean.initial_type.units                      = '';
         
         voltage.rescale.recursive_mean.final.value                             = 0;
         voltage.rescale.recursive_mean.final.name                              = 'Stop rescaling at';
         voltage.rescale.recursive_mean.final.descript                          = 'Final rescale point';
         voltage.rescale.recursive_mean.final.type                              = 'positive float';
         voltage.rescale.recursive_mean.final.units                             = 'seconds or samples';
         
         voltage.rescale.recursive_mean.final_type.value                        = 'seconds';
         voltage.rescale.recursive_mean.final_type.name                         = 'type';
         voltage.rescale.recursive_mean.final_type.descript                     = 'Final rescale point in seconds or samples';
         voltage.rescale.recursive_mean.final_type.type                         = 'list';
         voltage.rescale.recursive_mean.final_type.list                         = {'seconds', 'samples'};
         voltage.rescale.recursive_mean.final_type.units                        = '';

         voltage.rescale.recursive_mean.jump_ahead.value                        = 10;
         voltage.rescale.recursive_mean.jump_ahead.name                         = 'jump ahead';
         voltage.rescale.recursive_mean.jump_ahead.descript                     = 'Get largest spike in this amount of time';
         voltage.rescale.recursive_mean.jump_ahead.type                         = 'positive float';
         voltage.rescale.recursive_mean.jump_ahead.units                        = 'seconds';
         
      % Voltage utilities - downsample
         voltage.utilities.downsample.method.value                              = 'ditch';
         voltage.utilities.downsample.method.name                               = 'method';
         voltage.utilities.downsample.method.descript                           = 'downsample by ditching samples or by averaging across samples';
         voltage.utilities.downsample.method.type                               = 'list';
         voltage.utilities.downsample.method.list                               = {'ditch', 'average'};
         voltage.utilities.downsample.method.units                              = '';

         voltage.utilities.downsample.factor.value                              = 3;
         voltage.utilities.downsample.factor.name                               = 'factor';
         voltage.utilities.downsample.factor.descript                           = 'downsample by a factor, e.g. for 3 we only keep a third of samples';
         voltage.utilities.downsample.factor.type                               = 'positive integer';
         voltage.utilities.downsample.factor.units                              = 'integer';

      % Voltage utilities - truncate
         voltage.utilities.truncate.method.value                                = 'start';
         voltage.utilities.truncate.method.name                                 = 'method';
         voltage.utilities.truncate.method.descript                             = 'truncate from the start or end';
         voltage.utilities.truncate.method.type                                 = 'list';
         voltage.utilities.truncate.method.list                                 = {'start', 'end'};
         voltage.utilities.truncate.method.units                                = '';
         
         voltage.utilities.truncate.type.value                                  = 'seconds';
         voltage.utilities.truncate.type.name                                   = 'type';
         voltage.utilities.truncate.type.descript                               = 'truncate in seconds or samples';
         voltage.utilities.truncate.type.type                                   = 'list';
         voltage.utilities.truncate.type.list                                   = {'seconds', 'samples'};
         voltage.utilities.truncate.type.units                                  = '';
         
         voltage.utilities.truncate.amount.value                                = 1000;
         voltage.utilities.truncate.amount.name                                 = 'amount';
         voltage.utilities.truncate.amount.descript                             = 'amount to truncate by';
         voltage.utilities.truncate.amount.type                                 = 'positive integer';
         voltage.utilities.truncate.amount.units                                = '';
         
         % Denoise - remove values under a certain size        
         voltage.denoise.threshold.min_positive_duration.value                  = 5;
         voltage.denoise.threshold.min_positive_duration.name                   = 'min positive duration';
         voltage.denoise.threshold.min_positive_duration.descript               = 'min consecutive duration of positive voltages to be retained';
         voltage.denoise.threshold.min_positive_duration.type                   = 'positive float';
         voltage.denoise.threshold.min_positive_duration.units                  = 'ms'; % was \mu s for us coz dialogue has tex interpreter

         voltage.denoise.threshold.min_negative_duration.value                  = 5;
         voltage.denoise.threshold.min_negative_duration.name                   = 'min negative duration';
         voltage.denoise.threshold.min_negative_duration.descript               = 'min consecutive duration of negative voltages to be retained';
         voltage.denoise.threshold.min_negative_duration.type                   = 'positive float';
         voltage.denoise.threshold.min_negative_duration.units                  = 'ms'; % was \mu s for us coz dialogue has tex interpreter

         voltage.denoise.threshold.positive_threshold.value                     = 5;
         voltage.denoise.threshold.positive_threshold.name                      = 'positive threshold';
         voltage.denoise.threshold.positive_threshold.descript                  = 'positive duration must have peak at least this big';
         voltage.denoise.threshold.positive_threshold.type                      = 'positive float';
         voltage.denoise.threshold.positive_threshold.units                     = 'std dev';

         voltage.denoise.threshold.negative_threshold.value                     = 1;
         voltage.denoise.threshold.negative_threshold.name                      = 'negative threshold';
         voltage.denoise.threshold.negative_threshold.descript                  = 'negative duration must have peak magnitude at least this big';
         voltage.denoise.threshold.negative_threshold.type                      = 'positive float';
         voltage.denoise.threshold.negative_threshold.units                     = 'std dev';

      % Denoise - bandpass filter
         voltage.denoise.filter.filter_order.value                              = 10;
         voltage.denoise.filter.filter_order.name                               = 'filter_order';
         voltage.denoise.filter.filter_order.descript                           = 'number of filter coefficients to use';
         voltage.denoise.filter.filter_order.type                               = 'positive integer';
         voltage.denoise.filter.filter_order.units                              = 'integer';

         voltage.denoise.filter.high_cutoff_frequency.value                     = 100;
         voltage.denoise.filter.high_cutoff_frequency.name                      = 'high cutoff frequency';
         voltage.denoise.filter.high_cutoff_frequency.descript                  = 'remove fast changing signal by filtering out frequencies above this';
         voltage.denoise.filter.high_cutoff_frequency.type                      = 'positive float';
         voltage.denoise.filter.high_cutoff_frequency.units                     = 'Hz';
           
         voltage.denoise.filter.low_cutoff_frequency.value                      = 1;
         voltage.denoise.filter.low_cutoff_frequency.name                       = 'low cutoff frequency';
         voltage.denoise.filter.low_cutoff_frequency.descript                   = 'remove drift and slow changing signal by filtering out frequencies below this';
         voltage.denoise.filter.low_cutoff_frequency.type                       = 'positive float';
         voltage.denoise.filter.low_cutoff_frequency.units                      = 'Hz';

      % Denoise - use a wavelet denoising technique
         voltage.denoise.wavelets.mother_wavelet.value                          = 'symlets'; % name of mother wavelet to use
         voltage.denoise.wavelets.mother_wavelet.name                           = 'mother wavelet';
         voltage.denoise.wavelets.mother_wavelet.descript                       = 'name of mother wavelet to use';
         voltage.denoise.wavelets.mother_wavelet.type                           = 'list';
         voltage.denoise.wavelets.mother_wavelet.list                           = {'symlets','daubechies','biorthogonal','coiflets','fejer-korovkin','discrete meyer','reverse biorthogonal'};
         voltage.denoise.wavelets.mother_wavelet.units                          = 'wavelet';
         
         voltage.denoise.wavelets.wavelet_number.value                          = 7; % number wavelet to use in family
         voltage.denoise.wavelets.wavelet_number.name                           = 'wavelet number';
         voltage.denoise.wavelets.wavelet_number.descript                       = 'index of wavelet to use in wavelet family';
         voltage.denoise.wavelets.wavelet_number.type                           = 'positive integer';
         voltage.denoise.wavelets.wavelet_number.units                          = 'integer';
                    
         voltage.denoise.wavelets.wavelet_level.value                           = 5; % number of wavelets scales/level to use 
         voltage.denoise.wavelets.wavelet_level.name                            = 'wavelet_level';
         voltage.denoise.wavelets.wavelet_level.descript                        = 'number of wavelet levels to use';
         voltage.denoise.wavelets.wavelet_level.type                            = 'positive integer';
         voltage.denoise.wavelets.wavelet_level.units                           = 'integer';
   
      % AP template - by thresholding
         voltage.identify_ap_templates.threshold.positive_threshold.value       = 5;
         voltage.identify_ap_templates.threshold.positive_threshold.name        = 'positive threshold';
         voltage.identify_ap_templates.threshold.positive_threshold.descript    = 'only spikes with peak voltage larger than this many std devs considered';
         voltage.identify_ap_templates.threshold.positive_threshold.type        = 'positive float';
         voltage.identify_ap_templates.threshold.positive_threshold.units       = 'std dev';
         
         voltage.identify_ap_templates.threshold.negative_threshold.value       = 3;
         voltage.identify_ap_templates.threshold.negative_threshold.name        = 'negative threshold';
         voltage.identify_ap_templates.threshold.negative_threshold.descript    = 'only spikes with a negative voltage of magnitude larger than this considered';
         voltage.identify_ap_templates.threshold.negative_threshold.type        = 'positive float';
         voltage.identify_ap_templates.threshold.negative_threshold.units       = 'std dev';
         
         voltage.identify_ap_templates.threshold.glitch_threshold.value       = 10;
         voltage.identify_ap_templates.threshold.glitch_threshold.name        = 'glitch threshold';
         voltage.identify_ap_templates.threshold.glitch_threshold.descript    = 'spikes larger than this are considered a giltch';
         voltage.identify_ap_templates.threshold.glitch_threshold.type        = 'positive float';
         voltage.identify_ap_templates.threshold.glitch_threshold.units       = 'avg spike amplitude';
        
         voltage.identify_ap_templates.threshold.min_positive_duration.value    = 1;
         voltage.identify_ap_templates.threshold.min_positive_duration.name     = 'min positive duration';
         voltage.identify_ap_templates.threshold.min_positive_duration.descript = 'min duration of positive component of spike (0 to ignore)';
         voltage.identify_ap_templates.threshold.min_positive_duration.type     = 'positive float';
         voltage.identify_ap_templates.threshold.min_positive_duration.units    = 'ms';
        
         voltage.identify_ap_templates.threshold.min_negative_duration.value    = 1;
         voltage.identify_ap_templates.threshold.min_negative_duration.name     = 'min negative duration';
         voltage.identify_ap_templates.threshold.min_negative_duration.descript = 'minimum duration for negative component of spike (0 to ignore)';
         voltage.identify_ap_templates.threshold.min_negative_duration.type     = 'positive float';
         voltage.identify_ap_templates.threshold.min_negative_duration.units    = 'ms';
   
         voltage.identify_ap_templates.threshold.avg_window.value               = 10; % MA window in seconds
         voltage.identify_ap_templates.threshold.avg_window.name                = 'averaging window';
         voltage.identify_ap_templates.threshold.avg_window.descript            = 'length of moving average window for estimating std dev to use in thresholding';
         voltage.identify_ap_templates.threshold.avg_window.type                = 'positive float';
         voltage.identify_ap_templates.threshold.avg_window.units               = 'seconds';
         
         voltage.identify_ap_templates.threshold.skip_window.value              = 5; % skip ahead window in seconds
         voltage.identify_ap_templates.threshold.skip_window.name               = 'skip ahead window';
         voltage.identify_ap_templates.threshold.skip_window.descript           = 'amount to skip ahead to next avg when estimating std dev to use in thresholding';
         voltage.identify_ap_templates.threshold.skip_window.type               = 'positive float';
         voltage.identify_ap_templates.threshold.skip_window.units              = 'seconds';

         voltage.identify_ap_templates.threshold.match_type.value               = 'corr';
         voltage.identify_ap_templates.threshold.match_type.name                = 'match type';
         voltage.identify_ap_templates.threshold.match_type.descript            = 'statistic to calculate similarity with';
         voltage.identify_ap_templates.threshold.match_type.type                = 'list';
         voltage.identify_ap_templates.threshold.match_type.units               = 'statistic';
         voltage.identify_ap_templates.threshold.match_type.list                = {'corr','cov'};
         
         voltage.identify_ap_templates.threshold.allow_diff_lengths.value       = true;
         voltage.identify_ap_templates.threshold.allow_diff_lengths.name        = 'allow diff lengths';
         voltage.identify_ap_templates.threshold.allow_diff_lengths.descript    = 'allow different length spikes and family templates if true';
         voltage.identify_ap_templates.threshold.allow_diff_lengths.type        = 'boolean';
         voltage.identify_ap_templates.threshold.allow_diff_lengths.units       = 'true or false';
         
         voltage.identify_ap_templates.threshold.match_similarity.value         = 0.95;
         voltage.identify_ap_templates.threshold.match_similarity.name          = 'match similarity';
         voltage.identify_ap_templates.threshold.match_similarity.descript      = 'similarity required for spike shapes to be considered from the same template';
         voltage.identify_ap_templates.threshold.match_similarity.type          = 'normalised float';
         voltage.identify_ap_templates.threshold.match_similarity.units         = '\in [0, 1]';
   
         voltage.identify_ap_templates.threshold.remove_small_templates.value   = 100;
         voltage.identify_ap_templates.threshold.remove_small_templates.name    = 'remove small templates';
         voltage.identify_ap_templates.threshold.remove_small_templates.descript= 'remove AP template estimates from this many or less spikes (0 to ignore)';
         voltage.identify_ap_templates.threshold.remove_small_templates.type    = 'positive integer';
         voltage.identify_ap_templates.threshold.remove_small_templates.units   = 'integer';
   
         voltage.identify_ap_templates.threshold.normalise_aps.value            = true;
         voltage.identify_ap_templates.threshold.normalise_aps.name             = 'normalise APs';
         voltage.identify_ap_templates.threshold.normalise_aps.descript         = 'normalise APs so all spikes have the same amplitude when averagig for the AP template';
         voltage.identify_ap_templates.threshold.normalise_aps.type             = 'boolean';
         voltage.identify_ap_templates.threshold.normalise_aps.units            = 'true or false';

%          voltage.identify_ap_templates.threshold.reset_time.value            = 0;
%          voltage.identify_ap_templates.threshold.reset_time.name             = 'reset time';
%          voltage.identify_ap_templates.threshold.reset_time.descript         = 'min time before this axon can spike again';
%          voltage.identify_ap_templates.threshold.reset_time.type             = 'positive float';
%          voltage.identify_ap_templates.threshold.reset_time.units            = 'ms, 0 to ignore';
   
         voltage.identify_ap_templates.wavelets.wavelet_threshold.value         = 3;
         voltage.identify_ap_templates.wavelets.wavelet_threshold.name          = 'wavelet threshold';
         voltage.identify_ap_templates.wavelets.wavelet_threshold.descript      = 'only wavelet coeffs with magnitude larger than this many std devs are retained';
         voltage.identify_ap_templates.wavelets.wavelet_threshold.type          = 'positive float';
         voltage.identify_ap_templates.wavelets.wavelet_threshold.units         = 'std dev';
         
         voltage.identify_ap_templates.wavelets.lower_frequency.value           = 60;
         voltage.identify_ap_templates.wavelets.lower_frequency.name            = 'lower frequency';
         voltage.identify_ap_templates.wavelets.lower_frequency.descript        = 'only wavelet coefficients for wavelets with frequencies higher than this are included (0 to automate)';
         voltage.identify_ap_templates.wavelets.lower_frequency.type            = 'positive float';
         voltage.identify_ap_templates.wavelets.lower_frequency.units           = 'Hz';
        
         voltage.identify_ap_templates.wavelets.upper_frequency.value           = 350;
         voltage.identify_ap_templates.wavelets.upper_frequency.name            = 'upper frequency';
         voltage.identify_ap_templates.wavelets.upper_frequency.descript        = 'only wavelet coefficients for wavelets with frequencies lower than this are included (0 to automate)';
         voltage.identify_ap_templates.wavelets.upper_frequency.type            = 'positive float';
         voltage.identify_ap_templates.wavelets.upper_frequency.units           = 'Hz';
        
         voltage.identify_ap_templates.wavelets.match_type.value                = 'corr';
         voltage.identify_ap_templates.wavelets.match_type.name                 = 'match type';
         voltage.identify_ap_templates.wavelets.match_type.descript             = 'statistic to calculate similarity with';
         voltage.identify_ap_templates.wavelets.match_type.type                 = 'list';
         voltage.identify_ap_templates.wavelets.match_type.units                = 'statistic';
         voltage.identify_ap_templates.wavelets.match_type.list                 = {'corr','cov'};
         
         voltage.identify_ap_templates.wavelets.match_similarity.value          = 0.9;
         voltage.identify_ap_templates.wavelets.match_similarity.name           = 'match similarity';
         voltage.identify_ap_templates.wavelets.match_similarity.descript       = 'similarity required for spike shapes to be considered from the same template';
         voltage.identify_ap_templates.wavelets.match_similarity.type           = 'normalised float';
         voltage.identify_ap_templates.wavelets.match_similarity.units          = '\in [0, 1]';

         voltage.identify_ap_templates.wavelets.normalise_aps.value             = true;
         voltage.identify_ap_templates.wavelets.normalise_aps.name              = 'normalise APs';
         voltage.identify_ap_templates.wavelets.normalise_aps.descript          = 'normalise APs so all spikes have the same amplitude when averagig for the AP template';
         voltage.identify_ap_templates.wavelets.normalise_aps.type              = 'boolean';
         voltage.identify_ap_templates.wavelets.normalise_aps.units             = 'true or false';

         voltage.identify_ap_templates.wavelets.remove_small_templates.value    = 2;
         voltage.identify_ap_templates.wavelets.remove_small_templates.name     = 'remove small templates';
         voltage.identify_ap_templates.wavelets.remove_small_templates.descript = 'remove AP template estimates from this many or less spikes (0 to ignore)';
         voltage.identify_ap_templates.wavelets.remove_small_templates.type     = 'positive integer';
         voltage.identify_ap_templates.wavelets.remove_small_templates.units    = 'integer';
            
         % Extract spikes by k-means clustering         
         voltage.identify_ap_templates.k_means.match_type.value                 = 'corr';
         voltage.identify_ap_templates.k_means.match_type.name                  = 'match type';
         voltage.identify_ap_templates.k_means.match_type.descript              = 'statistic to calculate similarity with';
         voltage.identify_ap_templates.k_means.match_type.type                  = 'list';
         voltage.identify_ap_templates.k_means.match_type.units                 = 'statistic';
         voltage.identify_ap_templates.k_means.match_type.list                  = {'corr','cov'};

         voltage.identify_ap_templates.k_means.positive_threshold.value         = 4;
         voltage.identify_ap_templates.k_means.positive_threshold.name          = 'positive threshold';
         voltage.identify_ap_templates.k_means.positive_threshold.descript      = 'only spikes with peak voltage larger than this many std devs considered';
         voltage.identify_ap_templates.k_means.positive_threshold.type          = 'positive float';
         voltage.identify_ap_templates.k_means.positive_threshold.units         = 'std dev';
         
         voltage.identify_ap_templates.k_means.negative_threshold.value         = 0;
         voltage.identify_ap_templates.k_means.negative_threshold.name          = 'negative threshold magnitude';
         voltage.identify_ap_templates.k_means.negative_threshold.descript      = 'only spikes with a negative voltage of magnitude larger than this considered (0 to ignore)';
         voltage.identify_ap_templates.k_means.negative_threshold.type          = 'positive float';
         voltage.identify_ap_templates.k_means.negative_threshold.units         = 'std dev';

         voltage.identify_ap_templates.k_means.min_positive_duration.value      = 1;
         voltage.identify_ap_templates.k_means.min_positive_duration.name       = 'min positive duration';
         voltage.identify_ap_templates.k_means.min_positive_duration.descript   = 'min duration of positive component of spike';
         voltage.identify_ap_templates.k_means.min_positive_duration.type       = 'positive float';
         voltage.identify_ap_templates.k_means.min_positive_duration.units      = 'ms';
         
         voltage.identify_ap_templates.k_means.min_negative_duration.value      = 0.1;
         voltage.identify_ap_templates.k_means.min_negative_duration.name       = 'min negative duration';
         voltage.identify_ap_templates.k_means.min_negative_duration.descript   = 'min duration of negative component of spike';
         voltage.identify_ap_templates.k_means.min_negative_duration.type       = 'positive float';
         voltage.identify_ap_templates.k_means.min_negative_duration.units      = 'ms';
         
         voltage.identify_ap_templates.k_means.avg_window.value                 = 10; % MA window in seconds
         voltage.identify_ap_templates.k_means.avg_window.name                  = 'averaging window';
         voltage.identify_ap_templates.k_means.avg_window.descript              = 'length of moving average window for estimating std dev to use in thresholding';
         voltage.identify_ap_templates.k_means.avg_window.type                  = 'positive float';
         voltage.identify_ap_templates.k_means.avg_window.units                 = 'seconds';
         
         voltage.identify_ap_templates.k_means.skip_window.value                = 5; % skip ahead window in seconds
         voltage.identify_ap_templates.k_means.skip_window.name                 = 'skip ahead window';
         voltage.identify_ap_templates.k_means.skip_window.descript             = 'amount to skip ahead to next avg when estimating std dev to use in thresholding';
         voltage.identify_ap_templates.k_means.skip_window.type                 = 'positive float';
         voltage.identify_ap_templates.k_means.skip_window.units                = 'seconds';
         
         voltage.identify_ap_templates.k_means.match_similarity.value           = 0.8;
         voltage.identify_ap_templates.k_means.match_similarity.name            = 'match similarity';
         voltage.identify_ap_templates.k_means.match_similarity.descript        = 'similarity required for spike shapes to be considered from the same template';
         voltage.identify_ap_templates.k_means.match_similarity.type            = 'normalised float';
         voltage.identify_ap_templates.k_means.match_similarity.units           = '\in [0, 1]';

         voltage.identify_ap_templates.k_means.normalise_aps.value              = true;
         voltage.identify_ap_templates.k_means.normalise_aps.name               = 'normalise APs';
         voltage.identify_ap_templates.k_means.normalise_aps.descript           = 'normalise APs so constituent spikes have the same weight when calculating avg';
         voltage.identify_ap_templates.k_means.normalise_aps.type               = 'boolean';
         voltage.identify_ap_templates.k_means.normalise_aps.units              = 'true or false';
         
         voltage.identify_ap_templates.k_means.remove_small_templates.value     = 2;
         voltage.identify_ap_templates.k_means.remove_small_templates.name      = 'remove small templates';
         voltage.identify_ap_templates.k_means.remove_small_templates.descript  = 'remove AP template estimates from this many or less spikes (0 to ignore)';
         voltage.identify_ap_templates.k_means.remove_small_templates.type      = 'positive integer';
         voltage.identify_ap_templates.k_means.remove_small_templates.units     = 'integer';
         
%          voltage.identify_ap_templates.wavelets.min_positive_duration.value     = 5;
%          voltage.identify_ap_templates.wavelets.min_positive_duration.name      = 'min positive duration';
%          voltage.identify_ap_templates.wavelets.min_positive_duration.descript  = 'min consecutive duration of positive voltages to be retained (0 to ignore)';
%          voltage.identify_ap_templates.wavelets.min_positive_duration.type      = 'positive float';
%          voltage.identify_ap_templates.wavelets.min_positive_duration.units     = 'ms';
%          
%          voltage.identify_ap_templates.wavelets.min_negative_duration.value     = 5;
%          voltage.identify_ap_templates.wavelets.min_negative_duration.name      = 'min negative duration';
%          voltage.identify_ap_templates.wavelets.min_negative_duration.descript  = 'min consecutive duration of positive voltages to be retained (0 to ignore)';
%          voltage.identify_ap_templates.wavelets.min_negative_duration.type      = 'positive float';
%          voltage.identify_ap_templates.wavelets.min_negative_duration.units     = 'ms';
%          
%          voltage.identify_ap_templates.wavelets.timeposneg.value                = 0;
%          voltage.identify_ap_templates.wavelets.timeposneg.name                 = 'min total duration';
%          voltage.identify_ap_templates.wavelets.timeposneg.descript             = 'negative coefficient must occur within this time to be considered part of same spike (0 to ignore)';
%          voltage.identify_ap_templates.wavelets.timeposneg.type                 = 'positive float';
%          voltage.identify_ap_templates.wavelets.timeposneg.units                = 'ms';

   	% Get spikes - by thresholding (get templates & extract spikes all in one go)
         voltage.extract_spikes.matched_filter.match_type.value                 = 'corr';
         voltage.extract_spikes.matched_filter.match_type.name                  = 'match type';
         voltage.extract_spikes.matched_filter.match_type.descript              = 'statistic to calculate similarity with';
         voltage.extract_spikes.matched_filter.match_type.type                  = 'list';
         voltage.extract_spikes.matched_filter.match_type.units                 = 'statistic';
         voltage.extract_spikes.matched_filter.match_type.list                  = {'corr','cov'};
         
         voltage.extract_spikes.matched_filter.match_similarity.value           = 0.95;
         voltage.extract_spikes.matched_filter.match_similarity.name            = 'match similarity';
         voltage.extract_spikes.matched_filter.match_similarity.descript        = 'similarity required for spike shapes to be considered from the same family';
         voltage.extract_spikes.matched_filter.match_similarity.type            = 'normalised float';
         voltage.extract_spikes.matched_filter.match_similarity.units           = '\in [0, 1]';

         voltage.extract_spikes.matched_filter.ap_peak_change.value             = 20;
         voltage.extract_spikes.matched_filter.ap_peak_change.name              = 'AP peak change';
         voltage.extract_spikes.matched_filter.ap_peak_change.descript          = 'how quickly can AP peak magnitude change with time';
         voltage.extract_spikes.matched_filter.ap_peak_change.type              = 'percentage';
         voltage.extract_spikes.matched_filter.ap_peak_change.units             = '% per second';
         
         voltage.extract_spikes.matched_filter.forgetting_factor.value          = 0.9;
         voltage.extract_spikes.matched_filter.forgetting_factor.name           = 'forgetting factor';
         voltage.extract_spikes.matched_filter.forgetting_factor.descript       = 'Controls how quickly we forget past samples (0 - forget fastest, 1 - forget slowest)';
         voltage.extract_spikes.matched_filter.forgetting_factor.type           = 'normalised float';
         voltage.extract_spikes.matched_filter.forgetting_factor.units          = '\in [0, 1]';

         voltage.extract_spikes.matched_filter.kappa_pos.value                  = 3;
         voltage.extract_spikes.matched_filter.kappa_pos.name                   = 'AP peak change';
         voltage.extract_spikes.matched_filter.kappa_pos.descript               = 'how much can a positive spike peak differ from the mean peak value';
         voltage.extract_spikes.matched_filter.kappa_pos.type                   = 'positive float';
         voltage.extract_spikes.matched_filter.kappa_pos.units                  = 'standard deviations';

         voltage.extract_spikes.matched_filter.kappa_neg.value                  = 3;
         voltage.extract_spikes.matched_filter.kappa_neg.name                   = 'AP peak change';
         voltage.extract_spikes.matched_filter.kappa_neg.descript               = 'how much can a negative spike peak differ from the mean negative peak value';
         voltage.extract_spikes.matched_filter.kappa_neg.type                   = 'positive float';
         voltage.extract_spikes.matched_filter.kappa_neg.units                  = 'standard deviations';

         voltage.extract_spikes.matched_filter.positive_threshold.value         = 5;
         voltage.extract_spikes.matched_filter.positive_threshold.name          = 'positive threshold';
         voltage.extract_spikes.matched_filter.positive_threshold.descript      = 'only spikes with peak voltage larger than this many std devs considered';
         voltage.extract_spikes.matched_filter.positive_threshold.type          = 'positive float';
         voltage.extract_spikes.matched_filter.positive_threshold.units         = 'std dev';
         
         voltage.extract_spikes.matched_filter.negative_threshold.value         = 5;
         voltage.extract_spikes.matched_filter.negative_threshold.name          = 'negative threshold magnitude';
         voltage.extract_spikes.matched_filter.negative_threshold.descript      = 'only spikes with a negative voltage of magnitude larger than this considered (0 to ignore)';
         voltage.extract_spikes.matched_filter.negative_threshold.type          = 'positive float';
         voltage.extract_spikes.matched_filter.negative_threshold.units         = 'std dev';
         
         voltage.extract_spikes.threshold.glitch_threshold.value                = 10;
         voltage.extract_spikes.threshold.glitch_threshold.name                 = 'glitch threshold';
         voltage.extract_spikes.threshold.glitch_threshold.descript             = 'spikes larger than this threshold times the average spike magnitude are considered a giltch';
         voltage.extract_spikes.threshold.glitch_threshold.type                 = 'positive float';
         voltage.extract_spikes.threshold.glitch_threshold.units                = 'avg spike amplitude';
         
         voltage.extract_spikes.matched_filter.avg_window.value                 = 10; % MA window in seconds
         voltage.extract_spikes.matched_filter.avg_window.name                  = 'averaging window';
         voltage.extract_spikes.matched_filter.avg_window.descript              = 'length of moving average window for estimating std dev to use in thresholding';
         voltage.extract_spikes.matched_filter.avg_window.type                  = 'positive float';
         voltage.extract_spikes.matched_filter.avg_window.units                 = 'seconds';
         
         voltage.extract_spikes.matched_filter.skip_window.value                = 5; % skip ahead window in seconds
         voltage.extract_spikes.matched_filter.skip_window.name                 = 'skip ahead window';
         voltage.extract_spikes.matched_filter.skip_window.descript             = 'amount to skip ahead to next avg when estimating std dev to use in thresholding';
         voltage.extract_spikes.matched_filter.skip_window.type                 = 'positive float';
         voltage.extract_spikes.matched_filter.skip_window.units                = 'seconds';

         voltage.extract_spikes.matched_filter.min_positive_duration.value      = 2;
         voltage.extract_spikes.matched_filter.min_positive_duration.name       = 'min positive duration';
         voltage.extract_spikes.matched_filter.min_positive_duration.descript   = 'min duration of positive component of spike';
         voltage.extract_spikes.matched_filter.min_positive_duration.type       = 'positive float';
         voltage.extract_spikes.matched_filter.min_positive_duration.units      = 'ms';
         
         voltage.extract_spikes.matched_filter.min_negative_duration.value      = 2;
         voltage.extract_spikes.matched_filter.min_negative_duration.name       = 'min negative duration';
         voltage.extract_spikes.matched_filter.min_negative_duration.descript   = 'min duration of negative component of spike';
         voltage.extract_spikes.matched_filter.min_negative_duration.type       = 'positive float';
         voltage.extract_spikes.matched_filter.min_negative_duration.units      = 'ms';
         
         voltage.extract_spikes.matched_filter.allow_new_aps.value              = false;
         voltage.extract_spikes.matched_filter.allow_new_aps.name               = 'allow new aps';
         voltage.extract_spikes.matched_filter.allow_new_aps.descript           = 'allow finding new AP templates';
         voltage.extract_spikes.matched_filter.allow_new_aps.type               = 'boolean';
         voltage.extract_spikes.matched_filter.allow_new_aps.units              = 'true or false';
         
         voltage.extract_spikes.matched_filter.min_spiking_threshold.value      = 100;
         voltage.extract_spikes.matched_filter.min_spiking_threshold.name       = 'min spiking threshold';
         voltage.extract_spikes.matched_filter.min_spiking_threshold.descript   = 'remove a family with number of spikes less than this threshold';
         voltage.extract_spikes.matched_filter.min_spiking_threshold.type       = 'positive integer';
         voltage.extract_spikes.matched_filter.min_spiking_threshold.units      = 'number of spikes';
            
         voltage.extract_spikes.matched_filter.allow_diff_lengths.value         = true;
         voltage.extract_spikes.matched_filter.allow_diff_lengths.name          = 'allow diff lengths';
         voltage.extract_spikes.matched_filter.allow_diff_lengths.descript      = 'allow different length spikes and family templates if true';
         voltage.extract_spikes.matched_filter.allow_diff_lengths.type          = 'boolean';
         voltage.extract_spikes.matched_filter.allow_diff_lengths.units         = 'true or false';
         
         voltage.extract_spikes.matched_filter.remove_small_templates.value     = 100;
         voltage.extract_spikes.matched_filter.remove_small_templates.name      = 'remove small templates';
         voltage.extract_spikes.matched_filter.remove_small_templates.descript  = 'remove AP template estimates from this many or less spikes (0 to ignore)';
         voltage.extract_spikes.matched_filter.remove_small_templates.type      = 'positive integer';
         voltage.extract_spikes.matched_filter.remove_small_templates.units     = 'integer';
   
         voltage.extract_spikes.matched_filter.normalise_aps.value              = true;
         voltage.extract_spikes.matched_filter.normalise_aps.name               = 'normalise APs';
         voltage.extract_spikes.matched_filter.normalise_aps.descript           = 'normalise APs so constituent spikes have the same weight when calculating avg';
         voltage.extract_spikes.matched_filter.normalise_aps.type               = 'boolean';
         voltage.extract_spikes.matched_filter.normalise_aps.units              = 'true or false';

         voltage.extract_spikes.matched_filter.plot_spikes_consecutively.value  = false;
         voltage.extract_spikes.matched_filter.plot_spikes_consecutively.name   = 'plot spikes consecutively';
         voltage.extract_spikes.matched_filter.plot_spikes_consecutively.descript= 'plot resulting spikes WITHOUT inactive time between them';
         voltage.extract_spikes.matched_filter.plot_spikes_consecutively.type   = 'boolean';
         voltage.extract_spikes.matched_filter.plot_spikes_consecutively.units  = 'true or false';

         % Extract spikes by k-means clustering
         voltage.extract_spikes.k_means.no_clusters.value                       = 0;
         voltage.extract_spikes.k_means.no_clusters.name                        = 'number of clusters';
         voltage.extract_spikes.k_means.no_clusters.descript                    = 'Clusters for the K-means, 0 for automatic';
         voltage.extract_spikes.k_means.no_clusters.type                        = 'positive integer';
         voltage.extract_spikes.k_means.no_clusters.units                       = 'integer';
         
         voltage.extract_spikes.k_means.match_type.value                        = 'corr';
         voltage.extract_spikes.k_means.match_type.name                        = 'match type';
         voltage.extract_spikes.k_means.match_type.descript                    = 'statistic to calculate similarity with';
         voltage.extract_spikes.k_means.match_type.type                        = 'list';
         voltage.extract_spikes.k_means.match_type.units                       = 'statistic';
         voltage.extract_spikes.k_means.match_type.list                        = {'corr','cov'};
         
         voltage.extract_spikes.k_means.match_similarity.value                 = 0.8;
         voltage.extract_spikes.k_means.match_similarity.name                  = 'match similarity';
         voltage.extract_spikes.k_means.match_similarity.descript              = 'similarity required for spike shapes to be considered from the same family';
         voltage.extract_spikes.k_means.match_similarity.type                  = 'normalised float';
         voltage.extract_spikes.k_means.match_similarity.units                 = '\in [0, 1]';

         voltage.extract_spikes.k_means.ap_peak_change.value                   = 20;
         voltage.extract_spikes.k_means.ap_peak_change.name                    = 'AP peak change';
         voltage.extract_spikes.k_means.ap_peak_change.descript                = 'how quickly can AP peak magnitude change with time';
         voltage.extract_spikes.k_means.ap_peak_change.type                    = 'percentage';
         voltage.extract_spikes.k_means.ap_peak_change.units                   = '% per second';

         voltage.extract_spikes.k_means.positive_threshold.value               = 4;
         voltage.extract_spikes.k_means.positive_threshold.name                = 'positive threshold';
         voltage.extract_spikes.k_means.positive_threshold.descript            = 'only spikes with peak voltage larger than this many std devs considered';
         voltage.extract_spikes.k_means.positive_threshold.type                = 'positive float';
         voltage.extract_spikes.k_means.positive_threshold.units               = 'std dev';
         
         voltage.extract_spikes.k_means.negative_threshold.value               = 0;
         voltage.extract_spikes.k_means.negative_threshold.name                = 'negative threshold magnitude';
         voltage.extract_spikes.k_means.negative_threshold.descript            = 'only spikes with a negative voltage of magnitude larger than this considered (0 to ignore)';
         voltage.extract_spikes.k_means.negative_threshold.type                = 'positive float';
         voltage.extract_spikes.k_means.negative_threshold.units               = 'std dev';
         
         voltage.extract_spikes.k_means.glitch_threshold.value                 = 10;
         voltage.extract_spikes.k_means.glitch_threshold.name                  = 'glitch threshold';
         voltage.extract_spikes.k_means.glitch_threshold.descript              = 'spikes larger than this are considered a giltch';
         voltage.extract_spikes.k_means.glitch_threshold.type                  = 'positive float';
         voltage.extract_spikes.k_means.glitch_threshold.units                 = 'avg spike amplitude';
         
         voltage.extract_spikes.k_means.avg_window.value                       = 10; % MA window in seconds
         voltage.extract_spikes.k_means.avg_window.name                        = 'averaging window';
         voltage.extract_spikes.k_means.avg_window.descript                    = 'length of moving average window for estimating std dev to use in thresholding';
         voltage.extract_spikes.k_means.avg_window.type                        = 'positive float';
         voltage.extract_spikes.k_means.avg_window.units                       = 'seconds';
         
         voltage.extract_spikes.k_means.skip_window.value                      = 5; % skip ahead window in seconds
         voltage.extract_spikes.k_means.skip_window.name                       = 'skip ahead window';
         voltage.extract_spikes.k_means.skip_window.descript                   = 'amount to skip ahead to next avg when estimating std dev to use in thresholding';
         voltage.extract_spikes.k_means.skip_window.type                       = 'positive float';
         voltage.extract_spikes.k_means.skip_window.units                      = 'seconds';

         voltage.extract_spikes.k_means.min_positive_duration.value            = 1;
         voltage.extract_spikes.k_means.min_positive_duration.name             = 'min positive duration';
         voltage.extract_spikes.k_means.min_positive_duration.descript         = 'min duration of positive component of spike';
         voltage.extract_spikes.k_means.min_positive_duration.type             = 'positive float';
         voltage.extract_spikes.k_means.min_positive_duration.units            = 'ms';
         
         voltage.extract_spikes.k_means.min_negative_duration.value            = 0.1;
         voltage.extract_spikes.k_means.min_negative_duration.name             = 'min negative duration';
         voltage.extract_spikes.k_means.min_negative_duration.descript         = 'min duration of negative component of spike';
         voltage.extract_spikes.k_means.min_negative_duration.type             = 'positive float';
         voltage.extract_spikes.k_means.min_negative_duration.units            = 'ms';
         
         voltage.extract_spikes.k_means.allow_new_aps.value                    = true;
         voltage.extract_spikes.k_means.allow_new_aps.name                     = 'allow new aps';
         voltage.extract_spikes.k_means.allow_new_aps.descript                 = 'allow finding new AP templates';
         voltage.extract_spikes.k_means.allow_new_aps.type                     = 'boolean';
         voltage.extract_spikes.k_means.allow_new_aps.units                    = 'true or false';
         
         voltage.extract_spikes.k_means.remove_small_templates.value           = 100;
         voltage.extract_spikes.k_means.remove_small_templates.name            = 'remove small templates';
         voltage.extract_spikes.k_means.remove_small_templates.descript        = 'remove AP template estimates from this many or less spikes (0 to ignore)';
         voltage.extract_spikes.k_means.remove_small_templates.type            = 'positive integer';
         voltage.extract_spikes.k_means.remove_small_templates.units           = 'integer';
   
         voltage.extract_spikes.k_means.normalise_aps.value                    = true;
         voltage.extract_spikes.k_means.normalise_aps.name                     = 'normalise APs';
         voltage.extract_spikes.k_means.normalise_aps.descript                 = 'normalise APs so constituent spikes have the same weight when calculating avg';
         voltage.extract_spikes.k_means.normalise_aps.type                     = 'boolean';
         voltage.extract_spikes.k_means.normalise_aps.units                    = 'true or false';

         voltage.extract_spikes.k_means.plot_spikes_consecutively.value        = false;
         voltage.extract_spikes.k_means.plot_spikes_consecutively.name         = 'plot spikes consecutively';
         voltage.extract_spikes.k_means.plot_spikes_consecutively.descript     = 'plot resulting spikes without inactive time between them';
         voltage.extract_spikes.k_means.plot_spikes_consecutively.type         = 'boolean';
         voltage.extract_spikes.k_means.plot_spikes_consecutively.units        = 'true or false';
   %% AP template data
   	% Get spikes - using matched filtering
         ap.extract_spikes.matched_filter.match_type.value                      = 'corr';
         ap.extract_spikes.matched_filter.match_type.name                       = 'match type';
         ap.extract_spikes.matched_filter.match_type.descript                   = 'statistic to calculate similarity with';
         ap.extract_spikes.matched_filter.match_type.type                       = 'list';
         ap.extract_spikes.matched_filter.match_type.units                      = 'statistic';
         ap.extract_spikes.matched_filter.match_type.list                       = {'corr','cov'};
         
         ap.extract_spikes.matched_filter.match_similarity.value                = 0.7;
         ap.extract_spikes.matched_filter.match_similarity.name                 = 'match similarity';
         ap.extract_spikes.matched_filter.match_similarity.descript             = 'similarity required for spike shapes to be considered from the same family';
         ap.extract_spikes.matched_filter.match_similarity.type                 = 'normalised float';
         ap.extract_spikes.matched_filter.match_similarity.units                = '\in [0, 1]';

         ap.extract_spikes.matched_filter.forgetting_factor.value               = 0.98;
         ap.extract_spikes.matched_filter.forgetting_factor.name                = 'forgetting factor';
         ap.extract_spikes.matched_filter.forgetting_factor.descript            = 'Controls how quickly we forget past samples (0 - forget fastest, 1 - forget slowest)';
         ap.extract_spikes.matched_filter.forgetting_factor.type                = 'normalised float';
         ap.extract_spikes.matched_filter.forgetting_factor.units               = '\in [0, 1]';
         
         ap.extract_spikes.matched_filter.kappa_pos.value                       = 4;
         ap.extract_spikes.matched_filter.kappa_pos.name                        = 'AP peak change';
         ap.extract_spikes.matched_filter.kappa_pos.descript                    = 'how much can a positive peak differ from the mean peak value';
         ap.extract_spikes.matched_filter.kappa_pos.type                        = 'positive float';
         ap.extract_spikes.matched_filter.kappa_pos.units                       = 'standard deviations';

         ap.extract_spikes.matched_filter.kappa_neg.value                       = 4;
         ap.extract_spikes.matched_filter.kappa_neg.name                        = 'AP peak change';
         ap.extract_spikes.matched_filter.kappa_neg.descript                    = 'how much can a negative peak differ from the mean negative peak value';
         ap.extract_spikes.matched_filter.kappa_neg.type                        = 'positive float';
         ap.extract_spikes.matched_filter.kappa_neg.units                       = 'standard deviations';

         ap.extract_spikes.matched_filter.positive_threshold.value              = 5;
         ap.extract_spikes.matched_filter.positive_threshold.name               = 'positive threshold';
         ap.extract_spikes.matched_filter.positive_threshold.descript           = 'only spikes with peak voltage larger than this many std devs considered';
         ap.extract_spikes.matched_filter.positive_threshold.type               = 'positive float';
         ap.extract_spikes.matched_filter.positive_threshold.units              = 'std dev';
         
         ap.extract_spikes.matched_filter.negative_threshold.value              = 4;
         ap.extract_spikes.matched_filter.negative_threshold.name               = 'negative threshold magnitude';
         ap.extract_spikes.matched_filter.negative_threshold.descript           = 'only spikes with a negative voltage of magnitude larger than this considered (0 to ignore)';
         ap.extract_spikes.matched_filter.negative_threshold.type               = 'positive float';
         ap.extract_spikes.matched_filter.negative_threshold.units              = 'std dev';
         
         ap.extract_spikes.matched_filter.glitch_threshold.value                = 15;
         ap.extract_spikes.matched_filter.glitch_threshold.name                 = 'glitch threshold';
         ap.extract_spikes.matched_filter.glitch_threshold.descript             = 'spikes larger than this are considered a giltch';
         ap.extract_spikes.matched_filter.glitch_threshold.type                 = 'positive float';
         ap.extract_spikes.matched_filter.glitch_threshold.units                = 'avg spike amplitude';
         
         ap.extract_spikes.matched_filter.avg_window.value                      = 6; % MA window in seconds
         ap.extract_spikes.matched_filter.avg_window.name                       = 'averaging window';
         ap.extract_spikes.matched_filter.avg_window.descript                   = 'length of moving average window for estimating std dev to use in thresholding';
         ap.extract_spikes.matched_filter.avg_window.type                       = 'positive float';
         ap.extract_spikes.matched_filter.avg_window.units                      = 'seconds';
         
         ap.extract_spikes.matched_filter.skip_window.value                     = 3; % skip ahead window in seconds
         ap.extract_spikes.matched_filter.skip_window.name                      = 'skip ahead window';
         ap.extract_spikes.matched_filter.skip_window.descript                  = 'amount to skip ahead to next avg when estimating std dev to use in thresholding';
         ap.extract_spikes.matched_filter.skip_window.type                      = 'positive float';
         ap.extract_spikes.matched_filter.skip_window.units                     = 'seconds';

         ap.extract_spikes.matched_filter.min_positive_duration.value           = 2;
         ap.extract_spikes.matched_filter.min_positive_duration.name            = 'min positive duration';
         ap.extract_spikes.matched_filter.min_positive_duration.descript        = 'min duration of positive component of spike';
         ap.extract_spikes.matched_filter.min_positive_duration.type            = 'positive float';
         ap.extract_spikes.matched_filter.min_positive_duration.units           = 'ms';
         
         ap.extract_spikes.matched_filter.min_negative_duration.value           = 2;
         ap.extract_spikes.matched_filter.min_negative_duration.name            = 'min negative duration';
         ap.extract_spikes.matched_filter.min_negative_duration.descript        = 'min duration of negative component of spike';
         ap.extract_spikes.matched_filter.min_negative_duration.type            = 'positive float';
         ap.extract_spikes.matched_filter.min_negative_duration.units           = 'ms';
         
         ap.extract_spikes.matched_filter.allow_new_aps.value                   = false;
         ap.extract_spikes.matched_filter.allow_new_aps.name                    = 'allow new aps';
         ap.extract_spikes.matched_filter.allow_new_aps.descript                = 'allow finding new AP templates';
         ap.extract_spikes.matched_filter.allow_new_aps.type                    = 'boolean';
         ap.extract_spikes.matched_filter.allow_new_aps.units                   = 'true or false';
         
         ap.extract_spikes.matched_filter.min_spiking_threshold.value           = 100;
         ap.extract_spikes.matched_filter.min_spiking_threshold.name            = 'min spiking threshold';
         ap.extract_spikes.matched_filter.min_spiking_threshold.descript        = 'remove a family with number of spikes less than this threshold';
         ap.extract_spikes.matched_filter.min_spiking_threshold.type            = 'positive integer';
         ap.extract_spikes.matched_filter.min_spiking_threshold.units           = 'number of spikes';
   
         ap.extract_spikes.matched_filter.plot_spikes_consecutively.value       = false;
         ap.extract_spikes.matched_filter.plot_spikes_consecutively.name        = 'plot spikes consecutively';
         ap.extract_spikes.matched_filter.plot_spikes_consecutively.descript    = 'plot resulting spikes without inactive time between them';
         ap.extract_spikes.matched_filter.plot_spikes_consecutively.type        = 'boolean';
         ap.extract_spikes.matched_filter.plot_spikes_consecutively.units       = 'true or false';
         
       % Merge spike templates
   
         ap.merge_templates.user_selection.template_to_merge_with.value         = 1;
         ap.merge_templates.user_selection.template_to_merge_with.name          = 'template to merge with';
         ap.merge_templates.user_selection.template_to_merge_with.descript      = 'template ID of template to merge into';
         ap.merge_templates.user_selection.template_to_merge_with.type          = 'positive integer';
         ap.merge_templates.user_selection.template_to_merge_with.units         = 'template ID';

         ap.merge_templates.user_selection.number_of_templates_to_merge.value   = 1;
         ap.merge_templates.user_selection.number_of_templates_to_merge.name    = 'number of templates to merge';
         ap.merge_templates.user_selection.number_of_templates_to_merge.descript= 'number of other templates to merge into chosen template';
         ap.merge_templates.user_selection.number_of_templates_to_merge.type    = 'positive integer';
         ap.merge_templates.user_selection.number_of_templates_to_merge.units   = 'templates';

         ap.delete_templates.user_selection.number_of_templates_to_remove.value = 1;
         ap.delete_templates.user_selection.number_of_templates_to_remove.name  = 'number of templates to remove';
         ap.delete_templates.user_selection.number_of_templates_to_remove.descript = 'number of templates to be removed';
         ap.delete_templates.user_selection.number_of_templates_to_remove.type  = 'positive integer';
         ap.delete_templates.user_selection.number_of_templates_to_remove.units = 'templates';

         
%       % Extract spikes by matching with AP templates as mother wavelets
%          ap.extract_spikes.wavelets.wavelet_coefficient_threshold.value         = 0.8;
%          ap.extract_spikes.wavelets.wavelet_coefficient_threshold.name          = 'wavelet coefficient threshold';
%          ap.extract_spikes.wavelets.wavelet_coefficient_threshold.descript      = 'only consider wavelett coeffs larger than this value';
%          ap.extract_spikes.wavelets.wavelet_coefficient_threshold.type          = 'normalised float';
%          ap.extract_spikes.wavelets.wavelet_coefficient_threshold.units         = 'wavelet coeffs';
%        
%          ap.extract_spikes.wavelets.ap_peak_change_rate.value                   = 20;
%          ap.extract_spikes.wavelets.ap_peak_change_rate.name                    = 'AP peak change rate';
%          ap.extract_spikes.wavelets.ap_peak_change_rate.descript                = 'how quickly can AP size change with time';
%          ap.extract_spikes.wavelets.ap_peak_change_rate.type                    = 'percentage';
%          ap.extract_spikes.wavelets.ap_peak_change_rate.units                   = '% change per second';
%        
%          ap.extract_spikes.wavelets.positive_threshold.value                    = 10;
%          ap.extract_spikes.wavelets.positive_threshold.name                     = 'positive threshold';
%          ap.extract_spikes.wavelets.positive_threshold.descript                 = 'only spikes with peak voltage larger than this many std devs considered';
%          ap.extract_spikes.wavelets.positive_threshold.type                     = 'positive float';
%          ap.extract_spikes.wavelets.positive_threshold.units                    = 'std dev';
%          
%          ap.extract_spikes.wavelets.negative_threshold.value                    = 0;
%          ap.extract_spikes.wavelets.negative_threshold.name                     = 'negative threshold magnitude';
%          ap.extract_spikes.wavelets.negative_threshold.descript                 = 'only spikes with a negative voltage of magnitude larger than this considered (0 to ignore)';
%          ap.extract_spikes.wavelets.negative_threshold.type                     = 'positive float';
%          ap.extract_spikes.wavelets.negative_threshold.units                    = 'std dev';
%         
%          ap.extract_spikes.wavelets.min_positive_duration.value                 = 40;
%          ap.extract_spikes.wavelets.min_positive_duration.name                  = 'min positive duration';
%          ap.extract_spikes.wavelets.min_positive_duration.descript              = 'min duration of positive component of spike';
%          ap.extract_spikes.wavelets.min_positive_duration.type                  = 'positive float';
%          ap.extract_spikes.wavelets.min_positive_duration.units                 = 'ms';
%          
%          ap.extract_spikes.wavelets.min_negative_duration.value                 = 0;
%          ap.extract_spikes.wavelets.min_negative_duration.name                  = 'min negative duration';
%          ap.extract_spikes.wavelets.min_negative_duration.descript              = 'negative voltage must occur within this time to be considered part of same spike (0 to ignore)';
%          ap.extract_spikes.wavelets.min_negative_duration.type                  = 'positive float';
%          ap.extract_spikes.wavelets.min_negative_duration.units                 = 'ms';

   %% Spike data
      % Firing rate
         spike.firing_rate.moving_average.moving_avg_window.value               = 10;
         spike.firing_rate.moving_average.moving_avg_window.name                = 'moving average window';
         spike.firing_rate.moving_average.moving_avg_window.descript            = 'duration of moving average window';
         spike.firing_rate.moving_average.moving_avg_window.type                = 'positive float';
         spike.firing_rate.moving_average.moving_avg_window.units               = 'seconds';
         
         spike.firing_rate.moving_average.skip_window.value                     = 2;
         spike.firing_rate.moving_average.skip_window.name                      = 'skip ahead window';
         spike.firing_rate.moving_average.skip_window.descript                  = 'skip ahead time to next averaging window';
         spike.firing_rate.moving_average.skip_window.type                      = 'positive float';
         spike.firing_rate.moving_average.skip_window.units                     = 'seconds';
         
         spike.firing_rate.moving_average.smoothing_kernel.value                = 'gaussian';
         spike.firing_rate.moving_average.smoothing_kernel.name                 = 'smooth spikes';
         spike.firing_rate.moving_average.smoothing_kernel.descript             = 'convolve spikes with kernel so rate estimate is smooth';
         spike.firing_rate.moving_average.smoothing_kernel.type                 = 'list';
         spike.firing_rate.moving_average.smoothing_kernel.list                 = {'none', 'gaussian', 'exponential', 'boxcar'};
         spike.firing_rate.moving_average.smoothing_kernel.units                = 'distribution';
         
         
      % Statistics
         spike.statistics.interspike_interval.num_bins.value                    = 20;
         spike.statistics.interspike_interval.num_bins.name                     = 'num bins';
         spike.statistics.interspike_interval.num_bins.descript                 = 'number of bins to use for ISI histogram'; 
         spike.statistics.interspike_interval.num_bins.type                     = 'positive integer';
         spike.statistics.interspike_interval.num_bins.units                    = 'integer';

         spike.statistics.interspike_interval.max_isi.value                     = 1;
         spike.statistics.interspike_interval.max_isi.name                      = 'max isi';
         spike.statistics.interspike_interval.max_isi.descript                  = 'max isi to consider (so tail doesn''t ruin plot)'; 
         spike.statistics.interspike_interval.max_isi.type                      = 'positive float';
         spike.statistics.interspike_interval.max_isi.units                     = 's';

         spike.statistics.amplitude_change.avg_window.value                     = 1;
         spike.statistics.amplitude_change.avg_window.name                      = 'moving average window';
         spike.statistics.amplitude_change.avg_window.descript                  = 'duration to average spike peaks over';
         spike.statistics.amplitude_change.avg_window.type                      = 'positive float';
         spike.statistics.amplitude_change.avg_window.units                     = 'seconds';

         spike.statistics.raster.no_params.value                                = 1;
         spike.statistics.raster.no_params.name                                = 'no params';
         spike.statistics.raster.no_params.descript                            = 'placeholder for future parameters';
         spike.statistics.raster.no_params.type                                = 'positive float';
         spike.statistics.raster.no_params.units                               = 'seconds';

         spike.statistics.psth.time_post_spike.value                            = 1;
         spike.statistics.psth.time_post_spike.name                             = 'time post spike';
         spike.statistics.psth.time_post_spike.descript                         = 'number of seconds postspike to build histogram from';
         spike.statistics.psth.time_post_spike.type                             = 'positive float';
         spike.statistics.psth.time_post_spike.units                            = 'seconds';

         spike.statistics.psth.num_bins.value                                   = 50;
         spike.statistics.psth.num_bins.name                                    = 'num bins';
         spike.statistics.psth.num_bins.descript                                = 'number of bins to build histogram from';
         spike.statistics.psth.num_bins.type                                    = 'positive integer';
         spike.statistics.psth.num_bins.units                                   = 'integer';

         spike.spike_operations.merge_spikes.within_template.value              = false;
         spike.spike_operations.merge_spikes.within_template.name               = 'within template';
         spike.spike_operations.merge_spikes.within_template.descript           = 'merge spike trains within each template if true, else merge all spike trains from all templates';
         spike.spike_operations.merge_spikes.within_template.type               = 'boolean';
         spike.spike_operations.merge_spikes.within_template.units              = 'true or false';
         
%          spike.export_to_excel.spike_rate_and_count.file_name.value              = false;
%          spike.export_to_excel.spike_rate_and_count.file_name.name               = 'within template';
%          spike.export_to_excel.spike_rate_and_count.file_name.descript           = 'merge spike trains within each template if true, else merge all spike trains from all templates';
%          spike.export_to_excel.spike_rate_and_count.file_name.type               = 'boolean';
%          spike.export_to_excel.spike_rate_and_count.file_name.units              = 'true or false';

   %% Rate data
      % Statistics
         rate.statistics.autocorrelation.max_lag.value                          = 10; % longer since spike rate has been smoother over time --> dt increased
         rate.statistics.autocorrelation.max_lag.name                           = 'mag lag';
         rate.statistics.autocorrelation.max_lag.descript                       = 'max lag considered between samples'; 
         rate.statistics.autocorrelation.max_lag.type                           = 'positive float';
         rate.statistics.autocorrelation.max_lag.units                          = 's';
      
   %% Compile params for each data type
   params.voltage = voltage;
   params.ap = ap;
   params.spike   = spike;
   params.rate    = rate;
   
