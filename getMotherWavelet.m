% [psi, xval] = getMotherWavelet(mother, params)
% Return the mother wavelet for chosen family & given parameters
% Orthog wavelets: input index number
% Continuous:      input time & freq spread
% Compression:     input vanishing moments for reconstruction & decomposition
% Outputs:
%  psi  - wavelet value
%  xval - 
function [filter, filtname] = getMotherWavelet(mother, params)
% Multires analysis
% - Orthog:
%   If preserving energy in the analysis stage is important, you must use an 
%   orthogonal wavelet. An orthogonal transform preserves energy. Consider 
%   using an orthogonal wavelet with compact support. Keep in mind that except 
%   for the Haar wavelet, orthogonal wavelets with compact support are not 
%   symmetric. The associated filters have nonlinear phase. Orthogonal
%   wavelets include: 
%           haar, db, sym, fk, coif
%
% - feature detection: identifying closely spaced features --> choose
%   wavelets with smaller support:
%           haar, db2, sym2
% 
% - for anal of discrete variance, use modwt
% - no overlap in wavelets within & across scales (min redundance) wavedec

% Time-freq anal: morse, amor, bump
% - use the continuous wavelet transform: cwt
% - defined by 2 params that control time and frequency spread
% - cwt has much more redundancy & is much more expensive to compute

% Compression:
% - bio-orthog spline 

% Denoising:
% - orthog wavelet is good. If using for images bio-orthog is good cuz have
%   linear phase

   % get short name of mother wavelet, & check wavelet number. If it's
   % not a valid index into the family, set number to empty & return empty
   % vectors
%    'Haar              		haar           '
%     'Daubechies        		db             '
%     'Symlets           		sym            '
%     'Coiflets          		coif           '
%     'BiorSplines       		bior           '
%     'ReverseBior       		rbio           '
%     'Meyer             		meyr           '
%     'DMeyer            		dmey           '
%     'Gaussian          		gaus           '
%     'Mexican_hat       		mexh           '
%     'Morlet            		morl           '
%     'Complex Gaussian  		cgau           '
%     'Shannon           		shan           '
%     'Frequency B-Spline		fbsp           '
%     'Complex Morlet    		cmor           '
%     'Fejer-Korovkin    		fk             '
   filter = []; filtname = [];
   switch lower(mother)
      case 'biorthogonal'
         filtname = 'bior';
         % Just use bior1.1, bior2.2, bior3.3, bior4.4 - what's the diff?!
         % 'bior1.1', 'bior1.3' , 'bior1.5'
         % 'bior2.2', 'bior2.4' , 'bior2.6', 'bior2.8'
         % 'bior3.1', 'bior3.3' , 'bior3.5', 'bior3.7'
         % 'bior3.9', 'bior4.4' , 'bior5.5', 'bior6.8'.
         if params<1 || params>6, params=[]; return; end
         if params<6
            filtname = sprintf('%s%.1f', filtname, params*1.1);
         else
            filtname = 'bior6.8';
         end
      case 'coiflets'
         filtname = 'coif';
         if params>5, params=[]; return; end
         filtname = sprintf('%s%d', filtname, params);
      case 'haar'
         filtname = 'haar';
         if params~=1, params=[]; return; end
      case 'daubechies'
         filtname = 'db';
         if params>45, params=[]; return; end
         filtname = sprintf('%s%d', filtname, params);
      case 'discrete meyer'
         filtname = 'dmey';
         if params~=1, params=[]; return; end
      case 'fejer-korovkin'
         filtname = 'fk';
         N = [4, 6, 8, 14, 18, 22];
         if ~any(params==N), params=[]; return; end
         filtname = sprintf('%s%d', filtname, params);
      case 'reverse biorthogonal'
         filtname = 'rbior';
         % Just use rbior1.1, rbior2.2, rbior3.3, rbior4.4 - what's the diff?!
         % 'rbior1.1', 'rbior1.3' , 'rbior1.5'
         % 'rbior2.2', 'rbior2.4' , 'rbior2.6', 'rbior2.8'
         % 'rbior3.1', 'rbior3.3' , 'rbior3.5', 'rbior3.7'
         % 'rbior3.9', 'rbior4.4' , 'rbior5.5', 'rbior6.8'.
         if params<1 || params>6, params=[]; return; end
         if params<6
            filtname = sprintf('%s%.1f', filtname, params*1.1);
         else
            filtname = 'rbior6.8';
         end
      case 'symlets'
         filtname = 'sym';
         if params<2 || params>45, params=[]; return; end
         filtname = sprintf('%s%d', filtname, params);
         
      % the following families require different parameterisation because
      % they're continuous (need time & freq spread params), or for
      % compression (need vanishing moments for reconstruction &
      % decomposition filters)
      % - I haven't implemented these yet
%       case 'biorsplines'
%          sname = 'bior';
%       case 'reversebior'
%          sname = 'rbio';
%       case 'meyer'
%          sname = 'meyr';
%       case 'dmeyer'
%       case 'gaussian'
%          sname = 'gaus';
%       case 'mexican_hat'
%          sname = 'mexh';
%       case 'morlet'
%          sname = 'morl';
%       case 'complex gaussian'
%          sname = 'cgau';
%       case 'shannon'
%          sname = 'shan';
%       case 'frequency b-spline'
%          sname = 'fbsp';
%       case 'complex morlet'
%          sname = 'cmor';
   end
   try
      filter  = wfilters(filtname);
   end
   % [psi,xval] = wavefun(sname, params);
end