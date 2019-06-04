% lhood = noiseLikelihood(pf, obs, Rsig)
%
% This gives the likelihood of each particle given the new observation. So, 
%   if the state is given by x, and the observation by y, such that 
%
%              x(t-1)  -->  x(t)  -->  x(t+1)
%                 |           |           |
%                 V           V           V
%              y(t-1)       y(t)       y(t+1)
%
%   then we want to calculate the posterior probability of each state (ie 
%   particle) given a prediction of the state and a new observation, y(t)
% - State transition function: transitions x(t) --> x(t+1)
% - Likelihood measurement: likelihood of state given new obs, which is
%        equiv to prob of new obs given state, P( y(t+1) | x(t+1) )
% - The particles represent different possibilities for the state, which
%   are transitioned and then the likelihood of each calculated, to max
%   the posterior probability of the state, or the likelihood of the obs. 
%
% - Our state is resistance, and our observations are voltage (spike) peaks,
%   and we transition the resistance linearly with time, 
%   so that our transition function is
%        transition:       R(t+1)  =  a * (t_sp - t_prevsp) + R(t)
%        lhood:       v_peak(t+1)  ~  N( R(t+1), Rsig)
%                  OR v_peak(t+1)  ~  N( R(t+1), vsig)
% - each particle is one possible value for the state, R(t+1).
%
% Note: if you have more than 1 state then pred_particles has format
%    num_particles x num_states if StateOrientation is row, & vice-versa if col
%
%    lhood = myMeasurementLikelihoodFcn(pred_particles, obs, varargin)
function lhood = noiseLikelihood( pf, particles, obs, Rsig )
   lhood = normpdf( obs, particles(:,1), Rsig );   
end