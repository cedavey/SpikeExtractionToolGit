% new_particles = evolveParticles(Rcoeff, Rsig, time) 
%
% This evolves current state particles at time t, to state particles at 
% time t+1. So, 
%   if the state is given by x, and the observation by y, such that 
%
%              x(t-1)  -->  x(t)  -->  x(t+1)
%                 |           |           |
%                 V           V           V
%              y(t-1)       y(t)       y(t+1)
%
%   then we want to calculate the prior probability of each state (ie 
%   particle) at the next timestep 
% - State transition function: transitions x(t) --> x(t+1)
% - Likelihood measurement: likelihood of state given new obs, which is
%        equiv to prob of new obs given state, P( y(t+1) | x(t+1) )
% - The particles represent different possibilities for the state, which
%   are transitioned and then the likelihood of each calculated, to max
%   the posterior probability of the state, or the likelihood of the obs. 
%
% - Our state is resistance, and our observations are voltage (spike) peaks,
%   and we transition the resistance linearly with time
%   so that our transition function is
%        transition:       R(t+1)  =  a * (t_sp - t_prevsp) + R(t)
%        lhood:       v_peak(t+1)  ~  N( R(t+1), Rsig)
% - each particle is one possible value for the state, R(t+1).
%
function new_particles = evolveResistance(pf, old_particles, Rcoeff, R, time, npoles, nzeros) 
   useVAR = true; % use full VAR model, or else simple y = mx + b wrt time

   if ~useVAR
      % evolve resistance according to linear regression wrt time: 
      %   R(t_sp) = a*(t_sp - t_prevsp) + R(t_sp)
      % OR if DC shift is NOT constant, then we'd require
      %   R(t_sp) = mu_new - mu_old + a*(t_sp - t_prevsp) + R(t_prevsp)

      % lambda = learning rate
   %    new_particles = Rcoeff(2) * (t_sp - t_prevsp) + old_particles;

%       new_particles = Rcoeff(end) * time(end) + old_particles;
      new_particles = Rcoeff(end) * (time(end)-time(end-1)) + old_particles;

    else
      % estimate current resistance based on linear regression model & current
      % states of lagged resistance

      % can make VAR matrix & calculate contribution from time input &
      % lagged resistance all in one go (uncomment line below to do this)
      % X = makeRecursiveLeastSquaresMatrix( time, R, min(length(R),length(time)), npoles, nzeros );

      %%% OR %%%
      
      % separate into regression matrix for input & autoregression matrix from
      % state variables (i.e. particles)

      % calculate contribution from time input
      T = [1; time(end:-1:end-nzeros+2)]'; % +2 coz nzeros also has dc component
      b = [Rcoeff(1); Rcoeff(npoles+2:end)];
      input = T * b;

      % calculate contribution from autoregression
      a = Rcoeff(2:npoles+1);
      new_particles = input + old_particles(:,2:end) * a(:);

   end

   % now shift old particles along so we're keeping current plus na past
   % timesteps, where 1st column is most recent timestep
   new_particles = [new_particles old_particles(:, 1:end-1)]; 
   
end