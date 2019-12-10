% MOUSEWAITINGFUNCTION runs a function and sets the mouse pointer to
% 'waiting' mode while running.
% 
% Syntax:
%    mouseWaitingFunction(handles.figure1,@fun,arguments)
%
% Description:
% Sets the mouse pointer to the 'loading', or 'waiting' symbol while it
% runs the input function @fun
%
% Inputs:
%   fig1 - it is the handle to the active figure (normally handles.figure1)
%   fun  - it is the function it runs
%
%   all other inputs are passed as arguments to fun.
%
% * When running step by step the mouse wheel might not be removed. Press
%  the scroll wheel or scroll on the mouse pad/track pad to get rid of it.
% 
% Artemio - 19/June/2019
function mouseWaitingFunction(fig1,fun,varargin)
   isApp = false;
   try
      % Set the mouse pointer to waiting to know the function is running.
      set(fig1, 'pointer', 'watch');
      drawnow;
   catch
      % If it enters here it means we are on App instead of GUI. Mouse
      % pointer won't be shown.
      isApp = true;
   end
   
   try
      % Call the function
      fun(varargin{1:end});
   catch E
      % If something went wrong
      % Mouse pointer back to normal.
      if ~isApp, set(fig1, 'pointer', 'arrow'); end
      % Propagate the error
      runtimeErrorHandler(E);
   end
   
   % Mouse pointer back to normal.
   if ~isApp, set(fig1, 'pointer', 'arrow'); end
end