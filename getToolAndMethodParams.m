% method_params = extractToolAndMethodParams(method_params, data_type, tool, method)
%
% Return the configuration parameters for a particular tool + method of
% implementation
function method_params = getToolAndMethodParams(params, data_type, tool, method)
   tool   = lower(title2Str(tool,  1,1,'_')); 
   method = lower(title2Str(method,1,1,'_')); 
   try
      method_params = params.(data_type).(tool).(method);
   catch ME
      str = sprintf('Something funky''s going on with tool list for %s data type, setting to something safe', data_type);
      displayErrorMsg(str);
      runtimeErrorHandler(ME,'message',str);
      % set tool & method list to something vanilla that will always work
      fnames = fieldnames(params.(data_type));
      tool   = fnames{1};
      fnames = fieldnames(params.(data_type).(tool));
      method = fnames{1};
      method_params = params.(data_type).(tool).(method);
   end

end