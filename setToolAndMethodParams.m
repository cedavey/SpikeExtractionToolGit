% method_params = extractToolAndMethodParams(method_params, data_type, tool, method)
%
% Return the configuration parameters for a particular tool + method of
% implementation
function params = setToolAndMethodParams(params, new_params, data_type, tool, method)
   tool   = lower(title2Str(tool,  1,1,'_')); 
   method = lower(title2Str(method,1,1,'_')); 
   params.(data_type).(tool).(method) = new_params;
   
end