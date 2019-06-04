% displayErrorMsg(msg, timeout)
%
% Display an error (or some other) message, not requiring user response.
% To ask the user a question and get their response as either 'yes' or 'no'
% user userConfirmation. To get different response options use questdlg, 
% or for many different questions and responses in the dialogue box use 
% inputdlg. For no responses but just to tell the user something, use 
% displayErrorMsg.
%
% See also:   questdlg   userConfirmation   inputsdlg    inputdlg
function displayErrorMsg(msg, timeout)
    if isempty(msg), return; end
    if nargin<2
        uiwait(msgbox(msg,'modal'));
    else
        uiwait(msgbox(msg,timeout));
    end
return