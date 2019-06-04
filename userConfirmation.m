% response = userConfirmation(question, title)
%
% Asks user a question and gets their response as either 'yes' or 'no'. To
% get different response options use questdlg, or for many different
% questions and responses in the dialogue box use inputdlg. For no
% responses but just to tell the user something, use displayErrorMsg.
%
% See also:   questdlg   displayErrorMsg   inputsdlg    inputdlg
function response = userConfirmation(question,title)
    %src is the handle of the object generating the callback (the source of the event)
    %evnt is the The event data structure (can be empty for some callbacks)
    response = questdlg(question,title,'Yes','No','Yes');
return