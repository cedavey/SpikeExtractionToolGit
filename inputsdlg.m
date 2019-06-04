function [Answer,Canceled] = inputsdlg(Prompt, Title, Formats, DefAns, Options)
%INPUTSDLG Enhanced input dialog box supporting multiple data types
% ANSWER = INPUTSDLG(PROMPT) creates a modal dialog box that returns user
% input for multiple prompts in the cell array ANSWER. PROMPT is a 1-D
% cell array containing the PROMPT strings.
%
% Alternatively, PROMPT can be a two-column cell array where the prompt
% string is supplied in the first column. Output ANSWER of the function is
% then in a structure with the field names defined in the second column of
% the PROMPT cell array.
%
% Moreover, PROMPT may have three columns, where the third column gives
% units (i.e., post-fix labels to the right of controls) to display.
%
% INPUTSDLG uses UIWAIT to suspend execution until the user responds.
%
% ANSWER = INPUTSDLG(PROMPT,NAME) specifies the title for the dialog.
%
% Note that INPUTSDLG(PROMPT) & INPUTSDLG(PROMPT,NAME) are similar to the
% standard INPUTDLG function, except for the dialog layout.
%
% ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS) can be used to specify the type
% of parameters to display with FORMATS matrix of structures. The
% dimension of FORMATS defines how PROMPT items are laid out in the dialog
% box. For example, if PROMPT has 6 elements and the size of FORMATS is
% 2x3 then, the items are shown in 2 rows, 3 columns format.
%
% The fields in FORMATS structure are:
%
%   type   - Type of control ['check',{'edit'},'list','range','text','none']
%   style  - UI control type used. One of:
%            [{'checkbox'},       for 'check' type
%             {'edit'}            for 'edit' type
%              'listbox',{'popupmenu'},'radiobutton','togglebutton'
%                                 for 'list' type
%             {'slider'}          for 'range' type
%             {'text'}]           for 'text' type
%   items  - Selection items for 'list' type (cell of strings)
%   format - Data format: ['string','float','integer','file','dir']
%   limits - [min max] (see below for details)
%   size   - [width height] in pixels. Alternatively, 0 to auto-size or -1
%            to auto-expand when figure is resized.
%
% FORMATS type field defines what type of prompt item to be shown.
%
%   type  Description
%   -------------------------------------------------------------------
%   edit  Standard edit box (single or multi-line mode)
%   check Check box for boolean item
%   list  Chose from a list of items ('listbox' style allows multiple item
%         selection)
%   range Use slider to chose a value over a range
%   text  Static text (e.g., for instructions)
%   none  A placeholder. May be used for its neighboring item to extend
%         over multiple columns or rows (i.e., "to merge cells")
%
% The allowed data format depends on the type of the field:
%
%   type    allowed format
%   --------------------------------------------
%   check   integer
%   edit    {text}, float, integer, file, dir
%   list    integer
%   range   float
%
% By leaving format field empty, a proper format is automatically chosen
% (or default to text format for edit type).
%
% Formats 'file' and 'dir' for 'edit' type uses the standard UIGETFILE,
% UIPUTFILE, and UIGETDIR functions to retrieve a file or directory name.
%
% The role of limits field varies depending on other parameters:
%
%   style         role of limits
%   ---------------------------------------------------
%   checkbox      limits(1) is the ANSWER value if the check box is not
%                 selected  box is not selected and limits(2) is the ANSWER
%                 if the check box is selected.
%   edit (text format)
%                 If diff(limits)>1, multi-line mode; else, single-line
%                 mode
%   edit (numeric format)
%                 This style defines the range of allowed values
%   edit (file format)
%                 If 0<=diff(limits)<=1 uses UIGETFILE in single select
%                 mode with single-line edit. If diff(limits)>1 uses
%                 UIGETFILE in multi-select mode with multi-line edit. If
%                 diff(limits)<0 usees UIPUTFILE with single- line edit
%   listbox       If diff(limits)>1, multiple items can be selected
%   slider        limits(1) defines the smallest value while
%                 limits(2) defines the largest value
%   none          If diff(limits)==0 space is left empty
%                 If diff(limits)>0 : lets the item from left to extend
%                 If diff(limits)<0 : lets the item from above to extend
%
% Similar to how PROMPT strings are laid out, when FORMATS.style is set to
% either 'radiobutton' or 'togglebutton', FORMATS.items are laid out
% according to the dimension of FORMATS.items.
%
% There are two quick format options as well:
%
%  Quick Format Option 1 (mimicing INPUTDLG behavior):
%   FORMATS can specify the number of lines for each edit-type prompt in
%   FORMATS. FORMATS may be a constant value or a column vector having
%   one element per PROMPT that specifies how many lines per input field.
%   FORMATS may also be a matrix where the first column specifies how
%   many rows for the input field and the second column specifies how
%   many columns wide the input field should be.
%
%  Quick Format Option 2:
%   FORMATS can specify the types of controls and use their default
%   configurations. This option, however, cannot be used to specify
%   'list' control as its items are not specified. To use this option,
%   provide a string (if only 1 control) or a cell array of strings. If
%   a cell array is given, its dimension is honored for the dialog
%   layout.
%
% ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS,DEFAULTANSWER) specifies the
% default answer to display for each PROMPT. DEFAULTANSWER must contain
% the same number of elements as PROMPT (that are not of 'none' style). If
% PROMPT does not provide ANSWER structure fields, DEFAULTANSWER should be
% a cell array with element type corresponding to FORMATS.format. Leave the
% cell element empty for a prompt with 'text' type. If ANSWER is a
% structure, DEFAULTANSWER must be a structure with the specified fields.
% (If additional fields are present in DEFAULTANSWER, they will be returned
% as parts of ANSWER.)
%
% ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS,DEFAULTANSWER,OPTIONS) specifies
% additional options. If OPTIONS is the string 'on', the dialog is made
% resizable. If OPTIONS is a structure, the fields recognized are:
%
%  Option Field Description {} indicates the default value
%  ----------------------------------------------------------------------
%  Resize        Make dialog resizable: 'on' | {'off'}
%  WindowStyle   Sets dialog window style: {'normal'} | 'modal'
%  Interpreter   Label text interpreter: 'latex' | {'tex'} | 'none'
%  ApplyButton   Adds Apply button: 'on' | {'off'}
%  Sep           Space b/w prompts in pixels: {10}
%
% [ANSWER,CANCELED] = INPUTSDLG(...) returns CANCELED = TRUE if user
% pressed Cancel button, closed the dialog, or pressed ESC. In such event,
% the content of ANSWER is set to the default values.
%
% Note on Apply Button feature. Pressing the Apply button makes the current
% change permanent. That is, pressing Cancel button after pressing Apply
% button only reverts ANSWER back to the states when the Apply button was
% pressed last. Also, if user pressed Apply button, CANCELED flag will not
% be set even if user canceled out of the dialog box.
%
% Examples:
%
% prompt={'Enter the matrix size for x^2:';'Enter the colormap name:'};
% name='Input for Peaks function';
% create a 1x2 struct array using:
% formats = struct('type',{'edit','edit'},'format',{'integer','text'},'limits',{[1 inf],[0,1]});
%   OR
% s1 = struct('type','edit','format','integer','limits',[0 inf]);
% s2 = struct('type','edit','format','text','limits',[0 1]);
% formats = [s1 s2];
% defaultanswer={20,'hsv'};
%
% [answer,canceled] = inputsdlg(prompt,name,formats,defaultanswer);
%
% formats(2).size = -1; % auto-expand width and auto-set height
% options.Resize='on';
% options.WindowStyle='normal';
% options.Interpreter='tex';
%
% answer = inputsdlg(prompt,name,formats,defaultanswer,options);
%
% prompt(:,2) = {'Ndim';'Cmap'};
% defaultanswer = struct(defaultanswer,prompt(:,2),1);
%
% answer = inputsdlg(prompt,name,formats,defaultanswer,options);
%
% See also INPUTDLG, DIALOG, ERRORDLG, HELPDLG, LISTDLG, MSGBOX,
%  QUESTDLG, TEXTWRAP, UIWAIT, WARNDLG, UIGETFILE, UIPUTFILE, UIGETDIR.

% Version 1.11 (Nov. 19, 2009)
% Written by: Takeshi Ikuma
% Contributer: Andreas Greuer
% Created: Nov. 16, 2009
% Revision History:
%  v.1.1 (Nov. 19, 2009)
%  * Fixed bugs (reported by AG):
%   - not returning Canceled output
%   - erroneous struct output behavior
%   - error if all row elements of a column are auto-expandable
%  * Added Apply button option
%  * Added support for Units (label to the right of controls)
%  * Updated the help text
%  v.1.11 (Nov. 20, 2009)
%  * Fixed bugs (reported by AG):
%   - incorrect Canceled output when Cancel button is pressed
%  v.1.12 (Nov. 20, 2009)
%  * Fixed bugs (reported by AG):
%   - again incorrect Canceled output behavior

%%%%%%%%%%%%%%%%%%%%
%%% Nargin Check %%%
%%%%%%%%%%%%%%%%%%%%
error(nargchk(0,5,nargin));
error(nargoutchk(0,2,nargout));

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle Input Args %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1, Prompt=''; end
if nargin<2, Title = ''; end
if nargin<3, Formats=struct([]); end
if nargin<4, DefAns = {}; end
if nargin<5, Options = struct([]); end

% Check Prompt input
[Prompt,FieldNames,Units,err] = checkprompt(Prompt);
if ~isempty(err), error(err{:}); end
NumQuest = numel(Prompt); % number of prompts

if isempty(Title)
   Title = ' ';
elseif iscellstr(Title)
   Title = Title{1}; % take the first entry
elseif ~ischar(Title)
   error('inputsdlg:InvalidInput','Title must be a string of cell string.');
end

% make sure that the Formats structure is valid & fill it in default
% values as needed
[Formats,err] = checkformats(Formats,NumQuest);
if ~isempty(err), error(err{:}); end

% make sure that the DefAns is valid & set Answer using DefAns and default
% values if DefAns not given
[Answer,AnsStr,err] = checkdefaults(DefAns,Formats,FieldNames);
if ~isempty(err), error(err{:}); end

% make sure that the Options is valid
[Options,err] = checkoptions(Options);
if ~isempty(err), error(err{:}); end

Applied = false; % set true by pressing Apply Button

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Dialog GUI %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% lay contents out on a dialog box
[handles,Formats,sinfo] = buildgui(Prompt,Units,Title,Formats,Options);

% set default values
for k = 1:NumQuest
   switch Formats(k).style
      case {'checkbox' 'listbox' 'popupmenu' 'slider'}
         set(handles.ctrl(k),'Value',Answer{k});
      case 'edit'
         if any(strcmp(Formats(k).format,{'integer','float'}))
            set(handles.ctrl(k),'String',num2str(Answer{k}));
         else
            set(handles.ctrl(k),'String',Answer{k});
         end
      case {'radiobutton' 'togglebutton'}
         h = get(handles.ctrl(k),'UserData');
         set(handles.ctrl(k),'SelectedObject',h(Answer{k}));
   end
end

% set callback functions
for k = 1:NumQuest % for all non-'text' controls
   if strcmp(Formats(k).style,'edit')
      switch Formats(k).format
         case {'float','integer'}
            % for numeric edit box, check for the range & set mouse down behavior
            set(handles.ctrl(k),'Callback',@(hObj,evd)checkRange(hObj,evd,k,Formats(k).limits));
         case 'file'
            mode = diff(Formats(k).limits);
            set(handles.ctrl(k),...
               'ButtonDownFcn',@(hObj,evd)openFilePrompt(hObj,evd,Formats(k).items,Prompt{k},mode),...
               'Enable','inactive');
         case 'dir'
            set(handles.ctrl(k),...
               'ButtonDownFcn',@(hObj,evd)openDirPrompt(hObj,evd,Prompt{k}),...
               'Enable','inactive');
      end
   end
end

set(handles.fig,'UserData','Cancel');
set(handles.btns(1), 'KeyPressFcn', @(hObj,evd)doControlKeyPress(hObj,evd,true), 'Callback' ,@(hObj,evd)doCallback(hObj,evd,true));
set(handles.btns(2), 'KeyPressFcn', @(hObj,evd)doControlKeyPress(hObj,evd,false), 'Callback' ,@(hObj,evd)doCallback(hObj,evd,false));
if numel(handles.btns)>2
   set(handles.btns(3), 'KeyPressFcn', @(hObj,evd)doControlKeyPress(hObj,evd), 'Callback' ,@doApply);
end
set(handles.fig, 'KeyPressFcn', @doFigureKeyPress,'ResizeFcn', @doResize);

% make sure we are on screen
movegui(handles.fig)

% if there is a figure out there and it's modal, we need to be modal too
if ~isempty(gcbf) && strcmp(get(gcbf,'WindowStyle'),'modal')
   set(handles.Figure,'WindowStyle','modal');
end

% ready to begin the show!
set(handles.fig,'Visible','on');
drawnow;

% set focus on the first uicontol
h = handles.ctrl(find(~strcmp('text',{Formats.type}),1,'first'));
if ~isempty(h)
   switch get(h,'type')
      case 'uicontrol', uicontrol(h);
      case 'uitoggletool', uicontrol(get(h,'SelectedObject'));
   end
end

% Go into uiwait if the figure handle is still valid.
% This is mostly the case during regular use.
if ishghandle(handles.fig), uiwait(handles.fig); end

% Check handle validity again since we may be out of uiwait because the
% figure was deleted.
Canceled = ~(Applied || (ishghandle(handles.fig) && strcmp(get(handles.fig,'UserData'),'OK')));
if Canceled % return the default answer
   if isempty(FieldNames), Answer = DefAns;
   else Answer = AnsStr; end % AnsStr contains the default value until getAnswer is called
else
   Answer = getAnswer(Answer,AnsStr); % get the final answers
end

% Close the figure if it's still open
if ishghandle(handles.fig), delete(handles.fig); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function answer = getAnswer(answer,ansstr)
      % retrieve answer from controls
      for i = 1:numel(answer)
         switch Formats(i).style
            case {'checkbox' 'listbox' 'popupmenu' 'slider'}
               answer{i} = get(handles.ctrl(i),'Value');
            case 'edit'
               switch Formats(i).format
                  case {'text' 'dir'}
                     answer{i} = get(handles.ctrl(i),'String');
                  case 'file'
                     data = get(handles.ctrl(i),'UserData');
                     if isempty(data)
                        answer{i} = get(handles.ctrl(i),'String');
                     else % multi-select
                        answer{i} = strcat(data{1},data{2});
                     end
                  otherwise
                     answer{i} = str2double(get(handles.ctrl(i),'String'));
               end
            case {'radiobutton' 'togglebutton'}
               answer{i} = find(get(handles.ctrl(i),'SelectedObject')==get(handles.ctrl(i),'UserData'));
         end
      end
      
      % if struct answer expected, copy
      if ~isempty(FieldNames)
         idx = find(~cellfun('isempty',FieldNames))';
         for i = idx
            ansstr.(FieldNames{i}) = answer{i};
         end
         answer = ansstr;
      end
   end

   function doApply(hObj,evd) %#ok
      Applied = true; % set the flag
      DefAns = getAnswer(Answer,AnsStr);
      if isstruct(DefAns)
         DefStr = DefAns;
         DefAns = cell(size(FieldNames));
         idx = find(~cellfun('isempty',FieldNames))';
         for i = idx
            DefAns{i} = DefStr.(FieldNames{i});
         end
      end
   end

   function checkRange(hObj,evd,k,lim) %#ok
      val = str2double(get(hObj,'String'));
      if ~isnan(val) && val>=lim(1) && val<=lim(2)
         Answer{k} = val;
      else
         if strcmp(Formats(k).format,'integer')
            msg = sprintf('%d, %d',lim(1),lim(2));
         else
            msg = sprintf('%g, %g',lim(1),lim(2));
         end
         h = errordlg(sprintf('This parameter must be within the range [%s].',msg),'Invalid Value','modal');
         uiwait(h);
         set(hObj,'String',num2str(Answer{k}));
      end
   end

   function doFigureKeyPress(obj, evd)
      switch(evd.Key)
         case {'return','space'}
            set(obj,'UserData','OK');
            uiresume(obj);
         case {'escape'}
            delete(obj);
      end
   end

   function doControlKeyPress(obj, evd, varargin)
      switch(evd.Key)
         case {'return'} % execute its callback function with varargin
            cbfcn = get(obj,'Callback');
            cbfcn(obj,evd);
         case 'escape'
            delete(gcbf)
      end
   end

   function doCallback(obj, evd, isok) %#ok
      if isok || Applied
         set(gcbf,'UserData','OK');
         uiresume(gcbf);
      else
         delete(gcbf)
      end
   end

   function openFilePrompt(hObj,evd,spec,prompt,mode) %#ok
      
      if mode<0 % uiputfile
         [f,p] = uiputfile(spec,prompt,get(hObj,'String'));
         if f~=0, set(hObj,'String',[p f]); end
      elseif mode<=1 % uigetfile
         [f,p] = uigetfile(spec,prompt,get(hObj,'String'));
         if f~=0, set(hObj,'String',[p f]); end
      else % uigetfile multi-select
         % previously chosen files are lost, but directly is kept
         data = get(hObj,'UserData');
         [f,p] = uigetfile(spec,prompt,data{1},'MultiSelect','on');
         if p~=0
            data = {p f};
            if ischar(f) % single file selected
               set(hObj,'String',[p f],'UserData',data);
            else % multiple files selected
               str = sprintf('"%s%s"',p,f{1});
               for n = 2:length(f)
                  str = sprintf('%s\n"%s%s"',str,p,f{n});
               end
               set(hObj,'String',str,'UserData',data);
            end
         end
      end
      uicontrol(hObj);
   end

   function openDirPrompt(hObj,evd,prompt) %#ok
      p = uigetdir(get(hObj,'String'),prompt);
      if p~=0, set(hObj,'String',p); end
      uicontrol(hObj);
   end

   function doResize(hObj,evd) %#ok
      % This function places all controls in proper place.
      % Must be called before the GUI is made visible as buildgui function
      % just creates uicontrols and do not place them in proper places.
      
      % get current figure size
      figPos = get(hObj,'Position');
      figSize = figPos(3:4);
      
      % determine the column width & row heights
      workarea = [Options.Sep Options.Sep figSize - 2*Options.Sep]; % Options.Sep margins around the figure
      btnarea = [0,0,workarea(3),sinfo.h_btns];
      ctrlarea = [0,btnarea(4)+Options.Sep,workarea(3),workarea(4)-sinfo.h_btns];
      
      dim = size(sinfo.map);
      num = numel(handles.ctrl);
      
      % determine the column widths & margin
      w = sinfo.w_ctrls+sinfo.w_labels+sinfo.w_units; % minimum widths of elements
      width = zeros(size(sinfo.map));
      cext = false(1,dim(2));
      for n = 1:dim(2)
         cmap = sinfo.map(:,n);
         cext(n) = any(sinfo.autoextend(cmap(cmap~=0),1));
      end
      if any(cext) % found auto-extendable element(s)
         m_col = Options.Sep; % column margin fixed
         w_total = ctrlarea(3)-(dim(2)+1)*m_col; % sum of control width
         
         % record the widths of non-expandable elements
         for n = find(~sinfo.autoextend(:,1))'
            idx = find(sinfo.map==n);
            [i,j] = ind2sub(dim,idx);
            J = unique(j);
            n_col = numel(J);
            width(idx) = (w(n)-(n_col-1)*m_col)/n_col;
         end
         w_col = max(width,[],1); % column width based on non-expandables
         
         % figure out how to distribute extra spaces among auto-expandable columns
         w_avail = w_total - sum(w_col(~cext)); % available width
         idx = w_col(cext)==0; % column where all elements are auto-expandable
         if all(idx) % all columns auto-expandable
            w_col(cext) = w_avail/sum(cext); % equally distributed
         else
            w_xcol = w_col(cext);
            if any(idx)
               w_xcol(idx) = mean(w_xcol(~idx));
            end % some columns are auto-expandable
            w_xcol = w_xcol + (w_avail-sum(w_xcol))/sum(cext); % equally distribute the excess
            w_col(cext) = w_xcol;
         end
         
         w_col = max(w_col,1); % make sure it's positive
         
         % set the expandable controls' width
         for n = find(sinfo.autoextend(:,1))'
            idx = find(sinfo.map==n);
            [i,j] = ind2sub(dim,idx);
            J = unique(j);
            if strcmp(Formats(n).type,'text')
               sinfo.w_labels(n) = sum(w_col(J));
            else
               sinfo.w_ctrls(n) = max(sum(w_col(J)) - sinfo.w_labels(n),1);
            end
         end
      else % no auto-extension, adjust column margin
         for n = 1:num
            idx = find(sinfo.map==n);
            [i,j] = ind2sub(dim,idx);
            J = unique(j);
            n_col = numel(J);
            width(idx) = w(n)/n_col;
         end
         w_col = max(width,[],1);
         m_col = (ctrlarea(3)-sum(w_col))/(dim(2)+1);
      end
      
      % resize static text with auto-height
      aex = sinfo.autoextend;
      for m = find(strcmp('text',{Formats.type}) & any(sinfo.autoextend(:,2)))
         idx = find(sinfo.map==m);
         [i,j] = ind2sub(dim,idx);
         J = unique(j);
         w = sum(w_col(J)) + (numel(J)-1)*m_col;
         if w<=100, w = 100; end
         
         % create dummy uicontrol to wrap text
         h = uicontrol('Parent',hObj,'Style','text','Position',[0 0 w 20],'Visible','off');
         msg = textwrap(h,Prompt(m));
         delete(h);
         str = msg{1};
         for n = 2:length(msg)
            str = sprintf('%s\n%s',str,msg{n});
         end
         
         % update the text wrapping
         set(handles.labels(m),'String',str);
         
         % get updated size
         pos = get(handles.labels(m),'Extent');
         sinfo.w_labels(m) = pos(3);
         sinfo.h_labels(m) = pos(4);
         sinfo.h_ctrls(m) = pos(4);
         
         aex(m,2) = false;
      end
      
      % determine the row heights & margin
      h = max([sinfo.h_ctrls,sinfo.h_labels,sinfo.h_units],[],2);
      height = zeros(size(sinfo.map));
      rext = false(dim(1),1);
      for n = 1:dim(1)
         cmap = sinfo.map(n,:);
         rext(n) = any(aex(cmap(cmap~=0),2));
      end
      if any(rext)
         m_row = Options.Sep; % column margin fixed
         h_total = ctrlarea(4)-(dim(1)+1)*m_row; % sum of control width
         
         % record the heightss of non-expandable elements
         for n = find(~aex(:,2))'
            idx = find(sinfo.map==n);
            [i,j] = ind2sub(dim,idx); %#ok
            I = unique(i);
            n_row = numel(I);
            height(idx) = (h(n)-(n_row-1)*m_row)/n_row;
         end
         h_row = max(height,[],2); % column width based on non-expandables
         
         % figure out how to distribute extra spaces among auto-expandable rows
         h_avail = h_total - sum(h_row(~rext)); % available width
         idx = h_row(rext)==0; % column where all elements are auto-expandable
         if all(idx) % all columns auto-expandable
            h_row(rext) = h_avail/sum(rext); % equally distributed
         else
            h_xrow = h_row(rext);
            if any(idx), h_xrow(idx) = mean(h_xrow(~idx)); end % some columns are auto-expandable
            h_xrow = h_xrow + (h_avail-sum(h_xrow))/sum(rext); % equally distribute the excess
            h_row(rext) = h_xrow;
         end
         
         h_row = max(h_row,1);
         
         % set the expandable controls width
         for n = find(aex(:,2))'
            idx = find(sinfo.map==n);
            [i,j] = ind2sub(dim,idx); %#ok
            I = unique(i);
            n_row = numel(I);
            newh = sum(h_row(I)) + (n_row-1)*m_row;
            change = newh - sinfo.h_ctrls(n);
            sinfo.h_ctrls(n) = newh;
            sinfo.CLoffset(n) = sinfo.CLoffset(n) + change;
            sinfo.CUoffset(n) = sinfo.CUoffset(n) + change;
         end
      else % no auto-extension in heights, adjust row margin
         for n = 1:num
            idx = find(sinfo.map==n);
            [i,j] = ind2sub(dim,idx); %#ok
            I = unique(i);
            n_row = numel(I);
            height(idx) = h(n)/n_row;
         end
         h_row = max(height,[],2);
         m_row = (ctrlarea(4)-sum(h_row))/(dim(1)+1);
      end
      
      % set control positions
      for m = 1:num
         [i,j] = ind2sub(dim,find(sinfo.map==m));
         x0 = sum(w_col(1:j(1)-1)) + m_col*j(1) + ctrlarea(1) + workarea(1);
         y0 = sum(h_row(i(1):end)) - sinfo.h_ctrls(m) + m_row*(dim(1)-i(1)+1) + ctrlarea(2) + workarea(2);
         
         posL = [x0,y0+sinfo.CLoffset(m)];
         posC = [x0+sinfo.w_labels(m),y0,sinfo.w_ctrls(m),sinfo.h_ctrls(m)];
         posU = [x0+sinfo.w_labels(m)+sinfo.w_ctrls(m)+2,y0+sinfo.CUoffset(m)];
         
         if handles.labels(m)~=0, set(handles.labels(m),'Position',posL); end
         if handles.ctrl(m)~=0, set(handles.ctrl(m),'Position',posC); end
         if handles.units(m)~=0, set(handles.units(m),'Position',posU); end
      end
      
      % set positions of buttons
      nbtns = numel(handles.btns);
      w_allbtns = 2*Options.Sep*(nbtns-1) + nbtns*sinfo.w_btns;
      pos = [workarea(1) + (btnarea(3) - w_allbtns)/2, workarea(2), sinfo.w_btns, sinfo.h_btns];
      for n = 1:nbtns
         set(handles.btns(n),'Position',pos);
         pos(1) = pos(1) + sinfo.w_btns + 2*Options.Sep;
      end
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILDGUI :: Builds the dialog box and returns handles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles,Formats,sinfo] = buildgui(Prompt,Unit,Title,Formats,Options)

DefaultSize.text = [0 0]; % auto size (max width: 250)
DefaultSize.edit = [165 0]; % auto size height (base:20)
DefaultSize.popupmenu = [100 20]; % auto size width
DefaultSize.listbox = [0 0];   % auto size
DefaultSize.togglebutton = [0 0]; % auto size
DefaultSize.radiobutton = [0 0]; % auto size
DefaultSize.checkbox = [0 18]; % auto size width
DefaultSize.slider = [75 15];
DefaultSize.pushbutton = [69 22];

% determine how to utilize 'none' space
% Place all the elements at (0,0)
free = reshape(strcmp('none',{Formats.type}),size(Formats)); % location of empty block to be occupied by neighbor entry
num = sum(~free(:)); % number of controls
dim = size(Formats); % display grid dimension
map = zeros(dim); % determine which control occupies which block(s)
order = zeros(1,num); % uicontrol placement order
n = 1;
for f = 1:prod(dim)
   % traverse row-first
   [j,i] = ind2sub(dim([2 1]),f);
   m = sub2ind(dim,i,j);
   
   if free(m)
      mode = diff(Formats(m).limits);
      [i,j] = ind2sub(dim,m);
      if mode>0 && j>1, map(m) = map(sub2ind(dim,i,j-1)); % left
      elseif mode<0 && i>1, map(m) = map(sub2ind(dim,i-1,j)); % above
      end % other wise, 0 (nothing occupying)
   else
      map(m) = n;
      order(n) = m;
      n = n + 1;
   end
end

% remove none's from Formats and order the rest in Prompt order
Formats = Formats(order);

% assign default size if Formats.size is non-positive
autosize = false(num,2);
autoextend = false(num,2);
for m = 1:num
   autoextend(m,:) = Formats(m).size<0;
   
   % get default size if size not specified
   if Formats(m).size(1) <=0
      Formats(m).size(1) = DefaultSize.(Formats(m).style)(1);
   end
   if Formats(m).size(2) <=0
      Formats(m).size(2) = DefaultSize.(Formats(m).style)(2);
   end
   
   % if fixed size, turn autosize off
   autosize(m,:) = Formats(m).size<=0;
end


FigColor=get(0,'DefaultUicontrolBackgroundcolor');

fig = dialog(           ...
   'Visible'     ,'off'   , ...
   'Name'        ,Title   , ...
   'Pointer'     ,'arrow'  , ...
   'Units'       ,'pixels'  , ...
   'UserData'    ,'Cancel'  , ...
   'Tag'         ,'Figure'   , ...
   'HandleVisibility' ,'on' , ...
   'Color'       ,FigColor  , ...
   'NextPlot'    ,'add'   , ...
   'WindowStyle' ,Options.WindowStyle, ...
   'DoubleBuffer','on'    , ...
   'Resize'      ,Options.Resize    ...
   );

figSize = get(fig,'Position');

%%%%%%%%%%%%%%%%%%%%%
%%% Set Positions %%%
%%%%%%%%%%%%%%%%%%%%%
CommonInfo = {'Units'  'pixels';
   'FontSize'      get(0,'FactoryUIControlFontSize');
   'HandleVisibility'  'callback'}';

props.edit = [CommonInfo {...
   'FontWeight'   get(fig,'DefaultUicontrolFontWeight');
   'Style'      'edit';
   'HorizontalAlignment' 'left';
   'BackgroundColor' 'white'}'];

props.checkbox = [CommonInfo {...
   'Style'      'checkbox';
   'FontWeight'   get(fig,'DefaultUicontrolFontWeight');
   'HorizontalAlignment' 'left';
   'BackgroundColor' FigColor}'];

props.popupmenu = [CommonInfo {...
   'FontWeight'   get(fig,'DefaultUicontrolFontWeight');
   'Style'      'popupmenu';
   'HorizontalAlignment' 'left';
   'BackgroundColor' 'white'}'];

props.listbox = [CommonInfo {...
   'FontWeight'   get(fig,'DefaultUicontrolFontWeight');
   'Style'      'listbox';
   'HorizontalAlignment' 'left';
   'BackgroundColor' 'white'}'];

props.slider = [CommonInfo {...
   'Style'      'slider';
   }'];

props.uibuttongroup = [CommonInfo {...
   'FontWeight'   get(fig,'DefaultUicontrolFontWeight');
   'BackgroundColor' FigColor;
   }'];

props.radiobutton = props.checkbox;
props.radiobutton{2,4} = 'radiobutton';

props.pushbutton = [CommonInfo {...
   'Style'        'pushbutton';
   'FontWeight'     get(fig,'DefaultUicontrolFontWeight');
   'HorizontalAlignment' 'center'}'];

props.togglebutton = props.pushbutton;
props.togglebutton{2,4} = 'togglebutton';

% Add VerticalAlignment here as it is not applicable to the above.
TextInfo = [CommonInfo {...
   'BackgroundColor'      FigColor;
   'HorizontalAlignment' 'left';
   'VerticalAlignment'   'bottom';
   'Color'                get(0,'FactoryUIControlForegroundColor');
   'Interpreter'          Options.Interpreter}'];

% Place the container (UIPANEL & AXES) for the elements (for the ease of
% resizing)
ax = axes('Parent',fig,'Units','pixels','Visible','off',...
   'Position',[0 0 figSize(3:4)],'XLim',[0 figSize(3)],'YLim',[0 figSize(4)]);

hPrompt = zeros(num,1);
hEdit = zeros(num,1);
hUnit = zeros(num,1);
Cwidth = zeros(num,1);
Cheight = zeros(num,1);
Lwidth = zeros(num,1);
Lheight = zeros(num,1);
Uwidth = zeros(num,1);
Uheight = zeros(num,1);
CLoffset = zeros(num,1);
CUoffset = zeros(num,1);
for m = 1:num
   Formats(m).size(autosize(m,:)) = 1; % temporary width
   
   idx = strcmp(Formats(m).style,{'radiobutton','togglebutton'});
   if any(idx)
      % create the UI Button Group object
      hEdit(m) = uibuttongroup('Parent',fig,props.uibuttongroup{:},...
         'Position',[0 0 Formats(m).size],'Title',Prompt{m});
      
      if idx(1),margin = [15 0]; % radiobutton
      else margin = [10 2]; end % togglebutton
      
      num_btns = numel(Formats(m).items);
      dim_btns = size(Formats(m).items);
      hButtons = zeros(dim_btns);
      btn_w = zeros(dim_btns);
      btn_h = zeros(dim_btns);
      for k = 1:num_btns
         [j,i] = ind2sub(dim_btns([2 1]),k);
         if isempty(Formats(m).items{i,j}), continue; end % if empty string, no button at this position
         hButtons(i,j) = uicontrol('Parent',hEdit(m),'Style',Formats(m).style,props.(Formats(m).style){:},...
            'String',Formats(m).items{i,j},'Min',0,'Max',1,'UserData',k);
         pos = get(hButtons(i,j),'Extent');
         btn_w(i,j) = pos(3);
         btn_h(i,j) = pos(4);
      end
      
      set(hEdit(m),'UserData',hButtons);
      
      btn_w = btn_w + margin(1);
      btn_h = btn_h + margin(2);
      
      col_w = max(btn_w,[],1);
      row_h = max(btn_h,[],2);
      
      % set positions of buttons
      kvalid = find(hButtons~=0);
      for k = reshape(kvalid,1,numel(kvalid))
         [i,j] = ind2sub(dim_btns,k); % i-col, j-row
         pos = [sum(col_w(1:j-1))+Options.Sep*j sum(row_h(i+1:end)) + Options.Sep*(dim_btns(1)-i+0.5) btn_w(k) btn_h(k)];
         set(hButtons(i,j),'Position',pos);
      end
      
      Cwidth(m) = sum(col_w)+Options.Sep*(dim_btns(2)+1);
      Cheight(m) = sum(row_h)+Options.Sep*(dim_btns(1)+1);
      
      % no auto-extension
      autoextend(m,:) = false;
   elseif strcmp(Formats(m).type,'text') % static text (no labels)
      % create label
      if autosize(m,1) % if not autosize, wrap as needed
         str = Prompt{m};
      else
         h = uicontrol('Parent',fig,'Style','text','Position',[0 0 Formats(m).size]);
         msg = textwrap(h,Prompt(m));
         delete(h);
         str = msg{1};
         for n = 2:length(msg)
            str = sprintf('%s\n%s',str,msg{n});
         end
      end
      
      hPrompt(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',str);
      pos = get(hPrompt(m),'Extent');
      Lwidth(m) = pos(3);
      Lheight(m) = pos(4);
      Cheight(m) = pos(4);
      CLoffset(m) = 1;
   else
      hEdit(m) = uicontrol('Parent',fig,'Style',Formats(m).style,props.(Formats(m).style){:},...
         'Position',[0 0 Formats(m).size]);
      
      Cwidth(m) = Formats(m).size(1);
      Cheight(m) = Formats(m).size(2);
      
      % set min and max except for edit-text box
      if ~strcmp(Formats(m).style,'edit') || ~any(strcmp(Formats(m).format,{'float','integer'}))
         set(hEdit(m),'Min',Formats(m).limits(1),'Max',Formats(m).limits(2));
      end
      
      switch lower(Formats(m).style)
         case 'edit'
            textmode = any(strcmp(Formats(m).format,{'text','file'}));
            if autosize(m,2) % auto-height
               if textmode
                  dlim = diff(Formats(m).limits);
                  if (strcmp(Formats(m).format,'file') && dlim<1), dlim = 1; end
                  Cheight(m) = 15*dlim + 5;
               else % numeric -> force single-line
                  Cheight(m) = 20;
               end
               set(hEdit(m),'Position',[0 0 Formats(m).size]);
            end
            
            % If format is not text, reset Min & Max to single-line mode
            if ~textmode, set(hEdit(m),'Min',0,'Max',1); end
            
            % create label
            hPrompt(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Prompt{m});
            pos = get(hPrompt(m),'Extent');
            Lwidth(m) = pos(3); % label element width
            Lheight(m) = pos(4); % label element height
            
            if ~isempty(Unit{m})
               hUnit(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Unit{m});
               pos = get(hUnit(m),'Extent');
               Uwidth(m) = pos(3)+5; % label element width
               Uheight(m) = pos(4); % label element height
            end
            
            % Align label with the editbox text location
            % editbox first-line text location varies depending on single- or multi-line mode
            multline = textmode && diff(Formats(m).limits)>1;
            if multline
               CLoffset(m) = 3 + Cheight(m) - Lheight(m);
               CUoffset(m) = 3 + Cheight(m) - Uheight(m);
            else
               CLoffset(m) = 3 + (Cheight(m) - Lheight(m))/2;
               CUoffset(m) = 3 + (Cheight(m) - Uheight(m))/2;
            end
            
            % height auto-extension posible only if multiline
            autoextend(m,2) = autoextend(m,2) && multline;
            
            if strcmp(Formats(m).format,'file') && diff(Formats(m).limits)>1
               set(hEdit(m),'UserData',{'' ''});
            end
         case 'checkbox' % no labels
            if autosize(m,1) % auto-width
               set(hEdit(m),'String',Prompt{m});
               pos = get(hEdit(m),'Extent');
               pos(3) = pos(3) + 15;
               Cwidth(m) = pos(3);%Formats(m).size(1) = pos(3);
               set(hEdit(m),'Position',pos);
            end
            
            % width auto-extension possible only if width too short
            autoextend(m,2) = false; % no height auto-extend
         case 'popupmenu'
            % force the height
            Cheight(m) = DefaultSize.popupmenu(2);
            
            % get the width of the widest entry
            if autosize(m,1) % auto-width
               for itm = 1:numel(Formats(m).items)
                  set(hEdit(m),'String',Formats(m).items{itm},'Value',1);
                  p = get(hEdit(m),'Extent');
                  if p(3)>Cwidth(m), Cwidth(m) = p(3); end
               end
               Cwidth(m) = Cwidth(m) + 25;
            end
            
            % re-set position
            set(hEdit(m),'Position',[0 0 Cwidth(m) Cheight(m)]);
            
            % set menu & choose the first entry
            set(hEdit(m),'String',Formats(m).items','Value',1);
            
            % create label
            hPrompt(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Prompt{m});
            pos = get(hPrompt(m),'Extent');
            Lwidth(m) = pos(3); % label element width
            Lheight(m) = pos(4); % label element height
            
            if ~isempty(Unit{m})
               hUnit(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Unit{m});
               pos = get(hUnit(m),'Extent');
               Uwidth(m) = pos(3)+5; % label element width
               Uheight(m) = pos(4); % label element height
            end
            
            % Align label with the editbox text location
            % editbox first-line text location varies depending on single- or multi-line mode
            CLoffset(m) = 1 + Cheight(m) - Lheight(m);
            CUoffset(m) = 1 + Cheight(m) - Uheight(m);
            
            % width auto-extension possible only if width too short
            autoextend(m,2) = false; % no height auto-extend
         case 'listbox'
            % get the max extent
            pos = [0 0 0 0];
            for itm = 1:numel(Formats(m).items)
               set(hEdit(m),'String',Formats(m).items{itm},'Value',1);
               p = get(hEdit(m),'Extent');
               if p(3)>pos(3), pos(3) = p(3); end
            end
            pos(3) = pos(3) + 19;
            pos(4) = 13*numel(Formats(m).items)+2;
            
            % auto-size & set autoextend
            if autosize(m,1), Cwidth(m) = pos(3); end
            if autosize(m,2), Cheight(m) = pos(4); end
            if any(autosize(m,:)), set(hEdit(m),'Position',pos); end
            
            % set table & leave it unselected
            if diff(Formats(m).limits)>1, val = [];
            else val = 1; end
            set(hEdit(m),'String',Formats(m).items','Value',val);
            
            % create label
            hPrompt(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Prompt{m});
            pos = get(hPrompt(m),'Extent');
            Lwidth(m) = pos(3); % label element width
            Lheight(m) = pos(4); % label element height
            
            if ~isempty(Unit{m})
               hUnit(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Unit{m});
               pos = get(hUnit(m),'Extent');
               Uwidth(m) = pos(3)+5; % label element width
               Uheight(m) = pos(4); % label element height
            end
            
            % Align label with the editbox text location
            % editbox first-line text location varies depending on single- or multi-line mode
            CLoffset(m) = Cheight(m)-12;
            CUoffset(m) = Cheight(m)-12;
            
         case 'slider'
            % set slider initial value
            val = mean(Formats(m).limits);
            set(hEdit(m),'Value',val,'TooltipString',num2str(val),'Callback',@(hObj,evt)set(hObj,'TooltipString',num2str(get(hObj,'Value'))));
            
            % create label
            hPrompt(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Prompt{m});
            pos = get(hPrompt(m),'Extent');
            Lwidth(m) = pos(3); % label element width
            Lheight(m) = pos(4); % label element height
            
            if ~isempty(Unit{m})
               hUnit(m) = text('Parent',ax,TextInfo{:},'Position',[0,0],'String',Unit{m});
               pos = get(hUnit(m),'Extent');
               Uwidth(m) = pos(3)+5; % label element width
               Uheight(m) = pos(4); % label element height
            end
            % Align label with the editbox text location
            % editbox first-line text location varies depending on single- or multi-line mode
            CLoffset(m) = (Formats(m).size(2) - Lheight(m))/2 + 4;
            CUoffset(m) = (Formats(m).size(2) - Uheight(m))/2 + 4;
            
            % only resizable in the direction of slider
            if diff(Formats(m).size)>=0, autoextend(m,2) = false;
            else autoextend(m,1) = false; end
      end
   end
end

hBtn(1) = uicontrol(fig,props.pushbutton{:},'String','OK',...
   'Position',[0 0 DefaultSize.pushbutton]);
hBtn(2) = uicontrol(fig,props.pushbutton{:},'String','Cancel',...
   'Position',[0 0 DefaultSize.pushbutton]);
if strcmpi(Options.ApplyButton,'on')
   hBtn(3) = uicontrol(fig,props.pushbutton{:},'String','Apply',...
      'Position',[0 0 DefaultSize.pushbutton]);
end

% output handle struct
handles.fig = fig;
handles.labels = hPrompt';
handles.ctrl = hEdit';
handles.units = hUnit';
handles.btns = hBtn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine the minimum figure size
w = Cwidth+Lwidth+Uwidth;
h = max([Cheight,Lheight,Uheight],[],2);

% distribute width & height according to map
width = zeros(size(map));
height = zeros(size(map));
for n = 1:num
   idx = find(map==n);
   [i,j] = ind2sub(dim,idx);
   I = unique(i);
   J = unique(j);
   n_row = numel(I);
   n_col = numel(J);
   
   width(idx) = (w(n)-(n_col-1)*Options.Sep)/n_col;
   height(idx) = (h(n)-(n_row-1)*Options.Sep)/n_row;
end

col_w = max(width,[],1);
row_h = max(height,[],2);

wmin_ctrls = sum(col_w) + Options.Sep*(dim(2)+1);
hmin_ctrls = sum(row_h) + Options.Sep*(dim(1)+1);
nbtns = numel(handles.btns);
btns_w = nbtns*DefaultSize.pushbutton(1) + 2*(nbtns-1)*Options.Sep;
btns_h = DefaultSize.pushbutton(2);

figSize(3:4) = [max(wmin_ctrls,btns_w)+3*Options.Sep hmin_ctrls+btns_h+3*Options.Sep];

set(fig,'Position',figSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sinfo.minfigsize = figSize(3:4);  % minimum figure size
sinfo.map = map;        % uicontrol map
sinfo.autoextend = autoextend; % auto-extend when resize (to use the full span of figure window)
sinfo.w_labels = Lwidth;    % width of labels
sinfo.w_ctrls = Cwidth;    % width of controls
sinfo.w_units = Uwidth;
sinfo.h_labels = Lheight;   % height of labels
sinfo.h_ctrls = Cheight;    % height of controls
sinfo.h_units = Uheight;
sinfo.CLoffset = CLoffset;   % label y offset w.r.t. control
sinfo.CUoffset = CUoffset;
sinfo.w_btns = DefaultSize.pushbutton(1);
sinfo.h_btns = DefaultSize.pushbutton(2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKPROMPT :: Check Prompt input is valid & fill default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Prompt,FieldNames,Units,err] = checkprompt(Prompt)

% default configuration
FieldNames = {}; % answer in a cell
Units = {}; % no units

% standard error
err = {'inputsdlg:InvalidInput','Prompt must be a cell string with up to three columns.'};

if isempty(Prompt), Prompt = {'Input:'};
elseif ~iscell(Prompt), Prompt = cellstr(Prompt);
end

[nrow,ncol] = size(Prompt);

% prompt given in a row -> transpose
if ncol>3
   if nrow<3, Prompt = Prompt'; [nrow,ncol] = size(Prompt);
   else return; end
end

% struct fields defined
if ncol>1 && ~all(cellfun('isempty',Prompt(:,2))), FieldNames = Prompt(:,2); end

% unit labels defined
if ncol>2, Units = Prompt(:,3);
else Units = repmat({''},nrow,1); end

% return only the labels
Prompt = Prompt(:,1);

err = {}; % all cleared

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKFORMATS :: Check Formats input is valid & fill default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Formats,err] = checkformats(Formats,NumQuest)

fields  = [{'type' 'style' 'items' 'format' 'limits'  'size'};...
           {'edit' 'edit'  {''}  'text'  [0 1]   [0 0]}]; % defaults
err = {};

if isempty(Formats)
   Formats = repmat(struct(fields{:}),NumQuest,1);
   return
end

% backward compatibility (NumLines)
if isnumeric(Formats)
   [rw,cl]=size(Formats);
   ok = rw==1;
   if ok
      OneVect = ones(NumQuest,1);
      if cl == 2, NumLines=Formats(OneVect,:);
      elseif cl == 1, NumLines=Formats(OneVect);
      elseif cl == NumQuest, NumLines = Formats';
      else ok = false;
      end
   end
   if rw == NumQuest && any(cl == [1 2]), NumLines = Formats;
   elseif ~ok
      err = {'MATLAB:inputdlg:IncorrectSize', 'NumLines size is incorrect.'};
      return;
   end
   
   % set to default edit control (column stacked)
   Formats = repmat(struct(fields{:}),NumQuest,1);
   
   % set limits according to NumLines(:,1)
   numlines = mat2cell([zeros(NumQuest,1) NumLines(:,1)],ones(NumQuest,1),2);
   [Formats.limits] = deal(numlines{:});
   
   % sets the width to be 10*NumLines(:,2)
   if (size(NumLines,2) == 2)
      sizes = mat2cell([zeros(NumQuest,1) NumLines(:,2)],ones(NumQuest,1),2);
      [Formats.size] = deal(sizes{:});
   end
   
   return;
elseif ischar(Formats) || iscellstr(Formats) % given type
   if ischar(Formats), Formats = cellstr(Formats); end
   Formats = cell2struct(Formats,'type',3);
elseif ~isstruct(Formats)
   err = {'inputsdlg:InvalidInput','FORMATS must be an array of structure.'};
   return
end

% check fields
idx = find(~isfield(Formats,fields(1,:)));
for n = idx % if does not exist, use default
   [Formats.(fields{1,n})] = deal([]);
end

% set string fields to lower case
c = lower(cellfun(@char,{Formats.type},'UniformOutput',false));
[Formats.type] = deal(c{:});
c = lower(cellfun(@char,{Formats.style},'UniformOutput',false));
[Formats.style] = deal(c{:});
c = lower(cellfun(@char,{Formats.format},'UniformOutput',false));
[Formats.format] = deal(c{:});

% check number of entries matching NumQuest (number of PROMPT elements)
if sum( ~strcmp('none',{Formats.type}) & ~strcmp('',{Formats.type}) )~=NumQuest
   err = {'inputsdlg:InvalidInput','FORMATS must have matching number of elements to PROMPT (exluding ''none'' type).'};
   return
end

% check type field contents
if ~isempty(setdiff({Formats.type},{'check','edit','list','range','text','none',''}))
   err = {'inputsdlg:InvalidInput','FORMATS.type must be one of {''check'',''edit'',''list'',''range'',''none''}.'};
   return
end

% check format
if ~isempty(setdiff({Formats.format},{'text','float','integer','file','dir',''}))
   err = {'inputsdlg:InvalidInput','FORMATS.format must be one of {''text'', ''float'', or ''integer''}.'};
   return
end

num = numel(Formats);
for n = 1:num
   if isempty(Formats(n).type), Formats(n).type = 'none'; end
   
   switch Formats(n).type
      case 'none'
         if isempty(Formats(n).limits), Formats(n).limits = [0 0]; end
      case 'text'
         if isempty(Formats(n).style)
            Formats(n).style = 'text';
         elseif ~strcmp(Formats(n).style,'text');
            err = {'inputsdlg:InvalidInput','FORMATS.style for ''text'' type must be ''text''.'};
            return
         end
      case 'check'
         % check style
         if isempty(Formats(n).style)
            Formats(n).style = 'checkbox';
         elseif ~strcmp(Formats(n).style,'checkbox')
            err = {'inputsdlg:InvalidInput','FORMATS.style for ''check'' type must be ''checkbox''.'};
            return
         end
         
         % check format
         if isempty(Formats(n).format)
            Formats(n).format = 'integer';
         elseif any(strcmp(Formats(n).format,{'text','float'}))
            err = {'inputsdlg:InvalidInput','FORMATS.format for ''check'' type must be either ''float'' or ''integer''.'};
            return
         end
      case 'edit'
         % check style
         if isempty(Formats(n).style)
            Formats(n).style = 'edit';
         elseif ~strcmp(Formats(n).style,'edit')
            err = {'inputsdlg:InvalidInput','FORMATS.style for ''edit'' type must be ''edit''.'};
            return
         end
         
         % check format
         if isempty(Formats(n).format)
            Formats(n).format = 'text';
         elseif strcmp(Formats(n).format,'file') && isempty(Formats(n).items)
            Formats(n).items = {'*.*' 'All Files'};
         end
      case 'list'
         % check style
         if isempty(Formats(n).style)
            Formats(n).style = 'popupmenu';
         elseif ~any(strcmp(Formats(n).style,{'listbox','popupmenu','radiobutton','togglebutton'}))
            err = {'inputsdlg:InvalidInput','FORMATS.style for ''list'' type must be either ''listbox'', ''popupmenu'', ''radiobutton'', or ''togglebutton''.'};
            return
         end
         
         % check items
         if isempty(Formats(n).items)
            err = {'inputsdlg:InvalidInput','FORMATS.items must contain strings for the ''list'' type.'};
            return
         end
         if ~iscell(Formats(n).items)
            if ischar(Formats(n).items)
               Formats(n).items = cellstr(Formats(n).items);
            elseif isnumeric(Formats(n).items)
               Formats(n).items = cellfun(@num2str,num2cell(Formats(n).items),'UniformOutput',false);
            else
               err = {'inputsdlg:InvalidInput','FORMATS.items must be either a cell of strings or of numbers.'};
               return
            end
         end
         
         % check format
         if isempty(Formats(n).format)
            Formats(n).format = 'integer';
         elseif any(strcmp(Formats(n).format,{'text','float'}))
            err = {'inputsdlg:InvalidInput','FORMATS.format must be ''integer'' for the ''list'' type.'};
            return
         end
      case 'range'
         % check style
         if isempty(Formats(n).style)
            Formats(n).style = 'slider';
         elseif ~strcmp(Formats(n).style,'slider')
            err = {'inputsdlg:InvalidInput','FORMATS.style for ''range'' type must be ''slider''.'};
            return
         end
         
         if isempty(Formats(n).format)
            Formats(n).format = 'float';
         elseif ~strcmp(Formats(n).format,'float')
            err = {'inputsdlg:InvalidInput','FORMATS.format for ''range'' type must be ''float''.'};
            return
         end
   end
   
   % check limits
   if isempty(Formats(n).limits)
      if strcmp(Formats(n).style,'edit') && any(strcmp(Formats(n).format,{'integer','float'}))
         Formats(n).limits = [-inf inf]; % default for numeric edit box
      else
         Formats(n).limits = [0 1]; % default for all other controls
      end
   else
      if ~isnumeric(Formats(n).limits) || ~any(numel(Formats(n).limits)==[1 2]) || any(isnan(Formats(n).limits)) ... % 2-element numeric vector and cannot be NaN
            || (~(strcmp(Formats(n).type,'check')||(strcmp(Formats(n).type,'edit')&&strcmp(Formats(n).format,'file'))||strcmp(Formats(n).type,'none')) && Formats(n).limits(1)>Formats(n).limits(2)) ... % if not check or none, has to be increasing
            || (~(strcmp(Formats(n).type,'edit') && any(strcmp(Formats(n).format,{'float','integer'}))) && any(isinf(Formats(n).limits))) % if not edit-text, has to be finite
         
         err = {'inputsdlg:InvalidInput','FORMATS.limits must be increasing and non-NaN.'};
         return
      end
      if numel(Formats(n).limits)==1, Formats(n).limits(2) = 0; end
   end
   
   % check size
   if isempty(Formats(n).size)
      Formats(n).size = [0 0]; % default to auto-size
   else
      if ~isnumeric(Formats(n).size) || ~any(numel(Formats(n).size)==[1 2]) || any(isnan(Formats(n).size))
         err = {'inputsdlg:InvalidInput','FORMATS.size must be non-NaN.'};
         return
      end
      if numel(Formats(n).size)==1
         Formats(n).size(2) = 0;
      end
   end
   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKDEFAULTS :: Check the specified default values are compatible
%%% with Formats and if one not given fill in an initial value
function [DefAns,DefStr,err] = checkdefaults(DefAns,Formats,FieldNames)

DefStr = struct([]); % struct not used

% trim Formats to only include relevant entries (non-'none' types)
Formats = Formats'; % go through row first
Formats = Formats(~strcmp('none',{Formats.type}));
len = numel(Formats);

if isempty(DefAns) % if DefAns not given
   DefAns = cell(len,1); % will set DefAns to default values
elseif isstruct(DefAns)
   if isempty(FieldNames) % struct return
      err = {'inputsdlg:InvalidInput','Default answer given in a structure but the prompts do not have associated answer field names in the second column.'};
      return;
   end
   if ~all(isfield(DefAns,FieldNames)|cellfun('isempty',FieldNames))
      err = {'inputsdlg:InvalidInput','Default answer structure is missing at least one answer field.'};
      return;
   end
   DefStr = DefAns;
   DefAns = cell(len,1);
   for k = 1:len
      if ~isempty(FieldNames{k}), DefAns{k} = DefStr.(FieldNames{k}); end
   end
elseif ~iscell(DefAns)
   err = {'inputsdlg:InvalidInput','Default answer must be given in a cell array or a structure.'};
   return;
elseif length(DefAns)~=len
   err = {'inputsdlg:InvalidINput','Default answer cell dimension disagrees with the number of prompt'};
   return;
end

err = {'inputsdlg:InvalidInput','Default value is not consistent with Formats.'};

% go through each default values
for k = 1:len
   if isempty(DefAns{k}) % set non-empty default values
      switch Formats(k).type
         case 'check' % off
            DefAns{k} = Formats(k).limits(1);
         case 'edit'
            switch Formats(k).format
               case {'text','file','dir'}
                  DefAns{k} = ''; % change to empty string
               case {'float','integer'}
                  liminf = isinf(Formats(k).limits);
                  if all(liminf) % both limits inf
                     DefAns{k} = 0;
                  elseif any(liminf) % 1 limit inf
                     DefAns{k} = Formats(k).limits(~liminf);
                  else
                     DefAns{k} = round(mean(Formats(k).limits));
                  end
            end
         case 'list' % first item
            DefAns{k} = 1;
         case 'range' % middle value
            DefAns{k} = mean(Formats(k).limits);
      end
   else % check given values
      msel = strcmp(Formats(k).style,'listbox') && diff(Formats(k).limits)>1;
      switch Formats(k).format
         case 'text'
            if ~(isempty(DefAns{k}) || ischar(DefAns{k})),
                1449
                k
                Formats(k)
                DefAns(k)
                return; 
            end
         case 'float'
            if ~isfloat(DefAns{k}) || numel(DefAns{k})~=1, 
                1457
                k
                Formats(k)
                DefAns(k)
                return; 
            end
         case 'integer' % can be multi-select if type=list
            if ~islogical(DefAns{k}) ...
                  && (~isnumeric(DefAns{k}) || any(DefAns{k}~=floor(DefAns{k})) || (~msel && numel(DefAns{k})~=1))
               1466
               k
               Formats(k)
               DefAns(k)
               return;
            end
         case 'file'
            dlim = diff(Formats(k).limits);
            if ~isempty(DefAns{k}) && dlim>=0 % for uisetfile
               if dlim<=1 % single-file
                  files = DefAns(k);
               else %dlim>1 % multiple-files
                  files = DefAns{k};
               end
               for n = 1:length(files)
                  d = dir(files{n});
                  if ~isempty(files{n}) && length(d)~=1
                    1483
                    k
                    n
                    files{n}
                     err{2} = 'Default file name does not resolve to a unique file.';
                     return;
                  end
               end
            end
         case 'dir'
            if isempty(DefAns{k}) || ~isdir(DefAns{k})
               err{2} = 'Default directory does not exist.';
                1495
                k
                Formats(k)
                DefAns(k)
               return;
            end
      end
      
      switch Formats(k).type
         case 'check' % must be one of the values in limits
            if all(DefAns{k} ~= Formats(k).limits), return; end
         case 'edit' % if numeric, must be within the limits
            % user input only known to 3dp, so our accuracy is +/- 1e-3
            eps=1e-3;
            if any(strcmp(Formats(k).format,{'float','integer'})) ...
                  && (DefAns{k}<Formats(k).limits(1)-eps || DefAns{k}>Formats(k).limits(2)+eps)
                1509
                k
                Formats(k)
                DefAns(k)
               return;
            end
            
         case 'list' % has to be valid index to the list
            if any(DefAns{k}<1) || any(DefAns{k}>numel(Formats(k).items)), 
                1518
                k
                Formats(k)
                DefAns(k)
                return; 
            end
            
         case 'range' % has to be within the limits
            if DefAns{k}<Formats(k).limits(1) || DefAns{k}>Formats(k).limits(2)
                1527
                k
                Formats(k)
                DefAns(k)
               return;
            end
      end
   end
end

% also initialize DefStr if FieldNames given
if isempty(DefStr) && ~isempty(FieldNames)
   idx = ~cellfun('isempty',FieldNames);
   StructArgs = [FieldNames(idx) DefAns(idx)]';
   DefStr = struct(StructArgs{:});
end

err = {}; % all good
end

function [Options,err] = checkoptions(Options)

err = {'inputsdlg:InvalidInput','Options must be ''on'', ''off'', or a struct.'};

Fields = {'Resize',   'off'
   'WindowStyle', 'normal'
   'Interpreter', 'tex'
   'ApplyButton', 'off'
   'Sep',     10}';

if isempty(Options)
   Options = struct(Fields{:});
elseif ~isstruct(Options)
   if ~(ischar(Options)||iscellstr(Options)), return; end % must be struct or a string
   Fields{2,1} = char(Options);
   Options = struct(Fields{:});
end

if ~isfield(Options,'Resize') || isempty(Options.Resize)
   Options.Resize = 'off';
elseif ~strcmpi(Options.Resize,{'on','off'})
   err{2} = 'Resize option must be ''on'' or ''off''.';
   return;
end
if ~isfield(Options,'WindowStyle') || isempty(Options.WindowStyle)
   Options.WindowStyle = 'normal';
elseif ~strcmpi(Options.WindowStyle,{'normal','modal','docked'})
   err{2} = 'WindowStyle option must be ''normal'' or ''modal''.';
   return;
end
if ~isfield(Options,'Interpreter') || isempty(Options.Interpreter)
   Options.Interpreter = 'tex';
elseif ~strcmpi(Options.Interpreter,{'latex','tex','none'})
   err{2} = 'Interpreter option must be ''latex'', ''tex'', or ''none''.';
   return;
end
if ~isfield(Options,'ApplyButton') || isempty(Options.ApplyButton)
   Options.ApplyButton = 'off';
elseif ~strcmpi(Options.ApplyButton,{'on','off'})
   err{2} = 'ApplyButton option must be ''on'' or ''off''.';
   return;
end
if ~isfield(Options,'Sep') || isempty(Options.Sep)
   Options.Sep = 10;
elseif ~isnumeric(Options.Sep) || Options.Sep<0
   err{2} = 'Sep option must be non-negative scalar value.';
   return;
end

err = {}; % all cleared
end

% Copyright (c)2009, Takeshi Ikuma
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   * Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer. *
%   Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%   documentation and/or other materials provided with the distribution.
%   * Neither the name of the <ORGANIZATION> nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
