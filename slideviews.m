function varargout = slideviews(varargin)
% slideviews Application M-file for slideviews.fig
%    FIG = slideviews launch slideviews GUI.
%    slideviews('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 31-Oct-2014 12:23:52

if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'new');
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end
    
elseif nargin==1
    im=varargin{1};
    
    fig = openfig(mfilename,'new');
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    handles.im=im;
    set(handles.slider1,'min',1,'max',size(handles.im,3),'sliderstep',[1/(size(handles.im,3)-1), .1])
    guidata(fig, handles);
    slider1_Callback('', '', handles, '')
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif nargin==2
    
    im=varargin{1};
    ftitle=varargin{2};
    
    fig = openfig(mfilename,'new');
    set(fig,'Name',ftitle)
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    handles.im=im;
    set(handles.slider1,'min',1,'max',size(handles.im,3),'sliderstep',[1/(size(handles.im,3)-1), .1])
    guidata(fig, handles);
    slider1_Callback('', '', handles, '')
    if nargout > 0
        varargout{1} = fig;
    end
    
end


% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)

% disp('callback')
val=get(handles.slider1,'value');
set(handles.slider1,'value',round(val));
set(handles.text1,'String',['Image-' num2str(round(val))]);
set(handles.figure1,'selected','on');
% guidata(slideviews, handles);
handles.figure1;
im=squeeze(handles.im(:,:,round(val)));

if isa(im,'uint8') || isa(im,'logical')
    imshow(im);
else
    imagesc(im)
    colormap('jet')
    colorbar
end
%
% imagesc(im)
%     colormap('jet')
%     colorbar

axis off

daspect([1,1,1])


% --------------------------------------------------------------------
function varargout = figure1_CreateFcn(h, eventdata, handles, varargin)


% % --------------------------------------------------------------------
function varargout = slider1_ButtonDownFcn(h, eventdata, handles, varargin)
%

% --------------------------------------------------------------------
function varargout = Pixval_On_Callback(h, eventdata, handles, varargin)

impixelinfo
imdistline
