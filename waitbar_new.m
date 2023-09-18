function fout = waitbar_new(x_val, parent, varargin)
% 	
% This m-file allows to specify a waitbar on a figure where the user wants to display it. The sample usage of waitbar_new is as given.  
%  
% % function test  
%         h = dialog('Units', 'Pixels', 'Position', [20 20 360 100],'visible','off');  
%         h1 = uicontrol('Parent',h,'Units', 'Pixels', 'Position', [10 70 75 25],'String','Test','Callback',@test_button);  
%         movegui(h,'center');  
%         set(h,'visible','on');  
%         
%         
%         function test_button(hObject,eventdata)  
%         wh=waitbar_new(0, gcf, 'r', [10 30 300 14]);  
%         %wh=waitbar_new(0, gcf, [0.8 .5 1]);  
%         %wh=waitbar_new(0, gcf);  
%         
%         for i=1:10000  
%             waitbar_new(i/10000);  
%         end  
%         delete(wh);
% %
%
% WAITBAR_NEW displays the waitbar on a figure specified.
% H = WAITBAR_NEW(0, figure_handle, color, position) displays a waitbar on
% the figure whose handle is figure_handle. color and position arguments
% are optional. If specified, color should be in one of Matlab ColorSpec format.
% position is a four element vector of the form [left, bottom, width, height]
% position vector should be in pixel values.
% H = WAITBAR_NEW(0, figure_handle, color) will create a waitbar at the
% left bottom of the figure specified.
% H = WAITBAR_NEW(0, figure_handle) will create a waitbar with default
% color and position settings.

if( nargin >= 2 )
    if isnumeric( parent )
        type = 2; 
        parent_handle = parent;
    else
        error('MATLAB:waitbar_new:InvalidInputs', ['Input arguments of type ' class(parent) ' not valid.'])
    end
elseif( nargin == 1 )
    axes_handle = findobj('Type', 'axes', 'Tag', 'waitbar_Axes');
    parent_handle = get(axes_handle, 'Parent');
    if isempty(parent_handle)
        return;
    else
        type = 1; % Update the waitbar
    end
else
    error('MATLAB:waitbar_new:InvalidArguments', 'Input arguments not valid.');
end
x_val = max(0, min(100*x_val, 100));
switch type
    case 1,  % waitbar(x_val)    update
        patch_handle = findobj(parent_handle, 'Type', 'patch');
        if ( isempty(parent_handle) || isempty(patch_handle) )
            error('MATLAB:waitbar_new:WaitbarHandlesNotFound', 'Couldn''t find waitbar handles.');
        end
        xpatch = [0 x_val x_val 0];
        set(patch_handle, 'XData', xpatch);
        
    case 2,  % waitbar(x_val,figure_handle, color, position)  initialize
        if nargin > 2,
            patch_color = varargin{1};
        else
            patch_color = 'b';
        end
        oldRootUnits = get(0, 'Units');
        set(0, 'Units', 'points');
        axFontSize = get(0, 'FactoryAxesFontSize');
        set(0, 'Units', oldRootUnits);
        colormap([]);
        if nargin > 3,
            axPos = varargin{2};
        else
            oldRootUnits = get(parent_handle, 'Units');
            set(parent_handle, 'Units', 'pixels');
            axPos = get(parent_handle, 'Position');
            axPos(1:2) = 4;
            axPos(3) = axPos(3)*.7;
            axPos(4) = 14;
            set(parent_handle, 'Units', oldRootUnits);
        end
        axes_handle = axes( 'Tag', 'waitbar_Axes',...
            'Parent', parent_handle,...
            'XLim',[0 100],...
            'YLim',[0 1],...
            'Box','on', ...
            'Units','pixels',...
            'FontSize', axFontSize,...
            'Position',axPos,...
            'XTickMode','manual',...
            'YTickMode','manual',...
            'XTick',[],...
            'YTick',[],...
            'XTickLabelMode','manual',...
            'XTickLabel',[],...
            'YTickLabelMode','manual',...
            'YTickLabel',[]);
        xpatch = [0 x_val x_val 0];
        ypatch = [0.05 0.05 0.9 0.9];
        patch_handle = patch(xpatch, ypatch, patch_color, 'EdgeColor', patch_color, 'EraseMode', 'none');
end  % case
drawnow;
if ( nargout == 1 )
    fout = axes_handle;
end
