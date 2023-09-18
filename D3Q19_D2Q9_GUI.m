function D3Q19_D2Q9_GUI

close all hidden
clc

%%


mainFIG = figure;
%
% data.hax0 = axes('position',[0.2, 0.1, 0.6, 0.5]);
% get(data.hax0 )




data.hax = axes('position',[0.45, 0.2, 0.5, 0.5],'color','none');

%% Create figure
% scs = get(0,'screensize');

% pos = [.1*scs(3), .1*scs(4), .5*scs(3), .5*scs(4)];
pos = [90 80 850 450];
set(mainFIG,'numbertitle','off',...
    'name','D2Q9 and D3Q19 Lattice Boltzmann Fluid Flow Modeling Software', ...
    'position',pos,...
    'menubar','none', ...
    'toolbar','figure', ...
    'resize','off',...
    'color','white',...
    'DockControls','off', ...
    'HandleVisibility','callback');

%% Change the figure manubar items

data.hToolbar = findall(mainFIG,'tag','FigureToolBar');
get(findall(data.hToolbar),'tag')

data.hNewButton = findall(data.hToolbar,'tag','Standard.NewFigure');
set(data.hNewButton, 'ClickedCallback',@NewProject, ...
    'TooltipString','New Project');

data.hOpenButton = findall(data.hToolbar,'tag','Standard.FileOpen');
set(data.hOpenButton, 'ClickedCallback',@OpenProject, ...
    'TooltipString','Open Project');

data.hSaveButton = findall(data.hToolbar,'tag','Standard.SaveFigure');
set(data.hSaveButton, 'ClickedCallback',@SaveProject, ...
    'TooltipString','Save  Project');

%% figure position



pos=get(mainFIG,'position');% Create button
W = pos(3) ;
H = pos(4) ;

data.W = W;
data.H = H;



%% File menu

f = uimenu(mainFIG,'Label','File');
uimenu(f,'Label','New Project','Callback',@NewProject);
uimenu(f,'Label','Open Project','Callback',@OpenProject,'Accelerator','O');
uimenu(f,'Label','Save Project','Callback',@SaveProject, 'Accelerator','S');
uimenu(f,'Label','Save Project As','Callback',@SaveProjectAs);

data.mnuFile = f;



%% edit

% f = uimenu(mainFIG,'Label','Edit');
%     uimenu(f,'Label','Delete Push Pull Data','Callback',@DeletePPdata)
%     uimenu(f,'Label','Edit Push Pull Data','Callback',@EditPPdata)


%% D2Q9

f = uimenu(mainFIG,'Label','D2Q9 Simulation Steps');
uimenu(f,'Label','Step-1: Load the image file','Callback',@LoadImage);
uimenu(f,'Label','Step-2: Load Input Parameters','Callback',@LoadInputParameters,'separator','on');
uimenu(f,'Label','Step-3: Initialize Indices ','Callback',@InitializeIndices_bounceback,'separator','on');
uimenu(f,'Label','Step-4a: Run Simulation no k','Callback',@RunSimulation_nok,'separator','on', 'enable','on');
% uimenu(f,'Label','Step-4b: Run Simulation with ksat','Callback',@RunSimulation_withk);
% uimenu(f,'Label','Step-4c: Run Simulation for asphalt foam','Callback',@RunSimulation_foam,'separator','on','enable','off');

data.mnuAnalMonDat = f;



%%
f = uimenu(mainFIG,'Label','D3Q19 Simulation Steps');
uimenu(f,'Label','Step-1: Load the image slices','Callback',@LoadImageSlices);
uimenu(f,'Label','Step-2: Load Input Parameters','Callback',@LoadInputParameters,'separator','on');
uimenu(f,'Label','Step-3: Initialize Indices ','Callback',@InitializeIndices_bounceback_D3Q9,'separator','on');
uimenu(f,'Label','Step-4b: Run D3Q9 Simulation with ksat','Callback',@RunD3Q9Simulation_withk);

data.mnuAnalMonDat = f;




%% view menu


f = uimenu(mainFIG,'Label','View & Plot');
uimenu(f,'Label','Show Variables in Workspace','Callback',@viewvariables, 'Accelerator','P')
uimenu(f,'Label','MakeMovie','Callback',@makemovie)
uimenu(f,'Label','Plot Saved Data','Callback',@plotsaveddata)


data.mnuView = f;


%% help
f = uimenu(mainFIG,'Label','Help');
uimenu(f,'Label','Tutorial ','Callback',@f_Documentation)

uimenu(f,'Label','Formulations used in the software','Callback',@f_Formulations,'separator','on')

uimenu(f,'Label','Developer','Callback',@f_logo,'separator','on')

data.mnuData = f;


colnames = {'Progress'};
pos = [W * (1-0.35), H * (1-0.20), W * 0.35, H * 0.20];
columneditable =  [false];
columnformat = {'char'};
RowName = { '1.Load Image','2.Load Parameters','3.Initialize Indices',  '4.Run Simulation'};

progdata = {'-----';'-----';'-----';'-----'};

data.tbl_Progress= uitable(mainFIG, 'ColumnName', colnames, ...
    'Position', pos,...
    'ColumnEditable', columneditable,'Enable','on',...
    'ColumnFormat', columnformat,...
    'RowName',RowName,...
    'data',progdata,...
    'ButtonDownFcn',@f_progtable_btndown, ...
    'HitTest','on');


mtable = get(data.tbl_Progress);


data.progdata = progdata;


data.figureWidth = W;
data.figureHeight = H;


guidata(mainFIG, data)


function plotsaveddata(mainFIG, HDL)

data   = guidata(mainFIG);

f_Plot_State


function makemovie(mainFIG, HDL)

data   = guidata(mainFIG);

f_makemov3



function viewvariables(mainFIG, HDL)

data   = guidata(mainFIG)




function RunSimulation_foam(mainFIG, HDL) %#ok<*INUSD>

data   = guidata(mainFIG);

f_run_D2Q9_simul_foamWA1(data)





function RunSimulation_nok(mainFIG, HDL) %#ok<*INUSD>

data   = guidata(mainFIG);

% for validation. no k
f_run_D3Q19_simul(data)



function RunD3Q9Simulation_withk(mainFIG, HDL)

data   = guidata(mainFIG);

f_run_D3Q19_simul(data)




function RunSimulation_withk(mainFIG, HDL)

data   = guidata(mainFIG);

% with ksat
f_run_D2Q9_simul_Korner3(data)


%%
function f_progtable_btndown(mainFIG, HDL)



data   = guidata(mainFIG);
disp('hello')

%get(data.tbl_Progress, 'selected')
mtable =   data.tbl_Progress;

jtable = mtable.getTable;
row = jtable.getSelectedRow + 1 % Java indexes start at 0
col = jtable.getSelectedColumn + 1





function f_Documentation(mainFIG, HDL)


web('./docs/tutorial.htm','-browser')


function f_Formulations(mainFIG, HDL)


web('./docs/Formulations.htm','-browser')



function LoadImageSlices(mainFIG, HDL)

data   = guidata(mainFIG);

%%
[im, ImgBaseName, FilePath]  = ImLoadStack;

if isempty(im)
    return
    
elseif size(im,3)==1
    figure,
    imshow(im)
    impixelinfo
    
    data.im3D=im;
    data.ImgBaseName = ImgBaseName ;
    data.FilePath = FilePath ;
    
    guidata(mainFIG, data)
    
else
    
    
    data.im3D=im;
    data.ImgBaseName = ImgBaseName ;
    data.FilePath = FilePath ;
    
    slideviews(data.im3D)
    guidata(mainFIG, data)
end



%
%         SI = 60; SJ=60; SK=60;
%         cx = 30; cy=30; cz = 30;
%         R = 15;
%         I = 0;
%
%         data.im3D= MakeCrack1(SI,SJ,SK,cx,cy,cz,R,I);
% 
% slideviews(data.im3D)
% guidata(mainFIG, data)


function LoadImage(mainFIG, HDL)

data   = guidata(mainFIG);

%--load image data
[KFileName, FilePath] = uigetfile2({'*.tif'});
data.im=imread([FilePath,KFileName]);
data.KFileName = KFileName;
data.FilePath = FilePath;

guidata(mainFIG, data)

showim(mainFIG, data)



function showim(mainFIG, HDL)

data   = guidata(mainFIG);
W= data.W;

cla(data.hax,'reset')
set(data.hax,'visible','on')


imshow(data.im)
axis on,
hp = impixelinfo;
set(hp,'Position',[W-300, 0,100,10]);

[SI,SJ]=size(data.im);
text(0.05*SJ, 0.2*SI,'255=W, 151=A, 0=S','color','blue')

data.progdata{1} =  'COMPLETE';
set(data.tbl_Progress, 'data',data.progdata);
guidata(mainFIG, data)



function InitializeIndices_bounceback_D3Q9(mainFIG, HDL)

data   = guidata(mainFIG);

s_initialize_indices_bounceback_D3Q19

guidata(mainFIG, data)



function InitializeIndices_bounceback(mainFIG, HDL)

data   = guidata(mainFIG);

s_initialize_indices_bounceback3
%  s_initialize_indices
guidata(mainFIG, data)





function LoadInputParameters(mainFIG, HDL)

data   = guidata(mainFIG);

set(data.hax,'visible','on')

W = data.figureWidth;
H = data.figureHeight;
SI = 0;

btnwdth = 0.10;
btnht = 0.045;



colnames = {'Value', 'Unit'};
rownames = {'dx = dy'; ...
    'dt'; ...
    'dm'; ...
    'Viscosity (W)';...
    'Viscosity (A)';...
    'Density (W)';...
    'Density (A)'; ...
    'omege_alpha';...
    'Surf. Tens. (sigma)';...
    'Grav. accel';...
    'Plot time step';...
    'Save time step'; ...
    'Perm. of solids'};


if ~isfield(data,'inputparam')
    
    try load('LE_inputparam.mat' )
    catch
        inputparam = {'1',             'mm/pixel';...
            '0.001',     'sec/time';...
            '1',            'g/mass';,
            '1',            'mm2/s'; ...
            '17',           'mm2/s'; ...
            '0.001',       'g/mm3';...
            '0.001',        'g/mm3';...
            '1.99999',    ' ';...
            '72',             'mN/m';...
            '9.81e3',       'mm/s2';...
            '10',              'timesteps ';...
            '1000',              'timesteps '; ...
            '0', ''};
    end
    
    
    data.inputparam = inputparam;
    
end

pos = [W * 0.01, H * .1, W * 0.37, H * 0.8];
columneditable =  [true false];
columnformat = {'numeric', 'char'};

data.tbl_D2Q9InputData= uitable(gcbf, 'ColumnName', colnames, ...
    'rowname', rownames, ...
    'Position', pos,...
    'ColumnEditable', columneditable,'Enable','on',...
    'ColumnFormat', columnformat,...
    'Data',data.inputparam,...
    'ColumnWidth',{60,60,60});


% Button OK
spos = [W * 0.01, pos(2)-24, W * btnwdth, H *btnht];
data.btnInputTableDone   = uicontrol(gcbf,'style','pushbutton',...
    'Position',spos,'string','Done', 'Enable','on', 'ForegroundColor','Blue',...
    'callback',@InputTableDone);



get(data.btnInputTableDone)



guidata(mainFIG, data)



function InputTableDone(mainFIG, HDL)

data   = guidata(mainFIG);

set(data.btnInputTableDone,'visible','off')
set(data.tbl_D2Q9InputData, 'ColumnEditable', [false, false])

data.inputparam = get(data.tbl_D2Q9InputData,'Data');
inputparam = data.inputparam ;
save LE_inputparam inputparam



data.dx = str2num(data.inputparam{1,1});
data.dt = str2num(data.inputparam{2,1});
data.dm = str2num(data.inputparam{3,1});
data.vw_r = str2num(data.inputparam{4,1});
data.va_r = str2num(data.inputparam{5,1});
data.rhow = str2num(data.inputparam{6,1});
data.rhoa =  str2num(data.inputparam{7,1});
data.omega_alpha =  str2num(data.inputparam{8,1});
data.sigma_r  =  str2num(data.inputparam{9,1});
data.grav_accel_r  =  str2num(data.inputparam{10,1});
data.tplotstep = str2num(data.inputparam{11,1});
data.tsavestep = str2num(data.inputparam{12,1});
data.Ns = str2num(data.inputparam{13,1});


data.co = 1/sqrt(3); % Lattice speed of sound
data.grav_accel_l=   data.grav_accel_r*(data.dt^2/data.dx) ;  %lb units
data.sigma_l = data.sigma_r * data.dt^2 / data.dm;
data.rhow_l = data.rhow * data.dx^3/data.dm; % mass / pixel3
data.rhoa_l = data.rhoa * data.dx^3/data.dm;
data.tao_w    =  0.5 + data.vw_r*data.dt/(data.co*data.dx)^2  ;   % 0.0126
data.tao_a     =  0.5 + data.va_r*data.dt/(data.co*data.dx)^2  ;   % 4.3e-5
data.vw_l   =   (1/3)*(data.tao_w-0.5);
data.va_l    =   (1/3)*(data.tao_a-0.5);
data.omega_w= 1/data.tao_w;
data.omega_a= 1/data.tao_a;



%data.ksat_l=data.ksat_w*data.dt/data.dx;



%--Lattice directional velocity vectors
e(1,:)=[ 1  0]; e(2,:)=[ 0  1]; e(3,:)=[-1  0]; e(4,:)=[ 0 -1];
e(5,:)=[ 1  1]; e(6,:)=[-1  1]; e(7,:)=[-1 -1]; e(8,:)=[ 1 -1];
e(9,:)=[ 0  0];

data.ao=[3,4,1,2,7,8,5,6,9]; %opposite directions

%--weight factors for each direction
WF(1:4)=1/9;      WF(5:8)=1/36;     WF(9)  =4/9;

data.e = e;
data.WF = WF;



data.progdata{2} =  'COMPLETE';
data.progdata{3} = '------';
data.progdata{4} = '------';

set(data.tbl_Progress, 'data',data.progdata);

guidata(mainFIG, data)












function f_logo(mainFIG, HDL)

data   = guidata(mainFIG) ;
mbox = msgbox({'D2Q9 and D3Q19 Lattice boltzmann software ';'_____________________' ;...
    'Developer:'; '  ';'M. Emin Kutay, Ph.D., P.E. '; ...
    'Professor';...
    'Department of Civil & Environmental Engineering';...
    '    ';...
    'Michigan State University';...
    '3554 Engineering Building';...
    'East Lansing, MI 48824-1226';...
    'Phone: (517) 353 9297';...
    'Fax : (517) 432 1827';...
    'E-mail: kutay@msu.edu';...
    'Web: http://www.egr.msu.edu/~kutay'      ; '-------------------';...
    'Please contact for bugs and questions (after trying to figure out the problem first ;))'},'About this software');

set(mbox','color','white')






%% New Project

function NewProject(mainFIG, HDL)

data            = guidata(mainFIG);

r               = questdlg ('All data in the workspace will be cleared, ok?', ...
    'New Project','Yes', 'No', 'No');
switch r
    
    case 'Yes'
        
        close(gcbf)
        clc
        D2Q9_GUI
        
    case 'No '
        
        return
        
        
end

%% Open Project

function OpenProject(mainFIG, HDL)


data            = guidata(mainFIG);
clc

[fname, fpath]  = uigetfile2('*.D2Q', 'Open Project');
if fname==0, return, end
data.fname      =   fname;
data.fpath      =   fpath;
try
    dataL = load([fpath, fname],'-mat');
catch
    errordlg('Not a *.D2Q file','File format error')
    return
end

dataL.data

if isfield(data,'im'), data = rmfield(data,'im'); end
if isfield(dataL.data,'im'),
    data.im = dataL.data.im;
    guidata(mainFIG, data)
    showim(mainFIG, data)
    
end


%if isfield(data,'progdata'), data = rmfield(data,'progdata'); end
if isfield(dataL.data,'progdata'),
    data.progdata = dataL.data.progdata ;
    set(data.tbl_Progress, 'data',data.progdata);
    
end

guidata(mainFIG, data)








%% Save Project As

function SaveProjectAs(mainFIG, HDL)

data            = guidata(mainFIG);
disp('Saving ...')
[fname, fpath]  = uiputfile2('*.D2Q', 'Save Project As');
if fname==0, return, end
data.fname      =   fname;
data.fpath      =   fpath;
save([fpath, fname],'data')
disp(['Saved to ' data.fpath, data.fname])
set(gcbf,'name',[fname, ' @ ' fpath] );
guidata(mainFIG, data)







%% Save Project

function SaveProject(mainFIG, HDL)

data            = guidata(mainFIG);

try
    disp('Saving ...')
    save([data.fpath, data.fname],'data')
    disp(['Saved to ' data.fpath, data.fname])
    
catch
    guidata(mainFIG, data)
    SaveProjectAs(mainFIG, data)
end







