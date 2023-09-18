im = data.im3D;
d= load( 'Data_Microscopic_Vels_D3Q9.mat','e');
e3D = d.e;
data.e3D = e3D;
data.ao3D=[2,1,4,3,6,5,9,10,7,8,13,14,11,12,17,18,15,16,19]; %opposite directions

d= load( 'Data_Weight_Factors_D3Q9.mat','Wght');
data.WF3D = d.Wght;


%--define coordinates
SI  =size(im,1);
SJ  =size(im,2);
SK = size(im,3);

[x,y,z]  = meshgrid(1:SJ,1:SI,1:SK);
x      = reshape(x, [SI*SJ*SK,1]);
y      = reshape(y, [SI*SJ*SK,1]);
z      = reshape(z, [SI*SJ*SK,1]);

Ncoors = length(x);

data.Ncoors              = Ncoors;
data.x                       =x;
data.y                       =y;
data.z                       =z;
data.SI = SI;
data.SJ = SJ;
data.SK = SK;

disp('Determining neighbor indices...')
pause(0.01)
tic;

W = data.W;
H = data.H;
Ineigbor = zeros(Ncoors,19);

for a = 1:19
    
    xn_e= x - e3D(a,1);
    yn_e = y - e3D(a,2);
    zn_e = z - e3D(a,3);
    
    % periodic boundary conditions
    xn_e(xn_e>SJ) = 1;
    xn_e(xn_e==0) = SJ;
    
    yn_e(yn_e>SI) = 1;
    yn_e(yn_e==0) = SI;
    
    zn_e(zn_e>SK) = 1;
    zn_e(zn_e==0) = SK;
    
    
    Ix_es = sub2ind([SI, SJ, SK], yn_e, xn_e, zn_e);
    
    Ix_e3D(:,a) = Ix_es; %#ok<SAGROW>
    
end%for a = 1:9

data.Ix_e3D =  Ix_e3D;

%Boundary nodes
yNorth = ones(SJ,1);
xNorth = (1:SJ)';
zNorth = (1:SK)';

Indtop = find(z==1);
Indbot = find(z==SK);
IndNorth = find(y==SI);
IndSouth = find(y==1);
IndEast = find(x==SJ);
IndWest = find(x==1);
Indmid= find(z==SK/2 & y==SI/2);
IndZmid= find(z==SK/2);


data.Indtop = Indtop;
data.Indbot = Indbot;
data.IndNorth = IndNorth;
data.IndSouth = IndSouth;
data.IndWest = IndWest;
data.IndEast = IndEast;
data.Indmid = Indmid;
data.IndZmid = IndZmid;

uiwait( msgbox('3D INDEXING DONE....'))

toc



data.progdata{3} =  'COMPLETE';
data.progdata{4} = '------';

