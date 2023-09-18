function [im, ImgBaseName, FilePath]  = ImLoadStack

[filename, FilePath] = uigetfile2('*.*', 'Select First Image to open (only black and white (binary) image where black = 0 = solid, white = 1 = pore ');
if filename==0
    im=[];
    ImgBaseName=[];
    FilePath=[];
    return
end

    
ImgBaseName=filename(1:size(filename,2)-7);
Ext=filename(size(filename,2)-3:size(filename,2));
tifs=dir([FilePath '*' Ext]);
k=0;
for i=1:size(tifs,1),
    SImg=size(ImgBaseName,2);
    if size(tifs(i).name,2)>=SImg
    if tifs(i).name(1:SImg)==ImgBaseName, 
        k=k+1;
    end 
end
end
try load def, catch result     = {ImgBaseName,Ext,'1',num2str(k),'1','[.5,.5,.5]'}; end
result=inputdlg({'Image Base Name: ', 'Extension', 'First image number to start: ',...
                 ['Number of Files (Total ' num2str(k) ' files)'],...
                 'Numbering style (type 1 for numbering 1, 2, 5...), (type 0 for 001, 002, 003 ...)'},...
                 'Input Image Information ',1,result);
ImgBaseName=result{1};
Ext=result{2};
StartImage=str2num(result{3});
NumFiles=str2num(result{4});
NumbStyle = str2num(result{5});

clear tifs

save def result

imp=imread([FilePath,filename]);
im(size(imp,1),size(imp,2), NumFiles) = logical(0);

EndImage=NumFiles+StartImage-1;
cancelButton=0;
h = waitbar(0,'Loading Images', ...
    'Position',[ 243   343   270    93],...
    'Name','Please wait ...',...
    'CreateCancelBtn','cancelButton = 1;closereq;' );
for k=StartImage:EndImage
    
    waitbar((k-StartImage+1)/NumFiles,h,['Reading Images, Current Image= ' num2str(k)])
    if cancelButton==1, break, end
    Z=k-(StartImage)+1;
    
%     try
    if NumbStyle ==1
    imm=imread([FilePath ImgBaseName num2str(k) Ext]);
    else
    imm=imread([FilePath ImgBaseName numtostr(k) Ext]);
    end

        %imm=imread([FilePath ImgBaseName numtostr(k) Ext]);
        if size(imm,3)>1,
            
            imm=rgb2gray(imm);
            
        end

        im(:,:,Z)=imm;
%     catch
%         errordlg(['File ' [FilePath ImgBaseName numtostr(k) Ext] 'does not exist.'])
%         return
%     end
    
end
close(h)



