function im= MakeCrack1(SI,SJ,SK,cx,cy,cz,R,I)
%
%   Program generates a 3D sphere image 
%
%   Usage:
%
%           im=MakeSphereImage(SI,SJ,SK,cx,cy,cz,R)
%
%           SI, SJ, SK : image size in I, J and K directions
%           cx, cy, cz : centroid of the sphere
%           R          : radius of the sphere
%           I          : color index of the sphere (1-255)
%   
%
% by:  Muhammed Emin Kutay


im(SI,SJ,SK)=uint8(0);
%im=im;
tic
disp('generating the image...')
for x=1:SJ %(SJ/2-2):(SJ/2+2)
    
    for y=1:SI %(SI/2-1):(SI/2+1)
    
        for z=1:SK
            
            r=sqrt((x-cx)^2 + (y-cy)^2); % + (z-cz)^2);
            
            if r<=R
                
                im(y,x,z)=255;
                
             end
            
        end
    end
end
toc

for i=1:size(im,3)
    
    imshow(im(:,:,i))
    title(['Slice- ' num2str(i)])
    pause(0.01)

end
