function f_makemov


try load lastfld, catch fld = pwd; end
fld = uigetdir(fld,'Pick a directory of mat files');
if fld ==0, return, end

save lastfld fld

r = inputdlg('Prefix ? ')

prefx = r{1};

fnames=dir([fld, '/' prefx '*.mat'])



    
    for i=3:length(fnames)
        fdates{i-2,2} = fnames(i).name;
        fdates{i-2,1} = fnames(i).date;
%         fdates{i,1} = getfield(dir([fld,fnames{i}]),'date');
%         fdates{i,2} = fnames{i};
    end
    
    fdates = sortrows(fdates,1);
    
    
    avipath = [fld,  '/mov' datestr(now,'mmddyy-HHMMSS') ' .avi'];
    

% try    
%     aviobj = avifile(avipath,'Compression','Cinepak','quality',50);
% catch
    aviobj = avifile(avipath,'Compression','none','quality',50);
% end


%%
 load([fld ,'/' fdates{1,2}], '-mat')
mainFIG = figure('name',data.KFileName,'numbertitle','off','position',[ 200    30   600   600],'color','white');
    
ax(1) = axes('position',[0.05,0.05,0.40,0.40]);
ax(2) = axes('position',[0.55,0.05,0.40,0.40]);
ax(3) = axes('position',[0.05,0.55,0.40,0.40]);
ax(4) = axes('position',[0.55,0.55,0.40,0.40]);

    % come up with a nice colormap
    d1=colormap('jet');
    d2=colormap('bone');
    d1=flipud(d1);
    d1(1:end/2+1,:)=d2(end/2:end,:);
    
    

for i=1:1:length(fdates)
    
    
    
    load([fld ,'/' fdates{i,2}], '-mat','rho','t','SI','SJ','vx','vy')
    
    U_anal = (grav_accel_l)/(2*vw_l)*(((SJ-1)/2)^2-(-(SJ-1)/2+((1:SJ)-1)).^2); 

           vxM   = reshape(vx,[SI,SJ]);
           vyM   = reshape(vy,[SI,SJ]);
           v1    = sqrt(vxM.^2+vyM.^2);
        
           UxProfile=v1(Indmid);
%            tend =20000;
%            Umax_LB=zeros(tend,1);
%            t1=1:1:tend;
            Umax_LB(t)=max(UxProfile);

                disp(['t = ' num2str(t)])

                raoM = reshape(rho,[SI,SJ]);


                subplot(ax(1))
                cla
                imagesc(raoM), daspect([1,1,1])
                title('Density')
                colormap(d1)
                hold on
                hpl = plot(xblack, yblack,'y.'); set(hpl,'markersize',3)
                colorbar


                subplot(ax(2))
                cla
                imagesc(v1), daspect([1,1,1])
                title('vy')
                hold on
                quiver(vxM, vyM,'color','red');

                hpl = plot(xblack, yblack,'y.'); set(hpl,'markersize',3)
                colorbar

                subplot(ax(3))
                plot(xmid,UxProfile,'b-',xmid,U_anal,'bo')
                legend('LBM','Anal',1)
                title('Velocity Profile')

                subplot(ax(4))
                plot(t1(1:t),Umax_LB(1:t),'b-')
                xlabel('time step, ts'); ylabel('max velocity, lu/ts')
                title('Transient to steady state')
                pause(0.5)
                   

    pause(0.01)
    frame = getframe(gcf);
    aviobj = addframe(aviobj,frame); 
   


end

aviobj = close(aviobj);
disp(['Movie saved to ' avipath])

