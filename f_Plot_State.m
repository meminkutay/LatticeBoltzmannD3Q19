function f_Plot_State


[fname,fld] = uigetfile2('pick a mat file')



mainFIG = figure('name',fname,'numbertitle','off','position',[ 200    30   600   600],'color','white');
    
ax(1) = axes('position',[0.05,0.05,0.40,0.40]);
ax(2) = axes('position',[0.55,0.05,0.40,0.40]);
ax(3) = axes('position',[0.05,0.55,0.40,0.40]);
ax(4) = axes('position',[0.55,0.55,0.40,0.40]);

    % come up with a nice colormap
    d1=colormap('jet');
    d2=colormap('bone');
    d1=flipud(d1);
    d1(1:end/2+1,:)=d2(end/2:end,:);
    
    
    
    load([fld ,'/' fname], '-mat','IndR','rho','waterarea','t','SI','SJ','vx','vy','xblack','yblack','data')
    
         vyM   = reshape(vy,[SI,SJ]);
                    %         v1    = sqrt(vx1.^2+vy1.^2);

                   % UyProfile=vy1(Indmid);

                    disp(['t = ' num2str(t)])

                    raoM = reshape(rho,[SI,SJ]);

                    IndRM  =  reshape(IndR,[SI,SJ]);


                    subplot(ax(1))
                    cla
                    hpl = plot(xblack, yblack,'y.'); set(hpl,'markersize',5)
                    title(['Color indices, t=' num2str(t*data.dt) 's'])
                    hold on
                              
                    hs = surf(IndRM); shading interp,  view(180,-89), axis ij, axis off
                    set(hs,'facealpha',1)
                    %imagesc(1-IndRM), 
                    daspect([1,1,1]), 
%                     colormap hot
                   colormap(d1)
          
%                     quiver(vx1, vy1,'color','black');
                 
                    colorbar

                    subplot(ax(2))
                    cla
                    imagesc(raoM), daspect([1,1,1])
                    title('Density')
                    colorbar

                    subplot(ax(3))
                       cla
                        imagesc(vyM), daspect([1,1,1])
                        title('vy')
                        colorbar
                   

                    subplot(ax(4))
                    plot(1:t, waterarea(1:t),'b.-')
                    title('Water area')
                    pause(0.05)


