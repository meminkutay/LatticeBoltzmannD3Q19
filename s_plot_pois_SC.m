if ~exist('rhoM','var')

    mainFIG = figure('name',data.KFileName,'numbertitle','off','position',[  50         50        1200         600],'color','white');

    ax(1) = axes('position',[0.05,0.05,0.45,0.40]);
    ax(2) = axes('position',[0.55,0.05,0.40,0.40]);
    ax(3) = axes('position',[0.05,0.55,0.40,0.40]);
    ax(4) = axes('position',[0.55,0.55,0.40,0.40]);

    % come up with a nice colormap
    d1=colormap('jet');
    d2=colormap('bone');
    d1=flipud(d1);
    d1(1:end/2+1,:)=d2(end/2:end,:);
    
    
end

 U_anal = (grav_accel_l)/(2*vw_l)*(((SJ-1)/2)^2-(-(SJ-1)/2+((1:SJ)-1)).^2); 

%            vxM   = reshape(vx,[SI,SJ]);
%            vyM   = reshape(vy,[SI,SJ]);
%            v1    = sqrt(vxM.^2+vyM.^2);
%         
%            UxProfile=v1(Indmid);
% 
%            t1=1:1:tend;
%            Umax_LB(t)=max(UxProfile);
           
                    disp(['t = ' num2str(t)])

                    rhoM = reshape(rho,[SI,SJ]);

                    subplot(ax(1))
                    cla
                    imagesc(rhoM), daspect([1,1,1]), , axis ij
                    hold on
%                     quiver(vxM, vyM,'color','red');
                    colormap(d1)                   
                    colorbar
                    hpl = plot(xblack, yblack,'y.'); set(hpl,'markersize',3)

                    subplot(ax(2))
                    plot(xmid,UxProfile,'b-',xmid,U_anal,'bo')
                    legend('LBM','Anal',1)
                    title('Velocity Profile')

                    subplot(ax(3))
                    cla
                    imagesc(vyM), daspect([1,1,1])
                    title('vy')
                    hold on
                    quiver(vxM, vyM,'color','red');
                    colorbar
                    hpl = plot(xblack, yblack,'y.'); set(hpl,'markersize',3)
                    
                    subplot(ax(4))
                    plot(t1(1:t),Umax_LB(1:t),'b-')
                    xlabel('time step, ts'); ylabel('max velocity, lu/ts')
                    title('Transient to steady state')
                    pause(0.5)
                    

                   
                    
