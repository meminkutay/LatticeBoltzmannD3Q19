if ~exist('rhoM','var')
    
    mainFIG = figure('name','D3Q9 simul','numbertitle','off','position',[  50         50        1200         600],'color','white');
    
    ax(1) = axes('position',[0.05,0.05,0.45,0.40]);
    ax(2) = axes('position',[0.55,0.05,0.40,0.40]);
    ax(3) = axes('position',[0.05,0.55,0.40,0.40]);
    ax(4) = axes('position',[0.55,0.55,0.40,0.40]);
    
    % come up with a nice colormap
    d1=colormap('jet');
    d2=colormap('bone');
    d1=flipud(d1);
    d1(1:end/2+1,:)=d2(end/2:end,:);
    
    
    % these are here because I want to define them only once, not at every
    % loop.
    [SPx, SPy, SPz] = meshgrid(1:0.1*SJ:SJ, 1:0.1*SI:SI, [1:0.1*SK:SK]*0+1);
    V = smooth3(im);
    [X, Y, Z] = meshgrid(1:size(V ,2), 1:size(V ,1), 1:size(V ,3));
    
    
end


disp(['t = ' num2str(t)])
rhoM = reshape(rho,[SI,SJ,SK]);
%%
subplot(ax(1))
cla
imagesc(squeeze(rhoM(SI/2,:,:))'), daspect([1,1,1]), , axis ij
hold on
colormap(d1)
colorbar
hpl = plot3(xblack, yblack, zblack, 'y.'); set(hpl,'markersize',3)
title('Density --> pressure')

%%
subplot(ax(2))
plot(t1(1:t), kzzg(1:t),'b-')
xlabel('time step, ts'); ylabel('kzzg, lu/ts')
title('Transient to steady state')



%%
subplot(ax(3))

cla
imagesc(squeeze(vzM(SI/2,:,:))'), daspect([1,1,1])
title('Magnitude of the velocity')
hold on
%quiver(vxM, vzM, 'color','red');
colorbar
hpl = plot3(xblack, yblack, zblack, 'y.'); set(hpl,'markersize',3)


%%
subplot(ax(4))

cla
fcolor =[0.75    0.8472    0.75];
% fcolor =[rand rand rand];
hiso = patch(isosurface(X,Y,Z,V),...
    'FaceColor',fcolor,...
    'EdgeColor','none', ...
    'FaceAlpha',0.2);


camlight right
lighting phong
material shiny
axis equal
view(3)

lightangle(45,30);

set(gcf,'Renderer','zbuffer'); lighting phong
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)
set(gca,'xlim',[1,SJ],'ylim',[1,SJ],'zlim',[1,SK]);
daspect([1,1,1])

xlabel('X (pixels)')
ylabel('Y (pixels)')
zlabel('Z (pixels)')

hold on

streamline(vxM, vyM, vzM, SPx, SPy, SPz)


pause(0.05)




