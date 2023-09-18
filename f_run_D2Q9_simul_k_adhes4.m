function f_run_D2Q9_simul_k_adhes4(data)
clc

Ncoors              = data.Ncoors;
rhow_l               = data.rhow_l;
rhoa_l                = data.rhoa_l;
omega_w           =  data.omega_w;
omega_a            =  data.omega_a;
im                     = data.im;
e                       = data.e;
ao                     = data.ao;
Ix_e                  = data. Ix_e;

             
                 
WF                     = data.WF;
sigma_l             =data.sigma_l;
co                      = data.co;
omega_alpha     =  data.omega_alpha;




vw_l                    = data.vw_l;
% va_l                    = data.va_l;
grav_accel_l       = data.grav_accel_l;
SI                      = data.SI;
SJ                      = data.SJ;



%% tags of the fluids 0=solid, 1=fluid one, 2=fluid two
    FlTag  = zeros(Ncoors,1); 

    
% -Indices of the fluids in the image
    IndFL = [255, 151];
    
    for i=1:length(IndFL)
        FlTag(im==IndFL(i)) = i;
    end
    
    
%% Define densities and velocities at the fluid locations

   rho   = zeros(Ncoors, 1); % densities
    vx    = zeros(Ncoors, 1);
    vy    = zeros(Ncoors, 1);

 rho(FlTag==1,1) = rhow_l;% fluid 1 =  water
 rho(FlTag==2,1) = rhoa_l;% fluid 2  = air   
 rho(FlTag==0,1) = rhoa_l;% Solid area = air  
 
    
%% Color indices   

alpha(FlTag==1,1) = 1; % fluid 1 =  water
alpha(FlTag==2,1) = 0;   % fluid 2  = air   
alpha(FlTag==0,1) = 0;   % Solid area = air   

  
 %% Color indicator function     

IndR = zeros(Ncoors,1);
IndR(alpha>=1) = 1 ;
omega=omega_w.*IndR+omega_a.*(1-IndR);


% 
% %% mid nodes
%         ymid = ones(SJ,1) * round(SI/2);
%         xmid = [1:SJ]';
%         Indmid = sub2ind([SI,SJ],ymid,xmid);
%         

%% Black - SOLID  nodes

        [yblack, xblack] = find(im==0);
        
        IndBlack = sub2ind([SI,SJ], yblack, xblack);

        
% %        [yvoids, xvoids] = find(im~=0);
% %         
% %         IndVoids = sub2ind([SI,SJ], yvoids, xvoids);
        
%%  Initial Components of the surface tension term
        
          IndRM   = reshape(IndR,[SI,SJ]);
          [GxM, GyM] =gradient(IndRM);
        
          Gx   = reshape(GxM,[SI*SJ,1]);
          Gy   = reshape(GyM,[SI*SJ,1]);
           GM=sqrt(Gx.^2+Gy.^2);
                
           Inds = find(GM~=0); % to prevent NaN in Si    
        
        
        
%% Initial equilibrium distribution functions

        Fieq=zeros(Ncoors,9);
        Si = Fieq;
        alphai_eq  = Fieq;
        for a = 1:9
            
            Si(Inds,a) = + WF(a)*(9/2)*sigma_l.* (((Gx(Inds) + Gy(Inds)).^2./GM(Inds) ...
                                                                        -2*GM(Inds)).*( (e(a,1) + e(a,2)).^2 - 2/3));

            
            Fieq(:,a) = WF(a)*rho(:).* (1   + 3   *  (e(a,1) *  vx(:)  + e(a,2) * vy(:)) ...
                                            + 9/2 *  (e(a,1) *  vx(:)  + e(a,2) * vy(:)).^2 ...
                                            - 3/2  *  (vx(:).^2 + vy(:).^2)) ...
                                            + Si(:,a);
          
          alphai_eq(:,a)=(alpha./rho).*Fieq(:,a);   
            
        end

        Fi = Fieq ;     %at t=0        
        alphai=alphai_eq;    %at t=0



%% INPUT Gad K_l

ns =zeros(Ncoors,1);
ns(FlTag==0,1) = 1; % solid area permeability: 1= completely solid, bounce back, 0 = no solid.


%% --------- BEGINNING OF THE MAIN LOOP-------
t = 0;                
while t<200000000

    t = t+1;
    
  
%% Propogate and calculate new distribution function
        Fi_x = Fi*0;
        alphai_x = alphai*0;
        BodyF = alphai(:,1)*0;
    
              for a = 1:9
                  
                    % applying body force only the regions (lattices) occupied by the
                    % water and water/air interface. 
                    BodyF (IndR>0) = (-WF(a).*rhow_l.*(e(a,2).*grav_accel_l))./co^2; 
                    
                     
                    Fi_x        (:, a)    = Fi(Ix_e(:,a), a)        - omega.*( Fi(Ix_e(:,a), a)                   -   Fieq(Ix_e(:,a), a) ) - BodyF ;
                    alphai_x (:, a)    = alphai(Ix_e(:,a), a) - omega_alpha *( alphai(Ix_e(:,a), a)  - alphai_eq(Ix_e(:,a), a) );
                    
                    % permeability of solid regions -  used the Method - 3
                    % in S.D.C Walsh et al. (2009), Computers and Geosiences Vol 35,
                    % pp1186-1193.
                    Fi_x        (:, a)   =  (1-ns).* Fi_x (:, a) + ns.* Fi( Ix_e(:,a), ao(a)) ;
                    alphai_x        (:, a)   =  (1-ns).* alphai_x (:, a) + ns.* alphai( Ix_e(:,a), ao(a)) ;
                    
                    
                   
              end
              
%% Calculate new density, velocity, and alpha 

        rho  = sum(Fi_x,2);
        vx      = ((Fi + Fi_x)/2*e(:,1))./rho ;
        vy      = ((Fi + Fi_x)/2*e(:,2))./rho ;
        alpha        = sum(alphai_x ,2);
 
 
              Fi = Fi_x;
        alphai = alphai_x;


    
%% check maximum velocity
   % if max(abs(vx))>1 || max(abs(vy))>1
     if max(abs(vx))>5 || max(abs(vy))>5   % by Saravana
        uiwait(errordlg({['Maximum v_x= ', num2str(max(abs(vx))), ', Maximum v_y= ', num2str(max(abs(vy)))], '-----------------'...
                       ' YOU MUST REDUCE THE TIME STEP'},'ERROR: TOO LARGE VELOCITY'))
                   t=t-1;
                   s_plot_stuff
                   
         return
    end
    
        

%% Solid and other boundary conditions 

% %  Boundary conditions in the north, south, east and west
%  
% freq = 1/50;
% 
% vynorth_d =vynorth* sin(2*pi*freq*t);
% 
% [vx, vy, Fi, rho] = f_North_velocity_y_BC(Fi, vx, vy, rho, INorth, vynorth_d);
%     
%% New color indicator founction (IndR) and omega     
    
        IndR(alpha>=1) = 1 ;
        IndR(alpha<=0) = 0 ;
        Ialpha = find(alpha>0 & alpha<1);
        IndR(Ialpha)=(sin((alpha(Ialpha)-0.5).*pi)+1)./2;
        
        

        waterarea(t) = sum(sum(IndR));
        
        omega=omega_w.*IndR+omega_a.*(1-IndR);

          IndRM   = reshape(IndR,[SI,SJ]);
          [GxM, GyM] =gradient(IndRM);
        
          Gx   = reshape(GxM,[SI*SJ,1]);
          Gy   = reshape(GyM,[SI*SJ,1]);
     

        GM=sqrt(Gx.^2+Gy.^2);
        Inds = find(GM~=0); % to prevent NaN in Si 
        
       %find(isnan(GM)==1);
       
        

    %calculate new eqm distribution functions
        Fieq=zeros(Ncoors,a);
        Si = Fieq;
        alphai_eq = Fieq;
        
        for a = 1:9
            
            Si(Inds,a) =  WF(a)*(9/2)*sigma_l.* (((Gx(Inds) + Gy(Inds)).^2./GM(Inds) ...
                                                                        -2*GM(Inds)).*( (e(a,1) + e(a,2)).^2 - 2/3));

            
            Fieq(:,a) = WF(a)*rho.* (1   + 3     *  (e(a,1) *  vx  + e(a,2) * vy) ...
                                                        + 9/2 *  (e(a,1) *  vx  + e(a,2) * vy).^2 ...
                                                        - 3/2  *  (vx.^2 + vy.^2)) ...
                                                        + Si(:,a);  
                                        
          alphai_eq(:,a)=(alpha./rho).*Fi(:,a);   
            
       
          
        end

        

 
        
                
 %%  Plot the velocity every n time steps
 
        if round(t/data.tplotstep) == t/data.tplotstep
           s_plot_stuff
  
        end


%% save all data at every ... cycles
if t/ data.tsavestep==round(t/ data.tsavestep)
    save([data.FilePath '/' data.KFileName(1:end-4) 't' num2str(t) '.mat'])
end
        
    
    
end

    
disp('Done')    
    
    
    
    
    
        