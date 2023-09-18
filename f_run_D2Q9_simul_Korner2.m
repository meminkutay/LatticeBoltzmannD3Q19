function f_run_D2Q9_simul_Korner2(data)
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
Indmid = data.Indmid;     
xmid = data.xmid;  

Ns = data.Ns;
IndNorth = data.IndNorth;
IndSouth = data.IndSouth;
IndWest = data.IndWest;
IndEast = data.IndEast;

             
                 
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
 
%% Black - SOLID  nodes (used in plot fuction too!)

        [yblack, xblack] = find(im==0);
        
        IndBlack = sub2ind([SI,SJ], yblack, xblack); 
 
%% Initial equilibrium distribution functions

        Fieq=zeros(Ncoors,9);
      
        for a = 1:9
                 
            Fieq(:,a) = WF(a)*rho(:).* (1   + 3   *  (e(a,1) *  vx(:)  + e(a,2) * vy(:)) ...
                                            + 9/2 *  (e(a,1) *  vx(:)  + e(a,2) * vy(:)).^2 ...
                                            - 3/2  *  (vx(:).^2 + vy(:).^2));
                                            
        
        end

        Fi = Fieq ;     %at t=0        
      
%% INPUT Gad K_l

ns =zeros(Ncoors,1);
ns(FlTag==0,1) = Ns; % solid area permeability: 1= completely solid, bounce back, 0 = no solid.


%% --------- BEGINNING OF THE MAIN LOOP-------
t = 0;  
tend = 2000000;

while t< tend

    t = t+1;
    
  
%% Propogate and calculate new distribution function
        Fi_x = Fi*0;
        
           
              for a = 1:9
                  
                    % applying body force only the regions (lattices) occupied by the
                    % water and water/air interface. 
                    BodyF = (-WF(a).*rhow_l.*(e(a,2).*grav_accel_l))./co^2; 
                    
                     
                    Fi_x (:, a)    =  Fi(Ix_e(:,a), a)  - omega_w.*( Fi(Ix_e(:,a), a)  -   Fieq(Ix_e(:,a), a) ) - BodyF ;
                    
                    % permeability of solid regions -  used the Method - 3
                    % in S.D.C Walsh et al. (2009), Computers and Geosiences Vol 35,
                    % pp1186-1193.
                    Fi_x (:, a)   =  (1-ns).* Fi_x (:, a) + ns.* Fi( Ix_e(:,a), ao(a)) ;
                                   
                   
              end
                    Fi = Fi_x;
               
                    Fi(IndWest,5) = Fi(IndWest,7) - (1/2)*(Fi(IndWest,2) - Fi(IndWest,4));
                    Fi(IndWest,1) = Fi(IndWest,3);
                    Fi(IndWest,8) = Fi(IndWest,6) + (1/2)*(Fi(IndWest,2) - Fi(IndWest,4));


                    Fi(IndEast,6)=Fi(IndEast,8) - (1/2)*(Fi(IndEast,2) - Fi(IndEast,4));
                    Fi(IndEast,3)=Fi(IndEast,1);
                    Fi(IndEast,7)=Fi(IndEast,5) + (1/2)*(Fi(IndEast,2) - Fi(IndEast,4));
              
%% Calculate new density, velocity, and alpha 

        rho  = sum(Fi,2);
%         vx      = ((Fi + Fi_x)/2*e(:,1))./rho ;
%         vy      = ((Fi + Fi_x)/2*e(:,2))./rho ;
        vx   = ((Fi)*e(:,1))./rho ;
        vy   = ((Fi)*e(:,2))./rho ;
 
       
 
%% check maximum velocity
   % if max(abs(vx))>1 || max(abs(vy))>1
%      if max(abs(vx))>5 || max(abs(vy))>5   % by Saravana
%         uiwait(errordlg({['Maximum v_x= ', num2str(max(abs(vx))), ', Maximum v_y= ', num2str(max(abs(vy)))], '-----------------'...
%                        ' YOU MUST REDUCE THE TIME STEP'},'ERROR: TOO LARGE VELOCITY'))
%                    t=t-1;
%                    s_plot_stuff
%                    
%          return
%     end
    
        

%% Solid and other boundary conditions 

% %  Boundary conditions in the north, south, east and west
%  
% freq = 1/50;
% 
% vynorth_d =vynorth* sin(2*pi*freq*t);
% 
% [vx, vy, Fi, rho] = f_North_velocity_y_BC(Fi, vx, vy, rho, INorth, vynorth_d);
%     
       

%% calculate new eqm distribution functions
        Fieq=zeros(Ncoors,a);
       
        
        for a = 1:9
            
                        
            Fieq(:,a) = WF(a)*rho.* (1   + 3     *  (e(a,1) *  vx  + e(a,2) * vy) ...
                                                        + 9/2 *  (e(a,1) *  vx  + e(a,2) * vy).^2 ...
                                                        - 3/2  *  (vx.^2 + vy.^2));                                                         
                                        
                   
        end

        

 
        
                
 %%  Plot the velocity every n time steps
 
           vxM   = reshape(vx,[SI,SJ]);
           vyM   = reshape(vy,[SI,SJ]);
           v1    = sqrt(vxM.^2+vyM.^2);
        
           UxProfile=v1(Indmid);
            t1=1:1:tend;
            Umax_LB(t)=max(UxProfile);
           
        if round(t/data.tplotstep) == t/data.tplotstep
           s_plot_pois_SC
  
        end



%% save all data at every ... cycles
if t/ data.tsavestep==round(t/ data.tsavestep)
    save([data.FilePath '/' data.KFileName(1:end-4) 't' num2str(t) '.mat'])
end
        
    
    
end

    
disp('Done')    
    
    
    
    
    
        