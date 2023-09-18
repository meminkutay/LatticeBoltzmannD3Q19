function f_run_D2Q9_simul_Korner3(data)
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

        Fi_out = Fieq ;     %at t=0        
      
%% INPUT Gad K_l

ns =zeros(Ncoors,1);
ns(FlTag==0,1) = 0.1 % solid area permeability: 1= completely solid, bounce back, 0 = no solid.


%% --------- BEGINNING OF THE MAIN LOOP-------
t = 0;  
tend = 2000000;

while t< tend

    t = t+1;
    
  
%% Propogate and calculate new distribution function
        Fi_in = Fi_out*0;
        
           %streaming
           for a = 1:9
               Fi_in (:, a)    =  Fi_out(Ix_e(:,a), a);
           
           end
           
                    Fi_in(IndWest,5) = Fi_in(IndWest,7) - (1/2)*(Fi_in(IndWest,2) - Fi_in(IndWest,4));
                    Fi_in(IndWest,1) = Fi_in(IndWest,3);
                    Fi_in(IndWest,8) = Fi_in(IndWest,6) + (1/2)*(Fi_in(IndWest,2) - Fi_in(IndWest,4));


                    Fi_in(IndEast,6)=Fi_in(IndEast,8) - (1/2)*(Fi_in(IndEast,2) - Fi_in(IndEast,4));
                    Fi_in(IndEast,3)=Fi_in(IndEast,1);
                    Fi_in(IndEast,7)=Fi_in(IndEast,5) + (1/2)*(Fi_in(IndEast,2) - Fi_in(IndEast,4));
                    
            % Calculate new density, velocity, and alpha 

            rho  = sum(Fi_in,2);
            vx   = ((Fi_in)*e(:,1))./rho ;
            vy   = ((Fi_in)*e(:,2))./rho ;
 
           % calculate new eqm distribution functions
            Fieq=zeros(Ncoors,a);
            
            for a = 1:9


                Fieq(:,a) = WF(a)*rho.* (1   + 3     *  (e(a,1) *  vx  + e(a,2) * vy) ...
                                                            + 9/2 *  (e(a,1) *  vx  + e(a,2) * vy).^2 ...
                                                            - 3/2  *  (vx.^2 + vy.^2));                                                         


            end

        

              for a = 1:9
                         
                    BodyF = (-WF(a).*rhow_l.*(e(a,2).*grav_accel_l))./co^2; 
                    
                     
                    Fi_out (:, a)    =  Fi_in(:, a)  - omega_w.*( Fi_in(:, a)  -   Fieq(:, a) ) - BodyF ;
                    
                    % permeability of solid regions -  used the Method - 3
                    % in S.D.C Walsh et al. (2009), Computers and Geosiences Vol 35,
                    % pp1186-1193.
                    Fi_out (:, a)   =  (1-ns).* Fi_out (:, a) + ns.* Fi_in( :, ao(a)) ;
                                   
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
% if t/ data.tsavestep==round(t/ data.tsavestep)
%     save([data.FilePath '/' data.KFileName(1:end-4) 't' num2str(t) '.mat'])
% end
        
    
    
end

    
disp('Done')    
    
    
    
    
    
        