function f_run_D2Q9_simul_k_SC_SCMP_BB3_ns(data)
clc

Ncoors              = data.Ncoors;
rhow_l               = data.rhow_l;
rhoa_l                = data.rhoa_l;
omega_w           =  data.omega_w;
omega_a            =  data.omega_a;
im                     = data.im;
e                       = data.e;
ao                     = data.ao;
%Ix_e                  = data. Ix_e;
IndNorth = data.IndNorth;
IndSouth = data.IndSouth;
IndWest = data.IndWest;
IndEast = data.IndEast;
Indmid = data.Indmid;             
I_xnyn =  data. I_xnyn; 
xmid = data.xmid;

IndNeigborSolid = data.IndNeigborSolid;

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
    
 %% Black - SOLID  nodes

        [yblack, xblack] = find(im==0);
        
        IndBlack = sub2ind([SI,SJ], yblack, xblack);   
%% Define densities and velocities at the fluid locations

   rho   = zeros(Ncoors, 1); % densities
    vx    = zeros(Ncoors, 1);
    vy    = zeros(Ncoors, 1);

 rho(FlTag==1,1) = rhow_l;% fluid 1 =  water
 rho(FlTag==2,1) = rhoa_l;% fluid 2  = air   
 rho(FlTag==0,1) = rhoa_l;% Solid area = air  
 
        
%% Initial equilibrium distribution functions

        Fieq=zeros(Ncoors,9);
        
      
       for a = 1:9
            
          Fieq(:,a) = WF(a)*rho(:).* (1   + 3   *  (e(a,1) *  vx(:)  + e(a,2) * vy(:)) ...
                                            + 9/2 *  (e(a,1) *  vx(:)  + e(a,2) * vy(:)).^2 ...
                                            - 3/2  *  (vx(:).^2 + vy(:).^2));
          
                    
        end

        Fi = Fieq ;     %at t=0   

%% Pearmeable Solids

ns =zeros(Ncoors,1);
ns(FlTag==0,1) =1; % solid area permeability: 1= completely solid, bounce back, 0 = no solid.

% %% Gravity Force
% Fg_x = 0;
% 
% Fg_y = zeros(Ncoors,1);
% Fg_y = rhow_l*grav_accel_l;
% Fg_y(FlTag==0,1) =0;
% 

%% --------- BEGINNING OF THE MAIN LOOP-------
t = 0; 
tend =20000;
Umax_LB=zeros(tend,1);

while t<tend

    t = t+1;
    
%% Interparticle Forces

        rho_o =200; 
        Fsai_o=4;
        Gs =-120;
        Fsai = Fsai_o*exp(-rho_o./rho);
        
        Fint_ex = 0;
        Fint_ey = 0;
        
        for a=1:8
                           
           Fint_ex = Fint_ex + WF(a) * Fsai(I_xnyn(:,a)) * e(a,1);
           Fint_ey = Fint_ey + WF(a) * Fsai(I_xnyn(:,a)) * e(a,2);
                              

        end
        
        
            Fint_ex = -Gs*Fsai.*Fint_ex;
            Fint_ey = -Gs*Fsai.*Fint_ey;
            
%% Surface Forces

        Ga = ones(Ncoors,1)*(0);
        Ga(FlTag==0,1) =0;
       
        Fadh_ex_dum = 0;
        Fadh_ey_dum = 0;
        
        IfNeigSolid = FlTag(I_xnyn)==0; % column mtx with 1 and 0 s
        
        for a=1:8
                           
           Fadh_ex_dum = Fadh_ex_dum + IfNeigSolid(:,a)*WF(a) * e(a,1);
           Fadh_ey_dum = Fadh_ey_dum + IfNeigSolid(:,a)*WF(a) * e(a,2);
                              

        end
        
        
            Fadh_ex = -Ga.*Fsai.*Fadh_ex_dum;
            Fadh_ey = -Ga.*Fsai.*Fadh_ey_dum;
            
%% Propogate and calculate new distribution function
        Fi_next = Fi*0;
       
        %BodyF = Fi(:,1)*0;
    
              for a = 1:9
              
                    BodyF = (-WF(a).*rho.*(e(a,2).*grav_accel_l))./co^2; 
                    Fi_next     (I_xnyn(:,a), a)    = Fi(:, a)- BodyF(:);%    - omega_w.*( Fi(:, a)    -   Fieq(:, a) ) - BodyF;
              end      
                    Fi = Fi_next;
                                        
                    Fi(IndWest,5) = Fi(IndWest,7) - (1/2)*(Fi(IndWest,2) - Fi(IndWest,4));
                    Fi(IndWest,1) = Fi(IndWest,3);
                    Fi(IndWest,8) = Fi(IndWest,6) + (1/2)*(Fi(IndWest,2) - Fi(IndWest,4));


                    Fi(IndEast,6)=Fi(IndEast,8) - (1/2)*(Fi(IndEast,2) - Fi(IndEast,4));
                    Fi(IndEast,3)=Fi(IndEast,1);
                    Fi(IndEast,7)=Fi(IndEast,5) + (1/2)*(Fi(IndEast,2) - Fi(IndEast,4));
                    
                    rho  = sum(Fi,2);
                    vx      = (Fi*e(:,1))./rho + (1/omega_w)*Fint_ex./rho + (1/omega_w)*Fadh_ex./rho; 
                    vy      = (Fi*e(:,2))./rho + (1/omega_w)*Fint_ey./rho + (1/omega_w)*Fadh_ey./rho; 
                    
             for a = 1:9
                 
                    Fieq(:,a) = WF(a)*rho.* (1   + 3     *  (e(a,1) *  vx  + e(a,2) * vy) ...
                                                        + 9/2 *  (e(a,1) *  vx  + e(a,2) * vy).^2 ...
                                                        - 3/2  *  (vx.^2 + vy.^2));
                                                    
                    Fi(:, a)    = Fi(:, a) - omega_w.*( Fi(:, a)    -   Fieq(:, a) ); 
                   % Fi(IndNeigborSolid(:,a), ao(a))    = Fi(IndNeigborSolid(:,a), a);   
                    Fi(:, a)    = Fi(:, a) + ns.*(Fi(I_xnyn(:,a), ao(a)) - Fi(:, a)) ;     
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
    
    
    
    
    
        