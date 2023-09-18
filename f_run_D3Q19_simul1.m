function f_run_D3Q19_simul(data)
clc

rhow_l               = data.rhow_l;
omega_w           =  data.omega_w;
im                     = data.im3D;
e                       = data.e3D;
ao                     = data.ao3D;
Ns = data.Ns;
Ix_e                  = data. Ix_e3D;
dx = data.dx;
dt = data.dt;
 
IndNorth = data.IndNorth;
IndSouth = data.IndSouth;
IndWest = data.IndWest;
IndEast = data.IndEast;
IndZmid = data.IndZmid;  
Indmid = data.Indmid;   

Indtop = data.Indtop;
Indbot = data.Indbot;
              
WF                     = data.WF3D;
co                      = data.co;



vw_l                    = data.vw_l;
%va_l                    = data.va_l;
grav_accel_l       = data.grav_accel_l;
SI                      = data.SI;
SJ                      = data.SJ;
SK                      = data.SK;

Ncoors              = SI*SJ*SK; 
x = data.x ;
y = data.y ;
z = data.z ;

    
%% Define densities and velocities at the fluid locations

   rho   = rhow_l+  zeros(Ncoors, 1); % densities
    vx    = zeros(Ncoors, 1);
    vy    = zeros(Ncoors, 1);
    vz    = zeros(Ncoors, 1);

 
%% Black - SOLID  nodes (used in plot fuction too!)

        IndBlack = find(im==0);
        
        [yblack, xblack, zblack] = ind2sub([SI,SJ, SK], IndBlack); 
 
%% Initial equilibrium distribution functions

        Fieq=zeros(Ncoors,19);
      
        for a = 1:19
                 
            Fieq(:,a) = WF(a)*rho.* (1   + 3   *  (e(a,1) *  vx  + e(a,2) * vy + e(a,3) * vz) ...
                                            + 9/2 *  (e(a,1) *  vx + e(a,2) * vy +  e(a,3) * vz).^2 ...
                                            - 3/2  *  (vx.^2 + vy.^2 + vz.^2));
      
        end

        Fi_out = Fieq ;     %at t=0        
      
%% INPUT Gad K_l

ns =zeros(Ncoors,1);
ns(im==0,1) = Ns; % solid area permeability: 1= completely solid, bounce back, 0 = no solid.
ns(im==255,1) = 0;
ns(im==151,1) = Ns;

%% --------- BEGINNING OF THE MAIN LOOP-------
t = 0;  
tend = 100000;

pin = 1.002;
pout = 1.0;
c1 = 1/4;
c2 = 1/6;

kzzg = zeros(tend,1);

while t< tend

    t = t+1;
    
  
%% Propogate and calculate new distribution function
        Fi_in = Fi_out*0;
        
           %streaming
           for a = 1:19
               Fi_in (:, a)    =  Fi_out(Ix_e(:,a), a);
           
           end
           
                    Fi_in(IndNorth,4) = Fi_in(IndNorth,3); 
                    Fi_in(IndNorth,9) = Fi_in(IndNorth,7);
                    Fi_in(IndNorth,10) = Fi_in(IndNorth,8); 
                    Fi_in(IndNorth,17) = Fi_in(IndNorth,15);
                    Fi_in(IndNorth,18) = Fi_in(IndNorth,16); 
                    
                    Fi_in(IndSouth,3) = Fi_in(IndSouth,4); 
                    Fi_in(IndSouth,7) = Fi_in(IndSouth,9);
                    Fi_in(IndSouth,8) = Fi_in(IndSouth,10); 
                    Fi_in(IndSouth,15) = Fi_in(IndSouth,17);
                    Fi_in(IndSouth,16) = Fi_in(IndSouth,18); 
           
           
                    Fi_in(IndWest,2) = Fi_in(IndWest,1); %- (1/2)*(Fi_in(IndWest,2) - Fi_in(IndWest,4));
                    Fi_in(IndWest,8) = Fi_in(IndWest,10);
                    Fi_in(IndWest,9) = Fi_in(IndWest,7); 
                    Fi_in(IndWest,13) = Fi_in(IndWest,11);
                    Fi_in(IndWest,14) = Fi_in(IndWest,12); 


                    Fi_in(IndEast,1)=Fi_in(IndEast,2);% - (1/2)*(Fi_in(IndEast,2) - Fi_in(IndEast,4));
                    Fi_in(IndEast,7)=Fi_in(IndEast,9);
                    Fi_in(IndEast,10)=Fi_in(IndEast,8);
                    Fi_in(IndEast,11)=Fi_in(IndEast,13);
                    Fi_in(IndEast,12)=Fi_in(IndEast,14);

% %--apply boundary conditions at inlet nodes
%     rho(Indtop)=pin;
%     
% 	vx(Indtop)=0; 
%     vy(Indtop)=0;
% 	vz(Indtop)= 1-(sum(Fi_in(Indtop,1:4),2)+sum(Fi_in(Indtop,7:10),2)+Fi_in(Indtop,19)+ 2*(Fi_in(Indtop,6)+sum(Fi_in(Indtop,12:13),2)+sum(Fi_in(Indtop,16:17),2)))./pin;
%     
%     Fi_in(Indtop,5)=   Fi_in(Indtop,6) + 1/3*pin*vz(Indtop);
%     Fi_in(Indtop,15)=  Fi_in(Indtop,17) - c1*(Fi_in(Indtop,3) - Fi_in(Indtop,4)) + c2*pin*vz(Indtop);
%     Fi_in(Indtop,18)=  Fi_in(Indtop,16) + c1*(Fi_in(Indtop,3) - Fi_in(Indtop,4)) + c2*pin*vz(Indtop);
%     Fi_in(Indtop,11)=  Fi_in(Indtop,13) - c1*(Fi_in(Indtop,1) - Fi_in(Indtop,2)) + c2*pin*vz(Indtop);
%     Fi_in(Indtop,14)=  Fi_in(Indtop,12) + c1*(Fi_in(Indtop,1) -Fi_in(Indtop,2)) + c2*pin*vz(Indtop);   
%  
    
% %--apply boundary conditions at outlet nodes
% 	rho(Indbot)=pout;
%     
% 	vx(Indbot)=0; 
%     vy(Indbot)=0; 
%  	vz(Indbot)= -1 + (sum(Fi_in(Indbot,1:4),2)+sum(Fi_in(Indbot,7:10),2)+Fi_in(Indbot,19)+2*(Fi_in(Indbot,5)+sum(Fi_in(Indbot,14:15),2)+Fi_in(Indbot,11)+Fi_in(Indbot,18)))./pout;
%    
%     Fi_in(Indbot,6)=   Fi_in(Indbot,5)-1/3*pout*vz(Indbot);
%     Fi_in(Indbot,17)=  Fi_in(Indbot,15) + c1*(Fi_in(Indbot,3)-Fi_in(Indbot,4))- c2*pout*vz(Indbot);
%     Fi_in(Indbot,16)=  Fi_in(Indbot,18) - c1*(Fi_in(Indbot,3)-Fi_in(Indbot,4))-c2*pout*vz(Indbot);
%     Fi_in(Indbot,13)=  Fi_in(Indbot,11) + c1*(Fi_in(Indbot,1)-Fi_in(Indbot,2))-c2*pout*vz(Indbot);
%     Fi_in(Indbot,12)=  Fi_in(Indbot,14) - c1*(Fi_in(Indbot,1)-Fi_in(Indbot,2))-c2*pout*vz(Indbot);  
%    
                    
            % Calculate new density, velocity, and alpha 

            rho  = sum(Fi_in,2);
            vx   = (((Fi_in)*e(:,1)).*(1-ns))./rho ;
            vy   = (((Fi_in)*e(:,2)).*(1-ns))./rho ;
            vz   = (((Fi_in)*e(:,3)).*(1-ns))./rho ;
            
            
  
 
           % calculate new eqm distribution functions

            
            for a = 1:19
                
              Fieq(:,a) = WF(a)*rho.* (1   + 3   *  (e(a,1) *  vx  + e(a,2) * vy + e(a,3) * vz) ...
                                            + 9/2 *  (e(a,1) *  vx + e(a,2) * vy +  e(a,3) * vz).^2 ...
                                            - 3/2  *  (vx.^2 + vy.^2 + vz.^2));

            end

        

              for a = 1:19
                         
                    BodyF = (-WF(a).*rhow_l.*(e(a,3).*grav_accel_l))./co^2; 
                    
                     
                    Fi_out (:, a)    =  Fi_in(:, a)  - omega_w.*( Fi_in(:, a)  -   Fieq(:, a) ) - BodyF ;
                    
                    % permeability of solid regions -  used the Method - 3
                    % in S.D.C Walsh et al. (2009), Computers and Geosiences Vol 35,
                    % pp1186-1193.
                    Fi_out (:, a)   =  (1-ns).* Fi_out (:, a) + ns.* Fi_in( :, ao(a)) ;
                                   
              end
              
    
 %%  Plot the velocity every n time steps
 
           vxM   = reshape(vx,[SI,SJ, SK]);
           vyM   = reshape(vy,[SI,SJ, SK]);
           vzM   = reshape(vz,[SI,SJ, SK]);
           v1    = sqrt(vxM.^2+vyM.^2+vzM.^2);
           vmax = max(v1);
         
           UxProfile=v1(Indmid);
           Uzz = v1(IndZmid);
           
           t1=1:1:tend;
           Umax_LB(t)=max(max(vmax));
           
           
            V1=sum(sum(Uzz));
            A = (SI-1)*(SJ-1);
            delP = (pin - pout)*co^2/(SK-1);
            gamma = grav_accel_l*rhow_l;
            
             kzzg(t)=V1/A;
             

        if round(t/data.tplotstep) == t/data.tplotstep
            
           s_plot_D3Q19
           kzzp=V1*gamma/(A*delP);
           kzzg1=V1/A
           kzzg_w=kzzg1*dx/dt
  
        end


       Cumkzzg = cumsum(kzzg);
       
%% save all data at every ... cycles
if t/ data.tsavestep==round(t/ data.tsavestep)
    save([data.FilePath '/' 't' num2str(t) '.mat'])
end
        
    
    
end

    
disp('Done')    
    
    
    
    
    
        