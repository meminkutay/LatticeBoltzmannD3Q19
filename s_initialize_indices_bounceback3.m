            im = data.im;
            e = data.e;


            %--define coordinates
            SI=size(im,1);
            SJ=size(im,2);
            [x,y]  = meshgrid(1:SJ,1:SI);
            x      = reshape(x, [SI*SJ,1]);
            y      = reshape(y, [SI*SJ,1]);


            Ncoors = length(x);

           data.Ncoors              = Ncoors;
           data.x                       =x;
           data.y                       =y;
           data.SI = SI;
           data.SJ = SJ;


                disp('Determining neighbor indices...')
                pause(0.1)
                tic;
   
               W = data.W;
                H = data.H;
                Ineigbor = zeros(Ncoors,9);

                                    for a = 1:9

                                    xn_e= x - e(a,1);
                                    yn_e = y - e(a,2);

                                    % periodic boundary conditions
                                    xn_e(xn_e>SJ) = 1;
                                    xn_e(xn_e==0) = SJ;
                              
                                    yn_e(yn_e>SI) = 1;
                                    yn_e(yn_e==0) = SI;
                           
                                    
                                    Ix_es = sub2ind([SI,SJ], yn_e, xn_e);
                                    Ix_e(:,a) = Ix_es; %#ok<SAGROW>
                                    

                                     
                                         
                                    end%for a = 1:9
                     data.Ix_e =  Ix_e;
                        
%Boundary nodes
        yNorth = ones(SJ,1); 
        xNorth = (1:SJ)';
        IndNorth = sub2ind([SI,SJ],yNorth,xNorth);
        
        xSouth=xNorth;
        ySouth = zeros(SJ,1)+SI;
        IndSouth = sub2ind([SI,SJ],ySouth,xSouth);
        
        yWest= y(1:SI);
        xWest= ones(SI,1);
        IndWest=sub2ind([SI,SJ],yWest,xWest);
        
        yEast= y(1:SI);
        xEast= x((SI*SJ-SI+1):SI*SJ);
        IndEast=sub2ind([SI,SJ],yEast,xEast);

% mid nodes
        ymid = ones(SJ,1) * round(SI/2);
        xmid = [1:SJ]';
        Indmid = sub2ind([SI,SJ],ymid,xmid);
        
        data.IndNorth = IndNorth;
        data.IndSouth = IndSouth;
        data.IndWest = IndWest;
        data.IndEast = IndEast;
        data.Indmid = Indmid;
      %  data.I_xnyn = I_xnyn;
        data.xmid = xmid;
      %  data.IndNeigborSolid = IndNeigborSolid;
        
                uiwait( msgbox('INDEXING DONE....INDEXING DONE....INDEXING DONE....INDEXING DONE....INDEXING DONE....INDEXING DONE....'))

                 toc

                        

                data.progdata{3} =  'COMPLETE';
                data.progdata{4} = '------';

                set(data.tbl_Progress, 'data',data.progdata);
                imshow(data.im), axis on,    
                
                hp = impixelinfo;
                haxpos = get(data.hax,'position');
                set(hp,'Position',[W-300, 0,100,10]);

                [SI,SJ]=size(im);
                text(0.05*SJ, 0.2*SI,'255=W, 151=A, 0=S','color','blue')