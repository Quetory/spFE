function animate_mode(elements,coordinates, shape,NDOF)
Nperiods = 0.5;
freq = 0.025;
sc   = ones(3,1);
gain = 0.2;

if size(coordinates,2)==2    %2D objects 
    if size(elements,2)==3   %triangle         
        X=reshape(coordinates(elements',1),size(elements,2),size(elements,1));
        Y=reshape(coordinates(elements',2),size(elements,2),size(elements,1));   
        fill(X,Y,[0.3 0.3 0.9]);
    end
   
elseif size(coordinates,2)==3   %3D object
    if size(elements,2)==8  || size(elements,2)==20   
     
       faces1=elements(:,1:4);
       faces2=elements(:,5:8);
       faces3=elements(:,[1 4 8 5]);
       faces4=elements(:,[2 3 7 6]);
       faces5=elements(:,[1 2 6 5]);
       faces6=elements(:,[3 4 8 7]);
       
       faces=[faces1; faces2; faces3; faces4; faces5; faces6];

       figure
        if ~isreal(shape)
                vmag = abs(shape);
                vph  = angle(shape);
                shape = vmag.*sin(vph);
        end

            
       vanim_o = reshape(shape,NDOF,size(coordinates,1)).';
       set(gcf,'position',[680   208   851   770]);
       if NDOF==3 
           [mval,midx]=max(max(abs(vanim_o)));
           vanim = vanim_o/mval*gain;
           alimits = [min(coordinates-abs(vanim)).' max(coordinates+abs(vanim)).'];
           sc(midx) = 1.1;
           alimits = [1./sc sc].*alimits;
           colormap jet
           
           Clim = sqrt(sum(vanim_o.^2,2));
           for t = 1: ceil(Nperiods/freq)
               clf
               mode = coordinates + vanim*cos(2*pi*freq*t);
               Col = vanim_o*cos(2*pi*freq*t);
               C = sqrt(sum(Col.^2,2));
               X=reshape(mode(faces',1),size(faces,2),size(faces,1));
               Y=reshape(mode(faces',2),size(faces,2),size(faces,1)); 
               Z=reshape(mode(faces',3),size(faces,2),size(faces,1)); 
               Col = reshape(C(faces'),size(faces,2),size(faces,1));
               patch(X,Y,Z,Col);%,'EdgeColor','interp' );%[0.3 0.8 0.9]);
               axis equal
               xlim(alimits(1,:));
               ylim(alimits(2,:));
               zlim(alimits(3,:));
               caxis([min(Clim) max(Clim)]);
               colorbar
               view(3)
               
               grid on
               drawnow;
           end
       elseif NDOF==1
            alimits = [1.05*min(coordinates).' 1.05*max(coordinates).'];
            if ~isreal(vanim_o)
                vmag = abs(vanim_o);
                vph  = angle(vanim_o);
                [s,l]=bounds(vmag.*sin(vph));
                vanim = vmag.*sin(vph);
            else                
                [s,l]=bounds(vanim_o);
                vanim = vanim_o;
            end
            
            X=reshape(coordinates(faces',1),size(faces,2),size(faces,1));
            Y=reshape(coordinates(faces',2),size(faces,2),size(faces,1)); 
            Z = reshape(coordinates(faces',3),size(faces,2),size(faces,1));
            P = zeros(size(Z));
            h = patch(X,Y,Z,P,'EdgeColor','interp' );
%             set(gcf,'position',[680   208   851   770]);
            axis equal
            xlim(alimits(1,:));
            ylim(alimits(2,:));
            zlim(alimits(3,:));
            caxis([s l]);   
            colormap jet
            colorbar
%             view(3)
           view(35,-24)
            for t = 1: ceil(Nperiods/freq)
               mode = vanim*sin(2*pi*freq*t);
               P = reshape(mode(faces'),size(faces,2),size(faces,1));
       
               h.CData=P;
               drawnow

            end
           
       end
           
    end
    
end


