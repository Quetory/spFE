function animate_mode(elements,coordinates, shape,NDOF)
Nperiods = 2;
freq = 0.025;
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
       vanim = reshape(shape,NDOF,size(coordinates,1)).';
      
       if NDOF==3 
           alimits = [min(coordinates-vanim).' max(coordinates+vanim).'];
           for t = 1: ceil(Nperiods/freq)
               clf
               mode = coordinates + vanim*sin(2*pi*freq*t);
               X=reshape(mode(faces',1),size(faces,2),size(faces,1));
               Y=reshape(mode(faces',2),size(faces,2),size(faces,1)); 
               Z=reshape(mode(faces',3),size(faces,2),size(faces,1)); 
               patch(X,Y,Z,[0.3 0.3 0.9]);
               axis equal
               xlim(alimits(1,:));
               ylim(alimits(2,:));
               zlim(alimits(3,:));
               view(3)
               drawnow;
           end
       elseif NDOF==1
            alimits = [1.05*min(coordinates).' 1.05*max(coordinates).'];
            [s,l]=bounds(vanim);
            
            X=reshape(coordinates(faces',1),size(faces,2),size(faces,1));
            Y=reshape(coordinates(faces',2),size(faces,2),size(faces,1)); 
            Z = reshape(coordinates(faces',3),size(faces,2),size(faces,1));
            P = zeros(size(Z));
            h = patch(X,Y,Z,P,'EdgeColor','interp' );
            set(gcf,'position',[680   208   851   770]);
            axis equal
            xlim(alimits(1,:));
            ylim(alimits(2,:));
            zlim(alimits(3,:));
            caxis([s l]);   
            colormap jet
            colorbar
            view(3)
           
            for t = 1: ceil(Nperiods/freq)
               mode = vanim*sin(2*pi*freq*t);
               P = reshape(mode(faces'),size(faces,2),size(faces,1));
       
               h.CData=P;
               drawnow

            end
           
       end
           
    end
    
end


