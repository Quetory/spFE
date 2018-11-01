function show_mesh(elements,coordinates)
if size(coordinates,2)==2    %2D objects 
%     if size(elements,2)==3   %triangle         
        X=reshape(coordinates(elements',1),size(elements,2),size(elements,1));
        Y=reshape(coordinates(elements',2),size(elements,2),size(elements,1));   
        patch(X,Y,[0.3 0.3 0.9]);
%     end
   
elseif size(coordinates,2)==3   %3D object
    if size(elements,2)==4      %tetrahedra
        tetramesh(elements,coordinates,'FaceAlpha',1);camorbit(20,0);
    end
    if size(elements,2)==8  || size(elements,2)==20    %tetrahedra
       
      %% face 1 (J-I-L-K), face 2 (I-J-N-M), face 3 (J-K-O-N),
      %% face 4 (K-L-P-O), face 5 (L-I-M-P), face 6 (M-N-O-P)
       faces1=elements(:,1:4);
       faces2=elements(:,5:8);
       faces3=elements(:,[1 4 8 5]);
       faces4=elements(:,[2 3 7 6]);
       faces5=elements(:,[1 2 6 5]);
       faces6=elements(:,[3 4 8 7]);
       
       faces=[faces1; faces2; faces3; faces4; faces5; faces6];
       
       X=reshape(coordinates(faces',1),size(faces,2),size(faces,1));
       Y=reshape(coordinates(faces',2),size(faces,2),size(faces,1)); 
       Z=reshape(coordinates(faces',3),size(faces,2),size(faces,1)); 
       
%        subplot(221)
       patch(X,Y,Z,[0.3 0.3 0.9],'EdgeColor',[0.6 0.6 0.6]);    
       view(3)
       axis equal
       xlabel('X');
       ylabel('Y');
       zlabel('Z');
       grid on
%        subplot(222)
%        patch(X,Y,Z,[0.3 0.3 0.9],'EdgeColor',[0.6 0.6 0.6]);    
%        view([0 0 ])
%        axis equal
%        xlabel('X');
%        ylabel('Y');
%        zlabel('Z');
%        grid on
%        subplot(223)
%        patch(X,Y,Z,[0.3 0.3 0.9],'EdgeColor',[0.6 0.6 0.6]);    
%        view([90 0])
%        axis equal
%        xlabel('X');
%        ylabel('Y');
%        zlabel('Z');
%        grid on
%        subplot(224)
%        patch(X,Y,Z,[0.3 0.3 0.9],'EdgeColor',[0.6 0.6 0.6]);    
%        view([0 90])
%        axis equal
%        xlabel('X');
%        ylabel('Y');
%        zlabel('Z');       
%        grid on
    end
    
end




