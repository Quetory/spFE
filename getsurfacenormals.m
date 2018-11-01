function [n, areas, faces] = getsurfacenormals(elements, coordinates) 

faces(:,:,1)=elements(:,1:4);
faces(:,:,2)=elements(:,5:8);
faces(:,:,3)=elements(:,[1 4 8 5]);
faces(:,:,4)=elements(:,[2 3 7 6]);
faces(:,:,5)=elements(:,[1 2 6 5]);
faces(:,:,6)=elements(:,[3 4 8 7]);

idx = arrayfun(@(i) sum( any( isnan( coordinates( faces(:,:,i),: ) ) ) ) == 0, 1:6 );
faces = faces(:,:,idx);       

%face corners index 
A = faces(:,1); 
B = faces(:,2); 
C = faces(:,3);

%face normals 
n = cross(coordinates(B,:)-coordinates(A,:),coordinates(C,:)-coordinates(A,:),2); %area weighted
areas = vecnorm(n,2,2);
n = diag(areas)\n;


%% DEBUG
% x = coordinates(A(1),1);
% y = coordinates(A(1),2);
% z = coordinates(A(1),3);
% 
% u = coordinates(B(1),:)-coordinates(A(1),:);
% v = coordinates(C(1),:)-coordinates(A(1),:);
% 
% hold all
% quiver3(x,y,z,u(1),0,0)
% quiver3(x,y,z,0,v(2),0,'r')
% quiver3(x,y,z,n(1,1),n(1,2),n(1,3),'k')
