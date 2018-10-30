function [n, areas] = getsurfacenormals(elements, coordinates) 
%Vertex normals of a triangulated mesh, area weighted, left-hand-rule 
% N = patchnormals(FV) -struct with fields, faces Nx3 and vertices Mx3 
%N: vertex normals as Mx3

FV.vertices = coordinates;

faces(:,:,1)=elements(:,1:4);
faces(:,:,2)=elements(:,5:8);
faces(:,:,3)=elements(:,[1 4 8 5]);
faces(:,:,4)=elements(:,[2 3 7 6]);
faces(:,:,5)=elements(:,[1 2 6 5]);
faces(:,:,6)=elements(:,[3 4 8 7]);

idx = arrayfun(@(i) sum(any(isnan(coordinates(faces(:,:,i)))))==0,1:6);
FV.faces = faces(:,:,idx);       

%face corners index 
A = FV.faces(:,1); 
B = FV.faces(:,2); 
C = FV.faces(:,3);

%face normals 
n = cross(FV.vertices(B,:)-FV.vertices(A,:),FV.vertices(C,:)-FV.vertices(A,:),2); %area weighted
areas = (sqrt(sum(n.^2, 2)));

%% DEBUG
% x = FV.vertices(A(1),1);
% y = FV.vertices(A(1),2);
% z = FV.vertices(A(1),3);
% 
% u = FV.vertices(B(1),:)-FV.vertices(A(1),:);
% v = FV.vertices(C(1),:)-FV.vertices(A(1),:);
% 
% hold all
% quiver3(x,y,z,u(1),0,0)
% quiver3(x,y,z,0,v(2),0,'r')
% quiver3(x,y,z,n(1,1),n(1,2),n(1,3),'k')
