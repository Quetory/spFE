clear,close all
clc
%%
numx = 2;
numy = 2;
numz = 2;

[xc,yc,zc]=meshgrid( 0:1/numx:1, 0:1/numy:1,0:1/numz:1);
nodes = [xc(:) yc(:) zc(:)];

coordinates = nodes;%([1 2 4 3 5 6 8 7],:);

plot3(coordinates(:,1),coordinates(:,2),coordinates(:,2),'ko');

return
%%

node_pattern=[ 1 2 nnx+1 nnx+2 nny+1 nny+2 nny+nnx+1 nny+nnx+2 ]; % Node pattern 1 2 3 4 5 6 7 8 % according to me notations
element_pattern=[1 2 numx+1 numx+2]; %Element Pattern 1 2 3 4 ; according to my notations
element=make_elem_hexa(node_pattern,numx,numy,numz,inc_u,inc_v,inc_w,nnx);


%%
elements = 1:8;
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


patch('Faces',faces,'Vertices',[X(:) Y(:) Z(:)], 'FaceAlpha',0.1,'FaceColor','r' );
view(3)