clear,close all
clc
%%
Lx = 1;
Ly = 0.2;
Lz = .01;

Nex = 20;
Ney = 5;
Nez = 1;

[XYZ, ELEM ] = hex_mesh_3D( [Lx Ly Lz], [Nex Ney Nez], 0);

[NN,NDOF] = size(XYZ);

%% Assemble system matrices
mat.rho = 2700;
mat.E = 70e9;
mat.nu = 0.33;

[M, K] = assemble_system_matrices(ELEM, XYZ, mat, 'STRUCT');

%%
tol = 1e-6;
DN = find(XYZ(:,1)==0);
nd = length(DN);

Di = (3*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';
Di = Di(:);

K2=K;
M2=M;
K2(Di,:)=[];
K2(:,Di)=[];

M2(Di,:)=[];
M2(:,Di)=[];

%%
tol = 1e-6;
DN = find(XYZ(:,1)==1 & XYZ(:,3)==Lz);
nd = length(DN);

Fi = (3*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';
Fi = Fi(:);

F2 = spalloc(size(K,1),1, nd);
Fidx = Fi(3:3:end);

F2(Fidx,1) = 1000*ones(nd,1)/nd;
F2(Di,:)=[];

%%
x = K2\F2;


X = zeros(NDOF*NN,1);
xi = zeros(NDOF*NN,1);
xi(Di)=1;
X(~xi,:)=x;

show_mesh(ELEM,XYZ+reshape(X,3,[]).')

Xd = reshape(X(Fi),3,[]).';

%%
I = Lz^3*Ly/12;
delta = 1000*Lx^3/3/mat.E/I;

disp([ delta, max(Xd(:,3))])

