clear,close all
clc

Fload = 1000;

%%
Lx = 1;
Ly = 0.2;
Lz = .01;

Nex = 20;
Ney = 6;
Nez = 1;

[XYZ, ELEM ] = hex_mesh_3D( [Lx Ly Lz], [Nex Ney Nez], 0);

[NN,NDOF] = size(XYZ);

%% Assemble system matrices
mat.rho = 2700;
mat.E = 70e9;
mat.nu = 0.33;

[M, K] = assemble_system_matrices(ELEM, XYZ, mat, 'STRUCT');

%% Fixed boundary nodes. XYZ DoFs constrained
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

%% Apply Fload*Ly*Lz pressure in vertical direction on element faces @ X=Lx
tol = 1e-6;
DN = find(XYZ(:,1)==1 );
nd = length(DN);

[Fe, nF, faces] = apply_pressure(ELEM, XYZ, DN, Fload);
Fe2 = Fe;
Fe2(Di,:)=[];

%%
x = K2\Fe2;

X = zeros(NDOF*NN,1);
xi = zeros(NDOF*NN,1);
xi(Di)=1;
X(~xi,:)=x;

show_mesh(ELEM,XYZ+reshape(X,3,[]).')
hold all

[NFE,NPF] = size(faces);
Fidx = NDOF*repmat(nF,1,3)-repmat((NDOF-1):-1:0,nd,1);

Xdisp = reshape(X,3,[]).';
for ii = 1:nd
    
    DX = XYZ(nF(ii),1)+Xdisp(nF(ii),1);
    DY = XYZ(nF(ii),2)+Xdisp(nF(ii),2);
    DZ = XYZ(nF(ii),3)+Xdisp(nF(ii),3);
    U = full(Fe(Fidx(ii,1)))/1e3;  
    V = full(Fe(Fidx(ii,2)))/1e3;
    W = full(Fe(Fidx(ii,3)))/1e3;
    
    quiver3(DX, DY, DZ, U, V, W);
end
    
Xd = reshape(X(Fi),3,[]).';

%%
I = Lz^3*Ly/12;
A = Ly*Lz;
G = mat.E/2/(1+mat.nu);
kappa = 5/6;

delta = Fload*Lx^3/3/mat.E/I;
dtimo = Fload*Lx/kappa/A/G+delta;
dZ = max(Xd(:,3));
disp([ delta,dtimo, dZ, (delta-dZ)/delta*100, (dtimo-dZ)/delta*100])

