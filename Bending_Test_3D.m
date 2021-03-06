clear,close all
clc

%% Define geometry
Lx = 10;
Ly = 1;
Lz = 20;

Nex = 10;
Ney = 1;
Nez = 30;

[XYZ, ELEM ] = hex_mesh_3D( [Lx Ly Lz], [Nex Ney Nez], 1);

[NN,NDOF] = size(XYZ);

%% Assemble system matrices
mat.rho = 2700;
mat.E = 1.0e3;
mat.nu = 0.0;

[M, K] = assemble_system_matrices(ELEM, XYZ, mat, 'STRUCT');

%% Fixed boundary nodes. XYZ DoFs constrained
DN = find(XYZ(:,3)==0);
nd = length(DN);

Di = (3*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';

% DN = find(XYZ(:,1)==0 & XYZ(:,2)==12);
% nd = length(DN);
% Di2 = (3*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';
% Di2(2,:)=[];

Di = Di(:);

K2=K;
M2=M;
K2(Di,:)=[];
K2(:,Di)=[];

M2(Di,:)=[];
M2(:,Di)=[];


%% Modal analysis
% [V,D ] = eigs(K2,M2,6,'sm');%,'IsCholesky',true,'CholeskyPermutation',s);
% N = vecnorm(V,2,1);
% Vn = V/diag(N); % Normalize shape functions
% 
% X = zeros(NDOF*NN,1);
% xi = zeros(NDOF*NN,1);
% xi(Di)=1;
% X(~xi,:)=Vn(:,2);
% l = sqrt(diag(D))/2/pi;
% disp(['Resonance Frequency: ' num2str(l(2),3) ' Hz'])
% 
% animate_mode(ELEM,XYZ,X,NDOF);
% 
% return
%% Apply Fload*Ly*Lz pressure in vertical direction on element faces @ X=Lx
tol = 1e-6;
DN = find(XYZ(:,3)== Lz & XYZ(:,1)==Lx );
nd = length(DN);

Fload = 10;
Pload = 2;

% [Fe, nF, faces] = apply_pressure(ELEM, XYZ, DN, Pload, [0;1;0]);
% Fe2 = Fe;
% Fe2(Di,:)=[];

Fe = sparse(3*DN-1,1,Fload/2,size(K,1),1);

DN = find(XYZ(:,3)== Lz & XYZ(:,1)==0 );
nd = length(DN);

Fe2 = sparse(3*DN-1,1,Fload/2,size(K,1),1);
Fe2 = Fe+Fe2;
% Fe2 = Fe;
Fe2(Di,:)=[];
%%
x = K2\Fe2;

X = zeros(NDOF*NN,1);
xi = zeros(NDOF*NN,1);
xi(Di)=1;
X(~xi,:)=x;

show_mesh(ELEM,XYZ+reshape(X,3,[]).')
% hold all
% 
% [NFE,NPF] = size(faces);
nF=DN.';
Fidx = NDOF*repmat(nF.',1,3)-repmat((NDOF-1):-1:0,nd,1);
% 
% Xdisp = reshape(X,3,[]).';
% for ii = 1:nd
%     
%     DX = XYZ(nF(ii),1)+Xdisp(nF(ii),1);
%     DY = XYZ(nF(ii),2)+Xdisp(nF(ii),2);
%     DZ = XYZ(nF(ii),3)+Xdisp(nF(ii),3);
%     U = full(Fe(Fidx(ii,1)))/1e3;  
%     V = full(Fe(Fidx(ii,2)))/1e3;
%     W = full(Fe(Fidx(ii,3)))/1e3;
%     
%     quiver3(DX, DY, DZ, U, V, W);
% end
%     
Xd = X(Fidx);
disp(Xd)

return
%%
I = Lz^3*Ly/12;
A = Ly*Lz;
G = mat.E/2/(1+mat.nu);
kappa = 5/6;

delta = Fload*Lx^3/3/mat.E/I;
dtimo = Fload*Lx/kappa/A/G+delta;

c1 = 0.246;
c2=0.229;
phi = 1*Lx/c2/Ly/Lz^3/G;
dZ_ana = Ly/2*tan(phi);

dZ = max(Xd(:,3));
clc
disp(['timoshenko    ', 'Actual    ',  'relative error %'])
disp([ dtimo, dZ, (dtimo-dZ)/dtimo*100])

