clear,close all
clc

%% Define geometry for structural plate
Lx_s = 0.5;
Ly_s = 0.5;
Lz_s = .01;

Nex = 20;
Ney = 20;
Nez = 1;

[XYZs, ELEMs ] = hex_mesh_3D( [Lx_s Ly_s Lz_s], [Nex Ney Nez], 0);

[NN,NDOF] = size(XYZs);

%% Assemble system matrices for structural plate
mat(1).rho = 1000;
mat(1).E = 800e6*7.4940;
mat(1).nu = 0.3;

[Ms, Ks] = assemble_system_matrices(ELEMs, XYZs, mat(1), 'STRUCT');

%% Define geometry for acoustic volume
Lz_a = 1;

Nez = 12;

[XYZa, ELEMa ] = hex_mesh_3D( [Lx_s Ly_s Lz_a], [Nex Ney Nez], 0);

[NNa,~] = size(XYZa);

%offset z-axis 
XYZa(:,3) = XYZa(:,3) + Lz_s;

%% Assemble system matrices for acoustic volume

mat(2).rho = 1.18;
mat(2).c = 343;
[Mf, Kf] = assemble_system_matrices(ELEMa, XYZa, mat(2), 'ACOU');


%%
ELEM = [ELEMs; ELEMa+NN];
XYZ  = [XYZs; XYZa];

show_mesh(ELEM, XYZ)

%%
M = blkdiag(Ms,Mf);
K = blkdiag(Ks,Kf);

s.tot = size(K,1);
s.off = length(XYZs)*2;
%%
DNs = find(XYZs(:,3)==Lz_s);
DNa = find(XYZa(:,3)==Lz_s)+NN;

[R] = FSI_coupling_matrix(ELEM,XYZ,DNs,DNa,s);

%%
M = M + mat(2).rho*R.';
K = K - R;

%%
DN = unique([ find(XYZs(:,2)==0) ; find(XYZs(:,2)==Ly_s) ; find(XYZs(:,1)==0) ; find(XYZs(:,1)==Lx_s)]);

nd = length(DN);

Di = (NDOF*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';
Di = Di(:);

K2=K;
M2=M;
K2(Di,:)=[];
K2(:,Di)=[];

M2(Di,:)=[];
M2(:,Di)=[];
%%
[V,D ] = eigs(K2,M2,6,'sm');%,'IsCholesky',true,'CholeskyPermutation',s);
N = vecnorm(V,2,1);
Vn = V/diag(N); % Normalize shape functions

X = zeros(NDOF*NN+NNa,1);
xi = zeros(NDOF*NN+NNa,1);
xi(Di)=1;

n = 3;
X(~xi,:)=V(:,n);
l = sqrt(diag(D))/2/pi;

clc
disp(real(l))
fres = CalcClampedPlateFres(Lx_s,Ly_s,Lz_s,mat(1));
disp(['First plate resonance (mode 11) analytic: ' num2str(fres,3) ' Hz'])
disp(['Resonance Frequency: ' num2str(l(n),3) ' Hz'])

close all
animate_mode(ELEMa,XYZa,X(NDOF*NN+1:end),1);

%%
animate_mode(ELEMs,XYZs,X(1:NDOF*NN),3);

%%
tol = 1e-6;
DN = find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & XYZ(:,3)==0 );
nd = length(DN);

Fload = 1e3;

Fe = sparse(size(K,1),1);
Fe(DN*3) = Fload;

Fe2 = Fe;
Fe2(Di) = [];



%%

