clear,close all
clc

%% Define geometry for structural plate
Lx_s = 0.5;
Ly_s = 0.5;
Lz_s = .01;

Nex = 16;
Ney = 16;
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

Nez = 18;

[XYZa, ELEMa ] = hex_mesh_3D( [Lx_s Ly_s Lz_a], [Nex Ney Nez], 0);

[NNa,~] = size(XYZa);

%offset z-axis 
XYZa(:,3) = XYZa(:,3) + Lz_s;

%% Assemble system matrices for acoustic volume

mat(2).rho   = 1.18; % air
mat(2).c     = 343;
mat(2).sigma = 9e3; % flow resistivity Rockwool Flexi A Plate,
mat(2).h     = Lz_a;
mat(2).cs = [];
mat(2).K_s = [];
mat(2).rho_s = [];

f = 1000;
mat(2) = APM(f,mat(2));

[Mf, Kf] = assemble_system_matrices(ELEMa, XYZa, mat(2), 'ACOU');
return
%% Define geometry for second structural plate

Nez = 1;

[XYZs2, ELEMs2 ] = hex_mesh_3D( [Lx_s Ly_s Lz_s], [Nex Ney Nez], 0);


%offset z-axis 
XYZs2(:,3) = XYZs2(:,3) + Lz_s + Lz_a;

%% Assemble system matrices for second structural plate

mat(3).rho = 1000;
mat(3).E = 800e6*7.4940;
mat(3).nu = 0.3;
[Ms2, Ks2] = assemble_system_matrices(ELEMs2, XYZs2, mat(3), 'STRUCT');

%%
ELEM = [ELEMs; ELEMa+NN; ELEMs2+NN+NNa];
XYZ  = [XYZs; XYZa; XYZs2];

show_mesh(ELEM, XYZ)

%% Assemble global system matrices
M = blkdiag(Ms,Mf,Ms2);
K = blkdiag(Ks,Kf,Ks2);

%% Create FSI coupling matrix
DNs = find(XYZs(:,3)==Lz_s);
DNa = find(XYZa(:,3)==Lz_s)+NN;

s.tot = size(K,1);
s.off = [0 length(XYZs)*(NDOF-1)]; % correction for DoFs counting
[R] = FSI_coupling_matrix(ELEM,XYZ,DNs,DNa,s);


DNs = find(XYZs2(:,3)==Lz_s+Lz_a)+NN+NNa;
DNa = find(XYZa(:,3)==Lz_s+Lz_a)+NN;

s.off = [-length(XYZa)*(NDOF-1) length(XYZs)*(NDOF-1)]; % correction for DoFs counting
[R2] = FSI_coupling_matrix(ELEM,XYZ,DNs,DNa,s);

%%
M = M + mat(2).rho*R.' + mat(2).rho*R2.';
K = K - R - R2 ;

%% Find nodes for fixed support constraint
% find edges of plate 1
DN = unique([ find(XYZs(:,2)==0) ; find(XYZs(:,2)==Ly_s) ; find(XYZs(:,1)==0) ; find(XYZs(:,1)==Lx_s)]);
% DN = unique([ find(XYZs(:,3)==0) ; find(XYZs(:,3)==Ly_s)]);
% find edges of plate 
DN2 = unique([ find(XYZs2(:,2)==0) ; find(XYZs2(:,2)==Ly_s) ; find(XYZs2(:,1)==0) ; find(XYZs2(:,1)==Lx_s)]) + NN + NNa;

nd = length(DN);

Di = (NDOF*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';
Di = Di(:);

nd2 = length(DN2);
Di2 = (NDOF*repmat(DN2,1,3)-repmat([2 1 0],nd2,1)).';
Di2 = Di2 - length(XYZa)*(NDOF-1);% correct for missing 2 DoF of the pressure nodes
Di2 = Di2(:);

%% Apply zero displacment constraints to system matrices

A=K;
B=M;
A([Di; Di2],:)=[];
A(:,[Di; Di2])=[];

B([Di; Di2],:)=[];
B(:,[Di; Di2])=[];


%%
% 
% % [ V,D ] = eigs(-Kf,Mf,10,'sm','Tolerance',1e-24,'MaxIterations',1000,'SubspaceDimension',30);%,'IsCholesky',true,'CholeskyPermutation',s);
% % [V,D] = eigs_unsym(M2,K2,(2*pi*166)^2,100,8);
% % [V,D] = eigs_unsym_v2(M2,K2,10,(2*pi*166)^2,50,8);
% % [V,D,W] = eigs_unsym_v3(B,A,(2*pi*166),80,8);
% neig = 10;
% [W,V,D,resl,resr] = eigs_unsym_v2(A,B,(2*pi*166)^2 , neig);
% 
% l = sqrt(D)/2/pi;
% disp(real(l))
% 
% return
% 
% %% Animate modes
% 
% n = 8; % Select mode number
% 
% X = zeros(NDOF*NN+NNa+NDOF*NN,1);
% xi = zeros(NDOF*NN+NNa+NDOF*NN,1);
% 
% xi([Di; Di2])=1;
% X(~xi,:)=V(:,n);
% 
% clc
% close all
% disp(real(l(n)))
% animate_mode(ELEMa,XYZa,X(NDOF*NN+1:NDOF*NN+NNa),1);
% 
% % Plot acoustic mode of cavity
% animate_mode(ELEMs,XYZs,X(1:NDOF*NN),3);
% disp(max(abs(X(1:NDOF*NN))))
% 
% animate_mode(ELEMs2,XYZs2,X(NDOF*NN+NNa+1:end),3);
% disp(max(abs(X(NDOF*NN+NNa+1:end))))
% 
% return
%% Static load 
tol = 1e-6;
DN = find(XYZs(:,1)==Lx_s/2 &  XYZs(:,2)==Ly_s/2  & XYZs(:,3)==0 );
nd = length(DN);

Fload = 1e3;

Fe = sparse(size(K,1),1);
Fe(DN*3) = Fload;


Fe([Di;Di2]) = [];

%% Harmonic response
Nf = 500;
f = logspace(log10(120),log10(400),Nf);

xi = zeros(2*NDOF*NN+NNa,1);
xi([Di;Di2])=1;

X = zeros(2*NDOF*NN+NNa,Nf); 
for ii = 1 : Nf
    w = 2*pi*f(ii);
    H = -w^2*B + A ;
    [L,U,P,Q,D] = lu(H) ;
    
    y = Q*(U\(L\(P*(D\Fe))));

    X(~xi,ii)=y;   
    clc
    disp(ii)
end
%%
DNo(1) = find(XYZs2(:,1)==Lx_s/2 &  XYZs2(:,2)==Ly_s/2  & XYZs2(:,3)==2*Lz_s+Lz_a )+NN+NNa;
DNo(2) = find(XYZs(:,1)==Lx_s/2 &  XYZs(:,2)==Ly_s/2  & XYZs(:,3)==0);

DNo(1) = 3*DNo(1)-NNa*(NDOF-1); % Z displacement
DNo(2) = 3*DNo(2); % Z displacement

figure(1)
subplot(211)
loglog(f, abs(X(DNo,:))/Fload)
grid 
subplot(212)
semilogx(f, angle(X(DNo,:))/pi*180)
grid 
return
%% 
Fidz = d(1)*3;
A = -(2*pi*f).^2.*X(Fidz,:);
dP = (X(DN(3),:) - X(DN(2),:))*Nez;
figure(2)
subplot(211)
loglog(f, abs(-dP)); % pressure
hold all
loglog(f, abs(A)*mat(2).rho,'.')
grid 
subplot(212)
semilogx(f, angle(-dP)/pi*180); % pressure
hold all
semilogx(f, angle(A)/pi*180,'.')
grid 

