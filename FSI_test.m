clear,close all
clc

%% Define geometry for structural plate
Lx_s = 0.5;
Ly_s = 0.5;
Lz_s = .01;

Nex = 15;
Ney = 15;
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
s.off = length(XYZs)*(NDOF-1);
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

A=K;
B=M;
A(Di,:)=[];
A(:,Di)=[];

B(Di,:)=[];
B(:,Di)=[];

%%

% [ V,D ] = eigs(-Kf,Mf,10,'sm','Tolerance',1e-24,'MaxIterations',1000,'SubspaceDimension',30);%,'IsCholesky',true,'CholeskyPermutation',s);
% [V,D] = eigs_unsym(M2,K2,(2*pi*166)^2,100,8);
% [V,D] = eigs_unsym_v2(M2,K2,10,(2*pi*166)^2,50,8);
% [V,D,W] = eigs_unsym_v3(B,A,(2*pi*166),80,8);
[W,V,D,resl,resr] = eigs_unsym_v2(-A,B,(2*pi*166)^2);

return
%%

%%

X = zeros(NDOF*NN+NNa,1);
xi = zeros(NDOF*NN+NNa,1);

% X = zeros(NN+NNa,1);
% xi = zeros(NN+NNa,1);
% Di=[];
xi(Di)=1;
%% Plot structural mode of plate
n = 4;
X(~xi,:)=V(:,n);
l = sqrt(-D)/2/pi;
disp(l)
% clc
% disp([ (1:8).' , real(l)])
% fres = CalcClampedPlateFres(Lx_s,Ly_s,Lz_s,mat(1));
% disp(['First plate resonance (mode 11) analytic: ' num2str(fres,3) ' Hz'])
% disp(['Resonance Frequency: ' num2str(l(3),3) ' Hz'])


close all
animate_mode(ELEMa,XYZa,X(NDOF*NN+1:end),1);
% 
% Plot acoustic mode of cavity
% animate_mode(ELEMs,XYZs,X(1:NDOF*NN),3);
return
%% Static load 
tol = 1e-6;
DN = find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & XYZ(:,3)==0 );
nd = length(DN);

Fload = 1e3;

Fe = sparse(size(K,1),1);
Fe(DN*3) = Fload;

Fe2 = Fe;
Fe2(Di) = [];

%% Harmonic response
Nf = 200;
f = logspace(log10(160),log10(180),Nf);

xi = zeros(NDOF*NN+NNa,1);
xi(Di)=1;

X = zeros(NDOF*NN+NNa,Nf); 
for ii = 1 : Nf
    w = 2*pi*f(ii);
    H = -w^2*M2 + K2 ;
    [L,U,P,Q,D] = lu(H) ;
    
    y = Q*(U\(L\(P*(D\Fe2))));

    X(~xi,ii)=y;   
    clc
    disp(ii)
end
%%
d = find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & XYZ(:,3)== 0.01 );
DN = d(2);

dn = arrayfun(@(x) find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & (XYZ(:,3)>=x-1e-4+0.01 & XYZ(:,3)<=x+1e-4+0.01) ),1/Nez :1/Nez :1);%,'UniformOutput',false)

DN = [DN;dn.']+(NDOF-1)*NN;

figure(1)
subplot(211)
loglog(f, abs(X(DN,:)))
grid 
subplot(212)
semilogx(f, angle(X(DN,:))/pi*180)
grid 

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

