clear all
clc
%%
%% Define geometry for structural plate
Lx_s = 0.5;
Ly_s = 0.5;
Lz_s = .01;

Lz_a = 1;

Nex = 16;
Ney = 16;
Nez = 1;
Nez_a = 18;

mat(1).rho = 1000;
mat(1).E = 800e6*7.4940;
mat(1).nu = 0.3;


mat(2).rho   = 1.18; % air
mat(2).c     = 343;
mat(2).sigma = 100; % flow resistivity Rockwool Flexi A Plate,
mat(2).h     = Lz_a;
mat(2).cs = [];
mat(2).K_s = [];
mat(2).rho_s = [];
% mat(2).a00=1.04;        %tortuosity
% mat(2).p=0.98;
% mat(2).L = 129e-6;
% mat(2).Lt= 198e-6;

%%
[XYZs, ELEMs ] = hex_mesh_3D( [Lx_s Ly_s Lz_s], [Nex Ney Nez], 0);

[NN,NDOF] = size(XYZs);

%% Assemble system matrices for structural plate

[Ms, Ks] = assemble_system_matrices(ELEMs, XYZs, mat(1), 'STRUCT');

%% Define geometry for acoustic volume
Nez = 18;

[XYZa, ELEMa ] = hex_mesh_3D( [Lx_s Ly_s Lz_a], [Nex Ney Nez_a], 0);

[NNa,~] = size(XYZa);

%offset z-axis 
XYZa(:,3) = XYZa(:,3) + Lz_s;

%%
ELEM = [ELEMs; ELEMa+NN];
XYZ  = [XYZs; XYZa];
%%
Nf = 130;
f = logspace(log10(120),log10(190),Nf);

xi = zeros(NDOF*NN+NNa,1);


X = zeros(NDOF*NN+NNa,Nf); 
% return
for ii = 1 : Nf
    w = 2*pi*f(ii);
    clear K M 
    % Assemble system matrices for acoustic volume
    mat_acou.rho   = 1.18; % air
    mat_acou.c     = 343;
    mat_acou.sigma = 100; % flow resistivity Rockwool Flexi A Plate,
    mat_acou.cs = [];
    mat_acou.K_s = [];
    mat_acou.rho_s = [];
    
    mat_acou = APM( f(ii), mat_acou );

%     mat(2) = JCA( f(ii), mat(2) );
%     mat_acou.cs = 343*(0.92+0.098i);
%     mat_acou.rho_s = 1.18*(1.124-0.203i);
%     mat_acou.K_s = mat_acou.rho_s*mat_acou.cs^2;
    [Mf, Kf] = assemble_system_matrices(ELEMa, XYZa, mat_acou, 'ACOU');

    % Assemble global system matrices
    M = blkdiag(Ms,Mf);%,Ms2);
    K = blkdiag(Ks,Kf);%,Ks2);
       
    % Create FSI coupling matrix
    DNs = find(XYZs(:,3)==Lz_s);
    DNa = find(XYZa(:,3)==Lz_s)+NN;

    s.tot = size(K,1);
    s.off = [0 length(XYZs)*(NDOF-1)]; % correction for DoFs counting
    [R] = FSI_coupling_matrix(ELEM,XYZ,DNs,DNa,s);
    
    %%
    M = M + R.'; %mat(2).rho_s
    K = K - R ;
%     r = R(1:size(Ks,2),size(Ks,2)+1:end);
%     Cf = -w*imag(Mf) + 1/w*imag(Kf);
%     C = [ sparse(size(Ks,1),size(Ks,2)) -r; -r.' -Cf] ;
    % Find nodes for fixed support constraint
    % find edges of plate 1
    DN = unique([ find(XYZs(:,2)==0) ; find(XYZs(:,2)==Ly_s) ; find(XYZs(:,1)==0) ; find(XYZs(:,1)==Lx_s)]);

    nd = length(DN);

    Di = (NDOF*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';
    Di = Di(:);

    % Apply zero displacment constraints to system matrices
    A=K;
    B=M;
    A(Di,:)=[];
    A(:,Di)=[];

    B(Di,:)=[];
    B(:,Di)=[];
% 
%     C(Di,:)=[];
%     C(:,Di)=[];
    % Static load 

    DN = find(XYZs(:,1)==Lx_s/2 &  XYZs(:,2)==Ly_s/2  & XYZs(:,3)==0 );
    nd = length(DN);

    Fload = 1e3;

    Fe = sparse(size(K,1),1);
    Fe(DN*3) = Fload;

    Fe(Di) = [];

    % Harmonic response
   
    H = -w^2*B + A ;
%     clear L U P Q D
%     [L,U,P,Q,D] = lu(H) ;

%     y = Q*(U\(L\(P*(D\Fe))));
    y =H\Fe;
    xi(Di)=1;
    X(~xi,ii)=y;   
    clc
    disp(ii)
    disp(mat_acou.rho_s)
end

%%
DNo(1) = find(XYZs(:,1)==Lx_s/2 &  XYZs(:,2)==Ly_s/2  & XYZs(:,3)==0);
DNo(1) = 3*DNo(1); % Z displacement

fabs = f;
Habs = X(DNo,:);
sigma=mat(2).sigma;
% save Transfer_abs_s100v4 Habs fabs X DNo sigma

% Hab = importdata('foo_abs.txt');
% Hab=importdata('acou_abs_const.txt');

%%
clc
Hab=importdata('acou_abs_miki.txt');

figure(1)
subplot(211)
loglog(f, abs(X(DNo,:)),'.-')
hold all
loglog(Hab(1:130,1), abs(Hab(1:130,2)+1i*Hab(1:130,3)),'.-')

grid on
subplot(212)
semilogx(f, angle(X(DNo,:))/pi*180)
hold all
semilogx(Hab(1:130,1), angle(Hab(1:130,2)+1i*Hab(1:130,3))/pi*180,'.-')
grid on
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

