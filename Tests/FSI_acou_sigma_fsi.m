clear all
clc
tic
%%
%% Define geometry for structural plate
Lx_s = 0.5;
Ly_s = 0.5;
Lz_s = .01;

Lz_a = 1.0;

Nex = 16;
Ney = 16;
Nez = 1;
Nez_a = 18;

mat(1).rho = 1000;
mat(1).E = 800e6*7.4940;
mat(1).nu = 0.3;

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
Nf = 400;
f = logspace(log10(120), 3,Nf);

xi = zeros(NDOF*NN+NNa,1);

X = zeros(NDOF*NN+NNa,Nf); 
%%

for ii = 1 : Nf
    w = 2*pi*f(ii);
    clear K M mat_acou
    % Assemble system matrices for acoustic volume
    mat_acou.rho   = 1.18; % air @ T=22 degC
    mat_acou.c     = 343;  % air @ T=22 degC
    mat_acou.sigma = 1000;
    mat_acou.cs = [];
    mat_acou.K_s = [];
    mat_acou.rho_s = [];
    
    mat_acou = APM( f(ii), mat_acou );

    [Mf, Kf] = assemble_system_matrices(ELEMa, XYZa, mat_acou, 'ACOU');

    % Assemble global system matrices
    M = blkdiag(Ms, -real(Mf));
    K = blkdiag(Ks, -real(Kf));
       
    % Create FSI coupling matrix
    DNs = find(XYZs(:,3)==Lz_s);
    DNa = find(XYZa(:,3)==Lz_s)+NN;

    s.tot = size(K,1);
    s.off = [0 length(XYZs)*(NDOF-1)]; % correction for DoFs counting
    [R] = FSI_coupling_matrix(ELEM,XYZ,DNs,DNa,s);
    
%     [Cfsi] = FSI_impedance_matrix(ELEM,XYZ,DNa,mat_acou,s);
     
    r = R(1:size(Ks,2),size(Ks,2)+1:end);
%     M = M + R.'; 
%     K = K - R ;
    
    Cf = -w*imag(Mf) + 1/w*imag(Kf);
    C = [ sparse(size(Ks,1),size(Ks,2)) -r ; -r.' -Cf ] ;

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

    C(Di,:)=[];
    C(:,Di)=[];
%     
    % Static load 
    DN = find(XYZs(:,1)==Lx_s/2 &  XYZs(:,2)==Ly_s/2  & XYZs(:,3)==0 );
    nd = length(DN);
    Fload = 1;

    Fe = sparse(size(K,1),1);
    Fe(DN*3) = Fload;

    Fe(Di) = [];

    % Harmonic response
    H = -w^2*B  + 1i*w*C + A ;% + 1i*w*C
    clear L U P Q D
    [L,U,P,Q,D] = lu(H) ;

    y = Q*(U\(L\(P*(D\Fe))));
%     y =H\Fe;
    xi(Di)=1;
    X(~xi,ii)=y;   
    clc
    disp(ii)
    disp(mat_acou.rho_s)
end

%% Define output node
DNo(1) = find(XYZs(:,1)==Lx_s/2 &  XYZs(:,2)==Ly_s/2  & XYZs(:,3)==0);
DNo(1) = 3*DNo(1); % Z displacement

DNo(2) = find(XYZa(:,1)==0.3125 &  XYZa(:,2)==0.93750e-1  & XYZa(:,3)==0.51);
DNo(2) = 3*NN+DNo(2); % Z displacement


Hn = X(DNo,:);
toc
return
save Transfer_sigma_10_fsi Hn DNo DN Fload f

% Import Ansys data
Hab=importdata('acou_fsi_s10_uz_120_1e3.txt');

fa = Hab(:,1);
Hans(1,:) = (Hab(1:end,2)+1i*Hab(1:end,3));

Hab=importdata('acou_fsi_s10_p_120_1e3.txt');
Hans(2,:)=(Hab(1:end,2)+1i*Hab(1:end,3));


%% Plot
close all
figure(1)
subplot(211)
loglog(f, abs(Hn(1,:))/Fload)
hold all
loglog(fa, abs(Hans(1,:))/Fload)
grid 
subplot(212)
semilogx(f, angle(Hn(1,:))/pi*180)
hold all
semilogx(fa, angle(Hans(1,:))/pi*180)
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

