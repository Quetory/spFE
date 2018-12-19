clear all
clc
%%
%% Define geometry for structural plate
Lx_s = 0.5;
Ly_s = 0.5;

Lz_a = 1;

Nex = 18;
Ney = 18;
Nez_a = 36;


mat.rho   = 1.18; % air
mat.c     = 343;
mat.sigma = 10; % flow resistivity Rockwool Flexi A Plate,
mat.h     = Lz_a;
mat.cs = [];
mat.K_s = [];
mat.rho_s = [];
% mesh
[XYZ, ELEM ] = hex_mesh_3D( [Lx_s Ly_s Lz_a], [Nex Ney Nez_a], 0);

[NNa,~] = size(XYZ);

%%
Nf = 100;
f = logspace(log10(120),log10(200),Nf);

xi = zeros(NNa,1);


X = zeros(NNa,Nf); 
% return
for ii = 1 : Nf
    w = 2*pi*f(ii);
    clear K M 
    % Assemble system matrices for acoustic volume
    mat_acou(ii).rho   = 1.18; % air
    mat_acou(ii).c     = 343;
    mat_acou(ii).sigma = 10; % flow resistivity Rockwool Flexi A Plate,
    mat_acou(ii).cs = [];
    mat_acou(ii).K_s = [];
    mat_acou(ii).rho_s = [];
    
    mat_acou(ii) = APM( f(ii), mat_acou(ii) );

    [M, K] = assemble_system_matrices(ELEM, XYZ, mat_acou(ii), 'ACOU');

    % Mass [kg/m^3/s] load 
    DN = find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & XYZ(:,3)==0 );
    nd = length(DN);

    Fload = 1;

    Fe = sparse(size(K,1),1);
    Fe(DN) = Fload;

    % Harmonic response
    H = -w^2*M + K ;
    clear L U P Q D
    [L,U,P,Q,D] = lu(H) ;
    y = Q*(U\(L\(P*(D\Fe))));

    X(:,ii)=y;   
    clc
    disp(ii)
    disp(mat_acou(ii).rho_s)
end

%%
DNo(1) = find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & XYZ(:,3)==1/9);
DNo(2) = find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & XYZ(:,3)==2/9);
DNo(3) = find(XYZ(:,1)==Lx_s/2 &  XYZ(:,2)==Ly_s/2  & XYZ(:,3)==3/9);



%%
close all
Hab=importdata('acou_abs_miki_no_fsi_large.txt');

fa = Hab(:,1);
Hans(1,:) = (Hab(1:end,2)+1i*Hab(1:end,3));

Hab=importdata('acou_abs_miki_no_fsi_large_2.txt');
Hans(2,:)=(Hab(1:end,2)+1i*Hab(1:end,3));

Hab=importdata('acou_abs_miki_no_fsi_large_3.txt');
Hans(3,:)=(Hab(1:end,2)+1i*Hab(1:end,3));



figure(1)
subplot(211)
loglog(f, abs(X(DNo,:).*(2*pi*1i*f)./abs([mat_acou.rho_s])))
hold all
loglog(fa, abs(Hans),'.')

grid on
subplot(212)
semilogx(f, angle(X(DNo,:).*(2*pi*1i*f)./abs([mat_acou.rho_s]))/pi*180)
hold all
semilogx(fa, angle(Hans)/pi*180,'.')
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
loglog(f, abs(A)*mat.rho,'.')
grid 
subplot(212)
semilogx(f, angle(-dP)/pi*180); % pressure
hold all
semilogx(f, angle(A)/pi*180,'.')
grid 

