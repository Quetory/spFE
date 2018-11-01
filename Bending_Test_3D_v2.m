clear,close all
clc

eas = 1;
%%
en = importdata('e_numbers.txt');
nc = csvread('node_coords.csv');

NN = length(nc);
XYZ = nc(:,2:end);

ELEM = importdata('elem_con.txt');


[NE, NPE ] = size(ELEM);
NDOF = size(XYZ,2);

show_mesh(ELEM,XYZ)

%% Assembly of global mass and stiffness matrix
[Kelem] = ElemStiffness(XYZ(ELEM(1,:),:), 20e9,0.3, eas);
[Melem] = ElemMass(XYZ(ELEM(1,:),:), 2700 );

Ke = repmat(Kelem,1,1,NE);
Me = repmat(Melem,1,1,NE);

% Row and column sparse index
Elements=NDOF*ELEM(:,kron( 1:NPE, ones(1,NDOF ) ) ) -kron(ones(NE,1),kron(ones(1,NPE), NDOF-1:-1:0 ) );

 
Y=reshape(repmat(Elements,1,NPE*NDOF)',NPE*NDOF,NPE*NDOF,NE);
X=permute(Y,[2 1 3]);

%assemble global matrices
K = sparse(X(:),Y(:),Ke(:));
M = sparse(X(:),Y(:),Me(:));

%remove small numerical round-off errors
Ku = triu(K,1);
K = Ku.'+Ku + diag(diag(K));

Mu = triu(M,1);
M = Mu.'+Mu + diag(diag(M));

%%
tol = 1e-6;
DN = find(XYZ(:,1)==0);

hold all

%Create nodal-> equation mapping
Nmap =ELEM(1,:).';
for e = 2:NE
    
    n = ELEM(e,:);
    [C,ia,ib] = union(Nmap,n);
    Nmap = [Nmap; n(ib).'];
end
nd = length(DN);



Di = (3*repmat(DN,1,3)-repmat([2 1 0],nd,1)).';
Di = Di(:);

%%
K2=K;
M2=M;
K2(Di,:)=[];
K2(:,Di)=[];

M2(Di,:)=[];
M2(:,Di)=[];

%%
[V, l] = eigs(K2,M2,12,'sm');

lam = diag(l);
% lam(1:6)=0;
f = sqrt(lam)/2/pi;
disp(f)

% T = V.'*M2*V;
% V = V/diag(sqrt(diag(T)));
% replace zero displacement in modeshape
v = zeros(NDOF*NN,length(lam));
xi = zeros(NDOF*NN,1);
xi(Di)=1;
v(~xi,:)=V;


%%
close all
animate_mode(ELEM,XYZ,v(:,2),3)
