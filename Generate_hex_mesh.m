clear,close all
clc
%%
Lx = 1;
Ly = 1;
Lz = .25;

Nex = 4;
Ney = 4;
Nez = 4;

NPE = 8;

%% Generate mesh + coordinates
dx = linspace(0,Lx,Nex+1);
dy= linspace(0,Ly,Ney+1);
dz = linspace(0,Lz,Nez+1);

[x,y,z]=meshgrid(dx,dy,dz);

xc = movsum(dx,2,'Endpoints','discard')/2;
yc = movsum(dy,2,'Endpoints','discard')/2;
zc = movsum(dz,2,'Endpoints','discard')/2;

[xc,yc,zc]=meshgrid(xc,yc,zc);

coord = [x(:) y(:) z(:)];
center=[xc(:) yc(:) zc(:)];

clear x y z xc yc zc

%% Generate element connectivity matrix
NE = size(center,1);


E_r2 =sum((coord(1,:)-center(1,:)).^2,2);

idx = arrayfun(@(i) find(sum(bsxfun(@minus,coord,center(i,:)).^2,2)<=E_r2*1.001),1:NE,'UniformOutput', false);
ELEM = cell2mat(idx).' ;
ELEM = ELEM(:,[1 3 4 2 5 7 8 6]);

show_mesh(ELEM,coord);

clear center

%%

rho = 1.184;
c = 346.13;
XYZ = coord(ELEM(1,:),:);
[Melem] = ElemAcouMass(XYZ,rho,c);
[Kelem] = ElemAcouStiffness(XYZ,rho);
NDOF = 1;


Ke = repmat(Kelem,1,1,NE);
Me = repmat(Melem,1,1,NE);

Elements=NDOF*ELEM(:,kron(1:NPE,ones(1,NDOF)))-kron(ones(NE,1),kron(ones(1,NPE),(NDOF-1):-1:0));

Y=reshape(repmat(Elements,1,NDOF*NPE)',NDOF*NPE,NDOF*NPE,NE);
X=permute(Y,[2 1 3]);

K=sparse(X(:),Y(:),Ke(:));
M=sparse(X(:),Y(:),Me(:));

uns = 1/2*(K.'-K);
K = K+uns;

uns = 1/2*(M.'-M);
M=M+uns;


clear Ke Me Kelem Melem Elements X Y 

%%
opts.spdB= 1;
opts.tol = 1e-12;
[V,L]=eigs(-K,M,20,'sm',opts);
V = real(V);
l = diag(L);
feig = real(sqrt(-l)/2/pi);


[~,idx]=min(abs(feig-(c/2/Lz)));
disp([c/2/Lz feig(idx)])

%%
% close all
% animate_mode(ELEM,coord,V(:,idx),1);

%%

BCnodes = find(coord(:,3)==Lz);

XYZ = NaN(size(coord));
XYZ((coord(:,3)==Lz),:)= coord((coord(:,3)==Lz),:);

nBC = length(BCnodes);

[idx,~] = (arrayfun(@(i)  ind2sub(size(ELEM),find((ELEM-BCnodes(i))==0)), 1:nBC,'UniformOutput', false));

BCelem = [];
for i=1:nBC
    BCelem = [BCelem ; idx{i}];
end
BCelem = unique(BCelem);
    
show_mesh(ELEM(BCelem,:),coord);

[n,areas ] = getsurfacenormals(ELEM(BCelem,:),XYZ);

% TO-DO add surface integral of FSI surfaces