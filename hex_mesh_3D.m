function [XYZ, ELEM ] = hex_mesh_3D( dim, num_elem, plt_key)

Lx = dim(1);
Ly = dim(2);
Lz = dim(3);

Nex = num_elem(1);
Ney = num_elem(2);
Nez = num_elem(3);

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

%% Generate element connectivity matrix
NE = size(center,1);

E_r2 =sum((coord(1,:)-center(1,:)).^2,2);

idx = arrayfun(@(i) find(sum(bsxfun(@minus,coord,center(i,:)).^2,2)<=E_r2*1.001),1:NE,'UniformOutput', false);
ELEM = cell2mat(idx).' ;
ELEM = ELEM(:,[1 3 4 2 5 7 8 6]); % reordering 

XYZ = coord;
%%
if plt_key==1
    figure(100)
    show_mesh(ELEM,coord);
end

