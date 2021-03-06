clear,close all
clc
%%
Lx = 0.5;
Ly = 1;
Lz = 3;

Nex = 1;
Ney = 1;
Nez = 1;

[XYZ, ELEM ] = hex_mesh_3D( [Lx Ly Lz], [Nex Ney Nez], 1);

R = @(a) [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a) ];

XYZ = (R(15).'*(XYZ*R(15)).').';
%% Assemble system matrices
mat.rho = 1.184;
mat.c = 346.13;

mat.E = 1;
mat.nu = 0.25;

[M, K] = assemble_system_matrices(ELEM, XYZ, mat, 'STRUCT');


%%
opts.spdB= 1;
opts.tol = 1e-16;


[V,L]=eigs(-K,M,40,'sm',opts);
V = real(V);
l = diag(L);
feig = real(sqrt(-l)/2/pi);


[~,idx]=min(abs(feig-(mat.c/2/Lz)));
disp([mat.c/2/Lz feig(idx)])

%%
animate_mode(ELEM,coord,V(:,idx),1);
return
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
S = XYZ(ELEM(BCelem(1),:),:);
S = S(5:8,:)
[Fout] = SurfaceInt (S,n(1,3)*1/areas(1))

