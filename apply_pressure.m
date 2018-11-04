function [Fe, node_num, faces] = apply_pressure(ELEM, XYZ, NN, Fload)

XYZf = NaN(size(XYZ));
XYZf(NN,:) = XYZ(NN,:);

NDOF = size(XYZ,2);
nd = length(NN);

en = arrayfun(@(i) find(ELEM - NN(i)==0),1:nd,'UniformOutput',false);
ind = [];
for ii = 1:length(en)
    [i,~]=ind2sub(size(ELEM),en{ii});  
    ind = [ind ;i];
end

ELEMb = ELEM(unique(ind),:);
% show_mesh(ELEMb,XYZf)

[n, areas, faces] = getsurfacenormals(ELEMb, XYZf);
A = sum(areas);
[NFE,NPF] = size(faces);

Fout = zeros(NDOF,NPF,NFE);

for ii = 1 : NFE
    Fout(:,:,ii) = SurfaceInt(XYZ(faces(ii,:),:),[0;0;Fload/A]).';
end

Fidx = NDOF*faces(:,kron(1:NPF,ones(1,NDOF)))-kron(ones(NFE,1),kron(ones(1,NPF),(NDOF-1):-1:0));

Fidx = reshape(Fidx.' ,3,[]).';

Fi = Fidx.';
node_num = unique(faces);

Fe = sparse(Fi(:),1,Fout(:));
