function [Rfsi] = FSI_coupling_matrix(ELEM, XYZ, NNs, NNf, s)

%
%          Y                    
%          |                    
%      4 ----- 3                
%      |   |   |               
%      |   +------ X          
%      |       |       
%      1 ----- 2            
%
% 4 DoF nodes , UX, UY, UZ, P

XYZs = NaN(size(XYZ));
XYZs(NNs,:) = XYZ(NNs,:);

XYZf = NaN(size(XYZ));
XYZf(NNf,:) = XYZ(NNf,:);

% structural part
NDOF = size(XYZs,2);
nd = length(NNs);

en = arrayfun(@(i) find(ELEM - NNs(i)==0),1:nd,'UniformOutput',false);
ind = [];
for ii = 1:length(en)
    [i,~]=ind2sub(size(ELEM),en{ii});  
    ind = [ind ;i];
end
struc_elem = unique(ind);

ELEMb = ELEM(struc_elem,:);

[n, ~, faces] = getsurfacenormals(ELEMb,XYZs);


% fluid part
en = arrayfun(@(i) find(ELEM - NNf(i)==0),1:nd,'UniformOutput',false);
ind = [];
for ii = 1:length(en)
    [i,~]=ind2sub(size(ELEM),en{ii});  
    ind = [ind ;i];
end

fluid_elem = unique(ind);
ELEMbf= ELEM(fluid_elem,:);


% Get correct normal for FSI, i.e. from structural side towards fluid side
n_dir = (mean(XYZ(ELEM(struc_elem,:),:))-mean(XYZ(ELEM(fluid_elem,:),:)))*n(1,:).';

if n_dir>0
    nf = n.';
else
    nf = -n.';
end

[NF,NPF] = size(faces);

if NPF == 4
    IP.XI = [-1  -1 ;...
              1  -1 ;...
              1   1 ;...
             -1   1  ]/sqrt(3);
    IP.WT = ones(1,4);
end

% R : Pressure from fluid nodes to structural 3-DoF forces

R = zeros(NPF*NDOF,NPF, NF);

for Nn = 1: NF
    for IN = 1:numel(IP.WT)
        XI = IP.XI(IN,:);
        WT = IP.WT(IN);
        if NPF == 4
            [N,DN] = shape2D4(XI);
            J      = [1 1 1;(DN*XYZs(faces(Nn,:),:))];
        end
        dJ   = det(J);
        R(:,:,Nn) = R(:,:,Nn) + kron(N.',eye(3))*(nf(:,Nn)*N)*dJ*WT;
    end
end

% Ordering for X rows, Y colums and R forces
% % Face 1 nodes: 
% 4x3 Uxyz - p1
% 4x3 Uxyz - p2
% 4x3 Uxyz - p3
% 4x3 Uxyz - p4
% Face 2 nodes: 
% 4x3 Uxyz - p1
% 4x3 Uxyz - p2
% 4x3 Uxyz - p3
% 4x3 Uxyz - p4

R_s_idx=repmat(NDOF*faces(:,kron(1:NPF,ones(1,NDOF)))-kron(ones(NF,1),kron(ones(1,NPF),(NDOF-1):-1:0)),1,4).';

X =R_s_idx(:) + s.off(1);


%%
[~, ~ , faces] = getsurfacenormals(ELEMbf,XYZf);

R_f_idx=faces + s.off(2);
Y = R_f_idx(:,kron(1:NPF,ones(1,NPF*NDOF))).';

Rfsi = sparse(X(:),Y(:),R(:),s.tot,s.tot);
