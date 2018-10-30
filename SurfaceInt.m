function [Fout] = SurfaceInt (XYZ,Fin)

%
%          Y                    
%          |                    
%      4 ----- 3                
%      |   |   |               
%      |   +------ X          
%      |       |       
%      1 ----- 2            
%

NPE = size(XYZ,1);


if NPE == 4
    IP.XI = [-1  -1 ;...
              1  -1 ;...
              1   1 ;...
             -1   1  ]/sqrt(3);
    IP.WT = ones(1,4);
end
    
Fout = zeros(NPE,1);
for IN = 1:numel(IP.WT)
    XI = IP.XI(IN,:);
    WT = IP.WT(IN);
    if NPE == 4
        [N,DN] = shape2D4(XI);
        J      = [1 1 1;(DN*XYZ)];
    end
    dJ   = det(J);
    Fout = Fout + N.'*N*Fin*dJ*WT;
end
