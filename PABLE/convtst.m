function [neig,ritz,ipritz,ritzvalue,VL,VR,resl,resr] = ...
           convtst(fulldual,semidual,nb,szoft,T,R,S,anorm,tolconv,tolclust)
%
% Purpose 
% =======
% 
% Compute the eigenvalue decomposition of T_j, and test for converengece
%
% Arguments
% =========
%
% fulldual  (input/output) logical 
%
% semidual  (input/output) logical 
%
% nb        (input/output) integer 
%
% szoft     (input/output) integer 
%
% T         (input/output) double precision array, dimension(szoft, szoft) 
%
% R, S      (input/output) double precision array, dimension(n, nb)
%
% anrom     (input/output) double precision 
%
% tolconv   (input/output) double precision 
%
% tolclust  (input/output) double precision 
%
% neig      (output) integer 
%
% ritz      (output) double precision array, dimension( neig )
%
% ipritz    (output) integer array, dimension( neig )
%
% ritzvalue (output) double precision array, dimension( szoft )  
%
% VL, VR    (output) double precision arrays, dimension( szoft,szoft ) 
%
% resl,resr (output) double precision array, dimension( szoft )  
% 
% ==========================================================================
%
% Compute the eigen-decomposition of T_j 
%
[VR,teig] = eig(T(1:szoft,1:szoft));  % T*VR = VR*teig, VR has been normalized
ritzvalue = diag( teig );
VL = inv(VR)';                        % left eigenvectors of T;
for i = 1:szoft                    
    VL(:,i) = VL(:,i)/norm(VL(:,i));  % normalize the left eigenvectors
end 
%
% Compute the residual norms
%     
for k=1:szoft
%  resl(k) = norm( conj(R)*conj(VL(szoft-nb+1:szoft, k)), 1)/anorm;% left resi 
%  resr(k) = norm( S*VR(szoft-nb+1:szoft, k),1)/anorm; % right resi
  resl(k) = norm( conj(R)*conj(VL(szoft-nb+1:szoft, k)), 1);  % left res 
  resr(k) = norm( S*VR(szoft-nb+1:szoft, k),1);               % right res
end
%
% Test for convergence 
%
ritz = [ ]; ipritz = [ ]; neig = 0; 
if fulldual == 1 || semidual == 1
%
%  if Level 2 and up, compute relative gaps of RR values 
%
   gap = ones(szoft,1)*max( abs( ritzvalue )  );
   for k = 1:szoft
       for i = k+1:szoft
           dist = abs( ritzvalue( i ) - ritzvalue( k ) )/...
                       abs( ritzvalue( k ) );
%
           if nb == 1  % for the unblocked case 
               if dist < gap( k ), gap( k ) = dist; end
               if dist < gap( i ), gap( i ) = dist; end
%
           else   % for the blocked case, gap between the clusters. 
%
               if dist > tolclust  
                  if dist < gap( k ), gap( k ) = dist; end
                  if dist < gap( i ), gap( i ) = dist; end
               end                  
           end 
%
       end
   end
%
   reslr = zeros(szoft,1);
   for k=1:szoft
       reslr(k)= resr(k)*resl(k)/gap(k); 
       reslr(k)= max( min([resl(k), resr(k), reslr(k)]), eps);
       if reslr(k) <= tolconv
          neig = neig+1;
          ritz = [ritz; ritzvalue(k)]; 
          ipritz = [ipritz; k];  
       end
   end
   
%
else
%
%  Level 1 ( i.e., fulldual=0 and semidual =0 )
%
   for k=1:szoft
       resl(k) = max( resl(k), eps); 
       resr(k) = max( resr(k), eps);
       minres = max( min([resr(k),resl(k)]), eps);
       if minres <= tolconv
          neig = neig+1;
          ritz = [ritz; ritzvalue(k)];
          ipritz = [ipritz; k];
       end
   end
%
end
