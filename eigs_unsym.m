function [V,D] = eigs_unsym(M,K,s0,k,n)

Nk = size(K,1);

q  = zeros(Nk,k);
qs = zeros(Nk,k);

%% Initialize first Krylov vector

qs(:,1) = -M*rand(Nk,1);
%% Generate k Krylov vectors
% Advance theKrylov subspace to get qs(:,i+1)
H = K-s0*M;
[L,U,P,Q,D] = lu(H);
     
for i = 1 : k
     qn = norm(qs(:,i),2);

     % Normalize qs(:,i) to obtain q(:,i)
     q(:,i) = qs(:,i)/qn;
     
     qs(:,i+1) = Q*(U\(L\(P*(D\q(:,i)))));
     qn = norm( qs(:,i+1) );
         
     qs(:,i+1) = qs(:,i+1)/qn;
     
     % Orthogonalize the new basis vector against the previous
     % Krylov basis vectors
     for j = 1 : i
        qs(:,i+1) = qs(:,i+1) - q(:, j)*( q(:,j)'*qs(:,i+1) );
     end

     %Apply reorthogonalization
     if norm( qs(:,i+1) ) < qn/sqrt(2)
         for j = 1 : i
            qs(:,i+1) = qs(:,i+1) - q(:, j)*( q(:,j)'*qs(:,i+1) );
        end
    end
     
     
end

qn = norm(qs(:,k),2);
q(:,k) = qs(:,k)/qn;

Ks = q.'*K*q;
Ms = q.'*M*q;

[Vs,D]=eigs(-Ks,Ms,n,'sm');
l = sqrt(-D)/2/pi;

disp(diag(l));
V = q*Vs;


%%
D = zeros(n,length(10:5:k));
l = zeros(n,length(10:5:k));
cnt = 1;
for j = 10:5:k
    Q = q(:,1:j);
    Kst = Q.'*K*Q;
    Mst = Q.'*M*Q; 
    D(:,cnt)=eigs(-Kst,Mst,n,'sm');
    l(:,cnt) = sqrt(-real(D(:,cnt)))/2/pi;
    
    cnt = cnt +1;
end

figure(100)
semilogy(15:5:k,abs(diff(real(l(2:6,:)),[],2)))
