function [powerAllocateXstar,xstar] = get_RR_procedure(powerAllocateX,Aql0,Aml,bm,K,MrPower,subbandNumber,Mt)

%cvx_begin sdp quiet
%    cvx_solver Mosek
%    variable  powerAllocateX(subbandNumber * Mt,subbandNumber * Mt) hermitian semidefinite
%    minimize trace(Aql0 * powerAllocateX)
%    subject to
%    for iUser = 1:K-1
%       real( trace(Aml(:,:,iUser) * powerAllocateX))  <= real(bm(iUser));
%    end
%    trace(powerAllocateX) <= MrPower;
%cvx_end





R = rank(powerAllocateX);
Rl = R;
mTerm = 0;
Loop =1;
while(Loop)
    mTerm = mTerm+1;
    
    [Vx, Dx] = eig(powerAllocateX);
    %diag(Dx)
    Vdecompose = Vx * Dx ^(1/2);
    %Dx
    %cvx_begin quiet
    %    cvx_solver Mosek
    %        variable deltaL(R,R) hermitian 
    %        variable gamma3 nonnegative
    %        minimize gamma3 
    %        subject to 
    %        gamma3 = norm(Vdecompose * deltaL(R,R),'fro');
    %        norm(Vdecompose * deltaL(R,R),'fro') = eps; 
    %    cvx_end
    %decomposeEle = zeros(R,R);
    %for iUser  = 1:K-1
    %    decomposeEle = decomposeEle + Vdecompose' *Aml(:,:,iUser) * Vdecompose;
    %end
    bvec = null(Vdecompose' *Aml(:,:,1) * Vdecompose);
    
    deltaL = bvec * bvec';
    %deltaL = linsolve(Vdecompose,ones(R,R) * eps * 1e-10);    
    %deltaL = linsolve(Vdecompose' * Aml(:,:,1) * Vdecompose,ones(R,R) /1e10);
    [~,Dd] = eig(deltaL);
    
    maxEig = max(diag(Dd));
    powerAllocateX= Vdecompose * (eye(R)-1/maxEig * deltaL ) * Vdecompose';
   % Rl = rank(powerAllocateX,1e-11);
   [~, Dx] = eig(powerAllocateX);
   %diag(Dx)
   
   R1 = length(find(diag(Dx)>1e-9));
        if R1 == 1
            break;
        end
  % Rl = rank(powerAllocateX,1e-5)
    
   if mTerm>1000
      break;
   end
end
powerAllocateXstar = powerAllocateX;

xstar = powerAllocateX(:,1) * sqrt(trace(powerAllocateX))/norm(powerAllocateX(:,1));
end