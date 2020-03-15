function  powerAllocatePuser = get_RR_MAX_MIN_procedure(powerAllocateX,Aql0,Aml,bm,K,MrPower,subbandNumber,Mt)

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

powerAllocatePuser = zeros(subbandNumber,K);

for iUser = 1:K

    R = rank(powerAllocateX(:,:,iUser));
    Rl = rank(powerAllocateX(:,:,iUser),1e-2);
    mTerm = 0;
    Loop =1;
    while(Rl>1)
        mTerm = mTerm+1;
   
        [Vx, Dx] = eig(powerAllocateX(:,:,iUser));
        %diag(Dx)
        Vdecompose = Vx * Dx ^(1/2);
        ax =zeros(subbandNumber,subbandNumber);
        for iK = 1:K-1
            
           ax = ax + Vdecompose' *Aml(:,:,iK) * Vdecompose;
        end
       % bvec = null(Vdecompose' *Aml(:,:,1) * Vdecompose);
       
         bvec = null(ax);
        deltaL = bvec * bvec';
        
        [~,Dd] = eig(deltaL);
    
        maxEig = max(diag(Dd));
        
        powerAllocateX(:,:,iUser)= Vdecompose * (eye(R)-1/maxEig * deltaL ) * Vdecompose';
         
        % Rl = rank(powerAllocateX,1e-11);
        [Vx, Dx] = eig(powerAllocateX(:,:,iUser));
        %diag(Dx)
   
        R1 = length(find(diag(Dx)>1e-2));
            if R1 == 1
                break;
            end
        % Rl = rank(powerAllocateX,1e-5)
    
        if mTerm>1000
            break;
        end
    end
    powerAllocatePuser(:,iUser) = powerAllocateX(:,1,iUser) * sqrt(trace(powerAllocateX(:,:,iUser)))/norm(powerAllocateX(:,1,iUser),'fro');
end
%powerAllocateXstar = powerAllocateX;

%xstar = powerAllocateX(:,1) * sqrt(trace(powerAllocateX))/norm(powerAllocateX(:,1));
end