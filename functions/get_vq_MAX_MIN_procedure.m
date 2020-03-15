function powerAllocatePuser = get_vq_MAX_MIN_procedure(powerAllocateX,B1Q,Aml,CHE,K,MrPower,subbandNumber,Mt)





powerAllocatePuser = zeros(subbandNumber,K);

for iUser = 1:K
    
        
    eta = (binornd(1,0.5,[subbandNumber,1]) - 0.5) * 2;
    ita = (binornd(1,0.5,[subbandNumber,1]) - 0.5) * 2;
    
    [Vq,Dq] = eig(powerAllocateX(:,:,iUser)^(1/2) * B1Q(:,:,iUser) * powerAllocateX(:,:,iUser)^(1/2));
    Qq = Vq' * powerAllocateX(:,:,iUser)^(1/2) * CHE(iUser) * eye(subbandNumber) * powerAllocateX(:,:,iUser)^(1/2) * Vq;
    dota2 = trace( Qq .* eye(subbandNumber) );
    Qq1 = real(Qq);
    Qq2 = imag(Qq);
    if abs(eta.' *  Qq1 * eta-dota2 ) <=1e-6      
        powerAllocatePuser(:,iUser) = powerAllocateX(:,:,iUser)^(1/2) * Vq * eta;
        
    elseif abs(ita.' *  Qq1 * ita-dota2 ) <=1e-6   
        powerAllocatePuser(:,iUser) = powerAllocateX(:,:,iUser)^(1/2) * Vq * ita;
        
    else
        gamma0 = (eta.' * Qq2 * ita + sqrt((ita.' * Qq2 * eta)^2 - (ita.' * Qq1 * ita-dota2) * (eta.' * Qq1 * eta - dota2)))/(eta.' * Qq1 * eta-dota2);
        powerAllocatePuser(:,iUser) = powerAllocateX(:,:,iUser)^(1/2) * Vq * (gamma0/sqrt(1+gamma0^2) * eta + 1/sqrt(1+gamma0^2) * ita * 1i);
        
    end
%powerAllocateXstar = powerAllocateX;

%xstar = powerAllocateX(:,1) * sqrt(trace(powerAllocateX))/norm(powerAllocateX(:,1));
end
k=1;
end