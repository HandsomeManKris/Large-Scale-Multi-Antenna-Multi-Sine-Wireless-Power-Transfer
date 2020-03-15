function [sprecoder,VoutUser,Vout] =  get_Vout_CHE_MAX_MIN_RR(ChannelInfo,TransceiverInfo,InitialM)









GainFre = ChannelInfo.subbandChannelGainFre;
Mt = TransceiverInfo.Mt;
Tolerance = TransceiverInfo.Tolerance;
MrPower = TransceiverInfo.MrPower;
K =  TransceiverInfo.K;
b2 = TransceiverInfo.b2;
b4 = TransceiverInfo.b4;
Mdiag = InitialM.Mdiag;
subbandNumber = ChannelInfo.subbandNumber;
powerAllocateP = InitialM.powerAllocateP;
powerAllocatePuser = InitialM.powerAllocatePuser ;
CHEMatrix = InitialM.CHEMatrix;
CHEdiag = InitialM.CHEdiag;
tk = InitialM.tk;
tkstar = InitialM.tkstar;
Loop = 1;
wq = 1;
A0 = diag([0.5 ones(1,subbandNumber-1)] * (-3 * b4));
CHE = zeros(1,K);
for iUser = 1:K
    CHE(iUser) = mean(CHEMatrix(:,iUser));
end

while(Loop)
   littleA = cvx(0);
    %computeCq1 Aq1 
   cqbar = zeros(K,1);
   %computeCpp1 App1 Ap1
    App1 = zeros(subbandNumber ,subbandNumber ,K);
    Cpp1 = zeros(subbandNumber ,subbandNumber ,K);
   % Ap1 = zeros(subbandNumber *K,subbandNumber * K);
    CHEbar = zeros(subbandNumber * K,subbandNumber * K);
    for iUser = 1: K 
        cqbar(iUser) = -tk(:,iUser)' * A0 * tk(:,iUser) * MrPower^2 * CHE(iUser)^4;
        Cpp1(:,:,iUser) = Cpp1(:,:,iUser)-1/2 * (b2  * MrPower * CHE(iUser)^2 + 3 *  MrPower^2 * CHE(iUser)^4 * b4 * tk(1,iUser)) * Mdiag(:,:,1);
        for iSubbandNumber = 2:subbandNumber
            Cpp1(:,:,iUser) = Cpp1(:,:,iUser)- 3* b4 *  MrPower^2 * CHE(iUser)^4 * tkstar(iSubbandNumber-1,iUser) * Mdiag(:,:,iSubbandNumber);
        end
        App1(:,:,iUser) = Cpp1(:,:,iUser) + Cpp1(:,:,iUser)';
       % Ap1((iUser-1)*subbandNumber +1 : (iUser-1)*subbandNumber +subbandNumber,(iUser-1)*subbandNumber +1 : (iUser-1)*subbandNumber +subbandNumber) = App1(:,:,iUser);
    end
   %Compute CHEbar
   cvx_begin sdp quiet
        cvx_solver Mosek
                  
            variable powerAllocateXstar(subbandNumber , subbandNumber,K) hermitian semidefinite
        
      %  variable totalPower nonnegative
        variable gamma2p nonnegative
        minimize -gamma2p
        subject to
            for iUser = 1:K
                real(  trace(App1(:,:,iUser) * powerAllocateXstar(:,:,iUser) ) + cqbar(iUser)) + gamma2p <= 0;
                littleA  = littleA + real(trace( CHE(iUser) * eye(subbandNumber) .* powerAllocateXstar(:,:,iUser)));
            end
            littleA == 1;
   cvx_end
    %Update X
    [V,D] = eig(CHEbar);
    [~,index] = min(diag(D));
     b = (1* inv(V(:,index)' * CHEdiag * V(:,index)))^(0.5);
    pstar = b * V(:,index);
   % powerAllocateXstar = xstar * xstar';
    %Update tk tkstar
    tk = zeros(subbandNumber,K);
    tkstar = zeros(subbandNumber-1,K);
    for iUser = 1:K
 %   for iSubbandNumber = 1:subbandNumber
        powerAllocatePuser(: ,iUser) =  pstar((iUser-1)* subbandNumber+1:(iUser-1)* subbandNumber+subbandNumber,1);
    
    end
    for iUser = 1:K
        for iSubbandNumber = 1:subbandNumber
           % Mdiag(:,:, iSubbandNumber) = diag(diag(Mpp,iSubbandNumber-1),iSubbandNumber-1);
            tk(iSubbandNumber,iUser) = powerAllocatePuser( : ,iUser)' * Mdiag(:,:, iSubbandNumber) * powerAllocatePuser( : ,iUser);
            if iSubbandNumber ~= 1
                tkstar(iSubbandNumber-1,iUser) = powerAllocatePuser( : ,iUser)' * Mdiag(:,:, iSubbandNumber)' * powerAllocatePuser( : ,iUser);
            end
        end
    
    end
    if (norm((pstar-powerAllocateP),'fro')/norm(pstar,'fro')) <= Tolerance/10 || Loop>=200
        
        break;
    else
         Loop = Loop+1;
        
        powerAllocateP = pstar;
    end
end
    
    Vout = 0;
    for iUser = 1:K
        VoutUser(iUser) = b2 *  MrPower * CHE(iUser)^2 * powerAllocatePuser(: ,iUser)' *  Mdiag(:,:,1) *  powerAllocatePuser(: ,iUser) + 1.5 * b4 * MrPower^2 * CHE(iUser)^4 * norm( powerAllocatePuser(: ,iUser)' * Mdiag(:,:,1) *  powerAllocatePuser(: ,iUser))^2;
        Vout = Vout+b2 *  MrPower * CHE(iUser)^2 * powerAllocatePuser(: ,iUser)' *  Mdiag(:,:,1) *  powerAllocatePuser(: ,iUser) + 1.5 * b4 * MrPower^2 * CHE(iUser)^4 * norm( powerAllocatePuser(: ,iUser)' * Mdiag(:,:,1) *  powerAllocatePuser(: ,iUser))^2;
        for iSubbandNumber = 1:subbandNumber
           % sprecoder((iSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+Mt,1) = pstar( iSubbandNumber) * GainFre(:,iSubbandNumber,1).'/norm(GainFre(:,iSubbandNumber,1));
           % sprecoder = xstar;
           
        if iSubbandNumber ~=1
            Vout = Vout + 3 * b4 * MrPower^2 * CHE(iUser)^4 * norm( powerAllocatePuser(: ,iUser)' * Mdiag(:,:,iSubbandNumber) *  powerAllocatePuser(: ,iUser))^2;
             VoutUser(iUser) = VoutUser(iUser) +  3 * b4 * MrPower^2 * CHE(iUser)^4 * norm( powerAllocatePuser(: ,iUser)' * Mdiag(:,:,iSubbandNumber) *  powerAllocatePuser(: ,iUser))^2;
        end
        end
    end
     sprecoder = zeros(subbandNumber *Mt ,1 );
    for iSubbandNumber = 1:subbandNumber
        for iUser = 1:K
            sprecoder((iSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+Mt,1) = sprecoder((iSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+Mt,1) + ...
            powerAllocatePuser(iSubbandNumber ,iUser)* conj(GainFre(:,iSubbandNumber,iUser))/sqrt(Mt);
        end
    end
    sprecoder = sqrt(MrPower/Mt) * sprecoder/norm(sprecoder,'fro'); 
        
        
        
end