function [sprecoder,VoutUser,Vout] =  get_Vout_MAX_MIN_RR(ChannelInfo,TransceiverInfo,InitialM)






GainFre = ChannelInfo.subbandChannelGainFre;
Mt = TransceiverInfo.Mt;
Tolerance = TransceiverInfo.Tolerance;
MrPower = TransceiverInfo.MrPower;
K =  TransceiverInfo.K;
b2 = TransceiverInfo.b2;
b4 = TransceiverInfo.b4;
Mdiag = InitialM.Mdiag;
subbandNumber = ChannelInfo.subbandNumber;
powerAllocateXlast = InitialM.powerAllocateX;
tk = InitialM.tk;
tkstar = InitialM.tkstar;
Loop = 1;
wq = 1;

A0 = diag([0.5 ones(1,subbandNumber-1)] * (-3 * b4));

while(Loop) 
    %update cqbar Cq1,Aq1
    gammaGroup = zeros(K,1);
    cqbar = zeros(K,1);
    Cq1 = zeros(subbandNumber * Mt,subbandNumber *  Mt,K);
    Aq1 = zeros(subbandNumber * Mt,subbandNumber *  Mt,K);
    for iUser  = 1:K
        cqbar(iUser) = -tk(:,iUser)' * A0 * tk(:,iUser);
         Cq1(:,:,iUser) =  -1/2 * (b2 + 3 * b4 * tk(1,iUser)) * Mdiag(:,:,1,iUser);
         for iSubbandNumber = 2:subbandNumber
             Cq1(:,:,iUser) = Cq1(:,:,iUser) - 3* b4 * tkstar(iSubbandNumber-1,iUser) * Mdiag(:,:,iSubbandNumber,iUser);
         end
          Aq1(:,:,iUser) = Cq1(:,:,iUser) + Cq1(:,:,iUser)';
    end
       
     
    cvx_begin sdp quiet
        cvx_solver Mosek
        variable  powerAllocateX(subbandNumber * Mt,subbandNumber * Mt)  hermitian semidefinite
        %variable  powerAllocateX(subbandNumber * Mt,subbandNumber * Mt) hermitian 
        variable gamma0  nonnegative
        minimize -gamma0
        subject to 
            for iUser = 1:K
               real(trace(Aq1(:,:,iUser) * powerAllocateX) + cqbar(iUser)) + gamma0 <= 0;
                
            end
           % real(trace(Aq1(:,:,1) * powerAllocateX) + cqbar(1)) + gamma0 <= 0;
           % real(trace(Aq1(:,:,2) * powerAllocateX) + cqbar(2)) + gamma0 <= 0;
          %  real(trace(Aq1(:,:,3) * powerAllocateX) + cqbar(3) + gamma0 <= 0;
            trace( powerAllocateX)<= MrPower;
    cvx_end
    
    for iUser = 1:K
        gammaGroup(iUser)  =   trace(Aq1(:,:,iUser) * powerAllocateX) + cqbar(iUser);
    end
   
    [~,userIndex] = max(gammaGroup);
    nTerm = 0;
    
   
    Aml = zeros(subbandNumber * Mt,subbandNumber *  Mt,K-1);
    bm = zeros(K-1);
    for iUser = 1:K
        if iUser ~= userIndex
            nTerm =nTerm +1;
            Aml(:,:,nTerm) = Aq1(:,:,iUser) - Aq1(:,:,userIndex);
            bm (nTerm) = cqbar(iUser) - cqbar(userIndex) ; 
        end
    end
        %RR algoritm
     % powerAllocateXstar =  powerAllocateX ;
     % xstar = powerAllocateX(:,1) * sqrt(trace(powerAllocateX))/norm(powerAllocateX(:,1));
    [powerAllocateXstar,xstar] = get_RR_procedure(powerAllocateX,Aq1(:,:,userIndex),Aml,bm,K,MrPower,subbandNumber,Mt);  
       
     
    %Update tk tkstar
    tk = zeros(subbandNumber,K);
    tkstar = zeros(subbandNumber-1,K);
    for iUser = 1:K
        for iSubbandNumber = 1:subbandNumber
            %Mdiag(:,:, iSubbandNumber,iUser) = diag(diag(Mpp(:,:,iUser),iSubbandNumber-1),iSubbandNumber-1);
            tk(iSubbandNumber,iUser) = trace(Mdiag(:,:, iSubbandNumber,iUser) * powerAllocateXstar);
            if iSubbandNumber ~= 1
                tkstar(iSubbandNumber-1,iUser) = trace( Mdiag(:,:,iSubbandNumber,iUser)' * powerAllocateXstar);
            end
        end
    
    end
    %  norm((powerAllocateXstar-powerAllocateXlast),'fro')/norm(powerAllocateXstar,'fro')
     if (norm((powerAllocateXstar-powerAllocateXlast),'fro')/norm(powerAllocateXstar,'fro')) <= Tolerance/10 || Loop>=50
        
        break;
    else
         Loop = Loop+1;
        
        powerAllocateXlast = powerAllocateXstar;
    end
end
    pstar = xstar;
    Vout = 0;
    for iUser = 1:K
        Vout = Vout+b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        VoutUser(iUser) = b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        for iSubbandNumber = 1:subbandNumber
           % sprecoder((iSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+Mt,1) = pstar( iSubbandNumber) * GainFre(:,iSubbandNumber,1).'/norm(GainFre(:,iSubbandNumber,1));
            sprecoder = xstar;
        if iSubbandNumber ~=1
            Vout = Vout + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
            VoutUser(iUser) =  VoutUser(iUser) + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
        end
        end
    end
   
end





