function [sprecoder,VoutUser,Vout] =  get_Vout_WSUM_RateRegion(ChannelInfo,TransceiverInfo,InitialM)





w1 = 0:0.1:2;
w2 = 2-w1;
wq(:,1) = w1;
wq(:,2) = w2;
for iCase = 1: length(w1)
    
  

    GainFre = ChannelInfo.subbandChannelGainFre;
    Mt = TransceiverInfo.Mt;
    Tolerance = TransceiverInfo.Tolerance;
    MrPower = TransceiverInfo.MrPower;
    K =  TransceiverInfo.K;
    b2 = TransceiverInfo.b2;
    b4 = TransceiverInfo.b4;
    Mdiag = InitialM.Mdiag;
    subbandNumber = ChannelInfo.subbandNumber;
    powerAllocateX = InitialM.powerAllocateX;
    tk = InitialM.tk;
    tkstar = InitialM.tkstar;
    Loop = 1;



    while(Loop)
   
    %computeCpp1
    App1 = zeros(subbandNumber * Mt,subbandNumber * Mt);
    Cpp1 = zeros(subbandNumber * Mt,subbandNumber * Mt,K);
    for iUser = 1: K 
        Cpp1(:,:,iUser) = Cpp1(:,:,iUser)-1/2 * (b2 + 3 * b4 * tk(1,iUser)) * Mdiag(:,:,1,iUser);
        for iSubbandNumber = 2:subbandNumber
            Cpp1(:,:,iUser) = Cpp1(:,:,iUser)- 3* b4 * tkstar(iSubbandNumber-1,iUser) * Mdiag(:,:,iSubbandNumber,iUser);
        end
        Cpp1(:,:,iUser) = Cpp1(:,:,iUser) * wq(iCase,iUser);
    end
    Cpp1 = Cpp1(:,:,1) + Cpp1(:,:,2);
    App1 = Cpp1 +Cpp1';
    %Update X
    [V,D] = eig(App1);
    [~,index] = min(diag(D));
   
    xstar = sqrt(MrPower ) * V(:,index);
    powerAllocateXstar = xstar * xstar';
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
    if (norm((powerAllocateXstar-powerAllocateX),'fro')/norm(powerAllocateXstar,'fro')) <= Tolerance/1000 || Loop>=200
        
        break;
    else
         Loop = Loop+1;
        
        powerAllocateX = powerAllocateXstar;
    end
    end
        pstar = xstar;
        Vout = 0;
        for iUser = 1:K
        Vout = Vout+b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        VoutUser(iCase,iUser) = b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        for iSubbandNumber = 1:subbandNumber
           % sprecoder((iSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+Mt,1) = pstar( iSubbandNumber) * GainFre(:,iSubbandNumber,1).'/norm(GainFre(:,iSubbandNumber,1));
            sprecoder = xstar;
        if iSubbandNumber ~=1
            Vout = Vout + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
            VoutUser(iCase,iUser) =  VoutUser(iCase,iUser) + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
        end
        end
        end
   
end
end