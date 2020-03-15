function [sprecoder,pstar,Vout] =  get_Vout_singleuser(ChannelInfo,TransceiverInfo,InitialM)









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
wq = 1/K;

while(Loop)
   
    %computeCpp1
    Cpp1 = zeros(subbandNumber,subbandNumber);
    for iUser = 1: K 
        Cpp1 = Cpp1-1/2 * (b2 + 3 * b4 * tk(1,iUser)) * Mdiag(:,:,1,iUser);
        for iSubbandNumber = 2:subbandNumber
            Cpp1 = Cpp1- 3* b4 * tkstar(iSubbandNumber-1,iUser) * Mdiag(:,:,iSubbandNumber,iUser);
        end
    end
    Cpp1 = Cpp1 * wq;
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
    if (norm((powerAllocateXstar-powerAllocateX),'fro')/norm(powerAllocateXstar,'fro')) <= Tolerance/10 || Loop>=200
        
        break;
    else
         Loop = Loop+1;
        
        powerAllocateX = powerAllocateXstar;
    end
end
    pstar = xstar;
    Vout = 0;
    Vout = b2 * pstar' * Mdiag(:,:,1,1) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,1) * pstar)^2;
for iSubbandNumber = 1:subbandNumber
    sprecoder((iSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+Mt,1) = pstar( iSubbandNumber) * GainFre(:,iSubbandNumber,1).'/norm(GainFre(:,iSubbandNumber,1));
    
    if iSubbandNumber ~=1
        Vout = Vout + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,1) * pstar)^2;
     end
end
end
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        