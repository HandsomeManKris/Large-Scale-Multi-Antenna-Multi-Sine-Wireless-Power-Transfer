function [sprecoder,VoutUser,Vout] =  get_Vout_UP_MU(ChannelInfo,TransceiverInfo)

%%
Mt = TransceiverInfo.Mt;
b2 = TransceiverInfo.b2;
b4 = TransceiverInfo.b4;
K = TransceiverInfo.K;
MrPower = TransceiverInfo.MrPower;
subbandNumber = ChannelInfo.subbandNumber;
GainFre = ChannelInfo.subbandChannelGainFre;
sprecoder = zeros(subbandNumber * Mt,1);
Mpp = zeros(subbandNumber * Mt,subbandNumber *Mt ,K);
Mdiag = zeros(subbandNumber * Mt,subbandNumber *Mt,subbandNumber,K);
wFre = zeros(Mt,subbandNumber);
wFreNorm = zeros(subbandNumber,1);
htotal = zeros(subbandNumber * Mt,K);
%%
for iSubbandNumber = 1:subbandNumber
    for iUser = 1:K
        wFre(:,iSubbandNumber) = wFre(:,iSubbandNumber)+  conj(GainFre(:,iSubbandNumber,iUser)) / norm(GainFre(:,iSubbandNumber,iUser));
        
    end   
    wFreNorm(iSubbandNumber) = norm(wFre(:,iSubbandNumber),'fro')^2;
end

for iSubbandNumber = 1: subbandNumber
    sprecoder((iSubbandNumber-1) * Mt + 1:(iSubbandNumber-1) * Mt + Mt,1) = sqrt(MrPower) * wFre(:,iSubbandNumber)/sqrt(sum(wFreNorm));
end

for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
       
       htotal((iSubbandNumber-1) * Mt +1 : (iSubbandNumber-1) * Mt +Mt,iUser) = GainFre(:,iSubbandNumber,iUser);
    end
    Mpp(:,:,iUser) = conj(htotal(:,iUser)) *  htotal(:,iUser).';
end

%get M''
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
         Mdiag(:,:, iSubbandNumber,iUser) = zeros(subbandNumber *Mt,subbandNumber *Mt);
         for jSubbandNumber = 1:(subbandNumber + 1-iSubbandNumber)        
         Mdiag((jSubbandNumber-1)*Mt+1:(jSubbandNumber-1)*Mt+Mt,(iSubbandNumber-1)*Mt+(jSubbandNumber-1)*  Mt+1:(iSubbandNumber-1)*Mt+(jSubbandNumber-1)*Mt+Mt, iSubbandNumber,iUser)= ...
         Mpp((jSubbandNumber-1)*Mt+1:(jSubbandNumber-1)*Mt+Mt,(iSubbandNumber-1)*Mt+(jSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+(jSubbandNumber-1)*Mt+Mt,iUser);    
        end
    end
    
end
%powerAllocateP = ones(subbandNumber,1) * sqrt(MrPower/subbandNumber); 


    pstar = sprecoder;
    Vout = 0;
    VoutUser = zeros(K,1);
    for iUser = 1:K
        Vout = Vout+b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        VoutUser(iUser) = b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        for iSubbandNumber = 1:subbandNumber
           
        if iSubbandNumber ~=1
            Vout = Vout + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
            VoutUser(iUser) =  VoutUser(iUser) + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
        end
        end
    end
end