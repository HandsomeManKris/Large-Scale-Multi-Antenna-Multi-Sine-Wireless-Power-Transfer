function Vout =  get_Vout_UP(ChannelInfo,TransceiverInfo)




b2 = TransceiverInfo.b2;
b4 = TransceiverInfo.b4;
K = TransceiverInfo.K;
MrPower = TransceiverInfo.MrPower;
subbandNumber = ChannelInfo.subbandNumber;
GainFre = ChannelInfo.subbandChannelGainFre;
hnorm = zeros(subbandNumber,K);
Mpp = zeros(subbandNumber,subbandNumber,K);
Mdiag = zeros(subbandNumber,subbandNumber,subbandNumber,K);
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        hnorm(iSubbandNumber,iUser) = norm(GainFre(:,iSubbandNumber,iUser),'fro');
    end
    Mpp(:,:,iUser) = hnorm(:,iUser) *  hnorm(:,iUser).';
end

powerAllocateP = ones(subbandNumber,1) * sqrt(MrPower/subbandNumber); 


for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        Mdiag(:,:, iSubbandNumber,iUser) = diag(diag(Mpp(:,:,iUser),iSubbandNumber-1),iSubbandNumber-1);
        
    end
    
end
pstar = powerAllocateP;
Vout = 0;
Vout = b2 * pstar' * Mdiag(:,:,1,1) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,1) * pstar)^2;
for iSubbandNumber = 1:subbandNumber   
    if iSubbandNumber ~=1
       Vout = Vout + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,1) * pstar)^2;
    end
end
