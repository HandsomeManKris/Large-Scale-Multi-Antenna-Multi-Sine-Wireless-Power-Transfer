function InitialM = get_compact_CHE_WSUM(TransceiverInfo, ChannelInfo)






wq = 1;
Mt = TransceiverInfo.Mt;
K = TransceiverInfo.K;
MrPower = TransceiverInfo.MrPower;
subbandNumber = ChannelInfo.subbandNumber;
GainFre = ChannelInfo.subbandChannelGainFre;
Mpp = zeros(subbandNumber,subbandNumber);
Mdiag = zeros(subbandNumber,subbandNumber,subbandNumber);
CHEMatrix = zeros(subbandNumber, K); 

%get channel hardening coefficient
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        CHEMatrix(iSubbandNumber,iUser) = GainFre(:,iSubbandNumber,iUser).' * conj(GainFre(:,iSubbandNumber,iUser))/Mt;
       
    end
end
%CHE = mean(mean(CHEMatrix)) ;%?????
%CHE = min(mean(CHEMatrix)) ; %??
CHE = mean(min(CHEMatrix)) ;  %??
CHEdiag = diag(ones(K * subbandNumber,1)) * CHE ;
powerAllocateP = ones(subbandNumber * K,1)  * sqrt(1 / (CHE * subbandNumber * K)); 
powerAllocatePuser = zeros(subbandNumber,K);
for iUser = 1:K
 %   for iSubbandNumber = 1:subbandNumber
        powerAllocatePuser(1 :subbandNumber ,iUser) = powerAllocateP((iUser-1)* subbandNumber+1:(iUser-1)* subbandNumber+subbandNumber,1);
    
end
Mpp = ones(subbandNumber , subbandNumber);

%powerAllocateX = powerAllocateP * powerAllocateP';
tk = zeros(subbandNumber,K);
tkstar = zeros(subbandNumber-1,K);
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        Mdiag(:,:, iSubbandNumber) = diag(diag(Mpp,iSubbandNumber-1),iSubbandNumber-1);
        tk(iSubbandNumber,iUser) = powerAllocatePuser( : ,iUser)' * Mdiag(:,:, iSubbandNumber) * powerAllocatePuser( : ,iUser);
        if iSubbandNumber ~= 1
            tkstar(iSubbandNumber-1,iUser) = powerAllocatePuser( : ,iUser)' * Mdiag(:,:, iSubbandNumber)' * powerAllocatePuser( : ,iUser);
        end
    end
    
end
InitialM.CHE = CHE;
InitialM.CHEdiag = CHEdiag;
%InitialM.heq = heq;
InitialM.Mpp = Mpp;
InitialM.Mdiag = Mdiag;
InitialM.powerAllocateP = powerAllocateP;
InitialM.tk = tk;
InitialM.tkstar = tkstar;
InitialM.powerAllocatePuser = powerAllocatePuser;  