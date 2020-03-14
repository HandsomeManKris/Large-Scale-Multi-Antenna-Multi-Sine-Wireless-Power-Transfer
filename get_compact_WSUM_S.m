function InitialM = get_compact_WSUM_S(TransceiverInfo, ChannelInfo)

wq = 1;
Mt = TransceiverInfo.Mt;
K = TransceiverInfo.K;
MrPower = TransceiverInfo.MrPower;
subbandNumber = ChannelInfo.subbandNumber;
GainFre = ChannelInfo.subbandChannelGainFre;
hFre = zeros(Mt,Mt,subbandNumber);
Mpp = zeros(subbandNumber,subbandNumber,K);
Mdiag = zeros(subbandNumber,subbandNumber,subbandNumber,K);
wbeamforming = zeros(Mt,subbandNumber);
heq = zeros(subbandNumber,K);
for iSubbandNumber = 1:subbandNumber
    for iUser = 1:K
        hFre(:,:,iSubbandNumber) = hFre(:,:,iSubbandNumber)+ wq * conj(GainFre(:,iSubbandNumber,iUser)) * GainFre(:,iSubbandNumber,iUser).';
    end
    [V1,D1] = eig( hFre(:,:,iSubbandNumber));
    [~,index] = max(diag(D1)); 
    wbeamforming(:,iSubbandNumber) = V1(:,index);
    
end
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        heq(iSubbandNumber,iUser) = wbeamforming(:,iSubbandNumber)' * conj(GainFre(:,iSubbandNumber,iUser));
    end
    Mpp(:,:,iUser) = heq(:,iUser) *  heq(:,iUser)';
end

    
powerAllocateP = ones(subbandNumber,1) * sqrt(MrPower/subbandNumber); 
powerAllocateX = powerAllocateP * powerAllocateP';
tk = zeros(subbandNumber,K);
tkstar = zeros(subbandNumber-1,K);
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        Mdiag(:,:, iSubbandNumber,iUser) = diag(diag(Mpp(:,:,iUser),iSubbandNumber-1),iSubbandNumber-1);
        tk(iSubbandNumber,iUser) = trace(Mdiag(:,:, iSubbandNumber,iUser) * powerAllocateX);
        if iSubbandNumber ~= 1
            tkstar(iSubbandNumber-1,iUser) = trace( Mdiag(:,:,iSubbandNumber,iUser)' * powerAllocateX);
        end
    end
    
end
InitialM.wbeamforming = wbeamforming;
InitialM.hFre = hFre;
InitialM.heq = heq;
InitialM.Mpp = Mpp;
InitialM.Mdiag = Mdiag;
InitialM.powerAllocateP = powerAllocateP;
InitialM.tk = tk;
InitialM.tkstar = tkstar;
InitialM.powerAllocateX = powerAllocateX;    