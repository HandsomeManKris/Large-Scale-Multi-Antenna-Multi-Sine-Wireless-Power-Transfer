function InitialM = get_compact_singleuser(TransceiverInfo, ChannelInfo)




Mt = TransceiverInfo.Mt;
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

InitialM.hnorm = hnorm;
InitialM.Mpp = Mpp;
InitialM.Mdiag = Mdiag;
InitialM.powerAllocateP = powerAllocateP;
InitialM.tk = tk;
InitialM.tkstar = tkstar;
InitialM.powerAllocateX = powerAllocateX;    
        
%PowerWaveformAmplitude = sqrt( MrPower/subbandNumber) * Amplitude /norm(Amplitude,'fro');
%InitialWaveform.PowerWaveformAmplitude = PowerWaveformAmplitude;
%InitialWaveform.PowerRatio = PowerRatio;
%InitialWaveform.InforRatio = InforRatio;
%PowerRatio = 0.5;
%InforRatio = 1 - PowerRatio;
%get superposed waveform: \sum\limits_{n=0}^{N-1}s_{p,n}^2+s_{I,n}^2=2P
%For vector Frobenius Norm is just 2-norm, for matrix, Frobenius Norm is
%for matrix
%Phase can be 0 because of decoupling
%InforWaveformAmplitude = sqrt(2 * MtPower * InforRatio) * Amplitude /norm(Amplitude,'fro');
%InitialWaveform.InforWaveformAmplitude = InforWaveformAmplitude;