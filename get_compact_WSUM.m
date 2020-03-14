function InitialM = get_compact_WSUM(TransceiverInfo, ChannelInfo)

Mt = TransceiverInfo.Mt;
K = TransceiverInfo.K;
MrPower = TransceiverInfo.MrPower;
subbandNumber = ChannelInfo.subbandNumber;
GainFre = ChannelInfo.subbandChannelGainFre;
htotal = zeros(subbandNumber * Mt,K);
Mpp = zeros(subbandNumber *Mt,subbandNumber *Mt,K);
Mdiag = zeros(subbandNumber*Mt,subbandNumber*Mt,subbandNumber,K);
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
       % hnorm(iSubbandNumber,iUser) = norm(GainFre(:,iSubbandNumber,iUser),'fro');
       htotal((iSubbandNumber-1) * Mt +1 : (iSubbandNumber-1) * Mt +Mt,iUser) = GainFre(:,iSubbandNumber,iUser);
    end
    Mpp(:,:,iUser) = conj(htotal(:,iUser)) *  htotal(:,iUser).';
end

powerAllocateP = ones(subbandNumber *Mt,1) * sqrt(MrPower/subbandNumber/Mt); 
powerAllocateX = powerAllocateP * powerAllocateP';
tk = zeros(subbandNumber,K);
tkstar = zeros(subbandNumber-1,K);
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
for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        %for iMt = 1:Mt
       % Mdiag(:,:, iSubbandNumber,iUser) = diag(diag(Mpp(:,:,iUser),(iSubbandNumber-1) * Mt + Mt-1),iSubbandNumber-1);
        tk(iSubbandNumber,iUser) = trace(Mdiag(:,:, iSubbandNumber,iUser) * powerAllocateX);
        if iSubbandNumber ~= 1
        tkstar(iSubbandNumber-1,iUser) = trace( Mdiag(:,:,iSubbandNumber,iUser)' * powerAllocateX);
        end
    end
    
end

InitialM.htotal = htotal;
InitialM.Mpp = Mpp;
InitialM.Mdiag = Mdiag;
InitialM.powerAllocateP = powerAllocateP;
InitialM.tk = tk;
InitialM.tkstar = tkstar;
InitialM.powerAllocateX = powerAllocateX;  