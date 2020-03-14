function Vout =  get_Vout_ASS(ChannelInfo,TransceiverInfo)





K = TransceiverInfo.K;
MrPower = TransceiverInfo.MrPower;
subbandNumber = ChannelInfo.subbandNumber;
GainFre = ChannelInfo.subbandChannelGainFre;
hnorm = zeros(subbandNumber ,K);
b2 = TransceiverInfo.b2;
b4 = TransceiverInfo.b4;

for iUser = 1:K
    for iSubbandNumber = 1:subbandNumber
        hnorm(iSubbandNumber,iUser) = norm(GainFre(:,iSubbandNumber,iUser),'fro');
      % htotal((iSubbandNumber-1) * Mt +1 : (iSubbandNumber-1) * Mt +Mt,iUser) = GainFre(:,iSubbandNumber,iUser);
    end
    [hmax(iUser),index] = max(hnorm(:,iUser));
end
s = conj(GainFre(:,index,1))/norm(GainFre(:,index,1),'fro') * sqrt(MrPower);

conj(GainFre(:,index,1))
Vout = b2 * s' * conj(GainFre(:,index,1)) * GainFre(:,index,1).' * s + 1.5 * b4 * s' * conj(GainFre(:,index,1)) * GainFre(:,index,1).' * s * s' * conj(GainFre(:,index,1)) * GainFre(:,index,1).' * s;
%Vout = b2 * s' * conj(GainFre(:,index,1)) * GainFre(:,index,1).' * s + 1.5 * b4 * s' * conj(GainFre(:,index,1)) * ...
%    GainFre(:,index,1).'*s' * conj(GainFre(:,index,1)) * GainFre(:,index,1).';
%
    
end

