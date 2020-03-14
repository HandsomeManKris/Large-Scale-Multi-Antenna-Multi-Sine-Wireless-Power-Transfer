function [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo)
% Function
% - get the channel gain for each subband including the amplitude and phase.
%
% - InputArg(s):TransceiverInfo:Mt,Mr ChannelInfo(s):centerFrequency,
% - subbandFrequency,tapDelay,tapGain,subbandNumber,subbandChannelGain
% -
% - OutputArgs(s):ChannelInfo:subbandChannelGain,subbandChannelGainAmplitude,subbandChannelGainPhase
% - Written by Kris 2020.2.12



K = TransceiverInfo.K;
Mt = TransceiverInfo.Mt;
%Mr = TransceiverInfo.Mr;
centerFrequency = ChannelInfo.centerFrequency;
subbandFrequency = ChannelInfo.subbandFrequency;
tapDelay = ChannelInfo.tapDelay;
tapGain = ChannelInfo.tapGain;
fadingType = ChannelInfo.fadingType;
subbandNumber = ChannelInfo.subbandNumber;
subbandChannelGain = zeros(subbandNumber*Mt,K); %h_{n,m}
subbandChannelGainFre = zeros(Mt,K);
%get fading type
%if strcmp(fadingType,'flat')
%    for iMt = 1: Mt
%        for iMr = 1:Mr
            %All subbands have same channel gain
%            subbandChannelGain(:,iMt,iMr) = repmat(sum(tapGain(:,iMt,iMr) .* exp(-1i * 2 * pi * centerFrequency * tapDelay)),subbandNumber,1);
%        end
%    end
%elseif strcmp(fadingType,'selective')
%    for iMt = 1:Mt
%        for iMr =1:Mr
%            for iSubbandNumber = 1:subbandNumber
%                subbandChannelGain(iSubbandNumber,iMt,iMr) = sum(tapGain(:,iMt,iMr) .* exp(-1i * 2 * pi * subbandFrequency(iSubbandNumber) * tapDelay));
%            end
%        end
%    end
%else
%    disp('No fadingType!!')
%end
%get fading type
if strcmp(fadingType,'flat')
    for iUser = 1:K
        for iSubbandNumber = 1:subbandNumber
            for iMt = 1: Mt                    
               %All subbands have same channel gain
                subbandChannelGain((iSubbandNumber-1) * Mt + iMt,iUser) = sum(tapGain(:,iMt,iUser) .* exp(-1i * 2 * pi * centerFrequency * tapDelay));
                
            end
                subbandChannelGainFre(:,iSubbandNumber,iUser) = subbandChannelGain((iSubbandNumber-1) * Mt + 1:(iSubbandNumber-1) * Mt+Mt,iUser);
        end
    end
elseif strcmp(fadingType,'selective')
    for iUser = 1:K
        for iSubbandNumber = 1:subbandNumber
            for iMt = 1:Mt
            subbandChannelGain((iSubbandNumber-1) * Mt + iMt,iUser) = sum(tapGain(:,iMt,iUser) .* exp(-1i * 2 * pi * subbandFrequency(iSubbandNumber) * tapDelay));
            end
            subbandChannelGainFre(:,iSubbandNumber,iUser) = subbandChannelGain((iSubbandNumber-1) * Mt + 1:(iSubbandNumber-1) * Mt+Mt,iUser);
        end
    end
else
    disp('wrong channel type')
end
ChannelInfo.subbandChannelGainFre = subbandChannelGainFre;            
ChannelInfo.subbandChannelGainAmplitudeFre = abs(subbandChannelGainFre);   
ChannelInfo.subbandChannelGainPhaseFre = angle(subbandChannelGainFre);            
ChannelInfo.subbandChannelGain = subbandChannelGain;
ChannelInfo.subbandChannelGainAmplitude = abs(subbandChannelGain);   
ChannelInfo.subbandChannelGainPhase = angle(subbandChannelGain);