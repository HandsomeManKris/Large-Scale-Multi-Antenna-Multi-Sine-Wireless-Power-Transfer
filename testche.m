initialization;
realization = 100;
for iRealization = 1:realization
    for iSubbandCase = 1:ChannelInfoGroup.nSubbandNumber 
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
        ChannelInfo.subbandNumber = ChannelInfoGroup.subbandNumberGroup(iSubbandCase);
        ChannelInfo.subbandFrequency = ChannelInfoGroup.subbandFrequencyGroup{iSubbandCase};
        ChannelInfo.frequencyGap = ChannelInfoGroup.frequencyGapGroup(iSubbandCase );
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
        [sprecoder,pstar,Vout(iRealization,iSubbandCase)] = get_Vout_CHE_WSUM1(ChannelInfo,TransceiverInfo,InitialM);
    end
end
VoutAve = mean(Vout,1) *1000;