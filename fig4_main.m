initialization;
realization = 1;
Vout = zeros(realization,length(subbandNumberGroup ));
Vout1 = zeros(realization,length(subbandNumberGroup ));
for iRealization = 1:realization
    for iSubbandCase = 1:ChannelInfoGroup.nSubbandNumber 
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP;
        ChannelInfo.subbandNumber = ChannelInfoGroup.subbandNumberGroup(iSubbandCase);
        
        ChannelInfo.subbandFrequency = ChannelInfoGroup.subbandFrequencyGroup{iSubbandCase};
        ChannelInfo.frequencyGap = ChannelInfoGroup.frequencyGapGroup(iSubbandCase );
        ChannelInfo.subbandNumber
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_singleuser(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_WSUM(TransceiverInfo, ChannelInfo);
        [sprecoder1,pstar1,Vout1(iRealization,iSubbandCase)] = get_Vout_singleuser(ChannelInfo,TransceiverInfo,InitialM1);
        [sprecoder,pstar,Vout(iRealization,iSubbandCase)] = get_Vout_WSUM(ChannelInfo,TransceiverInfo,InitialM);
    end
end
VoutAve = mean(Vout,1);
VoutAve1 = mean(Vout1,1);
%figure('name','figure4')


plot(ChannelInfoGroup.subbandNumberGroup, VoutAve * 1000,'ro-');hold on;
plot(ChannelInfoGroup.subbandNumberGroup, VoutAve1 * 1000,'r*-');hold on;
xlabel('Number of tones');
ylabel('Average v_{out}[mV]');
grid on;
set(gca,'xtick',[1 2 4 8 16 32 64]);
%set(gca,'xticklabel',{'1','2','','4','','','','8','','','','','','','','16','','','','','','','','','','','','','','','','32','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','64'}); 
xlim([1 64])