initialization;
realization = 100;
Vout = zeros(realization,ChannelInfoGroup.nSubbandNumber );
Vout1 = zeros(realization,ChannelInfoGroup.nSubbandNumber );
Vout2 = zeros(realization,ChannelInfoGroup.nSubbandNumber );
Vout3 =  zeros(realization,ChannelInfoGroup.nSubbandNumber );
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
        
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP;
         InitialM1 = get_compact_singleuser(TransceiverInfo, ChannelInfo);    
        [sprecoder1,pstar1,Vout1(iRealization,iSubbandCase)] = get_Vout_singleuser(ChannelInfo,TransceiverInfo,InitialM1);
        Vout3(iRealization,iSubbandCase) = get_Vout_ASS(ChannelInfo,TransceiverInfo);
        Vout2(iRealization,iSubbandCase) = get_Vout_UP(ChannelInfo,TransceiverInfo);
        
    end
end
VoutAve = mean(Vout,1) *1000;
VoutAve1 = mean(Vout1,1) *1000;
VoutAve2 = mean(Vout2,1) *1000;
VoutAve3 = mean(Vout3,1) *1000;


for iSubbandCase = 1:ChannelInfoGroup.nSubbandNumber
    Vouttotal(:,1) = VoutAve1;
    Vouttotal(:,2) = VoutAve;
    Vouttotal(:,3) = VoutAve2;
    Vouttotal(:,4) = VoutAve3;
end
        
        
bar((1:1:ChannelInfoGroup.nSubbandNumber),abs(Vouttotal))     ; 

legend('SU WPT','CHE WSUM','UP','ASS','location','northwest')
xlabel('Number of tones');
ylabel('Average v_{out}[mV]');
grid on;
%set(gca,'xtick',[1 2 4 8 16 32 64]);
%set(gca,'xticklabel',{'(1,1)','(1,2)','(1,4)','(1,8)','(1,16)','(1,32)','(1,64)'});
%set(gca,'xticklabel',{'(4,1)','(4,2)','(4,4)','(4,8)','(4,16)','(4,32)','(4,64)'});
set(gca,'xticklabel',{'(20,1)','(20,2)','(20,4)','(20,8)','(20,16)','(20,32)','(20,64)'});
%xlim([1 64])