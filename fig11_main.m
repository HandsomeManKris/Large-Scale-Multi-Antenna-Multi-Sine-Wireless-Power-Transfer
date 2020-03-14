tic
initialization;
realization =10;
Vout1 = zeros(realization,ChannelInfoGroup.nKcase);
Vout2 = zeros(realization,ChannelInfoGroup.nKcase);
Vout3 = zeros(realization,ChannelInfoGroup.nKcase);
Vout4 = zeros(realization,ChannelInfoGroup.nKcase);
VoutUser1Min = zeros(realization,nKcase);
VoutUser2Min = zeros(realization,nKcase);
VoutUser3Min = zeros(realization,nKcase);
VoutUser4Min = zeros(realization,nKcase);
for iRealization = 1:realization    
    for iKcase = 1:ChannelInfoGroup.nKcase
        VoutUser1 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser2 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
       
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);       
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        [sprecoder1,VoutUser1,Vout1(iRealization,iKcase)] = get_Vout_UP_MU(ChannelInfo,TransceiverInfo);        
        VoutUser1Min(iRealization,iKcase) = min(VoutUser1);
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
        InitialM = get_compact_CHE_WSUM(TransceiverInfo, ChannelInfo);
        [sprecoder2,VoutUser2,Vout2(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
        VoutUser2Min(iRealization,iKcase) = min(VoutUser2);
    end
end

initialization3;

for iRealization = 1:realization    
    for iKcase = 1:ChannelInfoGroup.nKcase
        VoutUser3 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser4 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
       
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);       
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        [sprecoder3,VoutUser3,Vout3(iRealization,iKcase)] = get_Vout_UP_MU(ChannelInfo,TransceiverInfo);        
        VoutUser3Min(iRealization,iKcase) = min(VoutUser3);
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
        InitialM = get_compact_CHE_WSUM(TransceiverInfo, ChannelInfo);
        [sprecoder4,VoutUser4,Vout4(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
        VoutUser4Min(iRealization,iKcase) = min(VoutUser4);
    end
end



VoutAve1 = mean(VoutUser1Min,1) * 1000;
VoutAve2 = mean(VoutUser2Min,1) * 1000;
VoutAve3 = mean(VoutUser3Min,1) * 1000;
VoutAve4 = mean(VoutUser4Min,1) * 1000;
Vouttotal2(1:4,2) = VoutAve1;
Vouttotal2(1:4,1) = VoutAve2;
Vouttotal2(5:8,2) = VoutAve3;
Vouttotal2(5:8,1) = VoutAve4;
bar((1:1:8),Vouttotal2)
xlabel('Average v_{out}[mV] as a function of (M,N,K)');
ylabel('Average v_{out}[mV]');
set(gca,'xticklabel',{'(8,3)','(8,5)','(8,7)','(8,9)','(16,3)','(16,5)','(16,7)','(16,9)'});
grid on;
legend('Randmized CHE Max-Min','MU UP','location','northwest')
toc