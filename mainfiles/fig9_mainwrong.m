tic
initialization;
realization =100;
Vout1 = zeros(realization,ChannelInfoGroup.nKcase);
Vout2 = zeros(realization,ChannelInfoGroup.nKcase);
Vout3 = zeros(realization,ChannelInfoGroup.nKcase);
for iRealization = 1:realization
    
    for iKcase = 1:ChannelInfoGroup.nKcase
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);
        
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);    
        [sprecoder1,VoutUser1,Vout1(iRealization,iKcase)] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
        
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
       % [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
       % [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
        [sprecoder2,VoutUser2,Vout2(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_RR1(ChannelInfo,TransceiverInfo,InitialM);
        [sprecoder3,VoutUser3,Vout3(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
    end
end
VoutAve1 = mean(Vout1,1)/K * 1000;
VoutAve2 = mean(Vout2,1)/K * 1000;
VoutAve3 = mean(Vout3,1)/K * 1000;


initialization1;
realization =1;
Vout4 = zeros(realization,ChannelInfoGroup.nKcase);
Vout5 = zeros(realization,ChannelInfoGroup.nKcase);
Vout6 = zeros(realization,ChannelInfoGroup.nKcase);
for iRealization = 1:realization
    
    for iKcase = 1:ChannelInfoGroup.nKcase
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);
        
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);    
        [sprecoder1,VoutUser1,Vout4(iRealization,iKcase)] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
        
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
       % [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
       % [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
        [sprecoder2,VoutUser2,Vout5(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_RR1(ChannelInfo,TransceiverInfo,InitialM);
        [sprecoder3,VoutUser3,Vout6(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
    end
end
VoutAve4 = mean(Vout4,1)/K * 1000;
VoutAve5 = mean(Vout5,1)/K * 1000;
VoutAve6 = mean(Vout6,1)/K * 1000;


initialization2;
realization =1;
Vout7 = zeros(realization,ChannelInfoGroup.nKcase);
Vout8 = zeros(realization,ChannelInfoGroup.nKcase);
Vout9 = zeros(realization,ChannelInfoGroup.nKcase);
for iRealization = 1:realization
    
    for iKcase = 1:ChannelInfoGroup.nKcase
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);
        
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);    
        [sprecoder1,VoutUser1,Vout7(iRealization,iKcase)] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
        
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
       % [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
       % [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
        [sprecoder2,VoutUser2,Vout8(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_RR1(ChannelInfo,TransceiverInfo,InitialM);
        [sprecoder3,VoutUser3,Vout9(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
    end
end
VoutAve7 = mean(Vout7,1)/K * 1000;
VoutAve8 = mean(Vout8,1)/K * 1000;
VoutAve9 = mean(Vout9,1)/K * 1000;






Vouttotal2(1:4,1) = VoutAve1;
Vouttotal2(1:4,2) = VoutAve2;
Vouttotal2(1:4,3) = VoutAve3;
Vouttotal2(5:8,1) = VoutAve4;
Vouttotal2(5:8,2) = VoutAve5;
Vouttotal2(5:8,3) = VoutAve6;
Vouttotal2(9:12,1) = VoutAve7;
Vouttotal2(9:12,2) = VoutAve8;
Vouttotal2(9:12,3) = VoutAve9;


bar((1:1:12),Vouttotal2)
xlabel('Average v_{out}[mV] as a function of (M,N,K)');
ylabel('Average v_{out}[mV]');
legend('Max-Min-Rand','CHE Max-Min-RR','Randmized CHE Max-Min')
set(gca,'xticklabel',{'(4,4,2)','(4,4,3)','(4,4,4)','(4,4,5)','(4,8,2)','(4,8,3)','(4,8,4)','(4,8,5)','(20,8,2)','(20,8,3)','(20,8,4)','(20,8,5)'});
grid on;
toc