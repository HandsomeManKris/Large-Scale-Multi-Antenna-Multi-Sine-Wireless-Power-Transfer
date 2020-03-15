tic
initialization;
realization =1;

Vout1 = zeros(realization,ChannelInfoGroup.nKcase);
Vout2 = zeros(realization,ChannelInfoGroup.nKcase);
Vout3 = zeros(realization,ChannelInfoGroup.nKcase);
Vout4 = zeros(realization,ChannelInfoGroup.nKcase);
Vout5 = zeros(realization,ChannelInfoGroup.nKcase);
Vout6 = zeros(realization,ChannelInfoGroup.nKcase);
Vout7 = zeros(realization,ChannelInfoGroup.nKcase);
Vout8 = zeros(realization,ChannelInfoGroup.nKcase);
Vout9 = zeros(realization,ChannelInfoGroup.nKcase);
VoutUser1Min = zeros(realization,nKcase);
VoutUser2Min = zeros(realization,nKcase);
VoutUser3Min = zeros(realization,nKcase);
VoutUser4Min = zeros(realization,nKcase);
VoutUser5Min = zeros(realization,nKcase);
VoutUser6Min = zeros(realization,nKcase);
VoutUser7Min = zeros(realization,nKcase);
VoutUser8Min = zeros(realization,nKcase);
VoutUser9Min = zeros(realization,nKcase);
for iRealization = 1:realization
    
    
    for iKcase = 1:ChannelInfoGroup.nKcase
        VoutUser1 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser2 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser3 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);
        
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);    
        [sprecoder1,VoutUser1,Vout1(iRealization,iKcase)] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
        
        VoutUser1Min(iRealization,iKcase) = min(VoutUser1);
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
       % [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
       % [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
        [sprecoder2,VoutUser2,Vout2(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_RR1(ChannelInfo,TransceiverInfo,InitialM);
        [sprecoder3,VoutUser3,Vout3(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
        
        VoutUser2Min(iRealization,iKcase) = min(VoutUser2);
          
        VoutUser3Min(iRealization,iKcase) = min(VoutUser3);
    end
end


initialization1;



for iRealization = 1:realization
    
    for iKcase = 1:ChannelInfoGroup.nKcase
        VoutUser4 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser5 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser6 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);
        
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);    
        [sprecoder4,VoutUser4,Vout4(iRealization,iKcase)] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
        
        VoutUser4Min(iRealization,iKcase) = min(VoutUser4);
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
       % [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
       % [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
        [sprecoder5,VoutUser5,Vout5(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_RR1(ChannelInfo,TransceiverInfo,InitialM);
        [sprecoder6,VoutUser6,Vout6(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
        
        VoutUser5Min(iRealization,iKcase) = min(VoutUser5);
          
        VoutUser6Min(iRealization,iKcase) = min(VoutUser6);
    end
end

initialization2;



for iRealization = 1:realization
    
    for iKcase = 1:ChannelInfoGroup.nKcase
        VoutUser7 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser8 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser9 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);
        
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);    
        [sprecoder7,VoutUser7,Vout7(iRealization,iKcase)] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
        
        VoutUser7Min(iRealization,iKcase) = min(VoutUser7);
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP1;
       % [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
       % [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
        [sprecoder8,VoutUser8,Vout8(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_RR1(ChannelInfo,TransceiverInfo,InitialM);
        [sprecoder9,VoutUser9,Vout9(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);       
        VoutUser8Min(iRealization,iKcase) = min(VoutUser8);          
        VoutUser9Min(iRealization,iKcase) = min(VoutUser9);
    end
end





VoutAve1 = mean(VoutUser1Min,1) * 1000;
VoutAve2 = mean(VoutUser2Min,1) * 1000;
VoutAve3 = mean(VoutUser3Min,1) * 1000;
VoutAve4 = mean(VoutUser4Min,1) * 1000;
VoutAve5 = mean(VoutUser5Min,1) * 1000;
VoutAve6 = mean(VoutUser6Min,1) * 1000;
VoutAve7 = mean(VoutUser7Min,1) * 1000;
VoutAve8 = mean(VoutUser8Min,1) * 1000;
VoutAve9 = mean(VoutUser9Min,1) * 1000;

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