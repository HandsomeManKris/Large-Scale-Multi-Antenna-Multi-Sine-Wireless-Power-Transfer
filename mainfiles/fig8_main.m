tic
initialization;
realization =10;
VoutuSER = zeros(K,realization);
VoutuSER1 = zeros(K,realization);
Vout = zeros(realization,1);
Vout1 = zeros(realization,1);
for iRealization = 1:realization
    TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP;
    [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
    [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
    InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo); 
    [sprecoder,VoutuSER(:,iRealization),Vout(iRealization)] = get_Vout_MAX_MIN_RR(ChannelInfo,TransceiverInfo,InitialM1);
    for iTcase = 1:ChannelInfoGroup.nTcase
        TransceiverInfo.T = ChannelInfoGroup.Tcase(iTcase);
        [sprecoder1(:,iTcase),VoutuSER1(:,iTcase,iRealization),Vout1(iTcase,iRealization)] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
    end
end
VoutAve = mean(Vout)/K * 1000;
VoutTotal(1) = VoutAve ;
for iTcase = 1:ChannelInfoGroup.nTcase
    VoutAve1(iTcase) = mean(Vout1(iTcase,:))/K * 1000;
    VoutTotal(iTcase+1) = VoutAve1(iTcase);
end




bar((1:1:ChannelInfoGroup.nTcase+1),real(VoutTotal));
xlabel('K=2')
ylabel('Average v_{out}[mv]')

ylim([real(min(VoutTotal)-0.005 )  real(max(VoutTotal)) ]);
set(gca,'xticklabel',{'RR','T=500','T=50','T=5','T=1'});
toc
    
