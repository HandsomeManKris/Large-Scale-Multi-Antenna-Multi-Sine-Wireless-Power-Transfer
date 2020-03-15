tic
initialization;
realization =450;
Vout1 = zeros(realization,ChannelInfoGroup.nKcase);
Vout2 = zeros(realization,ChannelInfoGroup.nKcase);
Vout3 = zeros(realization,ChannelInfoGroup.nKcase);
Vout4 = zeros(realization,ChannelInfoGroup.nKcase);
VoutUser1Min1 = zeros(realization * 2,1);
VoutUser1Min2 = zeros(realization * 5,1);
VoutUser2Min1 = zeros(realization * 2,1);
VoutUser2Min2 = zeros(realization * 5,1);
VoutUser3Min1 = zeros(realization * 2,1);
VoutUser3Min2 = zeros(realization * 5,1);
VoutUser4Min1 = zeros(realization * 2,1);
VoutUser4Min2 = zeros(realization * 5,1);
for iRealization = 1:realization    
    for iKcase = 1:ChannelInfoGroup.nKcase
        VoutUser1 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser2 = zeros(ChannelInfoGroup.Kcase(iKcase),1); 
        VoutUser3 = zeros(ChannelInfoGroup.Kcase(iKcase),1);
        VoutUser4 = zeros(ChannelInfoGroup.Kcase(iKcase),1); 
        TransceiverInfo.K = ChannelInfoGroup.Kcase(iKcase);       
        TransceiverInfo.MrPower = TransceiverInfo.MrPowerEIRP;
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
        InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);
        [sprecoder1,VoutUser1,Vout1] = get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM1);
        %VoutUser1Min(iRealization,iKcase) = min(VoutUser1);
        [sprecoder2,VoutUser2,Vout2] =  get_Vout_WSUM(ChannelInfo,TransceiverInfo,InitialM1);
       
        [sprecoder3,VoutUser3,Vout3] =  get_Vout_WSUM_FA(ChannelInfo,TransceiverInfo,InitialM1);
         
        TransceiverInfo.MrPower =TransceiverInfo.MrPowerEIRP1;
        InitialM = get_compact_CHE_WSUM(TransceiverInfo, ChannelInfo);
        [sprecoder4,VoutUser4,Vout4(iRealization,iKcase)] = get_Vout_CHE_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM);
        %VoutUser4Min(iRealization,iKcase) = min(VoutUser4);
        if iKcase ==1
            VoutUser1Min1((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser1);
            VoutUser2Min1((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser2);
            VoutUser3Min1((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser3);
            VoutUser4Min1((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser4);
        else
            VoutUser1Min2((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser1);
            VoutUser2Min2((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser2); 
            VoutUser3Min2((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser3); 
            VoutUser4Min2((iRealization-1) *ChannelInfoGroup.Kcase(iKcase) +1 : (iRealization-1) *ChannelInfoGroup.Kcase(iKcase) + ChannelInfoGroup.Kcase(iKcase),1) = (VoutUser4);
        end
    end
end


h1 = cdfplot(log10(real(VoutUser1Min1(:,1)) * 1000));
set(h1,'color','b','LineStyle','-');  hold on;
h2 = cdfplot(log10(real(VoutUser4Min1(:,1)) * 1000));
set(h2,'color','b','LineStyle','--'); hold on;
h3 = cdfplot(log10(real(VoutUser3Min1(:,1)) * 1000));
set(h3,'color','b','LineStyle','-.'); hold on;
h4 = cdfplot(log10(real(VoutUser2Min1(:,1)) * 1000));
set(h4,'color','b','LineStyle',':'); hold on;
h5 = cdfplot(log10(real(VoutUser1Min2(:,1)) * 1000));
set(h5,'color','r','LineStyle','-'); hold on;
h6 = cdfplot(log10(real(VoutUser4Min2(:,1)) * 1000));
set(h6,'color','r','LineStyle','--'); hold on;
h7 = cdfplot(log10(real(VoutUser3Min2(:,1)) * 1000));
set(h7,'color','r','LineStyle','-.'); hold on;
h8 = cdfplot(log10(real(VoutUser2Min2(:,1)) * 1000));
set(h8,'color','r','LineStyle',':');  hold on;
legend('MAX-MIN-Rand','Randomized CHE Max-Min','FA WSum','WSUM w_{q} =1','location','best')
xlabel(' v_{out}[mV]');
ylabel('CDF');
grid on;
set(gca,'xtick',[-3 -2 -1 0 1 2]);
set(gca,'xticklabel',{'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'});
xlim([-3 2])






toc