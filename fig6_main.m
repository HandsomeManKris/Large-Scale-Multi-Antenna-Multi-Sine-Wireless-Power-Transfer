initialization;
realization = 1;
TransceiverInfo.MrPower =TransceiverInfo.MrPowerRegion1;
[ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
[ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
InitialM = get_compact_CHE_WSUM1(TransceiverInfo, ChannelInfo);
[sprecoder,VoutUser,Vout] = get_Vout_CHE_WSUM_RateRegion(ChannelInfo,TransceiverInfo,InitialM);
    
TransceiverInfo.MrPower =TransceiverInfo.MrPowerRegion;
InitialM1 = get_compact_WSUM(TransceiverInfo, ChannelInfo);    
[sprecoder1,VoutUser1,Vout1] = get_Vout_WSUM_RateRegion(ChannelInfo,TransceiverInfo,InitialM1);







VoutAveCHE1 = real(mean(VoutUser(find(VoutUser(:,1)>0),1))) * 1000;
VoutAveCHE2 = real(mean(VoutUser(find(VoutUser(:,2)>0),2))) * 1000;
VoutAveWSUM1 = max(VoutUser1(:,1)) * 1000;
VoutAveWSUM2 = max(VoutUser1(:,2)) * 1000;
%VoutUser2(:,1) = [  VoutUser1(:,1).' VoutAveWSUM1/1000 ].';
%VoutUser2(:,2) = [  VoutUser1(:,1).' 0 ].';
%[VoutUser1(:,1),index] = sort(VoutUser1(:,1),'ascend');
%VoutUser1(:,2) = VoutUser1(index,2);
plot(real(VoutUser1(:,2) * 1000),real(VoutUser1(:,1) * 1000),'b-');hold on;
%plot([0 real(VoutUser1(:,2).' * 1000), VoutAveWSUM2],[VoutAveWSUM1 flipud(real(VoutUser1(:,1).' * 1000)) 0],'b-') ;hold on;
%x =  real(VoutUser1(:,1)) * 1000;
%y =  real(VoutUser1(:,2)) * 1000;
%k=convhull(x,y);
%x1=x(k);
%y1=y(k);
%xx=floor(x1.*10^(5))./(10^(5));
%ind=find(xx==0);
%[~,ind_ini]=max(x1);
%plot(x1(ind_ini(1):end),y1(ind_ini(1):end),'+','LineWidth',2,'Color', [0.4660    0.6740    0.1880]) ;hold on;
%plot(x1(ind_ini(1):end),y1(ind_ini(1):end),'+','LineWidth',2,'Color', [0.4660    0.6740    0.1880]) ;hold on;
%plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'+','LineWidth',2,'Color', [0.4660    0.6740    0.1880])
%plot(x1,y1,'+','LineWidth',2,'Color', [0.4660    0.6740    0.1880]);hold on;
plot([VoutAveCHE1,0],[0,VoutAveCHE2],'r--s');hold on;
plot([VoutAveWSUM2,0],[0,VoutAveWSUM1],'b--*');hold on;
legend('WSUM',' CHE TDMA','WSUM TDMA')
xlabel('v_{out} of user 1[mv]');
ylabel('v_{out} of user 2[mv]');
grid on;
xlim([real(min(VoutUser1(:,2)) * 1000)  real(max(VoutUser1(:,2))) * 1000]);
ylim([real(min(VoutUser1(:,1)) * 1000)  real(max(VoutUser1(:,1))) * 1000]);