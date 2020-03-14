tic
initialization;
%VoutAve = zeros(1,length(distanceGroup));
%Vout = zeros(1,500);

for iRealization = 1:100
        for iDistance =1: length(distanceGroup)
            TransceiverInfo.MrPower =ChannelInfoGroup.MrPowerGroup(iDistance);
            [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
            [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);
            InitialM = get_compact_singleuser(TransceiverInfo, ChannelInfo);
            [sprecoder,pstar,Vout(iRealization,iDistance)] = get_Vout_singleuser(ChannelInfo,TransceiverInfo,InitialM);
        end
end
figure('name','fig3a')
 VoutAve = mean(Vout,1);
 xlabel('Distance(m)');
 ylabel('Average v_{out}[V]');
 legend('SU WPT');
semilogy(ChannelInfoGroup.distanceGroup, VoutAve);hold on;
grid on;
%ylim([1e-2 1e0])
    toc
