initialization;
 Vout = zeros(100,length(distanceGroup));
for iRealization = 1:100
    for iDistance =1: length(distanceGroup)
        TransceiverInfo.MrPower =ChannelInfoGroup.MrPowerGroup(iDistance);
        [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
        [ChannelInfo] = get_channel_response(TransceiverInfo, ChannelInfo);            
        Vout(iRealization,iDistance) =  get_Vout_ASS(ChannelInfo,TransceiverInfo);
    end
end
VoutAve = mean(Vout,1);
 xlabel('Distance(m)');
 ylabel('Average v_{out}[V]');
 legend('ASS');
semilogy(ChannelInfoGroup.distanceGroup, VoutAve);hold on;
grid on;