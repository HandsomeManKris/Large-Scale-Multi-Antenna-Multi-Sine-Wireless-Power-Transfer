%%plot SISO channel response
initialization;
[ChannelInfo] = modelE(ChannelInfo,TransceiverInfo );
ResponseInfo  = ChannelInfo;
subbandFrequency = -1.25 * bandwidth +  centerFrequency:2e4: 1.25 * bandwidth +  centerFrequency;
basebandFrequency = subbandFrequency - centerFrequency;
partIndex = find(basebandFrequency <= 0.1 * bandwidth & basebandFrequency >= -0.1 * bandwidth);
subbandNumber = length(subbandFrequency);
ResponseInfo.subbandNumber = subbandNumber;
ResponseInfo.subbandFrequency = subbandFrequency;
[ResponseInfo] = get_channel_response(TransceiverInfo, ResponseInfo);
figure(1) 
plot(basebandFrequency / 1e6, ResponseInfo.subbandChannelGainAmplitude); hold on
plot(basebandFrequency(partIndex) / 1e6, ResponseInfo.subbandChannelGainAmplitude(partIndex));
grid on; grid minor;
legend('10 MHz', '1 MHz','northeast');
xlabel('Frequency [MHz]');
ylabel('Frequency response');
xlim([-5, 5]);