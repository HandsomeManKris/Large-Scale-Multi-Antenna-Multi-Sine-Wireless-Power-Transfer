%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the initialization file for this paper
fadingType = 'selective';                   %The fading type, 'flat','selective'
bandwidth = 1e7;                       %Bandwidth of basedband
subbandNumber = 16;                     %The number of subband:1,2,4,8,16
frequencyGap = bandwidth/subbandNumber;%\delta_f = B/N
centerFrequency = 2.4e9;              %carrier frequency
% -(B-B/N)*0.5+fc:B/N:(B-B/N)*0.5+fc freuqncy of each carrier
subbandFrequency = -0.5*(bandwidth-bandwidth/subbandNumber) +  centerFrequency:bandwidth/subbandNumber:0.5 * (bandwidth-bandwidth/subbandNumber) +  centerFrequency;
ChannelInfo.fadingType = fadingType;
ChannelInfo.bandwidth = bandwidth ; 
ChannelInfo.subbandNumber = subbandNumber ;
ChannelInfo.frequencyGap = frequencyGap ; 
ChannelInfo.centerFrequency = centerFrequency ; 
ChannelInfo.subbandFrequency = subbandFrequency ;
%%

%TransceiverInfo
Mt = 50;               %Number of transmit antennas
Mr = 1;               %Number of receive antennas
K = 2; 
Resistance = 50;
ni = 1;
%Is = 5e-6;
Vt = 25.86e-3;
b2 = Resistance/(2 * ni * Vt);
b4 =  Resistance^2/(24 * ni^3 * Vt^3);
%k2 = b2 *Is /(ni * Vt)/ Resistance;
%k4 = b4 *Is /(ni * Vt)/Resistance^2;
%Number of users
%k2 = 0.0034;          %diode parameter 
%k4 = 0.3829;          %diode parameter
%MtPowerdBm = -20 ;    %average transmite pwoer
%MrPowerdBm = -20 ;    %average receive pwoer
%MtPower = 10 .^ ((MtPowerdBm - 30) / 10);
%MrPower = 10 .^ ((MrPowerdBm - 30) / 10);
EIRPdBm = 36;
EIRP = 10 .^ ((EIRPdBm - 30) / 10);
MrPowerdBmEIRP = EIRPdBm -60.046 ;
MrPowerdBmRegion = EIRPdBm - 66.07;
MrPowerEIRP =  10 .^ ((MrPowerdBmEIRP - 30) / 10)/Mt;
MrPowerEIRP1 = 10 .^ ((MrPowerdBmEIRP - 30) / 10);
MrPowerRegion = 10 .^ ((MrPowerdBmRegion - 30) / 10)/Mt;
MrPowerRegion1 = 10 .^ ((MrPowerdBmRegion - 30) / 10);
distance = 10;
MtPower = 0.5 ;
MrPowerdBm = log10(MtPower) * 10 + 30 - distance/10 * 60.046;
MrPower = 10 .^ ((MrPowerdBm - 30) / 10);
T = 50;
%SNRdB = 30;
%SNR = db2pow(SNRdB);
%noisePowerdBm = MrPowerdBm -  SNRdB;
%noisePower = 10 .^ ((noisePowerdBm - 30) / 10);
TransceiverInfo.T = T;
TransceiverInfo.MrPowerRegion = MrPowerRegion;
TransceiverInfo.MrPowerRegion1 = MrPowerRegion1;
TransceiverInfo.MrPowerEIRP1 = MrPowerEIRP1;
Tolerance = 1e-6;
TransceiverInfo.MrPowerEIRP = MrPowerEIRP;
TransceiverInfo.Mt =Mt;
TransceiverInfo.Mr =Mr;
TransceiverInfo.MtPower =MtPower;
TransceiverInfo.MrPower =MrPower;
TransceiverInfo.K = K;
%TransceiverInfo.k2 =k2;
%TransceiverInfo.k4 =k4;
%TransceiverInfo.SNR =SNR;
TransceiverInfo.b2 =b2;
TransceiverInfo.b4 =b4;
TransceiverInfo.Resistance = Resistance;
TransceiverInfo.Tolerance = Tolerance;
%TransceiverInfo.noisePower = noisePower;

%%
%ChannelInfoGroup
distanceGroup = [10,12,14,16,18,20];
MrPowerdBmGroup = log10(MtPower) * 10 + 30 -  60.046 -20*log10(distanceGroup/distance);
MrPowerGroup = 10 .^ ((MrPowerdBmGroup - 30) / 10);
subbandNumberGroup = [1,2,4,8,16,32,64];
nSubbandNumber = length(subbandNumberGroup);
frequencyGapGroup = bandwidth./subbandNumberGroup;
subbandFrequencyGroup = cell(1,nSubbandNumber);
for iGroup =1 : nSubbandNumber
    subbandFrequencyGroup{iGroup} = -0.5*(bandwidth-bandwidth/subbandNumberGroup(iGroup)) +  centerFrequency:bandwidth/subbandNumberGroup(iGroup):0.5 * (bandwidth-bandwidth/subbandNumberGroup(iGroup)) +  centerFrequency;
end
Tcase = [500,50,5,1];
nTcase = length(Tcase);
Kcase = [2,3,4,5,6,7,8,9];
nKcase = length(Kcase);
ChannelInfoGroup.Kcase = Kcase;
ChannelInfoGroup.nKcase = nKcase;
ChannelInfoGroup.Tcase = Tcase;
ChannelInfoGroup.nTcase = nTcase;
ChannelInfoGroup.nSubbandNumber = nSubbandNumber;
ChannelInfoGroup.frequencyGapGroup = frequencyGapGroup;
ChannelInfoGroup.subbandFrequencyGroup = subbandFrequencyGroup;
ChannelInfoGroup.subbandNumberGroup = subbandNumberGroup;
 
ChannelInfoGroup.distanceGroup = distanceGroup;
 ChannelInfoGroup.MrPowerGroup = MrPowerGroup;