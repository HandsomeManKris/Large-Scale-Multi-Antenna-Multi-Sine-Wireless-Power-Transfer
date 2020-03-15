function [ChannelInfo] = modelE(ChannelInfo,TransceiverInfo )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The explanation of the function
%Input:TranceiverInfo:Mr,Mt,ChannelInfo
%
%Output: The ChannelInfo:tapGain,tapDelay,tapAverageRelativePower,tap
%
%Explanation: corresponding to typical large open space and office environments for
%NLOS conditions and 100ns average rms delay spread..
%
%Source:Channel Models E for typical large open space (indoor and outdoor), NLOS conditions, and
%250 ns rms delay spread.,
%
%Written by Kris
K = TransceiverInfo.K;
%Mr = TransceiverInfo.Mr;
Mt = TransceiverInfo.Mt;

tapNumber = 18;
%tap = [1:1:tapNumber];
%ns = 10^(-9)s
tapDelay = [0,10,20,30,50,80,110,140,180,230,280,330,380,430,490,560,640,730]' *  1e-9;
%Power from the model
tapAverageRelativePowerdB1 = [-2.6,-3.0,-3.5,-3.9,-4.5,-5.6,-6.9,-8.2,-9.8,-11.7,-13.9,-16.1,-18.3,-20.5,-22.9,-inf,-inf,-inf]';
tapAverageRelativePowerdB2 = [-inf,-inf,-inf,-inf,-1.8,-3.2,-4.5,-5.8,-7.1,-9.9,-10.3,-14.3,-14.7,-18.7,-19.9,-22.4,-inf,-inf]';
tapAverageRelativePowerdB3 = [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-7.9,-9.6,-14.2,-13.8,-18.6,-18.1,-22.8,-inf,-inf,-inf]';
tapAverageRelativePowerdB4 = [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-20.6,-20.5,-20.7,-24.6,]';
%Transform db to power
tapAverageRelativePower1 = db2pow(tapAverageRelativePowerdB1);
tapAverageRelativePower2 = db2pow(tapAverageRelativePowerdB2);
tapAverageRelativePower3 = db2pow(tapAverageRelativePowerdB3);
tapAverageRelativePower4 = db2pow(tapAverageRelativePowerdB4);
tapAverageRelativePower = tapAverageRelativePower1+tapAverageRelativePower2+tapAverageRelativePower3+tapAverageRelativePower4;
%Taps are modeled as CSCG random variables.each with an average power beta
%sum beta =1
%h_{n,m} = \sum\limits{l=0}^{L-1}a_le^{-1i * 2\pi * f_nt},multipath channel
%gain
%tapGain = (randn(tapNumber,Mr,Mt)+1i * randn(tapNumber,Mr,Mt )).* repmat(sqrt(tapAverageRelativePower/2),[1,Mt,Mr]) ; 


tapGain = (randn(tapNumber,Mt,K)+1i * randn(tapNumber,Mt ,K)).* repmat(sqrt(tapAverageRelativePower/2),[1,Mt,K]) ; 
    

ChannelInfo.tapGain = tapGain;
ChannelInfo.tapDelay = tapDelay;
ChannelInfo.tapAverageRelativePower = tapAverageRelativePower;