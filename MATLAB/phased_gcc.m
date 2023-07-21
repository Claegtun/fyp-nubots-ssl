clear all

N = 2;
rxULA = phased.ULA('Element',...
    phased.OmnidirectionalMicrophoneElement,...
    'ElementSpacing',0.5,...
    'NumElements',N);

rxpos1 = [0;0;0];
rxvel1 = [0;0;0];
rxax1 = azelaxes(90,0);

srcpos = [30;100;0];
srcvel = [0;0;0];
srcax = azelaxes(-90,0);
srcULA = phased.OmnidirectionalMicrophoneElement;

fc = 300e3;             % 300 kHz
c = 1500;               % 1500 m/s
dmax = 150.0;             % 150 m
pri = (2*dmax)/c;
prf = 1/pri;
bw = 100.0e3;           % 100 kHz
fs = 2*bw;
waveform = phased.LinearFMWaveform('SampleRate',...
    fs,'SweepBandwidth',bw,...
    'PRF',prf,'PulseWidth',pri/10);

signal = waveform();

nfft = 128;

radiator = phased.WidebandRadiator('Sensor',srcULA,...
    'PropagationSpeed',c,'SampleRate',fs,...
    'CarrierFrequency',fc,'NumSubbands',nfft);
collector1 = phased.WidebandCollector('Sensor',rxULA,...
    'PropagationSpeed',c,'SampleRate',fs,...
    'CarrierFrequency',fc,'NumSubbands',nfft);

channel1 = phased.WidebandFreeSpace('PropagationSpeed',c,...
    'SampleRate',fs,'OperatingFrequency',fc,'NumSubbands',nfft);

[~,ang1t] = rangeangle(rxpos1,srcpos,srcax);

sigt = radiator(signal,ang1t);

sigp1 = channel1(sigt(:,1),srcpos,rxpos1,srcvel,rxvel1);

[~,ang1r] = rangeangle(srcpos,rxpos1,rxax1);

sigr1 = collector1(sigp1,ang1r);

doa1 = phased.GCCEstimator('SensorArray',rxULA,'SampleRate',fs,...
    'PropagationSpeed',c);

angest1 = doa1(sigr1)


