clear all; close all; clc
load handel.mat

T =size(y,1)/Fs;

T_frame = 0.5;
N_frame = round(T_frame*Fs);
T_delay = 0.5e-3;
N_delay = round(T_delay*Fs);

y_1 = y(1:N_frame);
y_2 = y((1+N_delay):(N_frame+N_delay));
t = linspace(0,2,size(y_1,1));
figure(1)
plot(t,[y_1,y_2])

[R,lags] = xcorr(y_1,y_2);
figure(2)
stem(lags/Fs,R)
