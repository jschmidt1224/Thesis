clear all; close all; clc;

Subcarriers = 8;
M = 2;

msgLen = 1000*Subcarriers;

txmsg = randi([0 M-1], [Subcarriers msgLen/Subcarriers]);
txmod = qammod(txmsg, M);
txmod(2,:) = 2*txmod(2,:);
txsig = ifft(txmod, Subcarriers, 1);
txsig = txsig(:);

txspec = fft(txsig);

plot(abs(txspec))