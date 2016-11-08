clear all; close all; clc;

N = 2;                  % Number of transmit antennas
M = 2;                  % Number of receive antennas
K = 64;                 % Number of OFDM channels
EbNoVec = 1:2:20;        % Eb/No in dB
modOrds = 2:1:4;             % constellation size = 2^modOrd
msgLen = 5000*K;

rng(12);    % seed 12 for really bad channel
% 457 is decent
chanlen = 3;
numchans = 3

u = 10;
for k = 1:numchans
    chan = reshape(chans(:,:,k,:), [2 2 chanlen]);
    changains = fft(chan, 64, 3);
    for modord = modOrds
        txmsg = randi([0 2^modord-1], [1 N*msgLen]);
        for i = 1:length(EbNoVec)
            modmsg = qammod(txmsg,2^modord,0,'gray');
            modmsg = reshape(modmsg, K, []);
            txsig = ifft(modmsg, K, 1);
            
            txsig = reshape(txsig, [K msgLen/K N]);
            
            sigPower = sum(abs(txsig(:)).^2)/length(txsig(:));
            reqSigPow =  10*log10(sigPower);
            reqSNR = EbNoVec(i)+10*log10(modord);
            
            txsigCyclic = [txsig(K-u+1:K,:,:); txsig(:,:,:)];
            txsigCyclic = reshape(txsigCyclic, [(K+u)*msgLen/K N]);
            rxsig = zeros((K+u)*msgLen/K+chanlen-1, 2);
            rxsig(:,1) = conv(txsigCyclic(:,1), reshape(chan(1,1,:), chanlen, 1)) + conv(txsigCyclic(:,2), reshape(chan(1,2,:), chanlen, 1));
            rxsig(:,2) = conv(txsigCyclic(:,1), reshape(chan(2,1,:), chanlen, 1)) + conv(txsigCyclic(:,2), reshape(chan(2,2,:), chanlen, 1));
            
            rxsig = rxsig(1:end-chanlen+1,:);
            rxsig = awgn(rxsig, reqSNR, reqSigPow);
            rxsig = reshape(rxsig, [K+u msgLen/K N]);
            rxsig = rxsig(u+1:end,:,:);
            rxmod = fft(rxsig, K, 1);
            rxmsg = zeros(K, msgLen/K, N);
            for k = 1:K
                rx = reshape(rxmod(k,:,:), msgLen/K, N).';
                rxeq = (pinv(changains(:,:,k))*rx).';
                rxmsg(k,:,:) = rxeq;
            end
            
            rxmsg = reshape(rxmsg, 1, []);
            rxmsg = qamdemod(rxmsg, 2^modord, 0, 'gray');
            
            [numerr, ber] = biterr(txmsg, rxmsg);
            BER(i) = ber;
        end
        figure
        semilogy(EbNoVec,BER,'-r');
        title(['64 channel OFDM with 2X2 MIMO ', num2str(2^modord), '-QAM ']);
        legend('Zero-Forcing', 'location', 'best');
        xlabel('SNR (dB)');
        ylabel('BER');
    end
end
