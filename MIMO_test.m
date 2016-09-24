clear all; close all; clc;

N = 2;                  % Number of transmit antennas
M = 2;                  % Number of receive antennas
K = 64;                 % Number of OFDM channels
EbNoVec = 1:2:20;        % Eb/No in dB
modOrds = 2:1:4;             % constellation size = 2^modOrd
msgLen = 5000*K;

rng(12);    % seed 12 for really bad channel
% 457 is decent


%% MIMO
seed = [12 457 11];
msg = {'Bad channel'; 'Medium channel'; 'Good channel'};
for k = 1:length(seed);
    rng(seed(k));
    rayleighChan = (randn(M, N) +  1i*randn(M, N))/sqrt(2);
    for modord = modOrds
        for i = 1:length(EbNoVec)
            txmsg = randi([0 2^modord-1], [N msgLen]);
            txsig = qammod(txmsg,2^modord,0,'gray');
            sigPower = sum(abs(txsig(:)).^2)/length(txsig(:));
            
            %chansig = awgn(txsig, EbNoVec(i)+10*log10(modord), 10*log10(sigPower));
            reqSigPow =  10*log10(sigPower);
            reqSNR = EbNoVec(i)+10*log10(modord);
            chansig = awgn(rayleighChan*txsig,reqSNR,reqSigPow);
            
            %Zero forcing!
            rxsigZF = pinv(rayleighChan)*chansig;
            
            %MMSE!
            No = reqSigPow - reqSNR;
            No = 10^(0.1*No);
            rxsigMMSE = (inv(rayleighChan'*rayleighChan + No*eye(size(rayleighChan,2)))*rayleighChan')*chansig;
            % where No is the noisepower
            
            %Precoding
            [U, S, V] = svd(rayleighChan);
            txsigPre = V*txsig;
            chansigpre = awgn(rayleighChan*txsigPre, reqSNR, reqSigPow);
            rxsigPre = U'*chansigpre;
            rxsigPre = S\rxsigPre;
            
            rxmsgPre = qamdemod(rxsigPre,2^modord,0,'gray');
            
            rxmsgZF = qamdemod(rxsigZF, 2^modord, 0, 'gray');
            rxmsgMMSE = qamdemod(rxsigMMSE,2^modord,0,'gray');
            
            [numerr, ber] = biterr(txmsg, rxmsgZF);
            BER(1,i) = ber;
            [numerr, ber] = biterr(txmsg, rxmsgMMSE);
            BER(2,i) = ber;
            [numerr, ber] = biterr(txmsg, rxmsgPre);
            BER(3,i) = ber;
        end
        figure
        semilogy(EbNoVec,BER(1,:),'r', ...
            EbNoVec,BER(2,:),'g', ...
            EbNoVec,BER(3,:),'k');
        title(['2X2 MIMO ', num2str(2^modord), '-QAM ' msg{k}]);
        legend('Zero-Forcing', 'MMSE', 'Precoding', 'location', 'best');
        xlabel('SNR (dB)');
        ylabel('BER');
    end
end

%% OFDM
chans = [1 .4 .6; 1 .9 .4; .5 .8 .6];
msg = {'Good channel'; 'Medium channel'; 'Bad channel'};
u = 10;
for k = 1:length(chans)
    chan = chans(k,:);
    for modord = modOrds
        for i = 1:length(EbNoVec)
            txmsg = randi([0 2^modord-1], [K msgLen/K]);
            modmsg = qammod(txmsg,2^modord,0,'gray');
            txsig = ifft(modmsg, K, 1);
            
            txsigCyclic = [txsig(K-u+1:K,:); txsig(:,:)];
            txsigCyclic = txsigCyclic(:);
            
            sigPower = sum(abs(txsigCyclic(:)).^2)/length(txsigCyclic(:));
            reqSigPow =  10*log10(sigPower);
            reqSNR = EbNoVec(i)+10*log10(modord);
            rxsig = conv(txsigCyclic, chan);%, reqSNR, reqSigPow);
            rxsig = rxsig(1:end-length(chan)+1);
            rxsig = awgn(rxsig, reqSNR, reqSigPow);
            
            rxsig = reshape(rxsig, K+u, msgLen/K);
            rxsig = rxsig(u+1:end,:);
            
            rxmod = fft(rxsig, K, 1);
            chaninv = repmat(fft(chan', K), 1,msgLen/K);
            rxmodZF = rxmod ./ chaninv;
            
            No = reqSigPow - reqSNR;
            No = 10^(.1*No);
            rxmodMMSE = rxmod ./ (chaninv + No);
            
            rxmsgZF = qamdemod(rxmodZF, 2^modord, 0, 'gray');
            rxmsgMMSE = qamdemod(rxmodMMSE, 2^modord, 0, 'gray');
            
            [numerr, ber] = biterr(txmsg, rxmsgMMSE);
            BER(2,i) = ber;
            [numerr, ber] = biterr(txmsg, rxmsgZF);
            BER(1,i) = ber;
        end
        figure
        semilogy(EbNoVec,BER(1,:),'r', ...
            EbNoVec,BER(2,:),'g')
        title(['64 channel OFDM ', num2str(2^modord), '-QAM ' msg{k}]);
        legend('Zero-Forcing', 'MMSE', 'location', 'best');
        xlabel('SNR (dB)');
        ylabel('BER');
    end
end
%% Both
clear chans; clear BER;
chanlen = 3;
numchans = 3;
%Assuming a system with 2 transmit and 2 receive antennas close to each
%other
chans(1,1,1,:) = [1      .45     .6      ];
chans(1,2,1,:) = [1      .4      .65     ];
chans(2,1,1,:) = [.95    .35     .61     ];
chans(2,2,1,:) = [1.5    .35     .55     ];

%assuming a system with 2 transmit and 2 receive antennas with only the
%transmit close to each other
chans(1,1,2,:) = [1      .45     .6      ];
chans(1,2,2,:) = [1      .4      .65     ];
chans(2,1,2,:) = [1      .9      .4      ];
chans(2,2,2,:) = [1.5    .8      .45     ];

%Assuming all antennas have different paths
chans(1,1,3,:) = [1      .6     .8      ];
chans(1,2,3,:) = [1      .4     .65     ];
chans(2,1,3,:) = [.8     .2     0       ];
chans(2,2,3,:) = [.5     .9     .1      ];

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

