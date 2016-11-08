clear all; close all; clc;

% Initializing Variables
PrimaryUsers = 0; %num primary users
SecondaryUsers = 1; %num secondary users
Users = PrimaryUsers + SecondaryUsers; %total users
Antenna = 2; %num rx/tx antenna per user
T = 65; %time to simulate
t = 1;
MsgLen = 10000000;
Subcarriers = 8; %num subcarriers
Powers = zeros(1, Users, Subcarriers);
Q = zeros(Antenna, Antenna, Users, Subcarriers);
noisePower = eye(Antenna); %noise power
H = zeros(Antenna, Antenna, Users, Users, Subcarriers, T);
nu = 1;
y = zeros(1, SecondaryUsers, Subcarriers);
Y = zeros(Antenna, Antenna, SecondaryUsers, Subcarriers);
P = 8;
P2 = 10000;
M = 2;

% Initializing Primary Powers
% Assuming each primary user only occupies one channel
for user = 1:PrimaryUsers
    primchan = randi([1 Subcarriers]);
    for subcarrier = 1:Subcarriers
        if (subcarrier == primchan)
            Powers(:,user, subcarrier) = 1;
        else
            Powers(:,user, subcarrier) = 0;
        end
    end
end

% Initializing Secondary Powers
% Each secondary user has access to all the channels but start with the
% same power on each channel
for user = PrimaryUsers+1:Users
    for subcarrier = 1:Subcarriers
        Powers(:,user, subcarrier) = P * exp(nu*t^(-1/2)*y(:,user-PrimaryUsers, subcarrier));
        divisor = 0;
        for chan = 1:Subcarriers
            divisor = divisor + exp(nu*t^(-1/2)*y(:,user-PrimaryUsers,chan));
        end
        Powers(:,user, subcarrier) = Powers(:,user, subcarrier) / divisor;
    end
end


% Initialize Primary Covariance Matrix
% Constant over time
for user = 1:PrimaryUsers
    for subcarrier = 1:Subcarriers
        Q(:,:,user, subcarrier) = eye(Antenna);
    end
end

% Initialize Secondary Covariance Matrix
% Starts as I but can change over time
for user = PrimaryUsers+1:Users
    for subcarrier = 1:Subcarriers
        Q(:,:,user, subcarrier) = exp(nu*t^(-1/2)*Y(:,:,user-PrimaryUsers, subcarrier))/trace(exp(nu*t^(-1/2)*Y(:,:,user-PrimaryUsers, subcarrier)));
    end
end

rng(18);
% Initialize Channel Matrices
for userFrom = 1:Users
    for userTo = 1:Users
        for subcarrier = 1:subcarrier
            tmpH = normrnd(0, 1, 2, 2)/sqrt(2) + 1i * normrnd(0, 1, 2, 2)/sqrt(2);
            for t = 1:T
                H(:,:,userFrom, userTo, subcarrier, t) = normrnd(0, 1, 2, 2)/sqrt(2) + 1i * normrnd(0, 1, 2, 2)/sqrt(2);
            end
        end
    end
end
clear primaryUser secondaryUser subcarrier userFrom userTo tmpH;
rng('shuffle')

BER = zeros(T, Users);
for t = 1:T
    for user = PrimaryUsers+1:Users
        for subcarrier = 1:Subcarriers
            Powers(:,user, subcarrier) = exp(nu*t^(-1/2)*y(:,user-PrimaryUsers, subcarrier));
            divisor = 0;
            for chan = 1:Subcarriers
                divisor = divisor + exp(nu*t^(-1/2)*y(:,user-PrimaryUsers,chan));
            end
            Powers(:,user, subcarrier) = Powers(:,user, subcarrier) / divisor;
            Q(:,:,user, subcarrier) = exp(nu*t^(-1/2)*Y(:,:,user-PrimaryUsers, subcarrier))/trace(exp(nu*t^(-1/2)*Y(:,:,user-PrimaryUsers, subcarrier)));
        end
    end
    
    txmsg = randi([0 2^M-1], 2, MsgLen, Users, Subcarriers);
    txmod = qammod(txmsg(:), 2^M, 0, 'gray');
    txmod = reshape(txmod, 2, MsgLen, Users, Subcarriers);
    for userTo = 1:Users
        for subcarrier = 1:Subcarriers
            rxsig = zeros(2, MsgLen);
            for userFrom = 1:Users
                rxsig = rxsig + P2*H(:,:,userFrom,userTo,subcarrier)*Powers(:,userFrom,subcarrier)*Q(:,:,userFrom,subcarrier)*txmod(:,:,userFrom,subcarrier);
            end
            rxsig = rxsig + noisePower * (1/sqrt(2)*randn(2,MsgLen) + 1/sqrt(2)*1i*randn(2,MsgLen));
            rxmsg = qamdemod(pinv(H(:,:,userTo,userTo,subcarrier))*1/P2*rxsig,2^M,0,'gray');
            [numerr, ber] = biterr(txmsg(:,:,userTo,subcarrier), rxmsg);
            BER(t,userTo) = BER(t, userTo) + ber;
        end
        for subcarrier = 1:Subcarriers
            W = zeros(Antenna);
            for userFrom = 1:Users
                if (userTo ~= userFrom)
                    tmpH = H(:,:,userFrom,userTo,subcarrier,t);
                    W = W + tmpH * (Powers(:,userFrom,subcarrier) .* Q(:,:,userFrom,subcarrier)) * tmpH';
                end
            end
            W = W + noisePower;
            tmpH = H(:,:,userTo,userTo,subcarrier,t);
            tmpP = Powers(:,userTo,subcarrier) * Q(:,:,userTo,subcarrier);
            %rateinc = log(det(W + tmpH*tmpP*tmpH')) - log(det(W));
            %rate = rate + rateinc;
            if (userTo > PrimaryUsers)
                h = W^(-1/2) * tmpH;
                n = userTo - PrimaryUsers;
                M1 = h'*(eye(Antenna) + h * tmpP * h')^-1*h;
                y(:,n,subcarrier) = y(:,n,subcarrier) + P *trace(M1*Q(:,:,userTo,subcarrier));
                Y(:,:,n,subcarrier) = Y(:,:,n, subcarrier) + Powers(:,userTo,subcarrier)*M1;
            end
        end
    end
end

BER = BER / Subcarriers;
figure
hold on;
for u = 1:Users
    semilogy(1:T,BER(:,u));
end
legend('primary', 'secondary');