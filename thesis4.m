clear all; clc;

% Initializing Variables
PrimaryUsers = 0; %num primary users
SecondaryUsers = 1; %num secondary users
Users = PrimaryUsers + SecondaryUsers; %total users
Antenna = 2; %num rx/tx antenna per user
t = 1;
BCHn = 511;
BCHk = 313;
BCHblocks = 10;
BCHPren = 31;
BCHPrek = 11;
Packets = 1;
Subcarriers = 8; %num subcarriers
T = Packets * (BCHn+BCHPren);
Iters = 1;

PowerBits = 7;
ModBits = 2;
CodeBits = 2;

Mords = ones(1, Users, Subcarriers);
Codes = zeros(CodeBits, Users, Subcarriers);
throughput = zeros(T, Users, Subcarriers);

Powers = zeros(1, Users, Subcarriers);
PowersRounded = zeros(1, Users, Subcarriers);
PowersRX = zeros(1, Users, Subcarriers);
MRX = zeros(1,Users,Subcarriers);
Q = zeros(Antenna, Antenna, Users, Subcarriers);
noisePower = eye(Antenna); %noise power
H = zeros(Antenna, Antenna, Users, Users, Subcarriers, T);
nu = 1;
y = zeros(1, SecondaryUsers, Subcarriers);
Y = zeros(Antenna, Antenna, SecondaryUsers, Subcarriers);
%P = 10;
P2 = 10;
t=1;
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
        Powers(:,user, subcarrier) = exp(nu*t^(-1/2)*y(:,user-PrimaryUsers, subcarrier));
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

rng(19);
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

%M = 1;
t = 1;
for p = 1:Packets
    PowersRounded = round(Powers*2^PowerBits);
    txpre = zeros(Antenna, 1, Users, Subcarriers,(PowerBits+ModBits+CodeBits)*(BCHPren/BCHPrek));
    for u=1:Users
        for s=1:Subcarriers
            BinaryPower = de2bi(PowersRounded(:,u,s), PowerBits);
            BinaryOrd = de2bi(Mords(:,u,s), ModBits);
            Binary = gf([BinaryPower, BinaryOrd, Codes(:,u,s)']);
            Binary = reshape(Binary, [], BCHPrek);
            pre = bchenc(Binary, BCHPren, BCHPrek);
            for a=1:Antenna
                txpre(a,1,u,s,:) = qammod(pre.x, 2^Mords(:,u,s), 'gray','InputType', 'bit', 'UnitAveragePower', true);
            end
        end
    end
    GFtxmsg = gf(randi([0 1], Antenna, 1, Users, Subcarriers, BCHk));
    txmsg = GFtxmsg.x;
    GFtxCoded = gf(zeros(Antenna, 1, Users, Subcarriers, BCHn));
    for a=1:Antenna
        for u=1:Users
            for s=1:Subcarriers
                reshaped = reshape(GFtxmsg(a,1,u,s,:), [], BCHk);
                GFtxCoded(a,1,u,s,:) = reshape(bchenc(reshaped, BCHn, BCHk), 1, 1, 1, 1, []);
            end
        end
    end
    
    txmsgCoded = GFtxCoded.x;
    txmod = zeros(Antenna, 1, Users, Subcarriers, BCHn);
    for u=1:Users
        for s=1:Subcarriers
            tmpmsg = txmsgCoded(:,1,u,s,:);
            tmpsig = qammod(tmpmsg(:), 2^Mords(:,u,s), 'gray','InputType', 'bit', 'UnitAveragePower', true);
            txmod(:,1,u,s,:) = reshape(tmpsig, Antenna, 1, 1, 1, []);
        end
    end
    
    txsigmaster = repmat(cat(5, txpre, txmod), 1, Iters, 1, 1, 1);
    Len = BCHn + BCHPren;
    rxMaster = zeros(Antenna, Iters, Users, Subcarriers, BCHn+BCHPren);
    %rxPre = zeros(1,1,Users,Subcarriers,BCHPren);
    rxPreDec = zeros(1,1,Users,Subcarriers,BCHPrek);
    rxDec = zeros(Antenna,1,Users,Subcarriers,BCHk);
    for b = 1:(BCHn+BCHPren)
        txsig = txsigmaster(:,:,:,:,b);
        t
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
        
        if (b == BCHPren + 1)
            for u=1:Users
                for s=1:Subcarriers
                    rxpre = gf(mode(mode(rxMaster(:,:,u,s,1:BCHPren),2)));
                    rxpredec = bchdec(reshape(rxpre, 1, BCHPren),BCHPren,BCHPrek);
                    rxPreDec(1,1,u,s,:) = rxpredec.x;
                    PowersRX(:,u,s) = double(bi2de(reshape(rxPreDec(1,1,u,s,1:PowerBits),1,[])))./2^PowerBits;
                    MRX(1,u,s) = double(bi2de(reshape(rxPreDec(1,1,u,s,PowerBits+1:BCHPrek-CodeBits),1,[])));
                    %CodeRX(1,u,s) = rxpre(end-CodeBits+1:end);
                end
            end
        end
        for userTo = 1:Users
            for subcarrier = 1:Subcarriers
                rxsig = zeros(Antenna, Iters);
                for userFrom = 1:Users
                    rxsig = rxsig + P2*H(:,:,userFrom,userTo,subcarrier)*PowersRounded(:,userFrom,subcarrier)*Q(:,:,userFrom,subcarrier)*txsig(:,:,userFrom,subcarrier);
                end
                rxsig = rxsig + noisePower * (1/sqrt(2)*randn(Antenna,Iters) + 1/sqrt(2)*1i*randn(Antenna,Iters));
                rxsig = pinv(H(:,:,userTo,userTo,subcarrier)).*1/P2*rxsig;
                if (b <= (BCHPren))
                    rxmsg = qamdemod(rxsig, 2, 'gray','OutputType','bit','UnitAveragePower',true);
                    rxMaster(:,:,userTo,subcarrier,b) = rxmsg;
                else
                    rxsig = rxsig ./ PowersRX(1,userTo,subcarrier);
                    rxmsg = qamdemod(rxsig, 2^MRX(1,userTo,subcarrier), 'gray','OutputType','bit','UnitAveragePower',true);
                    rxMaster(:,:,userTo,subcarrier,b) = rxmsg;
                end
                
%                 %rxpre = gf(qamdemod(rxsig(1,:), 2, 'gray','OutputType','bit','UnitAveragePower',true));
%                 [rxpredecoded, numerr] = bchdec(rxpre, BCHPren, BCHPrek);
%                 numerr
%                 if (numerr == -1)
%                     rxpre = gf(qamdemod(rxsig(2,1:BCHPren), 2, 'gray','OutputType','bit','UnitAveragePower',true));
%                     [rxpredecoded, numerr] = bchdec(rxpre, BCHPren, BCHPrek);
%                     if (numerr == -1)
%                         continue;
%                     end
%                 end
%                 
%                 rxpre = rxpredecoded.x;
%                 Power = double(bi2de(rxpre(1:PowerBits)))/2^PowerBits;
%                 M = double(bi2de(rxpre(PowerBits+1:BCHPrek-CodeBits)));
%                 Code = rxpre(end-CodeBits+1:end);
%                 
%                 rxmsg = qamdemod(rxsig(:,BCHPren+1:end),2^M,'gray','OutputType','bit','UnitAveragePower',true);
%                 %rxmsg = reshape(rxmsg, Antenna, []);
%                 
%                 for a=1:Antenna
%                     reshaped = gf(reshape(rxmsg, [], BCHn));
%                     [GFrxDecoded, numerr] = bchdec(reshaped, BCHn, BCHk);
%                     %GFrxDecoded = reshape(bchdec(reshaped, BCHn, BCHk), Antenna, []);
%                 end
%                 rxmsgdecoded = GFrxDecoded.x;
%                 numerr
%                 if (numerr ~= -1)
%                     throughput(t,userTo,subcarrier) = BCHk / (BCHn + BCHPren);
%                 end
                
                
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
                if (userTo > PrimaryUsers)
                    h = W^(-1/2) * tmpH;
                    n = userTo - PrimaryUsers;
                    M1 = h'*(eye(Antenna) + h * tmpP * h')^-1*h;
                    y(:,n,subcarrier) = y(:,n,subcarrier) + trace(M1*Q(:,:,userTo,subcarrier));
                    Y(:,:,n,subcarrier) = Y(:,:,n, subcarrier) + Powers(:,userTo,subcarrier)*M1;
                end
            end
            
            
        end
        if (b == BCHn+BCHPren)
            for u=1:Users
                for s=1:Subcarriers
                    rx = gf(mode(rxMaster(:,:,u,s,BCHPren+1:end),2));
                    rxdec = bchdec(reshape(rx, [], BCHn),BCHn,BCHk);
                    rxDec(:,1,u,s,:) = reshape(rxdec.x, Antenna,1,1,1,[]);
                end
            end
        end
        t = t+1;
    end
end


figure
hold on;
for u = 1:Users
    semilogy(1:T,sum(throughput(:,u,:),3));
end
%%
for u = 1:Users
    for m=1:Mords
        plot(1:T,throughput(:,u,m))
    end
end

%%
figure
hold on;
for u = 1:Users
    plot(1:T,BER(:,u,1))
end
%legend('primary','1','2','3','4','5','6')
%ylim([0 11])




