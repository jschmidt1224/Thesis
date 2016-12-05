clear all; clc;

% Initializing Variables
PrimaryUsers = 0; %num primary users
SecondaryUsers = 1; %num secondary users
Users = PrimaryUsers + SecondaryUsers; %total users
Antenna = 2; %num rx/tx antenna per user
T = 100; %time to simulate
t = 1;
BCHn = 31;
BCHk = 6;
BCHblocks = 10;
MsgLen = BCHblocks*BCHk;
Subcarriers = 8; %num subcarriers

Mords = 5;
BER = zeros(T, Users, Mords, Subcarriers);
for M = 1:Mords
    Powers = zeros(1, Users, Subcarriers);
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
    %rng('shuffle')
    
    
    for t = 1:T
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
        
        GFtxmsg = gf(randi([0 1], Antenna, MsgLen*M, Users, Subcarriers));
        txmsg = GFtxmsg.x;
%         GFtxCoded = gf(zeros(Antenna, BCHblocks*BCHn*M, Users, Subcarriers));
%         for a=1:Antenna
%             for u=1:Users
%                 for s=1:Subcarriers
%                     reshaped = reshape(GFtxmsg(a,:,u,s), [], BCHk);
%                     GFtxCoded(a,:,u,s) = reshape(bchenc(reshaped, BCHn, BCHk), 1, [], 1, 1);
%                 end
%             end
%         end
%         txmsgCoded = GFtxCoded.x;
        
        %to turn off ECC change txmsgCoded to txmsg and go to 
        txmod = qammod(txmsg(:), 2^M, 'gray','InputType', 'bit', 'UnitAveragePower', true);
       
        txmod = reshape(txmod, Antenna, [], Users, Subcarriers);
        
        Len = size(txmod, 2);
        
        for userTo = 1:Users
            for subcarrier = 1:Subcarriers
                rxsig = zeros(Antenna, Len);
                for userFrom = 1:Users
                    rxsig = rxsig + P2*H(:,:,userFrom,userTo,subcarrier)*Powers(:,userFrom,subcarrier)*Q(:,:,userFrom,subcarrier)*txmod(:,:,userFrom,subcarrier);
                end
                
                rxsig = rxsig + noisePower * (1/sqrt(2)*randn(Antenna,Len) + 1/sqrt(2)*1i*randn(Antenna,Len));
                
                rxsig = pinv(H(:,:,userTo,userTo,subcarrier)).*1/P2*rxsig;
                
                rxmsg = qamdemod(rxsig,2^M,'gray','OutputType','bit','UnitAveragePower',true);
                rxmsg = reshape(rxmsg, Antenna, []);
                
%                 for a=1:Antenna
%                     reshaped = gf(reshape(rxmsg, [], BCHn));
%                     GFrxDecoded = reshape(bchdec(reshaped, BCHn, BCHk), Antenna, []);
%                 end
%                 rxmsgdecoded = GFrxDecoded.x;
                
                [numerr, ber] = biterr(txmsg(:,:,userTo,subcarrier), rxmsg);
                
                BER(t,userTo,M, subcarrier) = BER(t, userTo,M, subcarrier) + ber;
                
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
    end
end


figure
hold on;
throughput = zeros(T, Users, Mords);
for u = 1:Users
    for m = 1:Mords
        for s = 1:Subcarriers
            throughput(:,u,m) = throughput(:,u,m) + ((1-2*BER(:,u,m,s))*m);
        end
    end
end
%throughput = throughput * BCHk / BCHn;
for u = 1:Users
    semilogy(1:T,max(throughput(:,u,:),[],3));
end
%%
figure
hold on
for u = 1:Users
    for m=1:Mords
        plot(1:T,throughput(:,u,m))
    end
end

%%
figure
hold on;
for u = 1:Users
    for m=1:Mords
        plot(1:T,BER(:,u,m,1))
    end
end
%legend('primary','1','2','3','4','5','6')
%ylim([0 11])




