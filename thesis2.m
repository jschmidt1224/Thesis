%% Recreation
% Recreating the resuslts from the paper Transmit without Regrets: Online
% Optimization in MIMO-OFDM Cognitive Radio Systems

clear all; clc;

% Initializing Variables
PrimaryUsers = 0; %num primary users
SecondaryUsers = 1; %num secondary users
Users = PrimaryUsers + SecondaryUsers; %total users
Antenna = 2; %num rx/tx antenna per user
T = 2000; %time to simulate
t = 1;
Subcarriers = 8; %num subcarriers
primaryP = cell(PrimaryUsers, Subcarriers); %power of primary
secondaryP = cell(SecondaryUsers, Subcarriers); %power of secondary
primaryQ = cell(PrimaryUsers, Subcarriers); %norm sig covar
secondaryQ = cell(SecondaryUsers, Subcarriers); %norm sig covar
noisePower = eye(Antenna); %noise power
H = cell(Users, Users, Subcarriers, T); %channel matrices
nu = 1;
y = cell(SecondaryUsers, Subcarriers);
Y = cell(SecondaryUsers, Subcarriers);
P = 10;
Iters = 100;
rates = zeros(Users,T);

for i = 1:Iters
    % Initialize y and Y
    for secondaryUser = 1:SecondaryUsers
        for subcarrier = 1:Subcarriers
            Y{secondaryUser, subcarrier} = zeros(Antenna);
            y{secondaryUser, subcarrier} = 0;
        end
    end
    
    % Initializing Primary Powers
    % Assuming each primary user only occupies one channel
    for primaryUser = 1:PrimaryUsers
        primchan = randi([1 Subcarriers]);
        for subcarrier = 1:Subcarriers
            if (subcarrier == primchan)
                primaryP{primaryUser, subcarrier} = 10;
            else
                primaryP{primaryUser, subcarrier} = 0;
            end
        end
    end
    
    % Initializing Secondary Powers
    % Each secondary user has access to all the channels but start with the
    % same power on each channel
    for secondaryUser = 1:SecondaryUsers
        for subcarrier = 1:Subcarriers
            secondaryP{secondaryUser, subcarrier} = P * exp(nu*t^(-1/2)*y{secondaryUser, subcarrier});
            divisor = 0;
            for chan = 1:Subcarriers
                divisor = divisor + exp(nu*t^(-1/2)*y{secondaryUser,chan});
            end
            secondaryP{secondaryUser, subcarrier} = secondaryP{secondaryUser, subcarrier} / divisor;
        end
    end
    
    
    
    % Initialize Primary Covariance Matrix
    % Constant over time
    for primaryUser = 1:PrimaryUsers
        for subcarrier = 1:Subcarriers
            primaryQ{primaryUser, subcarrier} = eye(Antenna);
        end
    end
    
    % Initialize Secondary Covariance Matrix
    % Starts as I but can change over time
    for secondaryUser = 1:SecondaryUsers
        for subcarrier = 1:Subcarriers
            secondaryQ{secondaryUser, subcarrier} = exp(nu*t^(-1/2)*Y{secondaryUser, subcarrier})/trace(exp(nu*t^(-1/2)*Y{secondaryUser, subcarrier}));
        end
    end
    
    % Initialize Channel Matrices
    % Assuming stationary for now
    % TODO: update to jakes when I understand that better
    %rng(19);
    for userFrom = 1:Users
        for userTo = 1:Users
            for subcarrier = 1:subcarrier
                tmpH = normrnd(0, 1, 2, 2)/sqrt(2) + 1i * normrnd(0, 1, 2, 2)/sqrt(2);
                for t = 1:T
                    H{userFrom, userTo, subcarrier, t} = normrnd(0, 1, 2, 2)/sqrt(2) + 1i * normrnd(0, 1, 2, 2)/sqrt(2);
                end
            end
        end
    end
    clear primaryUser secondaryUser subcarrier userFrom userTo tmpH;
    
    
    
    
    for t = 1:T
        for secondaryUser = 1:SecondaryUsers
            for subcarrier = 1:Subcarriers
                secondaryP{secondaryUser, subcarrier} = P * exp(nu*t^(-1/2)*y{secondaryUser, subcarrier});
                divisor = 0;
                for chan = 1:Subcarriers
                    divisor = divisor + exp(nu*t^(-1/2)*y{secondaryUser,chan});
                end
                secondaryP{secondaryUser, subcarrier} = secondaryP{secondaryUser, subcarrier} / divisor;
                secondaryQ{secondaryUser, subcarrier} = exp(nu*t^(-1/2)*Y{secondaryUser, subcarrier})/trace(exp(nu*t^(-1/2)*Y{secondaryUser, subcarrier}));
            end
        end
        userP = [primaryP; secondaryP];
        userQ = [primaryQ; secondaryQ];
        for userTo = 1:Users
            rate = 0;
            for subcarrier = 1:Subcarriers
                W = zeros(Antenna);
                for userFrom = 1:Users
                    if (userTo ~= userFrom)
                        tmpH = H{userFrom, userTo, subcarrier, t};
                        W = W + tmpH * (userP{userFrom, subcarrier} .* userQ{userFrom, subcarrier}) * tmpH';
                    end
                end
                W = W + noisePower;
                tmpH = H{userTo, userTo, subcarrier, t};
                tmpP = userP{userTo, subcarrier} * userQ{userTo, subcarrier};
                rateinc = log(det(W + tmpH*tmpP*tmpH')) - log(det(W));
                rate = rate + rateinc;
                if (userTo > PrimaryUsers)
                    h = W^(-1/2) * tmpH;
                    n = userTo - PrimaryUsers;
                    M = h'*(eye(Antenna) + h * tmpP * h')^-1*h;
                    y{n, subcarrier} = y{n,subcarrier} + P *trace(M*userQ{userTo, subcarrier});
                    Y{n,subcarrier} = Y{n, subcarrier} + userP{userTo, subcarrier}*M;
                end
            end
            rates(userTo, t) = rates(userTo, t) + rate;
            
        end
    end
    
    i
end
%%
figure
hold on
for u = 1:Users
    plot(200:T, abs(rates(u, 200:end))/i);
end
%ylim([0 11])

