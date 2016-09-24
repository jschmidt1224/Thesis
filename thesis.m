%% Recreation
% Recreating the resuslts from the paper Transmit without Regrets: Online 
% Optimization in MIMO-OFDM Cognitive Radio Systems

clear all; close all; clc;

% Number of primary users
Primary = 1;
% Number of transmit antennas for each primary user
PM = ones(1, Primary) * 2;
% Number of receive antennas for each primary user
PN = ones(1, Primary) * 2;
% Number of secondary users
Secondary = 3;
% Number of transmit antennas for each secondary user
SM = ones(1, Secondary) * 2;
% Number of receive antennas for each secondary user
SN = ones(1, Secondary) * 2;
% Total number of users in the system
Users = Primary + Secondary;
% Number of time frames to simulate
T = 10;
% Number of subcarriers in use
Subcarriers = 8;
% Average power matrices for the primary users
primaryPower = cell(Primary, Subcarriers);
% Average power matrices for the secondary users
secondaryPower = cell(Primary, Subcarriers);
% Average noise power for each user for each carrier
No = cell(Subcarriers, Users);
% Channel matrices in the form H{channel, from, to}
H = cell(Subcarriers, Users, Users);
% Initializing H
for from = 1 : Users
    for to = 1 : Users
        for sub = 1 : Subcarriers
            %TODO: Fix this to be affected by antenna variables above
            H{sub, from, to} = normrnd(0, 1, 2, 2)/sqrt(2) + ...
                               1i * normrnd(0, 1, 2, 2)/sqrt(2);
        end
    end
end


for k = 1:Subcarriers
    W = 0;
    for s = 2:Secondary
        W = W + 
    end
    
    
end
