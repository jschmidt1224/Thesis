%% Recreation
% Recreating the resuslts from the paper Transmit without Regrets: Online 
% Optimization in MIMO-OFDM Cognitive Radio Systems

clear all; close all; clc;

% Number of primary users
Primary = 1;
% Number of transmit antennas for each primary user
PM = ones(1, Primary) * 3;
% Number of receive antennas for each primary user
PN = ones(1, Primary) * 3;
% Number of secondary users
Secondary = 2;
% Number of transmit antennas for each secondary user
SM = ones(1, Secondary) * 3;
% Number of receive antennas for each secondary user
SN = ones(1, Secondary) * 3;
% Total number of users in the system
Users = Primary + Secondary;
% Number of subcarriers in use
Subcarriers = 8;
% Average power matrices for the primary users
primaryPower = cell(Primary, 1);
% Average noise power for each user for each carrier
No = cell(Subcarriers, Users);
% Transmit matrices in the form H{channel, from, to}
H = cell(Subcarriers, Users, Users);

for k = 1:Subcarriers
    W = 0;
    for s = 2:Secondary
        
    end
    
    
end
