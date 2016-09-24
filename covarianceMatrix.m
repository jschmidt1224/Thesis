
clear all; close all; clc;

M = 8;
msg = zeros(2, 64);%[0 0; 0 1; 0 2; 0 3; 1 0; 1 1; 1 2; 1 3; 2 0; 2 1; 2 2; 2 3; 3 0; 3 1; 3 2; 3 3]';
        
for i = 1:64
    msg(:,i) = [floor((i-1)/8); mod(i-1,8)];
end
tx = qammod(msg, M);
total = 0
totalinv = 0
for i = 1:64
    total = total + tx(:,i) * tx(:,i)';
    totalinv = totalinv + tx(:,i)' * tx(:,i);
end
total = total / 64
totalinv = totalinv / 64
%%
msg = [0 1 2 3]
tx = qammod(msg, 4)
total2 = 0;
for i=1:4
    total2 = total2 + tx(i)*tx(i)';
end
total2 = total2 / 4

%%
n = 100000;
total3 = 0;
for i=1:n
    x = [normrnd(0, 1)/sqrt(2) + 1i*normrnd(0,1)/sqrt(2); normrnd(0, 1)/sqrt(2) + 1i*normrnd(0,1)/sqrt(2)];
    total3 = total3 + x * x';
end
total3 ./ n
    
    
    
    
    