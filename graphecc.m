close all; clear all;
load('P100-noBCH.mat');
BERno = BER;
throughputno = throughput;
load('511-502BCH-P100.mat');
figure
hold on
for u = 1:Users
    semilogy(1:T,502/511*max(throughput(:,u,:),[],3));
    semilogy(1:T,max(throughputno(:,u,:),[],3));
end
legend('ecc', 'no-ecc');

figure
hold on;
for u = 1:Users
    for m=1:Mords
        plot(1:T,BER(:,u,m,1));
        plot(1:T,BERno(:,u,m,1));
    end
end
legend('ecc1', 'no-ecc1','ecc2','no-ecc2','ecc3','no-ecc3');
ylim([0 .01])