%%
figure
%hold on
for u = 1:Users
    for m = 1:Mords
        semilogy(SNRset, BERSNR(:,u,m));
        hold on
    end
end
