function plot_rmse_err( rErr, phiErr, avgSNR, algNames )
%PLOT_RMSE_ERR plots the RMSE of the different algorithms vs. SNR

[SNR,snrSrtI] = sort(mean(avgSNR,2));

figure
subplot(1,2,1); errorbar(mean(SNR,2),mean(rErr(snrSrtI,1,:),3),std(rErr(snrSrtI,1,:),[],3),':k','LineWidth',1); hold on
subplot(1,2,1); errorbar(mean(SNR,2),mean(rErr(snrSrtI,2,:),3),std(rErr(snrSrtI,2,:),[],3),'b','LineWidth',1.5);
subplot(1,2,1); errorbar(mean(SNR,2),mean(rErr(snrSrtI,3,:),3),std(rErr(snrSrtI,3,:),[],3),'g','LineWidth',1.5);
subplot(1,2,1); errorbar(mean(SNR,2),mean(rErr(snrSrtI,4,:),3),std(rErr(snrSrtI,4,:),[],3),'Color',[0.5,0.5,0.5],'LineWidth',1.5);
hold off;
set(gca,'yscale','log'); grid on;
xlabel('SNR','FontSize',14); title('$$r$$','Interpreter','Latex','FontSize',18); ylabel('RMSE (log scale)','FontSize',14)

subplot(1,2,2); errorbar(mean(SNR,2),mean(phiErr(snrSrtI,1,:),3),std(phiErr(snrSrtI,1,:),[],3),':k','LineWidth',1); hold on
subplot(1,2,2); errorbar(mean(SNR,2),mean(phiErr(snrSrtI,2,:),3),std(phiErr(snrSrtI,2,:),[],3),'b','LineWidth',1.5);
subplot(1,2,2); errorbar(mean(SNR,2),mean(phiErr(snrSrtI,3,:),3),std(phiErr(snrSrtI,3,:),[],3),'g','LineWidth',1.5);
subplot(1,2,2); errorbar(mean(SNR,2),mean(phiErr(snrSrtI,4,:),3),std(phiErr(snrSrtI,4,:),[],3),'Color',[0.5,0.5,0.5],'LineWidth',1.5);
hold off;
set(gca,'yscale','log'); grid on;
xlabel('SNR','FontSize',14); title('$$\phi$$','Interpreter','Latex','FontSize',18); ylabel('RMSE (log scale)','FontSize',14)

legend(algNames,'FontSize',12);


end

