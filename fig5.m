function fig5
%FIG5 creates the plot in Fig. 5 of the book chapter

load('data.mat')

DMdim  = 2;              % dimensions of the diffusion maps coordinates to use in the Kalman filter
deltaT = 0.01;           % time step
tt     = 100:size(yM,2); % samples to consider - ignoring the first samples in the error calculations due to initialization effect errors

% Plot state trajectory:
figure
plot(theta(:,1),theta(:,2),'k')
grid on; axis equal
xlabel('$$x_1$$','Interpreter','latex','FontSize',16)
ylabel('$$x_2$$','Interpreter','latex','FontSize',16)

%% Constrcut diffusion maps with the modified mahalanobis distance

% Compute the modified mahalanobis distance for the (noisy) measurements:
mahDist       = modified_mahalanobis(yM);

% Compute diffusion maps coordinates and eigenvalues:
[psi, lambda] = diffusion_maps(mahDist, DMdim);

%% DMK framework:

% Apply DMK:
[~, yDMK_est] = dmk(psi, lambda, yM, deltaT, tt);

%% Compute estimations of the particle filter:

DriftRate1 = @(t,X) -0.5*(X-1).^3+(X-1); % set drift parameters
DriftRate2 = @(t,X) -0.5*(X-6).^3+(X-6); % set drift parameters

% Particle Filter:
y_est_pf = particle_filter( yM, DriftRate1, DriftRate2, deltaT, noiseStd*std(yT,[],2), procStd, InitLoc );


%% Plot trajectory examples - for each SNR:
    figure
    subplot(2,1,1), scatter(tt*deltaT,yM(2,tt),20,[0.6,0.6,0.6],'x');
    hold on
    xlabel('t [sec]','FontSize',14); ylabel('$$r$$','Interpreter','latex','FontSize',16)
    subplot(2,1,1), plot(tt*deltaT,yT(2,tt),':k','LineWidth',2);
    subplot(2,1,1), plot(tt*deltaT,y_est_pf(2,tt),'Color',[0.5,0.5,0.5],'LineWidth',2);
    subplot(2,1,1), plot(tt*deltaT,yDMK_est(2,tt),'b','LineWidth',1);
    lgd = legend('Noisy meas.','Clean meas.','PF estimation','DMK estimation');
    lgd.FontSize = 12; 
    xlim([2.6,5.4]); ylim([-2,12]);
    hold off;
    
    subplot(2,1,2), scatter(tt*deltaT,yM(1,tt),20,[0.6,0.6,0.6],'x');
    xlabel('t [sec]','FontSize',14); ylabel('$$\phi$$','Interpreter','latex','FontSize',16)
    hold on
    subplot(2,1,2), plot(tt*deltaT,yT(1,tt),':k','LineWidth',2);
    subplot(2,1,2), plot(tt*deltaT,y_est_pf(1,tt),'Color',[0.5,0.5,0.5],'LineWidth',2);
    subplot(2,1,2), plot(tt*deltaT,yDMK_est(1,tt),'b','LineWidth',1);
    xlim([2.6,5.4]); ylim([-30,30]);
    hold off;
   
end

