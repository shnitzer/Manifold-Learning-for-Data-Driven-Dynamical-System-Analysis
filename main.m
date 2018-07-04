% MATLAB code implementation of the non-Linear object tracking example 
% from: "Manifold Learning for Data-Driven Dynamical System Analysis",
% Chapter in The Koopman Operator in Systems and Control, Springer.
% ***************************************************************@
% This implementation generates the underlying diffusion processes
% and the corresponding measurements of the non-linear tracking problem and 
% recovers the underlying processes using the DMK and observer methods.
% Author:   Tal Shnitzer.
% Created:  7/4/18.
% ***************************************************************@

%% Configuration
% Signal generation parameters:
len      = 1000;                % signal length
deltaT   = 0.01;                % time step
procStd  = sqrt(2);             % standard deviation of the process noise
noiseStd = [1,sqrt(1.5),sqrt(2),sqrt(2.5),sqrt(3)];       % relative standard deviation of noise
DMdim    = 2;                   % dimensions of the diffusion maps coordinates to use in the Kalman filter
iteNum   = 10;                  % number of process generation iterations for RMSE calulations
jumpErr  = 60;                  % processes with changes larger than jumpErr will be ignored

AlgNo    = 4;                   % number of compared algorithms (dmk, noisy measurements, EKF)

%% Initialization

rErr      = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = measurement error, (:,2,:) = dmk, (:,3,:) = observer, (:,4,:) = optimal particle filter
phiErr    = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = measurement error, (:,2,:) = dmk, (:,3,:) = observer, (:,4,:) = optimal particle filter
theta1Err = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = measurement error, (:,2,:) = dmk, (:,3,:) = observer, (:,4,:) = optimal particle filter
theta2Err = zeros(length(noiseStd),AlgNo,iteNum); % (:,1,:) = measurement error, (:,2,:) = dmk, (:,3,:) = observer, (:,4,:) = optimal particle filter

avgSNR = zeros(length(noiseStd),iteNum);

%% Generate Process
% save('data.mat','yM','yT','theta','procStd','noiseStd','deltaT','InitLoc');
h = waitbar(0, 'Please wait while the errors for different noise-STD are computed');
for nstd = 1:length(noiseStd)
    waitbar(nstd/length(noiseStd), h);
    noiseStd_curr = noiseStd(nstd);

    for iter = 1:iteNum
        
        goodProcess = 0;
        % Generating a process without discontinuities caused by arctan or
        % discontinuities larger than 'jumpErr':
        while ~goodProcess
            % Generate new underlying process:
            InitLoc    = 1*randn(2,1)+[1; 5];                                     % initial process location
            DriftRate1 = @(t,X) -0.5*(X-1).^3+(X-1);                         % set drift parameters
            DriftRate2 = @(t,X) -0.5*(X-6).^3+(X-6);                         % set drift parameters
            DiffRate1  = @(t,X) procStd;                                    % set diffusion parameters
            DiffRate2  = @(t,X) procStd;                                    % set diffusion parameters
            SDE1       = sde(DriftRate1, DiffRate1, 'StartState', InitLoc(1));    % define SDE
            [thet1, ~] = SDE1.simulate(len-1, 'DeltaTime', deltaT);           % simulate process
            SDE2       = sde(DriftRate2, DiffRate2, 'StartState', InitLoc(2));    % define SDE
            [thet2, t] = SDE2.simulate(len-1, 'DeltaTime', deltaT);           % simulate process
            
            theta = [thet1, thet2];
            % Generate measurements:
            phiT = atand(theta(:,1)./theta(:,2)).';       % clean angle values
            rT   = sqrt(theta(:,1).^2 + theta(:,2).^2).'; % clean radius values
            
            % Ignoring processes which contain significant discontinuities due
            % to arctan and drastic chanes
            if ~any(abs(diff(phiT))>60) && ~any(theta(:,2)<0)
                goodProcess = 1;
            end
        end
        phiM = phiT + noiseStd_curr*std(phiT) * randn(size(phiT)); % noisy angle
        rM   = rT   + noiseStd_curr*std(rT)   * randn(size(rT));   % noisy radius
        
        avgSNR(nstd,iter) = ( var(rT) + var(phiT) ) / ( noiseStd_curr.^2*var(phiT) + noiseStd_curr.^2*var(rT) );
        
        yT = [phiT; rT];
        yM = [phiM; rM];
        
        tt = 100:size(yM,2); % samples to consider - ignoring the first samples in the error calculations due to initialization effect errors
        
        %% Constrcut diffusion maps with the modified mahalanobis distance
        
        % Compute the modified mahalanobis distance for the (noisy) measurements:
        mahDist       = modified_mahalanobis(yM);
        
        % Compute diffusion maps coordinates and eigenvalues:
        [psi, lambda] = diffusion_maps(mahDist, DMdim);
        
        %% DMK and Observer frameworks:
        
        % Apply DMK:
        [~, yDMK_est] = dmk(psi, lambda, yM, deltaT, tt);
        
        % Computing the errors for DMK:
        dmk_data.Tr_y  = yT; 
        dmk_data.Est_y = yDMK_est;
        aI = 2;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter)] = error_calc(dmk_data, tt);
        
        % Apply Observer:
        yObs_est = observer(psi, lambda, yM, deltaT, tt);
        % Computing the errors for the observer:
        Obs_data.Tr_y  = yT; 
        Obs_data.Est_y = yObs_est;
        aI = 3;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter)] = error_calc(Obs_data, tt);
        
        
        %% Compute estimations of compared algorithms and noisy measurements
        
        % Measurement error:
        meas_data.Tr_y  = yT; 
        meas_data.Est_y = yM;
        aI = 1;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter)] = error_calc(meas_data, tt);
        
        % Particle Filter:
        y_est_pf = particle_filter( yM, DriftRate1, DriftRate2, deltaT, noiseStd_curr*std(yT,[],2), procStd, InitLoc );
        % Particle Filter error:
        pf_data.Tr_y  = yT;
        pf_data.Est_y = y_est_pf;
        aI = 4;
        [rErr(nstd,aI,iter), phiErr(nstd,aI,iter)] = error_calc(pf_data, tt);
        
    end
    
end
close(h);

%% Displaying the results and comparing to other algorithms:

plot_rmse_err( rErr, phiErr, avgSNR, {'Meas.','DMK','Observer','PF'} )


