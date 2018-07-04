function [ y_est ] = observer( psi, lambda, yM, deltaT, tt )
%OBSERVER Constructs the diffusion maps based observer framework

gamma  = 0.01;      % Gain parameter - empirically chosen
alpha  = yM*psi;    % Lift Function

% Calculate the new representation using the observer framework:
psi_hat      = zeros(size(psi,1), size(psi,2));
psi_hat(1,:) = psi(1,:);
for i = 1:length(psi)-1
    psi_hat(i+1,:) = (1-gamma)*(psi_hat(i,:) + (diag(lambda)*psi_hat(i,:).').') + (gamma*pinv(alpha)*yM(:,i+1)).';
end

% Lift the new representation back to the measurement space:
y_est = alpha * (psi_hat.');

% Scaling the data to fit the measurements - required due to gain issues:
y_est = bsxfun(@minus,y_est,mean(y_est,2));
y_est = bsxfun(@rdivide,y_est,std(y_est,0,2)./std(yM,0,2));
y_est = bsxfun(@plus,y_est, mean(yM,2));

end

