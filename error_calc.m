% Compute measurement error
% ***************************************************************@

function [rErr, phiErr] = error_calc(data, tt)
%MEAS_ERROR computes the measurement error
% 'tt' containes the samples which will be included in the error
% calculation.

Tr_y   = data.Tr_y;
Est_y  = data.Est_y;
rErr   = sqrt(mean((Est_y(2,tt)-Tr_y(2,tt)).^2))/std(Tr_y(2,tt));
phiErr = sqrt(mean((Est_y(1,tt)-Tr_y(1,tt)).^2))/std(Tr_y(1,tt));



