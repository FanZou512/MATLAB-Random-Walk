function [fitresult, gof] = DecayFit(t, FPTcump)
%CREATEFIT(T,FPTCUMP)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : t
%      Y Output: FPTcump
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, FPTcump );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-1.01 0 0.99];
opts.StartPoint = [-1 0.05 1];
opts.Upper = [-0.99 Inf 1.01];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
