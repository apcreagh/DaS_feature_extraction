function RMSE = eval_ArchimedeanSpiral( cen )

% ARGUMENT INPUTS
% cen(1) = x center
% cen(2) = y center
%
% GLOBAL INPUTS
% xdat = x data points
% ydat = y data points
%
% ARGUMENT OUTPUT
% rms = root-mean-square fit for data points
%
% GLOBAL OUTPUTS
% afit = radius coefficient

% global variables
global xdat ydat afit %bfit

% number of data points
ndat = length( xdat );

% center location
xc = cen( 1 );
yc = cen( 2 );

% convert to r,theta
x = xdat - xc;
y = ydat - yc;
r = sqrt( x.*x + y.*y );
theta = atan2( y, x );
theta = unwrap( theta );

%% 
% linearized LSQ fit for r = a * theta 
p = polyfit( theta, r, 1 );
afit = p(1);

% evaluate fit
rfit = (afit * theta );
dr = r - rfit;
rms = sqrt( dr'*dr / ndat );

% y=r;
% yhat=rfit;
% 
% er=(y - yhat);   % Errors
% er_sq=(y - yhat).^2;   % Squared Error
% MSE=mean((y - yhat).^2);   % Mean Squared Error
RMSE = sqrt(mean((r - rfit).^2));  % Root Mean Squared Error

end 
