function [best_center, RMSE]= fit_ArchimedeanSpiral(x1, y1, centre)
%% Draw-a-Shape: Function to Fit an Emperical Archimedean Spiral
% Emperically:  r = afit*theta  x = xfit+r*cos(theta) y = yfit+r*sin(theta)
%--------------------------------------------------------------------------
% Input:
%      xdat = column vector of x locations along spiral
%      ydat = column vector of y locations along spiral
% 
% =========================================================================
% Output: 
%        afit (removed, [08-05-18]): radius coefficient
%        xfit: x center
%        yfit: y center
%        rmse: root-mean-square of residuals for data points
%_________________________________________________________________________%
%% 
options = optimset('MaxFunEvals',1000, 'Display','off');
xdat = x1;
ydat = y1;
ndat = length( xdat );
if isnan(centre)
% geometric center based on constant heading change
% % http://mathworld.wolfram.com/ArchimedeanSpiral.html
dx = diff(xdat);
dy = diff(ydat);
s = sqrt( dx.*dx + dy.*dy );
heading = unwrap( atan2(dy,dx) );
dphi = diff( heading );
% bapprox = atan( dphi/2 );
r = s(1:ndat-2) .* tan( pi/2-dphi );
xc = xdat(1:ndat-2) + r .* cos( heading(1:ndat-2) + pi/2 );
yc = ydat(1:ndat-2) + r .* sin( heading(1:ndat-2) + pi/2 );
xc = mean( xc );
yc = mean( yc );
cen = [ xc  yc ]';

[ best_center, RMSE] = fminsearch( 'eval_ArchimedeanSpiral', cen, options);
%afit;
xfit = best_center(1);
yfit = best_center(2);
else 
%     
cen = centre;
[ best_center, RMSE] = fminsearch( 'eval_ArchimedeanSpiral', cen );
afit;
xfit=centre(1);
yfit=centre(2);
end 
% search for best fit
% plot results
x = xdat - xfit;
y = ydat - yfit;
theta = unwrap( atan2( y, x ) );
th_max = max( theta );
th_min = min( theta );
dth = ( th_max - th_min ) / 100;
th = (th_min : dth : th_max )';
%rfit = (afit * th );

% figure
%   clf
%   plot( xdat, ydat, 'ro', rfit.*cos(th)+xfit, rfit.*sin(th)+yfit, 'k', xfit, yfit, 'b+' )
%     axis equal tight
%     hold on
%     box off
%     hold on
%     set(gca,'XTick',[])
%     hold on
%     set(gca,'YTick',[])
%     title( 'Best Fit of Archimedean Spiral' )
end 
%EOF