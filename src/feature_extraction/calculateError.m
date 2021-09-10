function [AUC]=calculateError(T, R)
        
normR = R - min(R(:));
normR = normR ./ max(normR(:)); 

normTH = T - min(T(:));
normTH = normTH ./ max(normTH(:)); 

refx=linspace(0,1,length(normTH));
refy=linspace(0,1,length(normR));

dTH=trapz(normTH, refx);
dR=trapz(normR, refy);

% Q=[normTH, normR];
AUC=abs(dTH)+abs(dR);

% figure
% plot(refx, normTH, 'b')
% hold on
% plot(refy, normR, 'r')
% hold on
% legend('X-Coordinates', 'Y-Coordinates')
 end 