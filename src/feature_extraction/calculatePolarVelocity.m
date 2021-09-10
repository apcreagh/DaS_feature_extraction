function [AV, timeAV, RV, timeRV, mag_error]=calculatePolarVelocity(x, y, t, xref, yref)

[THETA,R] = cart2pol(x,y);
    
%     [idx0,~] = crossing(THETA);
%     THETA(idx0(1):idx0(end))=[];
%     R(idx0(1):idx0(end))=[];
    
THETA=unwrap(THETA);

% figure
[AUC1]=calculateError(THETA, R);
[Tref,Rref] = cart2pol(xref,yref);
THETA=unwrap(THETA);
[AUC2]=calculateError(Tref, Rref);
% hold on

% mag_error=abs(AUC2-AUC1);
mag_error=abs(AUC2+AUC1);

dR=diff(R);
dt=diff(t);
RV=dR./dt;

%ELIMINATE NaNs

timeRV=t;
NanIndex=find(isnan(RV));
RV(NanIndex)=[];
timeRV(NanIndex)=[];

INFIndex=find(isinf(RV));
RV(INFIndex)=[];
timeRV(INFIndex)=[];

dTheta=diff(THETA);
AV=dTheta./dt;

%ELIMINATE NaNs
timeAV=t;
NanIndex=find(isnan(RV));
AV(NanIndex)=[];
timeAV(NanIndex)=[];

INFIndex=find(isinf(AV));
AV(INFIndex)=[];
timeAV(INFIndex)=[];


end 