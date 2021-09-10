function [xnew, ynew, tnew, remove_index]=removeErroneousTouchPoints(x, y, t, N)
%Remove erronous touch points 

if isempty(x) && isempty(y) && isempty(t)
    [xnew, ynew, tnew, remove_index]=deal([]);
    return; 
end 
xnew=x; ynew=y; tnew=t;
if ~exist('N', 'var') || ~isempty(N)  
    N=0.5;  end 

fs=1/nanmedian(diff(tnew)); %now deals with nans

if isinf(fs)
  fprintf('emergency patch active\n') ;
  difftnew = diff(tnew) ;
  fs = median(difftnew(find(difftnew))) ;
end

dX=sqrt(xnew.^2+ ynew.^2);
window=round(fs*N); %every N seconds
[~,TF,~,~]=hampel(dX,window,5);
remove_index=find(TF);

xnew(remove_index)=[]; ynew(remove_index)=[];  tnew(remove_index)=[]; 

end 
%EOF