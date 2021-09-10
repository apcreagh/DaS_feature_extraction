function [signal, time, no_move_index]=removeNaNs_DS(signal, time)
%Function to remove NaNs & INFs
no_move_index = find(isnan(signal));
signal(no_move_index(1:end))=[];
time(no_move_index(1:end))=[];
inf_index = find(isinf(signal));
signal(inf_index(1:end))=[];
time(inf_index(1:end))=[];
end 
%EOF