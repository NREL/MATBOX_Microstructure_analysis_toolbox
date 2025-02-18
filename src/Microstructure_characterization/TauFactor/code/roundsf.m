function [y]=roundsf(number,sfs,method)
sfs = 5;
%opt = {'round','floor','ceil','fix'};
og = 10.^(floor(log10(abs(number)) - sfs + 1));
y = feval(method,number./og).*og;
y(find(number==0)) = 0;