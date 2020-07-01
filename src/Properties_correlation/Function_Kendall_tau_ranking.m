function [Kendall_tau] = Function_Kendall_tau_ranking(x,y)

% Initialisation
number_pair = 0;
concordant_pairs = 0;
discordant_pairs = 0;
tied_pair = 0;

% Remove all NaN
idx = find(isnan(x));
idy = find(isnan(y));
id_NaN = [idx; idy];
x(id_NaN)=[];
y(id_NaN)=[];
n_ = length(x);

for i=1:1:n_-1 % First pair (xi,yi)
    for j=i+1:1:n_ % Second pair (xj,yj)
        same_sign = sign(x(i)-x(j)) * sign(y(i)-y(j));
        if same_sign ==0 % We do not want tied pair (xi=xj or/and yi=yj)
            % We do not count it
            tied_pair = tied_pair+1;
        elseif same_sign==1
            concordant_pairs=concordant_pairs+1;
            number_pair=number_pair+1;
        else
            discordant_pairs=discordant_pairs+1;
            number_pair=number_pair+1;
        end
    end
end
a=1;b=-1;c=-2*number_pair;
n_correctd = max( (-b-sqrt(b^2-4*a*c))/(2*a), (-b+sqrt(b^2-4*a*c))/(2*a));
Kendall_tau = (concordant_pairs-discordant_pairs) / (n_correctd*(n_correctd-1)/2);

end

