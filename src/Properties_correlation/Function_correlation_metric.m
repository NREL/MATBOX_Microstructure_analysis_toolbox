function [CorrMatrix] = Function_correlation_metric(ValMatrix,Parnames,metric)

n_par = length(Parnames);
CorrMatrix = zeros(n_par,n_par);

for row=1:1:n_par-1
    for column=row+1:1:n_par
        par1 = ValMatrix(row,:);
        par2 = ValMatrix(column,:);
        if strcmp(metric,'Kendall')
            CorrMatrix(row,column) = Function_Kendall_tau_ranking(par1,par2);
        end
        CorrMatrix(column,row) = CorrMatrix(row,column); % Symmetric
    end
end

end