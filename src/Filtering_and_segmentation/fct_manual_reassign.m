function [newM,newtype,foo] = fct_manual_reassign(oldM,p)
foo=[];
newtype = 'same';

n_old_label = length(p.oldlabels);
newM = oldM;
for k=1:1:n_old_label
    if p.oldlabels(k)~=p.newlabels(k)
        newM(oldM==p.oldlabels(k)) = p.newlabels(k);
    end
end

% Convert to uint8 or uint16
[newM] = fct_intconvert(newM);

end