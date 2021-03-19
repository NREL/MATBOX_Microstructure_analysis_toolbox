function [mapout] = my_changem(mapout, newcode, oldcode)
    % https://www.mathworks.com/matlabcentral/answers/401395-function-changem-or-substitute-values-of-a-matrix
    assert(numel(newcode) == numel(oldcode), 'newcode and oldecode must have the same number of elements');
    [toreplace, bywhat] = ismember(mapout, oldcode);
    mapout(toreplace) = newcode(bywhat(toreplace));
end