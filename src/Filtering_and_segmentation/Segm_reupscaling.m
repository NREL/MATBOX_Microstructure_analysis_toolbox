function [Labels] = Segm_reupscaling(Labels,BW_initial,BW_beforeerosiondownscaling,p)

parameters_scaling.scaling_factor = 1/p.downscaling_factor;
parameters_scaling.label_or_greylevel = 'Label';
parameters_scaling.background = 0;
% Scale
% disp 'Upscaling...';
Labels = function_scaling(uint32(Labels),parameters_scaling); % Much slower but better
% if dimension==3
%     Labels_upscaled = imresize3(uint16(Labels),p.downscaling_factor,'nearest');
% else
%     Labels_upscaled = imresize(uint16(Labels),p.downscaling_factor,'nearest');
% end
%figure; imagesc(Labels_upscaled(:,:,zdisplay)); axis equal; axis tight; colormap(tmp); pause(0.1);

sz = size(Labels);

% Remove values outside BW_beforeerosiondownscaling
Labels(BW_beforeerosiondownscaling==0)=0;

% Is there some missing values ?
tmp = zeros(sz);
tmp(Labels==0)=1;
missingpoints = double(BW_beforeerosiondownscaling).*double(tmp);
id_missingpoints = find(missingpoints);
if ~isempty(id_missingpoints)>0 % Yes, assign to nearest label lake
    [~,idx] = bwdist(~tmp);
    Labels(id_missingpoints) = Labels(idx(id_missingpoints));
end

%chess_pattern_removal = false;
% Find index for the phase
%index_phase = find(BW_beforeerosiondownscaling==1);
% Get all coordinates
%[P1,P2,P3] = ind2sub(sz,index_phase);
% % SLOW
% [Labels] = Function_correct_DPSD_identification(BW_beforeerosiondownscaling, Labels, p.cpsd_refining, p.details_convergence, chess_pattern_removal, P1, P2, P3, p.visualize_2D);

% Extrapolate
if p.extrapolate
    sz_initial = size(BW_initial);
    sz = size(Labels);
    dimension = length(sz);
    if dimension == 2
        sz = [sz 1];
        sz_initial = [sz_initial 1];
    end
    if sum(sz_initial~=sz)>0
        new_Labels = zeros(sz_initial);
        new_Labels(1:sz(1),1:sz(2),1:sz(3)) = Labels;
        cond1 = BW_initial==1;
        cond2 = new_Labels==0;
        idmissing = find(cond1.*cond2==1);
        BW = ones(sz_initial);
        BW(idmissing)=0;
        [~,IDX] = bwdist(BW);
        new_Labels(idmissing) = new_Labels(IDX(idmissing));
        Labels = new_Labels;
    end
end


end