function [M] = fct_combination_merge(Mlow,Mhigh,p)

sz = size(Mlow);
dimension = length(sz);

unique_lowlevel = unique(Mlow);
unique_lowlevel(unique_lowlevel==0)=[];
n_lowlevel = length(unique_lowlevel);

newlowlevel_labels = zeros(sz);
newlowlevel_labels_shared = zeros(sz);

iter = 0;
for k=1:1:n_lowlevel
    % Old label
    oldlabel = unique_lowlevel(k);

    idx = find(Mlow==oldlabel);
    vals = Mhigh(idx);

    [counts, groupnames] = groupcounts(vals);

    tot = sum(counts);
    idmax = find(counts==max(counts));
    if max(counts) >= p.threshold*tot
        newlowlevel_labels(idx)=groupnames(idmax);
    else
        iter = iter+1;
        newlowlevel_labels_shared(idx)=iter;
    end
end

unique_lowlevelshared = unique(newlowlevel_labels_shared);
unique_lowlevelshared(unique_lowlevelshared==0)=[];
if ~isempty(unique_lowlevelshared)
    m = max(max(max(newlowlevel_labels)));
    n_lowlevelshared = length(unique_lowlevelshared);
    for k=1:1:n_lowlevelshared
        idx = find(newlowlevel_labels_shared==unique_lowlevelshared(k));
        m=m+1;
        newlowlevel_labels(idx)=m;
    end
end

M = newlowlevel_labels;


if p.check_circularity % BETA
    minsize_circularity = 0.05; % [0,1] Circularity is poorly calculated for small particles
    maxdifference_hull=0.15; % [0,1] Circularity is calculated on convex hull only if size of the convex hull is not that different from the size of the particle


    % PERMUTATION/CIRCULARITY analysis:
    % Low level labels are merged only if it increase circularity
    % All combinations are evaluated. Low level combinations of decreasing circularity are merged until there is no more labels remaining
    unique_combined = unique(M);
    unique_combined(unique_combined==0)=[];
    n_combined= length(unique_combined);
    Mb = M;
    hasbeenmodified = false;
    for k=1:1:n_combined
        % k/n_combined
        idx = find(M == unique_combined(k));
        [IX,IY,IZ]=ind2sub(sz,idx);
        sub_tmp = Mlow(min(IX):max(IX),min(IY):max(IY),min(IZ):max(IZ));

        sz_tmp = size(sub_tmp);
        sub_low = zeros(sz_tmp+2);
        if dimension==2
            sub_low(2:end-1,2:end-1)=sub_tmp; clear sub_tmp;
        else
            sub_low(2:end-1,2:end-1,2:end-1)=sub_tmp; clear sub_tmp;
        end
        sz_tmp = size(sub_low);
        nsub_low = sum(sum(sum( sub_low~=0 )));

        lowlabels = Mlow(idx);
        unique_lowlabels = unique(lowlabels);
        n_lowlabels = length(unique_lowlabels);

        if n_lowlabels>1
            % All cases
            % https://www.mathworks.com/matlabcentral/answers/525295-how-do-i-create-a-matrix-with-all-binary-combinations
            cases = [dec2bin(0:2^n_lowlabels-1)' - '0']';
            cases(1,:) = []; % Remove 0  0 0 ... case;
            [n,~]=size(cases);
            circularity = zeros(n,1);
            BWsize = zeros(n,1);
            for kc = 1:1:n
                idl = find(cases(kc,:));
                ls = unique_lowlabels(idl);
                BW = zeros(sz_tmp);
                for kl=1:1:length(ls)
                    BW(sub_low==ls(kl))=1;
                end
                nBW = sum(sum(sum(BW)));

                if nBW >= minsize_circularity * nsub_low % Circularity is poorly calculated for small regions
                    % Calculate sphericity (dim=3) or circularity (dim=2)
                    if dimension == 2
                        % Get convex hull (remove impact of surface roughness, focus only on shape)
                        hull = bwconvhull(BW);
                        nhull = sum(sum(sum(hull)));
                        nBW = sum(sum(sum(BW)));
                        if abs(nhull-nBW)/nBW > maxdifference_hull % Convex hull is too different to be relevant
                            hull = BW;
                        end

                        % Area numeric
                        An = sum(sum(hull));
                        % Perimeter numeric
                        diffrow = abs(diff(hull,1,1));
                        diffcol = abs(diff(hull,1,2));
                        Pn = sum(sum(diffrow)) + sum(sum(diffcol));
                        % Perimeter numeric, corrected
                        Pnc = Pn*(pi/4);
                        % Specific perimeter numeric
                        Spn = Pnc/An;
                        % Equivalent radius based on area
                        re = sqrt(An/pi);
                        % Specific perimeter equivalent disc
                        Spe = 2/re;
                        % Circularity
                        circularity(kc,1) = Spe/Spn;

                        % if circularity(kc,1)>=0.959 && n_lowlabels==12
                        %     figure; imagesc(BW); axis equal;
                        %     figure; imagesc(hull); axis equal; title([num2str(kc,'%i') ', c=' num2str(circularity(kc,1),'%1.3f')]); pause(.5);
                        % end

                    else
                        warning('not yet implemented');
                    end
                    BWsize(kc,1)= nBW;
                else
                    circularity(kc,1)=0; % Circularity is poorly calculated for small regions
                    BWsize(kc,1)=0;
                end
                %figure; imagesc(hull); axis equal; title([num2str(kc), ' c=' num2str(circularity(kc,1))]);
            end
            % Sort cases per circularity. In case of circularity equality, sort per per size
            % if n_lowlabels==12
            %    keyboard
            % end


            sav_cases = cases;
            sav_circularity = circularity;

            circularity = round(circularity,4);
            cases = [circularity BWsize cases];
            cases = sortrows(cases,1,"descend");
            uni_circularity = unique(circularity);
            uni_circularity = sort(uni_circularity,'descend');
            n_circularity = length(uni_circularity);
            new_cases = [];
            if n_circularity<length(circularity)
                for kcirc = 1:1:n_circularity
                    ids = find(cases(:,1)==uni_circularity(kcirc));
                    tmp_case = sortrows( cases(ids,:),2,"descend");
                    new_cases = [new_cases; tmp_case(1,:)];
                end
                cases = new_cases;
            end
            cases(:,1:2)=[];

            [n,~]=size(cases);
            if sum(cases(1,:))<n_lowlabels % We should not have merged all labels as done before!
                hasbeenmodified = true;
                % Suboptimal circularity, let's find best combinations of cases
                available_labels = unique_lowlabels;
                new_combined_labels = [];
                for kc=1:1:n
                    if isempty(available_labels)
                        break
                    end
                    idl = find(cases(kc,:));
                    tomerged = unique_lowlabels(idl);

                    % t=zeros(size(sub_low));
                    % for kk=1:1:length(tomerged)
                    %     t(sub_low==tomerged(kk))=1;
                    % end
                    % figure; imagesc(t); axis equal; title(num2str(circularity(kc,1)),'%1.3f')

                    tmp = ismember(available_labels,tomerged);
                    if sum(tmp)==length(tomerged)
                        m = max(max(max( Mb)));
                        for kl=1:1:length(tomerged)
                            Mb( Mlow==tomerged(kl) ) = m+1;
                        end
                        % Remove from available labels
                        available_labels(tmp)=[];
                    end
                end
            end
        end
    end

    if hasbeenmodified
        M = Mb;
    end

end

[M] = fct_intconvert(M);

end