function [Mnew,newtype,idx_separation] = fct_separate_labels(M,p)

newtype = 'same';
sz = size(M);
dimension = length(sz);

if p.thereisabackrgound
    M(M==p.background) = p.complimentary_label;
end

if strcmp(p.choice,'Separate identically all labels')
    Mnew = M;
    id_separation = max(max(max(M)))+1;
       
    n_iter = round(p.mean_width/2);
    for k_iter = 1:1:n_iter

        if dimension == 2
            for k=1:1:sz(1)-1
                cond1 = M(k,:) ~= p.complimentary_label;
                cond2 = M(k+1,:) ~= p.complimentary_label;
                cond3 = M(k,:) ~= M(k+1,:);
                idx = cond1.*cond2.*cond3 == 1;
                Mnew(k,idx) = id_separation;
            end
            for k=sz(1):-1:2
                cond1 = M(k,:) ~= p.complimentary_label;
                cond2 = M(k-1,:) ~= p.complimentary_label;
                cond3 = M(k,:) ~= M(k-1,:);
                idx = cond1.*cond2.*cond3 == 1;
                Mnew(k,idx) = id_separation;
            end
            for k=1:1:sz(2)-1
                cond1 = M(:,k) ~= p.complimentary_label;
                cond2 = M(:,k+1) ~= p.complimentary_label;
                cond3 = M(:,k) ~= M(:,k+1);
                idx = cond1.*cond2.*cond3 == 1;
                Mnew(idx,k) = id_separation;
            end
            for k=sz(2):-1:2
                cond1 = M(:,k) ~= p.complimentary_label;
                cond2 = M(:,k-1) ~= p.complimentary_label;
                cond3 = M(:,k) ~= M(:,k-1);
                idx = cond1.*cond2.*cond3 == 1;
                Mnew(idx,k) = id_separation;
            end

        else
            for k=1:1:sz(1)-1
                sl1 = squeeze(M(k,:,:));
                sl2 = squeeze(M(k+1,:,:));
                cond1 = sl1 ~= p.complimentary_label;
                cond2 = sl2 ~= p.complimentary_label;
                cond3 = sl1~=sl2;
                idx = cond1.*cond2.*cond3 == 1;
                sl = squeeze(Mnew(k,:,:));
                sl(idx) = id_separation;
                Mnew(k,:,:) = sl;
            end
            for k=sz(1):-1:2
                sl1 = squeeze(M(k,:,:));
                sl2 = squeeze(M(k-1,:,:));
                cond1 = sl1 ~= p.complimentary_label;
                cond2 = sl2 ~= p.complimentary_label;
                cond3 = sl1~=sl2;
                idx = cond1.*cond2.*cond3 == 1;
                sl = squeeze(Mnew(k,:,:));
                sl(idx) = id_separation;
                Mnew(k,:,:) = sl;
            end
            for k=1:1:sz(2)-1
                sl1 = squeeze(M(:,k,:));
                sl2 = squeeze(M(:,k+1,:));
                cond1 = sl1 ~= p.complimentary_label;
                cond2 = sl2 ~= p.complimentary_label;
                cond3 = sl1~=sl2;
                idx = cond1.*cond2.*cond3 == 1;
                sl = squeeze(Mnew(:,k,:));
                sl(idx) = id_separation;
                Mnew(:,k,:) = sl;
            end
            for k=sz(2):-1:2
                sl1 = squeeze(M(:,k,:));
                sl2 = squeeze(M(:,k-1,:));
                cond1 = sl1 ~= p.complimentary_label;
                cond2 = sl2 ~= p.complimentary_label;
                cond3 = sl1~=sl2;
                idx = cond1.*cond2.*cond3 == 1;
                sl = squeeze(Mnew(:,k,:));
                sl(idx) = id_separation;
                Mnew(:,k,:) = sl;
            end
            for k=1:1:sz(3)-1
                cond1 = M(:,:,k) ~= p.complimentary_label;
                cond2 = M(:,:,k+1) ~= p.complimentary_label;
                cond3 = M(:,:,k) ~= M(:,:,k+1);
                idx = cond1.*cond2.*cond3 == 1;
                sl = Mnew(:,:,k);
                sl(idx) = id_separation;
                Mnew(:,:,k) = sl;
            end
            for k=sz(3):-1:2
                cond1 = M(:,:,k) ~= p.complimentary_label;
                cond2 = M(:,:,k-1) ~= p.complimentary_label;
                cond3 = M(:,:,k) ~= M(:,:,k-1);
                idx = cond1.*cond2.*cond3 == 1;
                sl = Mnew(:,:,k);
                sl(idx) = id_separation;
                Mnew(:,:,k) = sl;
            end
        end
        M = Mnew;
    end
    idx_separation = find(Mnew==id_separation);
    Mnew(idx_separation)=0;

elseif strcmp(p.choice,'Separate specifically some labels') || strcmp(p.choice,'Separate randomly labels')
    M = double(M);
    if strcmp(p.choice,'Separate randomly labels')
        labels = unique(M);
        labels(labels==p.complimentary_label)=[];
        nlabels = length(labels);
        widths = rand(nlabels,nlabels);
        widths(widths<p.ratioseparated) = 1;
        widths(widths~=1)=0;

        [x_, ~, cdf_] = distributions(p.equiprobability, p.mean_width, p.sigma_width, p.maxdeviation_width);
        id0 = find(x_<0);
        if ~isempty(id0)
            x_(id0)=[];
            cdf_(id0)=[];
        end

        idx = find(widths);
        n = length(idx);
        vals = pickfromdistriubtion(p.maxdeviation_width,p.equiprobability,p.mean_width,cdf_,x_,n);
        vals = round(vals);
        idodd = find(mod(vals,2));
        vals(idodd) = vals(idodd) + 1;
        widths(idx)=vals;
    else
        labels = p.sep_specwidth(:,1);
        widths = p.sep_specwidth(:,2:end);
    end
    Mnew = M;

    sz_rc = size(widths);
    for r=2:1:sz_rc(1)
        for c = 1:1:r-1

            if widths(r,c)>0

                labelA = labels(r);
                labelB = labels(c);
                cantor_AB = - (1/2 * (labelA+labelB) * (labelA+labelB+1) + labelB); % Bijection N*N <-> N
                n_iter = round(widths(r,c)/2);

                for k_iter = 1:1:n_iter

                    if dimension == 2
                        for k=1:1:sz(1)-1
                            cond_kA = M(k,:) == labelA;
                            cond_k2B = M(k+1,:) == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            Mnew(k,idx) = cantor_AB;

                            %cond_kA = M(k,:) == labelA;
                            cond_kC = M(k+1,:) == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            Mnew(k,idx) = cantor_AB;

                            cond_kB = M(k,:) == labelB;
                            cond_k2A = M(k+1,:) == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            Mnew(k,idx) = cantor_AB;

                            %cond_kB = M(k,:) == labelB;
                            %cond_kC = M(k+1,:) == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            Mnew(k,idx) = cantor_AB;
                        end

                        for k=sz(1):-1:2
                            cond_kA = M(k,:) == labelA;
                            cond_k2B = M(k-1,:) == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            Mnew(k,idx) = cantor_AB;

                            %cond_kA = M(k,:) == labelA;
                            cond_kC = M(k-1,:) == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            Mnew(k,idx) = cantor_AB;

                            cond_kB = M(k,:) == labelB;
                            cond_k2A = M(k-1,:) == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            Mnew(k,idx) = cantor_AB;

                            %cond_kB = M(k,:) == labelB;
                            %cond_kC = M(k-1,:) == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            Mnew(k,idx) = cantor_AB;
                        end



                        for k=1:1:sz(2)-1
                            cond_kA = M(:,k) == labelA;
                            cond_k2B = M(:,k+1) == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            Mnew(idx,k) = cantor_AB;

                            %cond_kA = M(:,k) == labelA;
                            cond_kC = M(:,k+1) == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            Mnew(idx,k) = cantor_AB;

                            cond_kB = M(:,k) == labelB;
                            cond_k2A = M(:,k+1) == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            Mnew(idx,k) = cantor_AB;

                            %cond_kB = M(:,k) == labelB;
                            %cond_kC = M(:,k+1) == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            Mnew(idx,k) = cantor_AB;
                        end

                        for k=sz(2):-1:2
                            cond_kA = M(:,k) == labelA;
                            cond_k2B = M(:,k-1) == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            Mnew(idx,k) = cantor_AB;

                            %cond_kA = M(:,k) == labelA;
                            cond_kC = M(:,k-1) == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            Mnew(idx,k) = cantor_AB;

                            cond_kB = M(:,k) == labelB;
                            cond_k2A = M(:,k-1) == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            Mnew(idx,k) = cantor_AB;

                            %cond_kB = M(:,k) == labelB;
                            %cond_kC = M(:,k-1) == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            Mnew(idx,k) = cantor_AB;
                        end

                    else
                        for k=1:1:sz(1)-1
                            sl1 = squeeze(M(k,:,:));
                            sl2 = squeeze(M(k+1,:,:));

                            cond_kA = sl1 == labelA;
                            cond_k2B = sl2 == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;

                            %cond_kA = sl1 == labelA;
                            cond_kC = sl2 == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;

                            cond_kB = sl1 == labelB;
                            cond_k2A = sl2 == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;

                            %cond_kB = sl1 == labelB;
                            %cond_kC = sl2 == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;
                        end

                        for k=sz(1):-1:2
                            sl1 = squeeze(M(k,:,:));
                            sl2 = squeeze(M(k-1,:,:));

                            cond_kA = sl1 == labelA;
                            cond_k2B = sl2 == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;

                            %cond_kA = sl1 == labelA;
                            cond_kC = sl2 == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;

                            cond_kB = sl1 == labelB;
                            cond_k2A = sl2 == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;

                            %cond_kB = sl1 == labelB;
                            %cond_kC = sl2 == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            sl = squeeze(Mnew(k,:,:));
                            sl(idx) = cantor_AB;
                            Mnew(k,:,:) = sl;
                        end

                        for k=1:1:sz(2)-1
                            sl1 = squeeze(M(:,k,:));
                            sl2 = squeeze(M(:,k+1,:));

                            cond_kA = sl1 == labelA;
                            cond_k2B = sl2 == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;

                            %cond_kA = sl1 == labelA;
                            cond_kC = sl2 == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;

                            cond_kB = sl1 == labelB;
                            cond_k2A = sl2 == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;

                            %cond_kB = sl1 == labelB;
                            %cond_kC = sl2 == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;
                        end

                        for k=sz(2):-1:2
                            sl1 = squeeze(M(:,k,:));
                            sl2 = squeeze(M(:,k-1,:));

                            cond_kA = sl1 == labelA;
                            cond_k2B = sl2 == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;

                            %cond_kA = sl1 == labelA;
                            cond_kC = sl2 == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;

                            cond_kB = sl1 == labelB;
                            cond_k2A = sl2 == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;

                            %cond_kB = sl1 == labelB;
                            %cond_kC = sl2 == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            sl = squeeze(Mnew(:,k,:));
                            sl(idx) = cantor_AB;
                            Mnew(:,k,:) = sl;
                        end

                        for k=1:1:sz(3)-1
                            cond_kA = M(:,:,k) == labelA;
                            cond_k2B = M(:,:,k+1) == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;

                            %cond_kA = M(:,:,k) == labelA;
                            cond_kC = M(:,:,k+1) == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;

                            cond_kB = M(:,:,k) == labelB;
                            cond_k2A = M(:,:,k+1) == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;

                            %cond_kB = M(:,:,k) == labelB;
                            %cond_kC = M(:,:,k+1) == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;
                        end

                        for k=sz(3):-1:2
                            cond_kA = M(:,:,k) == labelA;
                            cond_k2B = M(:,:,k-1) == labelB;
                            idx = cond_kA.*cond_k2B == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;

                            %cond_kA = M(:,:,k) == labelA;
                            cond_kC = M(:,:,k-1) == cantor_AB;
                            idx = cond_kA.*cond_kC == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;

                            cond_kB = M(:,:,k) == labelB;
                            cond_k2A = M(:,:,k-1) == labelA;
                            idx = cond_kB.*cond_k2A == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;

                            %cond_kB = M(:,:,k) == labelB;
                            %cond_kC = M(:,:,k-1) == cantor_AB;
                            idx = cond_kB.*cond_kC == 1;
                            sl = Mnew(:,:,k);
                            sl(idx) = cantor_AB;
                            Mnew(:,:,k) = sl;
                        end
                    end

                    M = Mnew;
                end
            end
        end
    end
    idx_separation = find(Mnew<0);
    Mnew(idx_separation)=0;

end

[Mnew] = fct_intconvert(Mnew);

%% FUNCTIONS
    function [x, pdf_, cdf_] = distributions(equiprobability, mean_, sigma_, maxdeviation)
        x = []; pdf_ = []; cdf_ =[];
        % Univariate Gaussian Distribution
        % Probability density function of a gaussian distribution
        pdf_GD = @(x,mu,sigma) 1./(2*pi*sigma.^2).^(0.5).*exp(-(x-mu).^2 ./ (2*sigma.^2));

        if ~equiprobability && maxdeviation~=0
            x=linspace(mean_-maxdeviation, mean_+maxdeviation, 1000);
            pdf_=pdf_GD(x,mean_,sigma_);
            cdf_ = pdf2cdf(x,pdf_);
        end

        function c = pdf2cdf(x,pdf)
            n=length(x);
            c=zeros(1,n);
            for i=2:1:n
                c(i) = trapz(x(1:i),pdf(1:i));
            end
            if c(end)<1
                c(end)=1;
            end
        end
    end

    function [val] = pickfromdistriubtion(maxdeviation,equiprobability,mean_val,cdf_,x_val,n)
        if maxdeviation==0
            val = ones(n,1)*mean_val;
        else
            if equiprobability
                val = rand(n,1)*(2*maxdeviation) + (mean_val-maxdeviation);
            else
                val = interp1(cdf_,x_val,rand(n,1));
            end
        end
    end


end