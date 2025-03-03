function [Deff, Tau, Bruggeman, eps] = Tortuosity_algorithm(M,p)
% Arguments:
% - M: segmented 3D volume or 2D image
% - p: parameter structure
%   * Background parameters
%       p.isbackground: true or false
%           if true:
%               p.background_label : integer
%               p.background_method: 'Symmetry' or 'Iterative' (recommended)
%                   if 'Iterative'
%                       p.iterate_tol     : float (e.g. 0.0005)
%                       p.iterate_maxiter : integer (e.g. 10)
%                       p.iterate_figure  : true or false
%   * Scale parameter
%       p.scale: 'One scale' or 'Dual scale'
%           if 'One scale'
%               p.onescale_label: integer
%           if 'Dual scale'
%               p.dualscale_input: 'Lower scale porosity and Bruggeman exponent' (recommended) or 'Lower scale diffusivity'
%               if 'Lower scale porosity and Bruggeman exponent'
%                   p.lowerscale_porosity   : 2D array: first row: labels, second row: porosities
%                   p.lowerscale_Bruggeman  : 2D array: first row: labels, second row: Bruggeman exponent
%               if 'Lower scale diffusivity'
%                   p.lowerscale_diffusivity : 2D array: first row: labels, second row: diffusivity
% pyrunfile (for TauFactor2_binary.py and TauFactor2_multiphase.py) in this code call TauFactor2:
% Kench et al. (2023). TauFactor 2: A GPU accelerated python tool for microstructural analysis. Journal of Open Source Software, 8(88), 5358.
% https://doi.org/10.21105/joss.05358.

% TauFactor 2 paramters
tau_verbose = true;
tau_convcritstd = p.tolerance;

%% DIRECTION
% TauFactor 2 solves along first axis
if p.direction == 2
    pd.action = 'Swap axis 1 with axis 2';
    [M,~,~] = fct_flipswap(M,pd);
elseif p.direction == 3
    pd.action = 'Swap axis 1 with axis 3';
    [M,~,~] = fct_flipswap(M,pd);
end
sz = size(M);

%% SOLVING
if p.isbackground

    unis = unique(M);
    if isscalar(unis) && unis==p.background_label
        Deff = NaN;
        Tau = NaN;
        Bruggeman = NaN;
        eps = 0;
    else
        if strcmp(p.background_method,'Symmetry')
            ps.Background = zeros(sz);
            ps.Background(M==p.background_label)=1;
            ps.perslice = false;
            [M] = fct_symmetry(M,ps);

            if strcmp(p.scale,'One scale')
                BW = zeros(sz);
                BW(M==p.onescale_label) = 1;
                eps = sum(sum(sum(BW))) / numel(BW);
                if eps>0
                    BWpython = py.numpy.array(BW);
                    res = pyrunfile("TauFactor2_binary.py", "res", img=BWpython, verbose=tau_verbose, conv_crit_std=tau_convcritstd);
                    Deff = double(res(1));
                    if Deff<=0
                        Deff = 0;
                        Tau = NaN;
                        Bruggeman = NaN;
                    else
                        Tau = double(res(2));
                        Bruggeman = 1 - log(Tau)/log(eps);
                    end
                else
                    Deff = 0;
                    Tau = NaN;
                    Bruggeman = NaN;
                end

            elseif strcmp(p.scale,'Dual scale')
                % Calculate bulk diffusivity
                if strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
                    labels_tmp = p.lowerscale_porosity(1,:);
                    dbulk_tmp = p.lowerscale_porosity(2,:).^p.lowerscale_Bruggeman(2,:);
                    eps = 0;
                    for ke = 1:length(labels_tmp)
                        eps = eps + sum(sum(sum( M==p.lowerscale_porosity(1,ke) ))) * p.lowerscale_porosity(2,ke);
                    end
                    eps = eps / numel(M);
                elseif strcmp(p.dualscale_input,'Lower scale diffusivity')
                    labels_tmp = p.lowerscale_diffusivity(1,:);
                    dbulk_tmp = p.lowerscale_diffusivity(2,:);
                    eps = NaN; % Unknown
                    Bruggeman = NaN; % Unknown
                end

                % Assign zero diffuvisity to label 0
                BW = zeros(sz);
                labels = []; dbulk = [];
                l = 0;
                %dmean = 0; % Uncomment to check tortuosity factor calculation
                for k=1:1:length(labels_tmp)
                    if dbulk_tmp(k)>0
                        l=l+1;
                        labels(l)=l;
                        dbulk(l)=dbulk_tmp(k);
                        %idx = find(M==labels_tmp(k)); % Uncomment to check tortuosity factor calculation
                        %dmean = dmean + length(idx)*dbulk(l);
                        %BW(idx)=l;
                        BW(M==labels_tmp(k))=l;
                    end
                end
                %dmean = dmean/numel(BW); % Uncomment to check tortuosity factor calculation

                BWpython = py.numpy.array(BW);
                res = pyrunfile("TauFactor2_multiphase.py", "res", img=BWpython, labels=labels, dbulk=dbulk, verbose=tau_verbose);
                Deff = double(res(1));
                if Deff<=0
                    Deff = 0;
                    Tau = NaN;
                    Bruggeman = NaN;
                else
                    Tau = double(res(2)); % Tau = dmean*1/Deff
                    if strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
                        Tau = eps/Deff;
                        Bruggeman = 1 - log(Tau)/log(eps);
                    end
                end
            end

        elseif strcmp(p.background_method,'Iterative')
            if strcmp(p.scale,'One scale')
                BW = zeros(sz);
                BW(M==p.onescale_label) = 1;
                labels = [1];
                dbulk = [1];
                eps = sum(sum(sum(BW))) / numel(BW);
                vf_background = sum(sum(sum(M==p.background_label)))/ numel(BW);
                eps = eps/(1-vf_background);

            elseif strcmp(p.scale,'Dual scale')
                % Calculate bulk diffusivity
                if strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
                    labels_tmp = p.lowerscale_porosity(1,:);
                    dbulk_tmp = p.lowerscale_porosity(2,:).^p.lowerscale_Bruggeman(2,:);
                    eps = 0;
                    for ke = 1:length(labels_tmp)
                        if p.lowerscale_porosity(1,ke)~=double(p.background_label)
                            eps = eps + sum(sum(sum( M==p.lowerscale_porosity(1,ke) ))) * p.lowerscale_porosity(2,ke);
                        end
                    end
                    eps = eps / numel(M);
                    vf_background = sum(sum(sum(M==p.background_label)))/ numel(M);
                    eps = eps/(1-vf_background);
                elseif strcmp(p.dualscale_input,'Lower scale diffusivity')
                    labels_tmp = p.lowerscale_diffusivity(1,:);
                    dbulk_tmp = p.lowerscale_diffusivity(2,:);
                    eps = NaN;
                    Bruggeman = NaN;
                end
                % Assign zero diffuvisity to label 0
                BW = zeros(sz);
                labels = []; dbulk = [];
                l = 0;
                %dmean = 0; % Uncomment to check tortuosity factor calculation
                for k=1:1:length(labels_tmp)
                    if dbulk_tmp(k)>0 && labels_tmp(k)~=double(p.background_label)
                        l=l+1;
                        labels(l)=l;
                        dbulk(l)=dbulk_tmp(k);
                        %idx = find(M==labels_tmp(k)); % Uncomment to check tortuosity factor calculation
                        %dmean = dmean + length(idx)*dbulk(l);
                        %BW(idx)=l;
                        BW(M==labels_tmp(k))=l;
                    end
                end
                %dmean = dmean/numel(BW); % Uncomment to check tortuosity factor calculation
            end

            idx_background = find(M==p.background_label);
            vals = BW;
            vals(idx_background)=[];
            dmean = 0;
            for k=1:1:length(labels)
                dmean = dmean + sum(vals==labels(k))*dbulk(k);
            end
            dmean = dmean/numel(vals);

            max_label = max(max(max(BW)));
            labels = [labels max_label+1];

            BW(idx_background) = max_label+1;
            BWpython = py.numpy.array(BW);

            dbulk = [dbulk dmean];
            diff = 1e9;
            iter = 1;
            while diff>p.iterate_tol && iter<=p.iterate_maxiter

                dbulk(end) = dmean;

                res = pyrunfile("TauFactor2_multiphase.py", "res", img=BWpython, labels=labels, dbulk=dbulk, verbose=tau_verbose);
                Deff = double(res(1));
                Tau = double(res(2)); % Tau = dmean*1/Deff

                diff = abs(dmean-Deff);
                Diter(iter) = Deff;

                dmean = Deff;
                iter = iter + 1;
            end
            Deff = Diter(end);

            Bruggeman = NaN;
            if Deff>0
                if strcmp(p.scale,'One scale')
                    Tau = eps/Deff; % Overwritte Tau
                    Bruggeman = 1 - log(Tau)/log(eps);
                else
                    if strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
                        Tau = eps/Deff; % Overwritte Tau
                        Bruggeman = 1 - log(Tau)/log(eps);
                    end
                end
            end

            if p.iterate_figure
                fig = figure;
                fig.Color = 'w';
                ax_ = axes(fig);
                hold(ax_,'on');
                plot([1:1:length(Diter)],Diter,'LineWidth',2,'LineStyle','--','Marker','o','MarkerSize',12);
                xlabel('Iterations')
                ylabel({'Effective diffusivity','(and background bulk diffusivity)'})
                set(ax_,'FontName','Times new roman','FontSize',12)
                grid on;
                fullpath = fullfile(p.folder,[p.name '_iteration']);
                savefig(fig,[fullpath '.fig']);
                saveas(fig,[fullpath '.png']);
            end
        end
    end

else % No background
    if strcmp(p.scale,'One scale')
        BW = zeros(sz);
        BW(M==p.onescale_label) = 1;
        eps = sum(sum(sum(BW))) / numel(BW);

        % if p.direction==3
        %     keyboard
        % end

        if sum(sum(sum(BW)))>0
            BWpython = py.numpy.array(BW);
            res = pyrunfile("TauFactor2_binary.py", "res", img=BWpython, verbose=tau_verbose, conv_crit_std=tau_convcritstd);
            Deff = double(res(1));
            if Deff<=0
                Deff = 0;
                Tau = NaN;
                Bruggeman = NaN;
            else
                Tau = double(res(2));
                Bruggeman = 1 - log(Tau)/log(eps);
            end
        else
            Deff = 0;
            Tau = NaN;
            Bruggeman = NaN;
        end

    elseif strcmp(p.scale,'Dual scale')
        % Calculate bulk diffusivity
        if strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
            labels_tmp = p.lowerscale_porosity(1,:);
            dbulk_tmp = p.lowerscale_porosity(2,:).^p.lowerscale_Bruggeman(2,:);
            eps = 0;
            for ke = 1:length(labels_tmp)
                eps = eps + sum(sum(sum( M==p.lowerscale_porosity(1,ke) ))) * p.lowerscale_porosity(2,ke);
            end
            eps = eps / numel(M);
        elseif strcmp(p.dualscale_input,'Lower scale diffusivity')
            labels_tmp = p.lowerscale_diffusivity(1,:);
            dbulk_tmp = p.lowerscale_diffusivity(2,:);
            eps = NaN;
            Bruggeman = NaN;
        end

        % Assign zero diffuvisity to label 0
        BW = zeros(sz);
        labels = []; dbulk = [];
        l = 0;
        %dmean = 0; % Uncomment to check tortuosity factor calculation
        for k=1:1:length(labels_tmp)
            if dbulk_tmp(k)>0
                l=l+1;
                labels(l)=l;
                dbulk(l)=dbulk_tmp(k);
                %idx = find(M==labels_tmp(k)); % Uncomment to check tortuosity factor calculation
                %dmean = dmean + length(idx)*dbulk(l);
                %BW(idx)=l;
                BW(M==labels_tmp(k))=l;
            end
        end
        %dmean = dmean/numel(BW); % Uncomment to check tortuosity factor calculation

        BWpython = py.numpy.array(BW);
        res = pyrunfile("TauFactor2_multiphase.py", "res", img=BWpython, labels=labels, dbulk=dbulk, verbose=tau_verbose);
        Deff = double(res(1));
        Tau = double(res(2)); % Tau = dmean*1/Deff
        if strcmp(p.dualscale_input,'Lower scale porosity and Bruggeman exponent')
            % Overwritte tortuosity
            Tau = eps/Deff;
            Bruggeman = 1 - log(Tau)/log(eps);
        end

    end
end

end
