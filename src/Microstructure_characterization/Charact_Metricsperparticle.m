function [] = Charact_Metricsperparticle(Mins, Msem, infovol, opt, p)


%% FOLDER
Current_folder = fullfile(infovol.volpath,'Metrics_perparticle',infovol.sub);
if ~exist(Current_folder,'dir') % Folder existence is checked, and created if necessary
    mkdir(Current_folder);
end

sz = size(Mins);
dimension = length(sz);

% Semantic label
if ~isempty(Msem)
    if infovol.isbackground
        p.cracksId = 2;
        p.ParticlesId = 3;
    else
        p.cracksId = 1;
        p.ParticlesId = 2;
    end
end

% Remove truncated particles
if p.removetruncatedparticles
    if infovol.isbackground
        BW = zeros(sz);
        BW(Mins==0)=1;
        tmp = ones(sz+2);
        if dimension == 2
            tmp(2:end-1,2:end-1) = BW;
        else
            tmp(2:end-1,2:end-1,2:end-1) = BW;
        end
        dmap = bwdist(tmp,'chessboard');
        if dimension == 2
            attheborders = unique(Mins(dmap(2:end-1,2:end-1)==1));
        else
            attheborders = unique(Mins( dmap(2:end-1,2:end-1,2:end-1)==1 ));
        end
    else
        x0 = unique(Mins(1,:,:));
        x1 = unique(Mins(end,:,:));
        y0 = unique(Mins(:,1,:));
        y1 = unique(Mins(:,end,:));
        z0 = unique(Mins(:,:,1));
        z1 = unique(Mins(:,:,end));
        attheborders = unique([x0;x1;y0;y1;z0;z1]);
    end

    if ~isempty(attheborders)
        for k=1:1:length(attheborders)
            if infovol.isbackground
                Mins(Mins==attheborders(k))=1;
            else
                Mins(Mins==attheborders(k))=0;
            end
        end
    end
end

% Instances
unis_ins = unique(Mins);
if infovol.isbackground
    unis_ins(unis_ins==0)=[]; % Background
    unis_ins(unis_ins==1)=[]; % Complementary volume
else
    unis_ins(unis_ins==0)=[]; % Complementary volume
end
n_ins = length(unis_ins);

% Initialize
instance_volume = zeros(n_ins,1);
eq_diameter = zeros(n_ins,1);
crack_volume = zeros(n_ins,1);
crack_ratio = zeros(n_ins,1);

% Loop over instances
if ~isempty(Msem)
    mapcrackratio = zeros(sz);
end

for ki=1:1:n_ins
    bool = Mins == unis_ins(ki);

    nvoxels = sum(sum(sum( bool )));
    instance_volume(ki) = nvoxels * infovol.voxelsize^dimension;

    % Equivalent diameter
    if dimension == 2
        eq_diameter(ki) = ((3*instance_volume(ki))/(4*pi))^(1/3);
    else
        eq_diameter(ki) = sqrt(instance_volume(ki)/pi);
    end

    % Crack volume and ratio
    if ~isempty(Msem)
        labelsInInstance = Msem(bool);
        n_cracks = sum(sum(sum(labelsInInstance==p.cracksId)));
        n_particle = sum(sum(sum(labelsInInstance==p.ParticlesId)));
        crack_volume(ki) = n_cracks * infovol.voxelsize^dimension;
        crack_ratio(ki) = n_cracks/(n_cracks+n_particle);
        mapcrackratio(bool) = crack_ratio(ki)*1000; %0.1 percent step
    end
    
end

if ~isempty(Msem)
    mapcrackratio = round(mapcrackratio);
    [mapcrackratio] = fct_intconvert(mapcrackratio);
    function_save_tif(mapcrackratio,fullfile(Current_folder,'crackratio_map.tif'))
end


%% PLOT DISTRIBUTIONS
if strcmp(infovol.unit,'um') || strcmp(infovol.unit,'micrometer') || strcmp(infovol.unit,'Micrometer') || strcmp(infovol.unit,'micrometers') || strcmp(infovol.unit,'Micrometers')
    axisunit = '(\mum)';
    Dunit = '\mum';
else
    axisunit = ['(' infovol.unit ')'];
    Dunit = infovol.unit;
end

% Equivalent diameter
density_fct_parameters.round_value = 2;
density_fct_parameters.smooth_cumulative_fct = true;
if p.weightvolume
    w = instance_volume;
    str_sub = 'Weighted with particle volume';
else
    w = ones(size(instance_volume));
    str_sub = 'No weight';
end
[eq_diameter_distributions, ~] = Function_probability_density(eq_diameter,w,density_fct_parameters);

parameters_distributionfigure.figureposition = [100 100 1500 800];
parameters_distributionfigure.fontname = opt.format.fontname;
parameters_distributionfigure.grid = opt.format.grid;
parameters_distributionfigure.minorgrid = opt.format.minorgrid;
parameters_distributionfigure.fullpath = Current_folder;
parameters_distributionfigure.save = opt.save.savefig;
parameters_distributionfigure.subaxe1_title = 'Cumulative function';
parameters_distributionfigure.subaxe2_title = 'Distribution function';
parameters_distributionfigure.xlabel = ['Particle equivalent diameter ' axisunit];
parameters_distributionfigure.axefontsize = opt.format.axefontsize;
parameters_distributionfigure.legendfontsize = opt.format.legendfontsize;
parameters_distributionfigure.titlefontsize = opt.format.titlefontsize;
parameters_distributionfigure.sgtitlefontsize = opt.format.sgtitlefontsize;
parameters_distributionfigure.unit = Dunit;
parameters_distributionfigure.closefig = opt.format.autoclosefig;
parameters_distributionfigure.figurename =  'Particle equivalent diameter';
parameters_distributionfigure.filename = 'Equivalent diameter';
parameters_distributionfigure.title = {'Particle equivalent diameter',str_sub};
function_probability_distribution_figure(eq_diameter_distributions,parameters_distributionfigure);

if ~isempty(Msem)
    % Among cracked particles
    idnotcracked = find(crack_volume==0);
    eq_diameter_cracked = eq_diameter;
    crack_ratio_cracked = crack_ratio;
    crack_volume_cracked = crack_volume;
    w_cracked = w;
    eq_diameter_cracked(idnotcracked)=[];
    crack_ratio_cracked(idnotcracked)=[];
    crack_volume_cracked(idnotcracked)=[];
    w_cracked(idnotcracked)=[];    

    % Equivalent diameter of cracked particles
    [eq_diameter_cracked_distributions, ~] = Function_probability_density(eq_diameter_cracked,w_cracked,density_fct_parameters);
    parameters_distributionfigure.figurename =  'Equivalent diameter of cracked particles';
    parameters_distributionfigure.filename = 'Equivalent diameter only cracked instances';
    parameters_distributionfigure.title = {'Equivalent diameter of cracked particles',str_sub};
    function_probability_distribution_figure(eq_diameter_cracked_distributions,parameters_distributionfigure);

    % Crack ratio distribution
    density_fct_parameters.round_value = 2;
    [crack_ratio_distributions, ~] = Function_probability_density(crack_ratio*100,w,density_fct_parameters);
    parameters_distributionfigure.xlabel = 'Particle crack ratio (%)';
    parameters_distributionfigure.unit = '';
    parameters_distributionfigure.figurename =  'Particle crack ratio';
    parameters_distributionfigure.filename = 'Crack ratio';
    parameters_distributionfigure.title = {'Particle crack ratio',str_sub};
    function_probability_distribution_figure(crack_ratio_distributions,parameters_distributionfigure);

    [crack_ratio_cracked_distributions, ~] = Function_probability_density(crack_ratio_cracked*100,w_cracked,density_fct_parameters);
    parameters_distributionfigure.xlabel = 'Particle crack ratio (%)';
    parameters_distributionfigure.unit = '';
    parameters_distributionfigure.figurename =  'Particle crack ratio (only cracked particles)';
    parameters_distributionfigure.filename = 'Crack ratio only cracked instances';
    parameters_distributionfigure.title = {'Particle crack ratio (only cracked particles)',str_sub};
    function_probability_distribution_figure(crack_ratio_cracked_distributions,parameters_distributionfigure);    

    if p.weightvolume
        occurences = round(w_cracked./min(w_cracked));
        eq_diameter_cracked_weighted = array2warray(eq_diameter_cracked,occurences);
        crack_ratio_cracked_weighted = array2warray(crack_ratio_cracked,occurences);
        eq_diameter_weighted = array2warray(eq_diameter,occurences);
        crack_ratio_weighted = array2warray(crack_ratio,occurences);
    else
        eq_diameter_cracked_weighted = untitled(eq_diameter_cracked);
        crack_ratio_cracked_weighted = crack_ratio_cracked;
        eq_diameter_weighted = eq_diameter;
        crack_ratio_weighted = crack_ratio;
    end    

    % Bivariate distribution
    Fig = figure;
    Fig.Name= 'Bivariate distribution'; % Figure name
    Fig.Color='white'; % Background colour
    ax = axes('Parent',Fig);
    hold(ax,'on');
    h=histogram2(eq_diameter_cracked_weighted,crack_ratio_cracked_weighted*100,[50 50],'DisplayStyle','tile','ShowEmptyBins','off','Normalization','probability');
    h.EdgeColor ="none";
    hc=colorbar;
    ylabel(hc, 'Probability','FontSize',opt.format.axefontsize,'FontName',opt.format.fontname)
    xlabel(['Equivalent diameter ' axisunit]);
    ylabel('Particle crack ratio (%)');
    h_title=title({'Bivariate distribution (only cracked particles)',str_sub});
    set(ax,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
    h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
    grid(ax,'on');
    hold(ax,'off'); % Relase figure
    if opt.save.savefig % Save figure
        filename= 'Bivariate distribution diameter-crackratio only cracked instances';
        function_savefig(Fig, Current_folder, filename, opt.save); % Call function
    end
    if opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end

    Fig = figure;
    Fig.Name= 'Bivariate distribution'; % Figure name
    Fig.Color='white'; % Background colour
    ax = axes('Parent',Fig);
    hold(ax,'on');
    h=histogram2(eq_diameter_weighted,crack_ratio_weighted*100,[50 50],'DisplayStyle','tile','ShowEmptyBins','off','Normalization','probability');
    h.EdgeColor ="none";colorbar
    hc=colorbar;
    ylabel(hc, 'Probability','FontSize',opt.format.axefontsize,'FontName',opt.format.fontname)    
    xlabel(['Equivalent diameter ' axisunit]);
    ylabel('Particle crack ratio (%)');
    h_title=title({'Bivariate distribution (all particles)',str_sub});
    set(ax,'FontName',opt.format.fontname,'FontSize',opt.format.axefontsize); % Fontname and fontsize
    h_title.FontSize = opt.format.titlefontsize; % Set title fontsize
    grid(ax,'on');
    hold(ax,'off'); % Relase figure
    if opt.save.savefig % Save figure
        filename= 'Bivariate distribution diameter-crackratio all instances';
        function_savefig(Fig, Current_folder, filename, opt.save); % Call function
    end
    if opt.format.autoclosefig
        close(Fig); % Do not keep open figures
    end   
end

end

function warray = array2warray(array,occurence)
warray = [];
for k=1:1:length(occurence)
    warray = [warray; ones(occurence(k),1)*array(k)];
end
end