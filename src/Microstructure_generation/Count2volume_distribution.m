function [particle_volumedistribution_normalized,particle_volumedistribution_fulldomain] = Count2volume_distribution(particle_diameter,particle_countdistribution,complementary_volume)
% Convert from count distribution to volume distribution assuming spherical particles
% Syntax: [volumedistribution_normalized,particle_volumedistribution_fulldomain] = Count2volume_distribution(particle_diameter,particle_countdistribution,complementary_volume)
% Inputs:
% - particle_diameter: 1D array of particle diameters of length n
% - particle_countdistribution: associated count distribution of length n
% - complementary_volume: scalar bewteen 0 and 1. The volume ratio of the particle complementary volume
% Outputs:
% - particle_volumedistribution_normalized: normalized particle volume distribution
% - particle_volumedistribution_fulldomain: normalized volume distribution (including complementary volume).

particle_volume = 4/3*pi*(particle_diameter/2).^3;
particle_countdistribution_normalized = particle_countdistribution/sum(particle_countdistribution); % Normalize

sum_particlevolume = particle_countdistribution_normalized.*particle_volume;
particle_volumedistribution_normalized = sum_particlevolume/sum(sum_particlevolume); % Normalize

particle_volumedistribution_fulldomain = particle_volumedistribution_normalized*(1-complementary_volume)./sum(particle_volumedistribution_normalized);
particle_volumedistribution_fulldomain = [particle_volumedistribution_fulldomain complementary_volume];
combined_x = cellstr([string(particle_diameter) 'Complementary volume']);

Fig = figure; % Create figure
Fig.Name= 'Particle distribution';
Fig.Color='white'; % Background colour
scrsz = get(0,'ScreenSize'); % Screen resolution
set(Fig,'position',[scrsz(1) scrsz(2) 0.9*scrsz(3) 0.8*scrsz(4)]);

t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
title(t,'Particle distribution','FontWeight','normal','Fontname','Times new roman','Fontsize',16)
t.Subtitle.String = 'Assuming spherical particles';
t.Subtitle.FontAngle = 'italic';
t.Subtitle.FontName = 'Times new roman';
t.Subtitle.FontSize = 14;

for k=1:1:4
    nexttile   
    if k==1
        b=bar(particle_diameter,particle_countdistribution,1.0);
        xlabel('Particle diameter')
        ylabel('Particle count (input)')
    elseif k==2
        b=bar(particle_diameter,particle_countdistribution_normalized,1.0);
        xlabel('Particle diameter')
        ylabel('Particle count (normalized)')
    elseif k==3
        b=bar(particle_diameter,particle_volumedistribution_normalized,1.0);
        xlabel('Particle diameter')
        ylabel('Particle volume distribution (normalized)')
    elseif k==4
        X = categorical(combined_x);
        X = reordercats(X,combined_x);
        b=bar(X,particle_volumedistribution_fulldomain,1.0);
        b.FaceColor = 'flat';
        b.CData(end,:) = [0.5 0 0.5];
        xlabel('Particle diameter and complementary volume')
        ylabel('Volume distribution (normalized)')
    end
    set(gca,'Fontname','Times new roman','Fontsize',12)
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(b(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
end

end



