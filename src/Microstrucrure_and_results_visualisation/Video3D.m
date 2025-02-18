function Video3D(M,labels,colmap,alphamap,p,seq,sav)

%% INITIALIZE VIEW
% Size and quality
hFig = uifigure;
scrsz = get(0,'ScreenSize'); % Screen resolution
set(hFig,'position',[scrsz(1) scrsz(2) scrsz(3)*1 scrsz(4)*1]); % Full screen figuredrawnow;
viewer = viewer3d('Parent',hFig);
viewer.RenderingQuality = 'High';

% Box and scale bar
viewer.Box = p.box;
viewer.ScaleBar = p.scalebar;
viewer.ScaleBarUnits = 'Voxels';
viewer.OrientationAxes = p.orientationaxes;

% Background, gradient and lighting
viewer.BackgroundColor = p.backgroundcolor;
viewer.GradientColor = p.gradientcolor;
viewer.BackgroundGradient = p.gradient;
viewer.Lighting="on";

% Volume
Vol = volshow(M,Parent=viewer);
Vol.Colormap = gray;
Vol.Alphamap = zeros(1,256);
Vol.RenderingStyle="GradientOpacity";
Vol.OverlayData=M;
Vol.OverlayRenderingStyle="LabelOverlay";
Vol.OverlayThreshold = 0.001;

% Phase color and opacity
Vol.OverlayColormap = colmap;
Vol.OverlayAlphamap = alphamap;

% Camera position
if p.rotation==0
    p.rotation = 0.01;
end
if p.elevation==0
    p.elevation = 0.01;
end
sz = size(M);
center = sz/2 + 0.5;
viewer.CameraTarget = center;
dist = sqrt(sz(1)^2 + sz(2)^2 + sz(3)^2);
vec = pi/(p.rotation); 
myPosition_initial = center + ([cos(vec) sin(vec) 1/p.elevation]*dist);
viewer.CameraPosition = myPosition_initial;

% New with 2024b
min_ = min(min(min(M)));
max_ = max(max(max(M)));

%% VIDEO
video_format = 'mpeg-4';
video_handle = VideoWriter([sav.folder sav.filename],video_format);
set(video_handle,'Quality',100); % Set video quality
set(video_handle,'FrameRate',seq.fps);
open(video_handle);

frame_number = 0;
previous_M = zeros(size(M));

for kseq = 1:1:seq.nseq
    if kseq==1
        previous_vec = pi/(p.rotation);
    else
        previous_vec = vecs(end);
    end

    seq_time = seq.time(kseq);
    seq_turn = seq.nturn(kseq);
    seq_clipdir = seq.clipdir(kseq);
    seq_clip_r0 = seq.clip_r0(kseq);
    seq_clip_r1 = seq.clip_r1(kseq);

    tmp = seq.clipphases(kseq).vals;
    if strcmp(tmp,'all')
        complementary = true;
        seq_clipphases = 1:1:max(labels);
    else       
        complementary = false;
        idtmp = strfind(tmp,'-');
        if isempty(idtmp)
            seq_clipphases = str2num(tmp);
        else
            phaseA = str2num(tmp(1:idtmp-1));
            phaseB = str2num(tmp(idtmp+1:end));
            seq_clipphases = phaseA:1:phaseB;
        end
    end
    %seq_clipphases = seq.clipphases(kseq).vals;
    
    seq_novisiblephases = seq.notvisibleohases(kseq).vals;

    % Set number of frame
    current_total = round(seq.fps*seq_time);

    % Get all rotation position
    if seq_turn~=0
        sign(seq_turn)
        if kseq==1
            vecs =  previous_vec + linspace(0, seq_turn*2*pi, current_total)';
        else
            vecs =  previous_vec + linspace(0, seq_turn*2*pi, current_total)';
            vecs(1)=[];
            current_total = current_total-1;
        end
    else
        vecs = previous_vec*ones(current_total,1);
    end
    myPosition = center + ([cos(vecs) sin(vecs) ones(size(vecs))/p.elevation]*dist);

    % Remove not visible phases
    possible_phases = M;
    if ~isempty(seq_novisiblephases)
        for k=1:1:length(seq_novisiblephases)
            possible_phases(M==seq_novisiblephases(k))=0;
        end
    end    
    possible_phases(possible_phases~=0)=1;

    % Clipping
    if ~isempty(seq_clipdir) && seq_clipdir~=0 && ~isempty(seq_clipphases) && any(seq_clipphases~=0)
        if seq_clipdir>0
            dir = ceil(seq_clipdir);
        else
            dir = abs(floor(seq_clipdir));
        end
        x0 = 1;
        x1 = current_total;

        a = (sz(dir) - 1) / (1 - 0);
        b = sz(dir)-a*1;
        y0 = round(a*seq_clip_r0+b);
        y1 = round(a*seq_clip_r1+b);

        a_slice = ( y1 -y0) / (x1 -x0);
        b_slice = y1 - a_slice*x1;
        slice_idx0 = round(a_slice*1 + b_slice);

        if kseq==1
            % Find visible and unclipped phase
            m = ismember(labels,seq_clipphases) + ismember(labels,seq_novisiblephases);
            i = find(m==0);
            notclipped_and_isvisible=zeros(sz);
            if ~isempty(i)
                for k=1:1:length(i)
                    notclipped_and_isvisible( (M.*possible_phases)==labels(i(k))) = labels(i(k));
                end
            end
        end
        % Only clipped phase
        clipped_phases = zeros(size(M));
        if complementary
            clipped_phases(M~=0)=1;
        else
            for k=1:1:length(seq_clipphases)
                clipped_phases(M==seq_clipphases(k)) = 1;
            end
        end 
        clipped_M = double(M).*double(possible_phases).*clipped_phases;
    end
    
    for current_framenumber=1:1:current_total
        frame_number = frame_number+1;

        % Clipping
        if ~isempty(seq_clipdir) && seq_clipdir~=0 && ~isempty(seq_clipphases) && any(seq_clipphases~=0)
            if kseq==1
                current_M = notclipped_and_isvisible;
            else
                current_M = previous_M;
            end

            slice_idx1 = round(a_slice*current_framenumber + b_slice);
            if seq_clipdir>0 && seq_clipdir<=1
                current_M(slice_idx0:slice_idx1,:,:) = clipped_M(slice_idx0:slice_idx1,:,:);
            elseif seq_clipdir>1 && seq_clipdir<=2
                current_M(:,slice_idx0:slice_idx1,:) = clipped_M(:,slice_idx0:slice_idx1,:);
            elseif seq_clipdir>2 && seq_clipdir<=3
                current_M(:,:,slice_idx0:slice_idx1) = clipped_M(:,:,slice_idx0:slice_idx1);
            elseif seq_clipdir<0 && seq_clipdir>=-1
                current_M(end-slice_idx1+1:end-slice_idx0+1,:,:) = clipped_M(end-slice_idx1+1:end-slice_idx0+1,:,:);      
            elseif seq_clipdir<-1 && seq_clipdir>=-2
                current_M(:,end-slice_idx1+1:end-slice_idx0+1,:) = clipped_M(:,end-slice_idx1+1:end-slice_idx0+1,:);       
            elseif seq_clipdir<-2 && seq_clipdir>=-3
                current_M(:,:,end-slice_idx1+1:end-slice_idx0+1) = clipped_M(:,:,end-slice_idx1+1:end-slice_idx0+1);                       
            end
            if kseq>1
                current_M(idx) = previous_M(idx);
            end
        else
            if kseq==1
                current_M = M.*possible_phases;
            else
                current_M = previous_M;
            end
        end

        % Ensure color will not change
        % ids = unique(current_M);
        % if length(ids)==labels
        %     Vol.OverlayColormap = colmap;
        %     Vol.OverlayAlphamap = alphamap;
        % else
        %     m = ismember(labels,ids);
        %     r = find(m==0);
        %     current_colmap = colmap;
        %     current_alphamap = alphamap;
        %     current_colmap(r,:) = [];
        %     current_alphamap(r) = [];
        %     Vol.OverlayColormap = current_colmap;
        %     Vol.OverlayAlphamap = current_alphamap;
        % end
        current_M(1,1,1)=max_;
        Vol.Data = current_M;
        Vol.OverlayData=current_M;

        Vol.OverlayColormap = colmap;
        Vol.OverlayAlphamap = alphamap;

        Vol.OverlayDisplayRange=[min_ max_];

        % Rotation
        viewer.CameraPosition = myPosition(current_framenumber,:);

        % Save frame
        if frame_number==1
            pause(p.buffertime); % Buffering time. If first frame is all blue bakcground, increase this value.
        else
            pause(0.1);
        end
        stored_frame(frame_number) = getframe(hFig);
        writeVideo(video_handle,stored_frame(frame_number))
    end
    previous_M = current_M;
    idx = find(previous_M~=0);

end

close(video_handle);

end