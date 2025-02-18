close all
clearvars
clc

%% IMPORT
tmp = function_load_tif('C:\Users\fussegli\Desktop\T\nmc-1-cal-greyscale_step_12.tif');
M = uint8(tmp);
imagesc(M(:,:,50));

%% ALGORITHM
% Interface has a fractal dimension of 2
% Loop over threshold A - threshold B and find interval that is closer to 2

p.boxmaxlength = 3;
p.topology_dimension = 3;
p.plot = false;
n=15; % Number of search interval
tol = 0.1; % Tolerance on dimension

f=zeros(n-1,1);
threshold_bounds = round(linspace(0,255,n));
thresholds = zeros(n-1,1);
for k=1:1:n-1
    k
    thresholdA = threshold_bounds(k)+1;
    thresholdB = threshold_bounds(k+1);
    thresholds(k) = round((thresholdA+thresholdB)/2);
    BW = zeros(size(M));
    cond1 = double(M>=thresholdA);
    cond2 = double(M<=thresholdB);

    BW( (cond1+cond2) == 2 )=1;

    [~,~,fractal_dimension,~] = Function_fractaldimension_boxcounting(BW,p);
    f(k)=fractal_dimension(2,1);
end
fd=abs(f-2);
[f fd fd<=tol]

Fig = figure;
Fig.Color='white';
ax_=axes('Parent',Fig);
hold(ax_,'on');
h_title=title('Interface detection through fractal dimension'); % Set title font
plot(thresholds,f,'LineWidth',2,'Marker','o','MarkerSize',12);
xlabel('(Threshold B - Threshold A)/2');
ylabel('Fractal dimension');
grid(ax_,'on'); % Display grid
set(ax_,'FontName','Times new roman','FontSize',12); % Fontname and fontsize
hold(ax_,'off');

id=find(fd==min(fd))
interface = zeros(size(M));
cond1 = double(M>=threshold_bounds(id)+1);
cond2 = double(M<=threshold_bounds(id+1));
interface( (cond1+cond2) == 2 )=125;

Seg=zeros(size(M));
Seg(M>thresholds(id))=125;

Microstructure_comparison_visualization_interface(M,interface,Seg)




