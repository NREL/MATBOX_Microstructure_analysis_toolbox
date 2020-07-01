%%%%%%%%%%%%%%%%%%%%%%   TauFactor

% TauFactor is a MatLab application for quantifying the calculate factor
% from microstructure data stored in tiff files. It is suitable for 2D or
% 3D geometry data for materials containing up to 3 distinct phases. For
% support with this application please contact Sam Cooper at:
% camsooper@gmail.com

% Copyright (c) 2016, Samuel J Cooper
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% TauFactor makes use of freezeColors.m available at:
% https://www.mathworks.com/matlabcentral/fileexchange/7943-freezecolors---unfreezecolors
% See the code for its own BSD license details.


function varargout = TauFactor(varargin)
%% TAUFACTOR MATLAB code for TauFactor.fig
%      TAUFACTOR, by itself, creates a new TAUFACTOR or raises the existing
%      singleton*.
%
%      H = TAUFACTOR returns the handle to a new TAUFACTOR or the handle to
%      the existing singleton*.
%
%      TAUFACTOR('CALLBACK',hObject,eventData,hand,...) calls the local
%      function named CALLBACK in TAUFACTOR.M with the given input arguments.
%
%      TAUFACTOR('Property','Value',...) creates a new TAUFACTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TauFactor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TauFactor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TauFactor

% Last Modified by GUIDE v2.5 19-Jan-2017 15:29:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct(...
    'gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TauFactor_OpeningFcn, ...
    'gui_OutputFcn',  @TauFactor_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT
function TauFactor_OpeningFcn(hObject, eventdata, hand, varargin)
%% Initialise various state variable
hand.InLineMode=0;
hand.mexCtrl='mex';
try AA=[1:2]';BB=uint32(AA);Mex3DTauIso_mex(2,3,[AA;AA],[AA;AA],1,AA,AA,...
        BB,BB,BB,BB,BB,BB,BB,BB,BB,BB,BB,BB);
catch 'MEX-file not compatible with the platform';
    hand.mexCtrl='mat';
    set(hand.Check_mex,'Enable','off');
    set(hand.Check_mex,'Value',0);
end

hand.WoSpace='base';
hand.impCheck=0;
hand.Blocked=0;
hand.output = hObject;
hand.RandomSecrets=1;
hand.Check_FreqPlots=0;
hand.filename='Dummy';
hand.w=0;
hand.Aniso=0;
hand.ConvErr=0;
% Open window and
% set(hand.Window,'Position',[250,37,42,39]);
% Create tau logo
[hand.logo_A, hand.logo_map, hand.logo_alpha] = imread('TauFactor_icon2.png');
h = imshow(hand.logo_A, hand.logo_map);
set(h, 'AlphaData', hand.logo_alpha);
% save an anitialisation fielding the global handle
hand.reset=hand;
drawnow
guidata(hObject, hand);

% --- Outputs from this function are returned to the command line.
function varargout = TauFactor_OutputFcn(hObject, eventdata, hand)
varargout{1} = hand.output;

function [OutputResults]=InLine(FindTau, FindMetrics, RVAmode, Net_Or, PhaDir, VoxDims)
%% Assign input variables to the global handle
% Inline mode allows the use to operate the functionsality of TauFactor
% from within other scripts or directly from the command line
hand.InLineMode=1;
hand.mexCtrl='mex';
hand.Check_mex=uicontrol('Parent',[],'Value',1);
try AA=[1:2]';BB=uint32(AA);Mex3DTauIso_mex(2,3,[AA;AA],[AA;AA],1,AA,AA,...
        BB,BB,BB,BB,BB,BB,BB,BB,BB,BB,BB,BB);
catch
    %     disp('MEX-file not compatible with the platform');
    hand.mexCtrl='mat';
    hand.Check_mex=uicontrol('Parent',[],'Value',0);
end
hand.WoSpace='caller';
hand.PhaDir=PhaDir;
hand.FindMetrics=FindMetrics;

hand.L1box=uicontrol('Parent',[],'String',VoxDims(1));
hand.L2box=uicontrol('Parent',[],'String',VoxDims(2));
hand.L3box=uicontrol('Parent',[],'String',VoxDims(3));

hand.Edit_D_Black=uicontrol('Parent',[],'String',1);
hand.Edit_D_Green=uicontrol('Parent',[],'String',0);
hand.Edit_D_White=uicontrol('Parent',[],'String',0);

hand.Check_Impedance=uicontrol('Parent',[],'Value',0);
hand.Check_VaryD=uicontrol('Parent',[],'Value',0);
% hand.Check_pdfSave=uicontrol('Parent',[],'Value',SaveResults);
if RVAmode==0
    hand.Check_RVA=uicontrol('Parent',[],'Value',0);
    hand.Pop_RV=uicontrol('Parent',[],'Value',1);
else
    hand.Check_RVA=uicontrol('Parent',[],'Value',1);
    hand.Pop_RV=uicontrol('Parent',[],'Value',RVAmode);
end

hObject=figure('Visible','off');
if ~exist('Name')
    Name='OutputResults';
end
hand.filename=[Name,'.zzz'];
hand.impCheck=0;
hand.Blocked=0;
hand.Check_FreqPlots=0;
hand.MBytesA=0;
hand.MemLoc={''};
hand.Aniso=0;
hand.NetVals=unique(Net_Or);
if length(hand.NetVals)==1
    if hand.NetVals==0
        hand.Net_Or=false(size(Net_Or));
    elseif hand.NetVals==2
        hand.Net_Or=unit8(2*ones(size(Net_Or),'uint8'));
    else
        hand.Net_Or=ones(size(Net_Or),'logical');
    end
elseif length(hand.NetVals)==2
    hand.Net_Or=Net_Or==hand.NetVals(2);
elseif length(hand.NetVals)==3
    Net=uint8(Net_Or==hand.NetVals(2));
    Net(Net_Or==hand.NetVals(3))=2;
    hand.Net_Or=Net;
else
    hand.Net_Or=Net_Or<0;
    hand.Dmap=Net_Or;
    set(hand.Check_VaryD,'value',1)
    set(hand.Edit_D_Black,'String',1000)
    set(hand.Edit_D_Green,'String','Dmap')
    set(hand.Edit_D_White,'String',0)
end
Net=1;
Net_Or=1;
hand.Check_Reverse=uicontrol('Parent',[],'Visible','off');
hand.Check_Blocked=uicontrol('Parent',[],'Visible','off');
hand.Check_B1=uicontrol('Parent',[],'Value',PhaDir(1,1));
hand.Check_B2=uicontrol('Parent',[],'Value',PhaDir(1,2));
hand.Check_B3=uicontrol('Parent',[],'Value',PhaDir(1,3));
hand.Check_W1=uicontrol('Parent',[],'Value',PhaDir(3,1));
hand.Check_W2=uicontrol('Parent',[],'Value',PhaDir(3,2));
hand.Check_W3=uicontrol('Parent',[],'Value',PhaDir(3,3));
if length(hand.NetVals)>2
    hand.Check_G1=uicontrol('Parent',[],'Value',PhaDir(2,1));
    hand.Check_G2=uicontrol('Parent',[],'Value',PhaDir(2,2));
    hand.Check_G3=uicontrol('Parent',[],'Value',PhaDir(2,3));
else
    hand.Check_G1=uicontrol('Parent',[],'Value',0);
    hand.Check_G2=uicontrol('Parent',[],'Value',0);
    hand.Check_G3=uicontrol('Parent',[],'Value',0);
end
hand.VFs(1)=mean(mean(mean(hand.Net_Or==0)));
hand.VFs(2)=mean(mean(mean(hand.Net_Or==1)));
hand.VFs(3)=mean(mean(mean(hand.Net_Or==2)));
[hand]=ExpectedTime(hObject, 1, hand);
if FindTau==1
    [hand]=Tortuosity(hObject, 1, hand);
end
if FindMetrics==1
    [hand]=Metrics(hObject, 1, hand);
end
% If all the simulations were on non percoalting volumes, no results with
% be generated, so a dummy variable is required
if ~isfield(hand,'OutputResults')
    hand.OutputResults.Reply=0;
    OutputResults=hand.OutputResults;
else
    OutputResults=hand.OutputResults;
end

function [hand]=LoadData_Callback(hObject, eventdata, hand)
%% Load data from tif stacks
set(hand.LoadData,...
    'String','Loading...',...
    'ForegroundColor',[.9 .9 .9]);
set(hand.LoadData,'Enable','off')
set(hand.Radio_BatchStack,'String','Batch Loaded')
set(hand.Radio_BatchStack,'Value',0)
drawnow
% if isfield(hand,'anno')
%     delete(hand.hs)
%     delete(hand.anno)
%     delete(hand.im)
%     hand=hand.reset;
%     hand.reset=hand;
% end
[filename,pathname]=uigetfile('*.*','Select .tif file.','MultiSelect','on');
if iscell(filename)==1
    set(hand.Radio_BatchStack,'Visible','On');
    hand.batch_filenames=filename;
    hand.batch_length=length(filename);
    hand.filename=char(hand.batch_filenames(1));
    for i=1:hand.batch_length
        if length(filename{i})>3
            if strcmp(lower(filename{i}(end-2:end)),'tif')~=1 && strcmp(lower(filename{i}(end-2:end)),'png')~=1 && strcmp(lower(filename{i}(end-3:end)),'tiff')~=1 && strcmp(lower(filename{i}(end-2:end)),'mat')~=1
                disp('File type not suitable, load *.tif')
                set(hand.LoadData,...
                    'String','Load Data',...
                    'ForegroundColor',[1 1 1]);
                set(hand.LoadData,'Enable','On')
                return
            end
        end
    end
else
    set(hand.Radio_BatchStack,'Value',0);
    set(hand.Radio_BatchStack,'Visible','Off');
    if length(filename)>3
        if strcmp(lower(filename(end-2:end)),'tif')~=1 && strcmp(lower(filename(end-2:end)),'png')~=1 && strcmp(lower(filename(end-3:end)),'tiff')~=1 && strcmp(lower(filename(end-2:end)),'mat')~=1
            disp('File type not suitable, load *.tif')
            set(hand.LoadData,...
                'String','Load Data',...
                'ForegroundColor',[1 1 1]);
            set(hand.LoadData,'Enable','On')
            return
        end
    end
    hand.batch_length=1;
    hand.filename=filename;
end
hand.pathname=pathname;
[hand]=Radio_BatchStack_Callback(hObject, eventdata, hand);
set(hand.LoadData,...
    'String','Load Data',...
    'ForegroundColor',[1 1 1]);
set(hand.LoadData,'Enable','On')
guidata(hObject, hand);

function [hand]=Radio_BatchStack_Callback(hObject, eventdata, hand)
if get(hand.Radio_BatchStack,'Value')==0
    set(hand.Radio_BatchStack,'String','Batch Loaded')
else
    set(hand.Radio_BatchStack,'String','Stack Loaded')
end
if isfield(hand,'anno')
    delete(hand.hs)
    delete(hand.anno)
    if isfield(hand,'VFanno')
        delete(hand.VFanno)
    end
    delete(hand.im)
end
[hand]=PrepData_Callback(hObject, eventdata, hand);
[hand]=PlotData_Callback(hObject, eventdata, hand);
[hand]=PrepWindow_Callback(hObject, eventdata, hand);
guidata(hObject, hand);

function [hand]=PrepData_Callback(hObject, eventdata, hand)
%% Convert the input geometry data into the required format, which labels the phases with 0, 1 and 2.
if hand.pathname==0 % If no file is selected, load a test dataset
    hand.pathname  = [cd,'\'];
    hand.filename  = 'Random70.tif';
    [Net]=ExampleNet_70;
    [a,b,c]=size(Net);
    hand.NetVals=[0 1];
else
    % Build Network Map
    if strcmp(hand.filename(end-2:end),'mat')
        Net=struct2cell(load([hand.pathname,hand.filename]));
        Net=Net{1};
        [a,b,c]=size(Net);
    else
        imageInfo = imfinfo([hand.pathname,hand.filename]);
        a=imageInfo.Height;
        b=imageInfo.Width;
        c=numel(imageInfo);
        Net = zeros(a,b,c);      % Preallocate the cell array
        if c>1;
            for i = 1:c
                Net(:,:,i) = imread([hand.pathname,hand.filename],'Index',i,'Info',imageInfo);
            end
        else
            Net=imread([hand.pathname,hand.filename]);
        end
    end
    if get(hand.Radio_BatchStack,'Value')==1
        c=hand.batch_length;
        Net = zeros(a,b,c);
        for i = 1:c
            if strcmp(hand.batch_filenames{i}(end-2:end),'mat')
                zz=struct2cell(load([hand.pathname,hand.batch_filenames{i}]));
                zz=zz{1};
            else
                imageInfo = imfinfo([hand.pathname,hand.batch_filenames{i}]);
                zz =(imread([hand.pathname,hand.batch_filenames{i}],'Index',1,'Info',imageInfo));
            end
            if size(zz(:,:,1))~=[a b]
                disp('2D images must be of equal dimension to build stack!');
                set(hand.Radio_BatchStack,'String','Batch Loaded')
                set(hand.Radio_BatchStack,'Value',0)
                return
            end
            Net(:,:,i)=zz(:,:,1);
            zz=1;
        end
    end
    hand.NetVals=unique(Net);
    if length(hand.NetVals)==1
        if hand.NetVals==0
            Net=zeros(size(Net),'logical');
        elseif hand.NetVals==2
            Net=unit8(2*ones(size(Net),'uint8'));
        else
            Net=ones(size(Net),'logical');
        end
    elseif length(hand.NetVals)==2
        Net=Net==hand.NetVals(2);
    elseif length(hand.NetVals)==3
        temp=uint8(Net==hand.NetVals(2));
        temp(Net==hand.NetVals(3))=2;
        Net=temp;
        temp=1;
    else
        if length(hand.NetVals)>3
            disp('Volume has more than 3 phases; interpreting as diffusivity map.');
        end
    end
end
if c==1 && length(Net(1,1,:))>1
    Net=Net(:,:,1);
end
if length(hand.NetVals)<4
    hand.Net_Or=uint8(Net);
    hand.VFs(1)=mean(mean(mean(hand.Net_Or==0)));
    hand.VFs(2)=mean(mean(mean(hand.Net_Or==1)));
    hand.VFs(3)=mean(mean(mean(hand.Net_Or==2)));
else
    if strcmp(class(Net),'double') || strcmp(class(Net),'single')
        hand.Net_Or=single(Net);
    else
        hand.Net_Or=uint8(Net);
    end
    hand.VFs(1)=mean(hand.Net_Or(:));
end
[hand]=ExpectedTime(hObject, eventdata, hand);

if isfield(hand,'MBytesA')
    hand.MBytesA=0;
    hand.MemLoc={''};
end
MBytes=whos;hand.MBytesA=[sum([MBytes.bytes])/1024^2]; %Memory
hand.MemLoc={'start'};
drawnow
guidata(hObject, hand);

function [hand]=PlotData_Callback(hObject, eventdata, hand)
%% Plots the initaial sliced representation of the data into the window after loading
if isfield(hand,'anno')
    set(hand.reset.Window,'Position',get(hand.Window,'Position'));
end
if length(hand.NetVals)==3
    hand.myColorMap=zeros(64,3);%[0 0 0;1 1 1;0 1 0];
    hand.myColorMap(32:33,2)=1;
    hand.myColorMap(64,:)=[1 1 1];
    Cmapping=32;
elseif length(hand.NetVals)<3
    hand.myColorMap=gray;
    Cmapping=64;
else
    hand.myColorMap=gray;
    Cmapping=0.5;
end
set(hand.Text_FileName,'Visible','On')
set(hand.Text_FileName,'String',hand.filename);
[a b c]=size(hand.Net_Or);
if length(hand.Net_Or(:))<1e9
    spacer=' x ';
else
    spacer='x';
end
set(hand.Text_VoxelDims,'Visible','On',...
    'String',[num2str(a),spacer,num2str(b),spacer,num2str(c)]);
if c>1
    if a<b
        wh1=[0.5 0.215*a/b];
    else
        wh1=[0.5*b/a 0.215];
    end
    if c<20
        LayersNo=c;
    else
        LayersNo=17;
    end
    if max(size(hand.Net_Or))<300
        ss=1;
    else
        ss=2;
    end
    p.step=ceil(c/LayersNo);
    p.layers=floor(c/p.step);
    for i=p.layers:-1:1
        hand.hs(i)=subplot(1,p.layers,i);
        set(hand.hs(i),'visible','off')
        hand.im(i)=image(Cmapping*hand.Net_Or(1:ss:end,1:ss:end,1+(i-1)*p.step));
        set(gca,'xtick',[],'ytick',[]);
        %         p.pos(i,:)=[0.17+0.007*i 0.47+0.007*i 0.5 0.5];
        p.pos(i,:)=[0.17+0.007*i 0.615+0.007*i wh1];
        set(hand.hs(i),'position',[0 0 0 0]);
    end
    %     set(hand.axes1,'visible','on')
    for i=p.layers:-1:1
        set(hand.hs(i),'position',p.pos(i,:));
        
        if length(hand.NetVals)==2
            colormap(gray);
        else
            colormap(hand.myColorMap);
        end
        pause(0.03);
    end
    hand.anno(1,1)= annotation(gcf,'doublearrow',[0.12 0.12],[0.63 0.8]);
    hand.anno(1,2)= annotation(gcf,'doublearrow',[0.19 0.63],[0.6 0.6]);
    hand.anno(1,3)= annotation(gcf,'doublearrow',[0.75 0.88],[0.62 0.73]);
    hand.anno(2,1)= annotation(gcf,'textbox',[0. 0.65 0.1 0.1],'String','1');
    hand.anno(2,2)= annotation(gcf,'textbox',[0.34 0.51 0.1 0.1],'String','2');
    hand.anno(2,3)= annotation(gcf,'textbox',[0.84 0.6 0.1 0.1],'String','3');
    set(hand.anno(2,:),'LineStyle','none','FontUnits','Normalized','FontSize',0.05);
    set(hand.anno,'units','characters')
    set(hand.Check_W3,'Enable','on');
    set(hand.Check_B3,'Enable','on');
    set(hand.L3box,'Enable','on');
else
    hand.hs=subplot(1,1,1);
    hand.im=imagesc(hand.Net_Or(:,:,1));
    if a<b
        wh1=[28 10.5*a/b];
    else
        wh1=[28*b/a 10.5];
    end
    set(hand.hs,'units','characters');
    set(hand.hs,'Position',[0.5*(42-wh1(1)) 25 wh1]);
    hand.anno(1,1)=annotation(gcf,'doublearrow',[0.12 0.12],[0.65 0.9]);
    hand.anno(1,2)=annotation(gcf,'doublearrow',[0.21 0.79],[0.62 0.62]);
    hand.anno(2,1)= annotation(gcf,'textbox',[0. 0.72 0.1 0.1],'String','1');
    hand.anno(2,2)= annotation(gcf,'textbox',[0.43 0.53 0.1 0.1],'String','2');
    set(hand.anno(2,1:2),'LineStyle','none','FontUnits','Normalized','FontSize',0.05);
    set(hand.anno(:,1:2),'units','characters');
    set(gca,'XColor',[0 0 1]);
    set(gca,'YColor',[0 0 1]);
    set(hand.Check_W3,'Enable','off');
    set(hand.Check_B3,'Enable','off');
    set(hand.Check_G3,'Enable','off');
    if get(hand.Check_Impedance, 'Value')==0;
        set(hand.L3box,'Enable','on');
    end
    set(hand.Check_W3,'Value',0);
    set(hand.Check_B3,'Value',0);
    set(hand.Check_G3,'Value',0);
end
%

if length(hand.NetVals)==2
    colormap(gray);
else
    colormap(hand.myColorMap);
end
set(gca,'xtick',[],'ytick',[]);
if length(hand.NetVals)>1 && length(hand.NetVals)<4
    hand.VFanno=1;
    hand.VFanno(1)= annotation(gcf,'textbox',[0.82 0.84 0.2 0.1],'String',[num2str(100*hand.VFs(1),2),'%']);
    hand.VFanno(2)= annotation(gcf,'textbox',[0.82 0.81 0.2 0.1],'String',[num2str(100*hand.VFs(2),2),'%']);
    set(hand.VFanno(2),'Color',[0.5 0.5 0.5]);
    if length(hand.NetVals)==3
        hand.VFanno(3)= annotation(gcf,'textbox',[0.82 0.78 0.2 0.1],'String',[num2str(100*hand.VFs(3),2),'%']);
        set(hand.VFanno(2),'Color',[0 1 0]);
        set(hand.VFanno(3),'Color',[0.5 0.5 0.5]);
    end
    set(hand.VFanno,'LineStyle','none');
else
    hand.VFanno=1;
    hand.VFanno(1)= annotation(gcf,'textbox',[1 1 0 0],'String','%');
end
guidata(hObject, hand);

function [hand]=PrepWindow_Callback(hObject, eventdata, hand)
%% Unhide all the various calcuation options in the window
[a,b,c]=size(hand.Net_Or);
set(hand.TextB_Green,'String','Green','ForegroundColor',[0 0 0])
set(hand.Calculate,'Enable','on');
set(hand.MetricsButton,'Enable','on');
set(hand.TimeBox,'Visible','on');
set(hand.RAMbox,'Visible','on');
set(hand.Check_mex,'Visible','on');
set(hand.Check_pdfSave,'Visible','on');
set(hand.Check_FluxMapSave,'Visible','on');
set(hand.Check_RVA,'Visible','on')
set(hand.Pop_RV,'Visible','on')
set(hand.Check_VaryD,'Visible','on')
set(hand.TextB_Direction,'Visible','on');
set(hand.TextB_Aniso,'Visible','on');
set(hand.TextB_White,'Visible','on');
set(hand.TextB_Black,'Visible','on');
set(hand.TextB_Green,'Visible','on');
set(hand.Check_B1,'Visible','on');
set(hand.Check_B2,'Visible','on');
set(hand.Check_B3,'Visible','on');
if get(hand.Check_VaryD,'value')==0
    set(hand.Check_W1,'Visible','on');
    set(hand.Check_W2,'Visible','on');
    set(hand.Check_W3,'Visible','on');
    set(hand.Check_G1,'Visible','on');
    set(hand.Check_G2,'Visible','on');
    set(hand.Check_G3,'Visible','on');
end
set(hand.Check_Impedance,'Visible','on')
if length(hand.NetVals)==3
    set(hand.Check_G1,'Enable','on');
    set(hand.Check_G2,'Enable','on');
    set(hand.Check_G3,'Enable','off');
    if c>1
        set(hand.Check_G3,'Enable','on');
    end
    if get(hand.Check_VaryD,'value')==1
        set(hand.Edit_D_Black,'visible','on')
        if length(hand.NetVals)>2
            set(hand.Edit_D_Green,'visible','on')
        end
        set(hand.Edit_D_White,'visible','on')
    end
    set(hand.Check_VaryD,'Enable','on')
    %     set(hand.Edit_D_Black,'String','1')
    set(hand.Edit_D_Black,'Enable','on')
    %     set(hand.Edit_D_Green,'String','0')
    set(hand.Edit_D_Green,'Enable','on')
    %     set(hand.Edit_D_White,'String','0')
    set(hand.Edit_D_White,'Enable','on')
elseif length(hand.NetVals)<3
    set(hand.Check_G1,'Enable','off');
    set(hand.Check_G2,'Enable','off');
    set(hand.Check_G3,'Enable','off');
    set(hand.Check_G1,'Value',0);
    set(hand.Check_G2,'Value',0);
    set(hand.Check_G3,'Value',0);
    set(hand.Edit_D_Green,'visible','off')
    set(hand.Check_VaryD,'Enable','on')
    %     set(hand.Edit_D_Black,'String','1')
    set(hand.Edit_D_Black,'Enable','on')
    %     set(hand.Edit_D_Green,'String','1')
    set(hand.Edit_D_Green,'Enable','on')
    %     set(hand.Edit_D_White,'String','1')
    set(hand.Edit_D_White,'Enable','on')
else
    set(hand.TextB_Green,'String','Mode:','ForegroundColor',[0.5 0.5 0.5])
    set(hand.Check_VaryD,'Value',1)
    set(hand.Check_VaryD,'Enable','off')
    Check_VaryD_Callback(hObject, eventdata, hand);
    set(hand.Edit_D_Black,'Enable','on')
    set(hand.Edit_D_Green,'String','Dmap')
    set(hand.Edit_D_Green,'Enable','on')
    set(hand.Edit_D_White,'Enable','on')
    set(hand.Check_Impedance,'Enable','off')
    if strcmp(class(hand.Net_Or),'double')
        set(hand.Edit_D_Black,'String',num2str(min(hand.NetVals)))
        set(hand.Edit_D_White,'String',num2str(max(hand.NetVals)))
    elseif strcmp(class(hand.Net_Or),'uint8')
        set(hand.Edit_D_Black,'String',num2str(min(hand.NetVals)/255))
        set(hand.Edit_D_White,'String',num2str(max(hand.NetVals)/255))
    end
end
set(hand.L1box,'Visible','on');
set(hand.L2box,'Visible','on');
set(hand.L3box,'Visible','on');
set(hand.TextB_Dir1,'Visible','on');
set(hand.TextB_Dir2,'Visible','on');
set(hand.TextB_Dir3,'Visible','on');
set(hand.MetricsButton,'Visible','on');
set(hand.Calculate,'Visible','on');
set(hand.Check_W1,'Enable','on');
set(hand.Check_W2,'Enable','on');
set(hand.Check_B1,'Enable','on');
set(hand.Check_B2,'Enable','on');
% set(hand.TextB_Aniso,'Enable','on');
% if license('test', 'Distrib_Computing_Toolbox') && 1==2
%     set(hand.Check_Impedance,'visible','on')
% else
%     set(hand.Check_Impedance,'visible','off')
% end
guidata(hObject, hand);

%% Whenever one the the checkboxes is selected, the time/memory approximations are updated
function Check_W1_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_W2_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_W3_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_B1_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_B2_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_B3_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_G1_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_G2_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);
function Check_G3_Callback(hObject, eventdata, hand)
ExpectedTime(hObject, eventdata, hand);

function TextB_Aniso_Callback(hObject, eventdata, hand)
%% Ghost of the Aniso button

%% Update expected time and check for anisotropy when a size box is changed
function L1box_Callback(hObject, eventdata, hand)
[hand]=VoxelSizeFun(hObject, eventdata, hand);
function L1box_CreateFcn(hObject, eventdata, hand)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function L2box_Callback(hObject, eventdata, hand)
[hand]=VoxelSizeFun(hObject, eventdata, hand);
function L2box_CreateFcn(hObject, eventdata, hand)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function L3box_Callback(hObject, eventdata, hand)
[hand]=VoxelSizeFun(hObject, eventdata, hand);
function L3box_CreateFcn(hObject, eventdata, hand)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [hand] = VoxelSizeFun(hObject, eventdata, hand)
%% Voxel dimension comparison function
[hand]=ExpectedTime(hObject, eventdata, hand);
if str2double(get(hand.L1box,'String'))~=str2double(get(hand.L2box,'String')) ||...
        str2double(get(hand.L2box,'String'))~=str2double(get(hand.L3box,'String'))
    hand.Aniso=1;
else
    hand.Aniso=0;
end

function [hand]=ExpectedTime(hObject, eventdata, hand)
%% Calculate the expected end time of the simulation
if get(hand.Check_Impedance, 'Value')==1;
    hand.impCheck=1;
    set(hand.Check_Reverse,'Visible','on')
    set(hand.Check_Blocked,'Visible','on')
    set(hand.L1box,'String',1);set(hand.L1box,'Enable','off');
    set(hand.L2box,'String',1);set(hand.L2box,'Enable','off');
    set(hand.L3box,'String',1);set(hand.L3box,'Enable','off');
else
    hand.impCheck=0;
    set(hand.L1box,'Enable','on');
    set(hand.L2box,'Enable','on');
    set(hand.L3box,'Enable','on');
    set(hand.Check_Reverse,'Visible','off')
    set(hand.Check_Blocked,'Visible','off')
    set(hand.Check_Blocked, 'Value',0)
    set(hand.Check_Reverse, 'Value',0)
end

[hand]=Check_Blocked_Callback(hObject, eventdata, hand);
SimsNo=...
    get(hand.Check_W1,'Value')+...
    get(hand.Check_W2,'Value')+...
    get(hand.Check_W3,'Value')+...
    get(hand.Check_G1,'Value')+...
    get(hand.Check_G2,'Value')+...
    get(hand.Check_G3,'Value')+...
    get(hand.Check_B1,'Value')+...
    get(hand.Check_B2,'Value')+...
    get(hand.Check_B3,'Value');
% if SimsNo==0
%     set(hand.Check_W1,'Value',strcmp(get(hand.Check_W1,'Enable'),'on'));
%     set(hand.Check_W2,'Value',strcmp(get(hand.Check_W2,'Enable'),'on'));
%     set(hand.Check_W3,'Value',strcmp(get(hand.Check_W3,'Enable'),'on'));
%     set(hand.Check_G1,'Value',strcmp(get(hand.Check_G1,'Enable'),'on'));
%     set(hand.Check_G2,'Value',strcmp(get(hand.Check_G2,'Enable'),'on'));
%     set(hand.Check_G3,'Value',strcmp(get(hand.Check_G3,'Enable'),'on'));
%     set(hand.Check_B1,'Value',strcmp(get(hand.Check_B1,'Enable'),'on'));
%     set(hand.Check_B2,'Value',strcmp(get(hand.Check_B2,'Enable'),'on'));
%     set(hand.Check_B3,'Value',strcmp(get(hand.Check_B3,'Enable'),'on'));
%     SimsNo=...
%         get(hand.Check_W1,'Value')+...
%         get(hand.Check_W2,'Value')+...
%         get(hand.Check_W3,'Value')+...
%         get(hand.Check_G1,'Value')+...
%         get(hand.Check_G2,'Value')+...
%         get(hand.Check_G3,'Value')+...
%         get(hand.Check_B1,'Value')+...
%         get(hand.Check_B2,'Value')+...
%         get(hand.Check_B3,'Value');
% end
if ~isfield(hand,'Net_Or')
    return
end
waitcoeff=1.5e-08;
if str2double(get(hand.L1box,'String'))~=str2double(get(hand.L2box,'String')) ||...
        str2double(get(hand.L2box,'String'))~=str2double(get(hand.L3box,'String'))
    waitcoeff=1.3*waitcoeff;
    hand.iter_max=max(size(hand.Net_Or))*40;
else
    hand.iter_max=max(size(hand.Net_Or))*30;
end

if get(hand.Check_Impedance, 'Value')==1
    hand.iter_max=hand.iter_max*14;
end
hand.time_approx=roundsf(SimsNo*numel(hand.Net_Or)*hand.iter_max/4*waitcoeff/86400,2,'round');
if get(hand.Check_Impedance,'Value')==1;
    hand.time_approx=hand.time_approx*8;
end

if get(hand.Check_VaryD,'Value')==1;
    hand.time_approx=hand.time_approx*1.1;
end

if get(hand.Check_RVA, 'Value')==1;
    hand.time_approx=hand.time_approx*5;
end

if get(hand.Check_mex,'Value')==1;
    hand.time_approx=hand.time_approx/3;
end

if isfield(hand,'batch_length')&& get(hand.Radio_BatchStack,'Value')==0
    hand.time_approx=hand.time_approx*hand.batch_length;
end
if hand.InLineMode==0
    [tocStr]=TimeString(hand.time_approx);
    set(hand.TimeBox,'String',tocStr);
    [hand]=ExpectedRAM(hObject, eventdata, hand);
end
guidata(hObject, hand);

function [hand]=ExpectedRAM(hObject, eventdata, hand)
%% Calculate the approximate required memory
memcoeff=35/1024^2;
if str2double(get(hand.L1box,'String'))~=str2double(get(hand.L2box,'String')) ||...
        str2double(get(hand.L2box,'String'))~=str2double(get(hand.L3box,'String'))
    memcoeff=1.3*memcoeff;
end
if get(hand.Check_Impedance,'Value')==1
    memcoeff=1.4*memcoeff;
end
if get(hand.Check_VaryD,'Value')==1;
    memcoeff=5*memcoeff;
end
[a b c]=size(hand.Net_Or);
hand.RAM_approx=memcoeff*(a+2)*(b+2)*(c+2);
if hand.RAM_approx<1024
    set(hand.RAMbox,'String',['RAM = ',num2str(roundsf(hand.RAM_approx,2,'round')),' MB']);
else
    set(hand.RAMbox,'String',['RAM = ',num2str(roundsf(hand.RAM_approx/1024,2,'round')),' GB']);
end
if ispc
    [~,sys] = memory;
    hand.RAM_availableMBs=sys.PhysicalMemory.Available/1024^2;
    if hand.RAM_availableMBs<hand.RAM_approx
        set(hand.RAMbox,'ForegroundColor',[1 0 0])
    else
        set(hand.RAMbox,'ForegroundColor',[0 0 1])
    end
end
guidata(hObject, hand);

function MetricsButton_Callback(hObject, eventdata, hand)
%% Call the Metrics function and indicated to the user to wait
set(hand.MetricsButton,...
    'ForegroundColor',[.2 .2 .2],...
    'String','wait...');
drawnow
if hand.batch_length~=1 && get(hand.Radio_BatchStack,'Value')==0
    for i=1:hand.batch_length
        hand.filename=char(hand.batch_filenames(i));
        [hand]=PrepData_Callback(hObject, eventdata, hand);
        [hand]=Metrics(hObject, eventdata, hand);
    end
else
    [hand]=Metrics(hObject, eventdata, hand);
end
set(hand.MetricsButton,...
    'ForegroundColor',[0 0 0],...
    'String','Metrics');
drawnow
guidata(hObject, hand);

function [hand]=Metrics(hObject, eventdata, hand)
if length(hand.filename)>53
    hand.fil=[hand.filename(1:50)];
elseif length(hand.filename)>4
    hand.fil=[hand.filename(1:end-4)];
else
    hand.fil=['OutputResults'];
end
hand.fil(~isstrprop(hand.fil,'alphanum'))='_';
if ~isstrprop(hand.fil(1), 'alpha')
    hand.fil=['X',hand.fil];
end
%% Call the Metrics RVA function (multiple times if required)
VoxDim(1)=str2double(get(hand.L1box,'String'));
VoxDim(2)=str2double(get(hand.L2box,'String'));
VoxDim(3)=str2double(get(hand.L3box,'String'));

% hand.Net_Or=flipud(rot90(rot90(rot90(hand.Net_Or))));
hand.Dir=0;
hand.RVArepeats=0;
if get(hand.Check_RVA, 'Value')==1 && get(hand.Pop_RV, 'Value')>1;
    hand.RVArepeats=...
        ((get(hand.Check_W1,'Value')+get(hand.Check_B1,'Value')+get(hand.Check_G1,'Value'))>=1)+...
        ((get(hand.Check_W2,'Value')+get(hand.Check_B2,'Value')+get(hand.Check_G2,'Value'))>=1)+...
        ((get(hand.Check_W3,'Value')+get(hand.Check_B3,'Value')+get(hand.Check_G3,'Value'))>=1);
    if (get(hand.Check_W1,'Value')+get(hand.Check_B1,'Value')+get(hand.Check_G1,'Value'))>=1
        hand.Dir=1;
        hand.RVArepeats=hand.RVArepeats-1;
        [hand]=MetricsRVA(VoxDim(1),VoxDim(2),VoxDim(3),get(hand.Check_RVA, 'Value'),hand);
        drawnow
    end
    if (get(hand.Check_W2,'Value')+get(hand.Check_B2,'Value')+get(hand.Check_G2,'Value'))>=1
        hand.Dir=2;
        hand.RVArepeats=hand.RVArepeats-1;
        [hand]=MetricsRVA(VoxDim(2),VoxDim(3),VoxDim(1),get(hand.Check_RVA, 'Value'),hand);
        drawnow
    end
    if (get(hand.Check_W3,'Value')+get(hand.Check_B3,'Value')+get(hand.Check_G3,'Value'))>=1
        hand.Dir=3;
        hand.RVArepeats=hand.RVArepeats-1;
        [hand]=MetricsRVA(VoxDim(3),VoxDim(1),VoxDim(2),get(hand.Check_RVA, 'Value'),hand);
    end
else
    [hand]=MetricsRVA(VoxDim(1),VoxDim(2),VoxDim(3),get(hand.Check_RVA, 'Value'),hand);
end
guidata(hObject, hand);

function Calculate_Callback(hObject, eventdata, hand)
%% Is the Tortuosity button is pressed, begin the initiation sequence
set(hand.Calculate,...
    'ForegroundColor',[.9 .9 .9],...
    'String','Calculating...');
[hand]=ExpectedTime(hObject, eventdata, hand);
set(hand.LoadData,'enable','off');
set(hand.MetricsButton,'enable','off');

if hand.time_approx<1
    set(hand.TimeBox,'String',['End = ',datestr(hand.time_approx+now, 'HH:MM PM')]);
else
    set(hand.TimeBox,'String',['End = ',datestr(hand.time_approx+now, 'ddd (HHPM)')]);
end
drawnow
if hand.batch_length~=1 && get(hand.Radio_BatchStack,'Value')==0
    %     fh=findall(0,'type','figure')';
    for i=1:hand.batch_length
        hand.filename=char(hand.batch_filenames(i));
        [hand]=PrepData_Callback(hObject, eventdata, hand);
        [hand]=Tortuosity(hObject, eventdata, hand);
    end
else
    [hand]=Tortuosity(hObject, eventdata, hand);
end

set(hand.Calculate,...
    'ForegroundColor',[1 1 1],...
    'String','Tortuosity');
set(hand.LoadData,'enable','on');
set(hand.MetricsButton,'enable','on');
drawnow
hand.MBytesA=0;
hand.MemLoc={''};
[num2cell(round(hand.MBytesA))', hand.MemLoc'];
guidata(hObject, hand);

function [hand]=Tortuosity(hObject, eventdata, hand)
xxx=tic;
hand.TauSet=zeros(3,3);
if length(hand.filename)>53
    hand.fil=[hand.filename(1:50)];
elseif length(hand.filename)>4
    hand.fil=[hand.filename(1:end-4)];
else
    hand.fil=['OutputResults'];
end
hand.fil(~isstrprop(hand.fil,'alphanum'))='_';
if ~isstrprop(hand.fil(1), 'alpha')
    hand.fil=['X',hand.fil];
end
%% Step through each phase/direction check box to run the tortuosity calculations
if str2double(get(hand.L1box,'String'))~=str2double(get(hand.L2box,'String')) ||...
        str2double(get(hand.L2box,'String'))~=str2double(get(hand.L3box,'String'))
    hand.Error_max=roundsf(1e-4/(max(size(hand.Net_Or)))^1.2,3,'round');
else
    hand.Error_max=roundsf(1e-4/(max(size(hand.Net_Or)))^1.1,3,'round');
end
if hand.iter_max<10
    hand.iter_max=10;
end
if get(hand.Check_Impedance,'Value')==0
    hand.check_f=round(hand.iter_max/200);
else
    hand.check_f=round(hand.iter_max/800);
end
hand.JacobiRat=200;
if hand.check_f<10
    hand.check_f=10;
end
if (get(hand.Check_B1,'Value')+get(hand.Check_B2,'Value')+get(hand.Check_B3,'Value')+...
        get(hand.Check_G1,'Value')+get(hand.Check_G2,'Value')+get(hand.Check_G3,'Value')+...
        get(hand.Check_W1,'Value')+get(hand.Check_W2,'Value')+get(hand.Check_W3,'Value'))==0
    return
end
if hand.InLineMode==0
    if isnan(str2double(get(hand.L1box,'String'))*str2double(get(hand.L2box,'String'))*str2double(get(hand.L3box,'String')))||...
            (str2double(get(hand.L1box,'String'))*str2double(get(hand.L2box,'String'))*str2double(get(hand.L3box,'String')))==0
        set(hand.L1box,'BackgroundColor',[1 0 0]);
        set(hand.L2box,'BackgroundColor',[1 0 0]);
        set(hand.L3box,'BackgroundColor',[1 0 0]);
        return
    else
        set(hand.L1box,'BackgroundColor',[1 1 1]);
        set(hand.L2box,'BackgroundColor',[1 1 1]);
        set(hand.L3box,'BackgroundColor',[1 1 1]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a,b,c]=size(hand.Net_Or);
    disp(['%%%%%%%%%%%% RESULTS %%%%%%%%%%%% - ',hand.filename,char(10),...
        'Voxel size = ',get(hand.L1box,'String'),' x ', get(hand.L2box,'String'),' x ', get(hand.L3box,'String'),' nm',char(10)...
        'Volume size = ',num2str(a),' x ', num2str(b),' x ', num2str(c),' voxels'...
        ]);
end
if get(hand.Check_VaryD,'value')==0
    if get(hand.Check_B1,'Value')==1
        hand.Net=hand.Net_Or==0;
        hand.Dir=1;
        hand.Pha='Black';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(1,1)=hand.Results.Tau;
    end
    if get(hand.Check_B2,'Value')==1
        [~,~,c]=size(hand.Net_Or);
        hand.Net=hand.Net_Or==0;
        if c>1
            hand.Net=logical(permute(hand.Net,[2 3 1]));
        else
            hand.Net=logical(rot90(hand.Net,3));
        end
        hand.Dir=2;
        hand.Pha='Black';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(1,2)=hand.Results.Tau;
    end
    if get(hand.Check_B3,'Value')==1
        hand.Net=hand.Net_Or==0;
        hand.Net=logical(flip(permute(hand.Net,[3 1 2]),3));
        hand.Dir=3;
        hand.Pha='Black';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(1,3)=hand.Results.Tau;
    end
    
    if get(hand.Check_G1,'Value')==1
        hand.Net=hand.Net_Or==1;
        hand.Dir=1;
        hand.Pha='Green';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(2,1)=hand.Results.Tau;
    end
    if get(hand.Check_G2,'Value')==1
        [~,~,c]=size(hand.Net_Or);
        hand.Net=hand.Net_Or==1;
        if c>1
            hand.Net=permute(hand.Net,[2 3 1]);
        else
            hand.Net=logical(rot90(hand.Net,3));
        end
        hand.Dir=2;
        hand.Pha='Green';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(2,2)=hand.Results.Tau;
    end
    if get(hand.Check_G3,'Value')==1
        hand.Net=hand.Net_Or==1;
        hand.Net=flip(permute(hand.Net,[3 1 2]),3);
        hand.Dir=3;
        hand.Pha='Green';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(2,3)=hand.Results.Tau;
    end
    
    if get(hand.Check_W1,'Value')==1
        if max(hand.Net_Or(:))>0
            hand.Net=hand.Net_Or==max(hand.Net_Or(:));
        else
            hand.Net=zeros(size(hand.Net_Or));
        end
        hand.Dir=1;
        hand.Pha='White';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(3,1)=hand.Results.Tau;
    end
    if get(hand.Check_W2,'Value')==1
        [~,~,c]=size(hand.Net_Or);
        if max(hand.Net_Or(:))>0
            hand.Net=hand.Net_Or==max(hand.Net_Or(:));
        else
            hand.Net=zeros(size(hand.Net_Or));
        end
        if c>1
            hand.Net=logical(permute(hand.Net,[2 3 1]));
        else
            hand.Net=logical(rot90(hand.Net,3));
        end
        hand.Dir=2;
        hand.Pha='White';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(3,2)=hand.Results.Tau;
    end
    if get(hand.Check_W3,'Value')==1
        if max(hand.Net_Or(:))>0
            hand.Net=hand.Net_Or==max(hand.Net_Or(:));
        else
            hand.Net=zeros(size(hand.Net_Or));
        end
        hand.Net=logical(flip(permute(hand.Net,[3 1 2]),3));
        hand.Dir=3;
        hand.Pha='White';
        [hand]=RVA_Tau(hand);
        [hand]=PrintTauResults(hObject, eventdata, hand);
        hand.TauSet(3,3)=hand.Results.Tau;
    end
    if hand.InLineMode==0
        if get(hand.Check_pdfSave,'Value')==1
            Save_XL(hand);
            if get(hand.Check_B1,'Value')==1 && get(hand.Check_B2,'Value')==1 && get(hand.Check_B3,'Value')==1 &&...
                    get(hand.Check_G1,'Value')==1 && get(hand.Check_G2,'Value')==1 && get(hand.Check_G3,'Value')==1 &&...
                    get(hand.Check_W1,'Value')==1 && get(hand.Check_W2,'Value')==1 && get(hand.Check_W3,'Value')==1
                if get(hand.Check_RVA, 'Value')==0;
                    set(hand.Check_RVA','Value',1)
                    [hand]=Metrics(hObject, eventdata, hand);
                    set(hand.Check_RVA','Value',0)
                else
                    [hand]=Metrics(hObject, eventdata, hand);
                end
            end
        end
    end
else % Variable D
    hand.Pha='Multi';
    hand.D_Strings{1}=(get(hand.Edit_D_Black,'String'));
    hand.D_Strings{2}=(get(hand.Edit_D_Green,'String'));
    hand.D_Strings{3}=(get(hand.Edit_D_White,'String'));
    if length(hand.NetVals)<4
        if length(hand.NetVals)<3 % If two phases
            hand.Ds=[str2double(hand.D_Strings(1)), str2double(hand.D_Strings(3))];
        else
            hand.Ds=[str2double(hand.D_Strings(1)), str2double(hand.D_Strings(2)), str2double(hand.D_Strings(3))];
        end
        if ~isnan(sum(hand.Ds)) && max(hand.Ds)==0 % All Ds zero
            disp('At least one phase needs to have a non-zero D value.')
            return
        elseif length(unique(hand.Ds))==1
            disp('If all phases have the same D value, the system becomes analytical.')
            return
        else
            disp([   'D_black = ',get(hand.Edit_D_Black,'String'),...
                '  || D_green = ', get(hand.Edit_D_Green,'String'),...
                '  || D_white = ', get(hand.Edit_D_White,'String')]);
            hand.Dmap=hand.Ds(1)*double(hand.Net_Or==0);
            hand.Dmap(hand.Net_Or==1)=hand.Ds(2);
            if length(hand.NetVals)==3
                hand.Dmap(hand.Net_Or==2)=hand.Ds(3);
            end
            if sum(isnan(hand.Ds))>0
                for i=1:length(hand.NetVals)
                    if isnan(hand.Ds(i))
                        Dist=1e-9*mean([str2double(get(hand.L1box,'String')),...
                            str2double(get(hand.L2box,'String')),...
                            str2double(get(hand.L3box,'String'))]);
                        if i==2 && length(hand.NetVals)==2
                            Str=hand.D_Strings{3};
                        else
                            Str=hand.D_Strings{i};
                        end
                        idx=find(Str=='v');
                        D=str2double(Str(1:idx-1));
                        v=str2double(Str(idx+1:end));
                        hand.Dmap(hand.Net_Or==(i-1))=0;
                        DM=bwdist(hand.Net_Or~=(i-1));
                        DM(DM>0)=1./(1/D+1./((DM(DM>0)-0.5)*Dist*v*sqrt(32/(3*pi))));
                        hand.Dmap=hand.Dmap+DM;
                    end
                end
            end
            hand.Dmap=padarray(hand.Dmap,[1,1,1],0);
            hand.D=mean(mean(mean(hand.Dmap(2:end-1,2:end-1,2:end-1))));
            hand.Dmap=1./hand.Dmap;
            if min(hand.Ds)>0
                hand.Net=ones(size(hand.Net_Or),'uint8');
                if get(hand.Check_B1,'Value')==1
                    hand.Dir=1;
                    [hand]=RVA_Tau(hand);
                end
                if get(hand.Check_B2,'Value')==1
                    hand.Dir=2;
                    if c>1
                        hand.Net=logical(permute(ones(size(hand.Net_Or),'uint8'),[2 3 1]));
                        hand.Dmap=(permute(hand.Dmap,[2 3 1]));
                    else
                        hand.Net=logical(rot90(ones(size(hand.Net_Or),'uint8'),3));
                        hand.Dmap=(rot90(hand.Dmap,3));
                    end
                    [hand]=RVA_Tau(hand);
                    if get(hand.Check_B3,'Value')==1
                        if c>1 %% undo Dmap adjustment
                            hand.Dmap=(permute(hand.Dmap,[3 1 2]));
                        else
                            hand.Dmap=rot90(hand.Dmap,1);
                        end
                    end
                end
                if get(hand.Check_B3,'Value')==1
                    hand.Dir=3;
                    hand.Net=flip(permute(true(size(hand.Net_Or)),[3 1 2]),3);
                    hand.Dmap=(flip(permute(hand.Dmap,[3 1 2]),3));
                    [hand]=RVA_Tau(hand);
                end
            else
                NZD=find(hand.Ds);
                hand.Net=zeros(size(hand.Net_Or),'uint8');
                for i=1:length(NZD)
                    hand.Net=hand.Net+uint8(hand.Net_Or==(NZD(i)-1));
                end
                if get(hand.Check_B1,'Value')==1
                    hand.Dir=1;
                    [hand]=RVA_Tau(hand);
                    [hand]=PrintTauResults(hObject, eventdata, hand);
                end
                if get(hand.Check_B2,'Value')==1
                    hand.Net=zeros(size(hand.Net_Or),'uint8');
                    for i=1:length(NZD)
                        hand.Net=hand.Net+uint8(hand.Net_Or==(NZD(i)-1));
                    end
                    [~,~,c]=size(hand.Net_Or);
                    if c>1
                        hand.Net=logical(permute(hand.Net,[2 3 1]));
                        hand.Dmap=(permute(hand.Dmap,[2 3 1]));
                    else
                        hand.Net=logical(rot90(hand.Net,3));
                        hand.Dmap=(rot90(hand.Dmap,3));
                    end
                    hand.Dir=2;
                    [hand]=RVA_Tau(hand);
                    [hand]=PrintTauResults(hObject, eventdata, hand);
                    if get(hand.Check_B3,'Value')==1
                        if c>1
                            hand.Dmap=(permute(hand.Dmap,[3 1 2]));
                        else
                            hand.Dmap=(rot90(hand.Dmap,1));
                        end
                    end
                end
                if get(hand.Check_B3,'Value')==1
                    hand.Net=zeros(size(hand.Net_Or),'uint8');
                    for i=1:length(NZD)
                        hand.Net=hand.Net+uint8(hand.Net_Or==(NZD(i)-1));
                    end
                    hand.Net=logical(flip(permute(hand.Net,[3 1 2]),3));
                    hand.Dmap=(flip(permute(hand.Dmap,[3 1 2]),3));
                    hand.Dir=3;
                    [hand]=RVA_Tau(hand);
                    [hand]=PrintTauResults(hObject, eventdata, hand);
                end
            end
        end
    else
        if hand.InLineMode==0
            if strcmp(hand.D_Strings{2},'Dmap')
                hand.Ds=[str2double(hand.D_Strings{1}),str2double(hand.D_Strings{3})];
                if strcmp(class(hand.Net_Or),'double') || strcmp(class(hand.Net_Or),'single')
                    hand.Dmap=hand.Net_Or;
                elseif strcmp(class(hand.Net_Or),'uint8')
                    hand.Dmap=hand.Ds(1)*ones(size(hand.Net_Or))+(hand.Ds(2)-hand.Ds(1))*double(hand.Net_Or)/255;
                end
            elseif strcmp(hand.D_Strings{2},'Dseg')
                ThreshFig=figure;
                set(ThreshFig,'Color',[1 1 1]);
                imagesc((hand.Net_Or(:,:,1)));
                title('Draw rectangle in CONDUCTIVE phase');
                rect_pB = uint32(round(getrect(ThreshFig)));
                rect_pB=mean(mean(hand.Net_Or(max(1,rect_pB(2)):min(b,sum(rect_pB([2,4]))),...
                    max(1,rect_pB(1)):min(a,sum(rect_pB([1,3]))),1)));
                title('Draw rectangle in INSULATING phase');
                rect_pW = uint32(round(getrect(ThreshFig)));
                rect_pW=mean(mean(hand.Net_Or(max(1,rect_pW(2)):min(b,sum(rect_pW([2,4]))),...
                    max(1,rect_pW(1)):min(a,sum(rect_pW([1,3]))),1)));
                Thresh=mean([rect_pB,rect_pW]);
                Alpha=(2/(rect_pB-rect_pW))*log(999);%99=(1/0.01)-1
                hand.Dmap=sigmf(hand.Net_Or,[Alpha Thresh]);
                imagesc((hand.Dmap(2:end-1,2:end-1,2)));
                pause(1);close(ThreshFig)
            end
        else
        end
        hand.Net=hand.Dmap>0;
        hand.D=mean(hand.Dmap(:));
        hand.Dmap=double(padarray(hand.Dmap,[1,1,1],0));
        hand.Dmap=1./hand.Dmap;
        if get(hand.Check_B1,'Value')==1
            hand.Dir=1;
            [hand]=RVA_Tau(hand);
            [hand]=PrintTauResults(hObject, eventdata, hand);
        end
        if get(hand.Check_B2,'Value')==1
            [~,~,c]=size(hand.Net_Or);
            if c>1
                hand.Net=logical(permute(hand.Dmap(2:end-1,2:end-1,2:end-1),[2 3 1]));
                hand.Dmap=(permute(hand.Dmap,[2 3 1]));
            else
                hand.Net=logical(rot90(hand.Dmap(2:end-1,2:end-1,2:end-1),3));
                hand.Dmap=(rot90(hand.Dmap,3));
            end
            hand.Dir=2;
            [hand]=RVA_Tau(hand);
            [hand]=PrintTauResults(hObject, eventdata, hand);
            if get(hand.Check_B3,'Value')==1
                if c>1
                    hand.Dmap=(permute(hand.Dmap,[3 1 2]));
                else
                    hand.Dmap=(rot90(hand.Dmap,1));
                end
            end
        end
        if get(hand.Check_B3,'Value')==1
            hand.Net=logical(flip(permute(hand.Dmap(2:end-1,2:end-1,2:end-1),[3 1 2]),3));
            hand.Dmap=(flip(permute(hand.Dmap,[3 1 2]),3));
            hand.Dir=3;
            [hand]=RVA_Tau(hand);
            [hand]=PrintTauResults(hObject, eventdata, hand);
        end
    end
end
if hand.InLineMode==0
    [a,b,c]=size(hand.Net_Or);
    [timeStr]=TimeString(toc(xxx)/86400);
    if isfield(hand.Results,'MBytes')
        if max(round(hand.Results.MBytes))<1024
            disp(['Calculation for ',num2str(a*b*c),' voxels (Memory < ',num2str(roundsf(max(hand.Results.MBytes),2,'ceil')),' MB, Time < ',timeStr,').',char(10)])
        else
            disp(['Calculation for ',num2str(round(a*b*c/1000^2)),' million voxels (Memory < ',num2str(roundsf(max(hand.Results.MBytes)/1024,2,'ceil')),' GB, Time < ',timeStr,').',char(10)])
        end
    end
    ExpectedTime(hObject, eventdata, hand);
    ExpectedRAM(hObject, eventdata, hand);
    if get(hand.Check_pdfSave,'Value')==1
        evalin(hand.WoSpace,['save(''',hand.pathname,hand.fil,''',''',hand.fil,''');']);
    end
end

function [hand]=PrintTauResults(hObject, eventdata, hand);
%% Print the results from each direction
if hand.InLineMode==0
    if hand.Results.Tau~=inf
        disp([...
            'Tau Factor = ',num2str(roundsf(hand.Results.Tau(end),4,'round')),' in direction ',num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.'... %,char(10),...
            ]);
    else
        disp([...
            'Tau Factor = inf in direction ',num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.'... %,char(10),...
            ]);
    end
end

function Check_Impedance_Callback(hObject, eventdata, hand)
%% Recalculate expected time if impedance is selected
[hand]=ExpectedTime(hObject, eventdata, hand);
guidata(hObject, hand);


function [hand]=Check_Blocked_Callback(hObject, eventdata, hand)
%% Recalculate expected time if impedance is selected
if get(hand.Check_Blocked, 'Value')==1;
    hand.Blocked=1;
else
    hand.Blocked=0;
end
guidata(hObject, hand);

function [hand]=RVA_Tau(hand)
%% Prepare subvolume of data for simulation is necessary
if get(hand.Check_RVA, 'Value')==1;
    hand.Net_full=hand.Net;
    if get(hand.Check_VaryD,'value')==1
        hand.Dmap_full=hand.Dmap(2:end-1,2:end-1,2:end-1);
    end
    steps=10;
    shrinkStep=1/steps;
    for i=steps:-1:1
        hand.shrink=(i-1)*shrinkStep;%volume percent of shrinkage
        switch get(hand.Pop_RV, 'Value')
            case 1
                [hand]=RVA_Net_Cube(hand);
            case 2
                [hand]=RVA_Net_Lconst(hand);
            case 3
                [hand]=RVA_Net_Aconst(hand);
            case 4
                [hand]=RVA_Net_Aconst(hand);
        end
        hand.Net=hand.Net_full(...
            hand.RVAdims(1,1):hand.RVAdims(1,2),...
            hand.RVAdims(2,1):hand.RVAdims(2,2),...
            hand.RVAdims(3,1):hand.RVAdims(3,2));
        if get(hand.Check_VaryD,'value')==1
            hand.Dmap=padarray(hand.Dmap_full(...
                hand.RVAdims(1,1):hand.RVAdims(1,2),...
                hand.RVAdims(2,1):hand.RVAdims(2,2),...
                hand.RVAdims(3,1):hand.RVAdims(3,2)),[1,1,1],inf);
        end
        hand.NetVol(steps+1-i)=numel(hand.Net)/numel(hand.Net_full);
        [hand]=Initialise(hand);
        if sum(hand.Net_Perc(:))~=0
            hand.Results.Tau=mean([hand.TauFacTop(end),hand.TauFacBot(end)]);
        else
            hand.Results.Tau=inf;
        end
        hand.Tau_RV(steps+1-i)=hand.Results.Tau;
        hand.Pore_RV(steps+1-i)=hand.VolFrac(end);
        
        if hand.InLineMode==0
            varname=[hand.fil,'.TauRVA',num2str(get(hand.Pop_RV, 'Value')),'_',num2str(hand.Pha(1)),num2str(num2str(hand.Dir))];
            assignin(hand.WoSpace,'temp',hand.NetVol);
            evalin(hand.WoSpace,[varname,'.NetVol = temp;']);
            
            assignin(hand.WoSpace,'temp',hand.Pore_RV);
            evalin(hand.WoSpace,[varname,'.VoleFrac = temp;']);
            
            assignin(hand.WoSpace,'temp',hand.Tau_RV);
            evalin(hand.WoSpace,[varname,'.Tau = temp;']);
            
            evalin(hand.WoSpace,'clear temp');
        else
            varname=['OutputResults.TauRVA',num2str(get(hand.Pop_RV, 'Value')),'_',num2str(hand.Pha(1)),num2str(num2str(hand.Dir))];
            eval(['hand.',varname,'.NetVol = hand.NetVol;']);
            eval(['hand.',varname,'.VoleFrac = hand.Pore_RV;']);
            eval(['hand.',varname,'.Tau = hand.Tau_RV;']);
        end
        %%%plotting
        if hand.InLineMode==0
            if i==steps
                hand.RVfig=figure(...
                    'Name',['TF_RVA: ','p',hand.Pha(1),'d',num2str(hand.Dir),'_',hand.fil],...
                    'Color',[1 1 1],...
                    'renderer','painters',...
                    'WindowStyle','normal',...
                    'PaperPositionMode','auto',...
                    'PaperOrientation','landscape');
            end
            if ~isnan(hand.VolFrac_Perc) && hand.VolFrac_Perc~=0
                
                figure(hand.RVfig);
                [hAx,hLine(1),hLine(2)] = plotyy(hand.NetVol,hand.Pore_RV,hand.NetVol,hand.Tau_RV); %use yyaxis in future
                legend('Volume Fraction','Tortuosity Factor')
                set(legend,'Interpreter','latex');
                xlabel('Fraction of original volume','Interpreter','Latex')
                ylabel(hAx(1),'Volume Fraction of Conductive Phase','Interpreter','Latex') % left y-axis
                ylabel(hAx(2),'Directional Tortuosity Factor','Interpreter','Latex') % right y-axis
                set(hAx(1),'ylim', [0 1]);
                set(hAx(1),'ytick',[0:0.2:1]);
                %                 set(hAx(2),'ylim', [1         roundsf(max(hand.Tau_RV(isfinite(hand.Tau_RV))+0.5),2,'ceil')]);
                %                 set(hAx(2),'ytick',[1 : 1 :   ceil(max(hand.Tau_RV(isfinite(hand.Tau_RV))+0.5))]);
                set(hAx,'TickLabelInterpreter','latex')
                set(hAx,'position',[0.1500 0.100 0.70 0.80])
                set(hAx,'XLim',[0 1])
                set(hAx,'LineWidth',1)
                set(hLine,'LineWidth',1.5)
                switch get(hand.Pop_RV, 'Value')
                    case 1
                        title({['\verb|',hand.filename,'|'];['RVA (Cubic) in direction ',...
                            num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.']},'Interpreter','latex');
                    case 2
                        title({['\verb|',hand.filename,'|'];['RVA (L=const.) in direction ',...
                            num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.']},'Interpreter','latex');
                    case 3
                        title({['\verb|',hand.filename,'|'];['RVA (A=const. from top) in direction ',...
                            num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.']},'Interpreter','latex');
                    case 4
                        title({['\verb|',hand.filename,'|'];['RVA (A=const. from base) in direction ',...
                            num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.']},'Interpreter','latex');
                        
                end
            end
            if get(hand.Check_pdfSave,'Value')==1 && i==1
                print(hand.RVfig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir),'_TauRVA'],'-dpdf');
            end
        end
    end
    hand.Pore_RV=1;
    hand.Tau_RV=1;
    hand.NetVol=1;
    hand.RVfig=1;
else
    [a,b,c]=size(hand.Net);
    hand.RVAdims=[...
        1,a;...
        1,b;...
        1,c];
    [hand]=Initialise(hand);
    if isfield(hand,'VolFrac')
        if hand.InLineMode==0
            varname=[hand.fil,'.Tau','_',num2str(hand.Pha(1)),num2str(num2str(hand.Dir))];
            
            assignin(hand.WoSpace,'temp',hand.VolFrac);
            evalin(hand.WoSpace,[varname,'.VoleFrac = temp;']);
            
            assignin(hand.WoSpace,'temp',hand.Results.Tau);
            evalin(hand.WoSpace,[varname,'.Tau = temp;']);
            
            evalin(hand.WoSpace,'clear temp');
        else
            varname=['OutputResults.Tau','_',num2str(hand.Pha(1)),num2str(num2str(hand.Dir))];
            
            eval(['hand.',varname,'.VoleFrac = hand.VolFrac;']);
            eval(['hand.',varname,'.Tau = hand.Results.Tau;']);
        end
        hand.Results.VolFrac=hand.VolFrac;
        if max(hand.Net_Perc(:))>0
            hand=rmfield(hand,{'Tconv','XFlux','T1Top','T2Top','T1Bot','T2Bot','Map',});
        end
    end
    hand.TauFacTop=0;
    hand.TauFacBot=0;
    hand.DeffTop=0;
    hand.DeffBot=0;
end

function [hand]=Initialise(hand)
%% First step of preperation for tortuosity calculation in which the ordering of calls is set
tic
hand.Results.SimTime=0;
if get(hand.Check_Impedance,'Value')==1
    hand.impCheck=1;
end
if get(hand.Check_Reverse, 'Value')==1 && hand.impCheck==1
    hand.Net=flipud(hand.Net);
end
[hand]=Percolation(hand);
hand.Net_Perc=logical(hand.Net_Perc);
[a,b,c]=size(hand.Net_Perc);
hand.VolFrac_Perc=sum(hand.Net_Perc(:))/sum(hand.Net(:));
if sum(sum(hand.Net_Perc(end,:,:)))==0% || Max_Path<5
    hand.Results.SimTime='Does not percolate in this phase+direction';
    hand.VolFrac=mean(hand.Net(:));
    if hand.VolFrac>0
        hand.Results.Tau=inf;
    else
        hand.Results.Tau=nan;
    end
    %     hand.Results.MBytes=nan;
else
    hand.impCheck=0;
    [hand]=Preparation1(hand);
    [hand]=Preparation2(hand);
    [hand]=Preparation3(hand);
    if get(hand.Check_Impedance, 'Value')==0
        eval(['[hand]=Iterate_',hand.mexCtrl,'(hand);'])
    else
        hand.impCheck=0;
        if sum(sum(hand.Map(2,2:end-1,2:end-1)))~=0 && hand.Blocked==0
            eval(['[hand]=Iterate_',hand.mexCtrl,'(hand);'])
            hand.Tau=mean([hand.TauFacBot(end),hand.TauFacTop(end)]);
        else
            hand.Results.Tau=inf;
            hand.Tau=1;
            hand.Tconv=single(hand.Map);
        end
        hand.impCheck=1;
        hand.D=complex(1);
        hand.delta_x=complex(hand.L_X);
        %         freqChar=2*pi*hand.D/(a*hand.Tau(end)*hand.delta_x)^2
        if hand.PercFlag==1
            hand.freqChar=hand.D/(a*hand.delta_x)^2;
        else
            hand.freqChar=hand.D/(hand.Max_Path*hand.delta_x)^2;
        end
        %         if sum(sum(sum(hand.Net_Perc)))==0
        %             freqChar=hand.D/(a*2*hand.delta_x)^2;
        %         else
        %             freqChar=hand.D/(a*hand.TauFacBot(end)*hand.delta_x)^2;
        %         end
        if sum(sum(hand.Map(2,2:end-1,2:end-1)))~=0 && hand.Blocked==0
            hand.y=-4:0.1:11;%-4:0.5:11;%-2:0.5:10
        else
            hand.y=-2:0.1:11;%-4:0.5:11;%-2:0.5:10
        end
        tic
        hand.impFig=figure(...
            'Name',['TF_Impedance: ','p',hand.Pha(1),'d',num2str(hand.Dir),'_',hand.fil],...
            'units','characters',...
            'position',[277 35 99 40],...
            'renderer','painters',...
            'Color',[1 1 1]);
        hand.freqSet=complex(hand.freqChar*2.^hand.y);
        hand.ImpedanceBotConv=ones(length(hand.y),1)*nan;
        for freqNo=1:length(hand.y)
            %             hand.y(freqNo)
            hand.freqNo=freqNo;
            hand.freqSet=complex(hand.freqChar*2.^hand.y);
            hand.freq=complex(hand.freqSet(freqNo));
            if hand.Check_FreqPlots==1
                [hand]=InitiatePlot1(hand);
            end
            [hand]=Preparation3imp(hand);
            if hand.Check_FreqPlots==1
                [hand]=InitiatePlot2(hand,(hand.Tconv(:,:,2)),hand.Map);
            end
            eval(['[hand]=Iterate_',hand.mexCtrl,'(hand);'])
            
            [hand]=ImpPlot(hand);
            if (abs(imag(hand.ImpedanceBotConv(freqNo))/max(real(hand.ImpedanceBotConv(:))))<5*hand.conTol &&...
                    abs(real(hand.ImpedanceBotConv(freqNo))/max(real(hand.ImpedanceBotConv(:))))<5*hand.conTol)...
                    || hand.ConvErr==2
                hand.ImpedanceBotConv=1; hand.ImpedanceBot=1;
                hand.ImpedanceTopConv=1; hand.ImpedanceTop=1;
                hand.TauFacBot=1; hand.TauFacTop=1;
                hand.DeffBot=1; hand.DeffTop=1;
                hand.freqSet=1;
                hand.ConvErr=0;
                % sound(1,1000);pause(0.02);sound(1,1000);
                hand.Results.SimTime=hand.Results.SimTime+toc;
                [hand]=Saving(hand);
                return
            end
        end
        [hand]=Saving(hand);
    end
    
    hand=rmfield(hand,'NN_aV');
    if isfield(hand,'R')
        hand=rmfield(hand,'R');
    end
    
end
[hand]=MBytesRecord(hand,whos); %Memory

function [hand]=Saving(hand)
if get(hand.Check_pdfSave,'Value')==1
    Name=[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir)];
    GifWait=2/length(hand.y);
    if get(hand.Check_RVA, 'Value')==0
        print(hand.impFig,[Name,'_Imp'],'-dpdf');
    else
        print(hand.impFig,[Name,'_Imp_RVA',num2str(round(100*hand.NetVol(end)))],'-dpdf');
    end
    print(hand.impFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir),'_Imp'],'-dpdf');
    hand.im=0;
    for i=1:hand.freqNo
        if i==1
            imwrite(hand.impFigGif1(:,:,i),hand.cccmap2,[Name,'plot.gif'],'gif','DelayTime',GifWait,'LoopCount',inf)
            imwrite(hand.impFigGif2(:,:,i),hand.cccmap2,[Name,'real.gif'],'gif','DelayTime',GifWait,'LoopCount',inf)
        else
            imwrite(hand.impFigGif1(:,:,i),[Name,'plot.gif'],'gif','WriteMode','append','DelayTime',GifWait)
            imwrite(hand.impFigGif2(:,:,i),[Name,'real.gif'],'gif','WriteMode','append','DelayTime',GifWait)
        end
    end
end

function [hand]=Percolation(hand)
Mask=false(size(hand.Net));
Mask(end,:,:)=hand.Net(end,:,:);
DM=bwdistgeodesic(logical(hand.Net),Mask,'cityblock');

hand.Min_Path=min(min(DM(1, :,:)))+1;
hand.Max_Path=max(DM(isfinite(DM)))+1;

hand.Net_Perc=isfinite(DM);

if hand.impCheck==0
    Mask(end,:,:)=0;
    Mask(1,:,:)=hand.Net_Perc(1,:,:);
    DM=bwdistgeodesic(logical(hand.Net_Perc),Mask,'cityblock');
    hand.Net_Perc=(isfinite(DM));
end

hand.MAA = uint16(sort(DM(isfinite(DM))));%
hand.MAA = uint16(find([numel(hand.MAA);diff(hand.MAA);numel(hand.MAA)]));
hand.MAA=mean(diff(hand.MAA));
if sum(sum(hand.Net_Perc(1,:,:)))==0
    hand.PercFlag=0;
else
    hand.PercFlag=1;
end
%% Plotting ditance maps with cut-off areas
% GrayPore=1-cat(3,Net_Perc,Net_Perc,Net_Perc);
% imagesc(DM);colormap(jet);
% hold on; h1=subimage(0*GrayPore);axis square;set(h1, 'AlphaData', 1-Net);hold off
% NonPerc=DM==inf;
% NonPerc=cat(3,NonPerc,NonPerc,NonPerc);
% hold on; h2=subimage(0.5*NonPerc);axis square;set(h2, 'AlphaData', NonPerc(:,:,1));hold off
% set(gca, 'XTick', [],'YTick',[]);

function [hand]=Preparation1(hand)
%% Second data preparation step where the checkerboarding and vectorisation of the is performed.
% This is also where the anisotropic weightings are accounted for
hand.VolFrac=sum(hand.Net(:))/numel(hand.Net);
[a,b,c]=size(hand.Net_Perc);
hand.Net=1;
switch hand.Dir
    case 1
        hand.L_X=1e-9*str2double(get(hand.L1box,'String'));
        hand.L_Y=1e-9*str2double(get(hand.L2box,'String'));
        hand.L_Z=1e-9*str2double(get(hand.L3box,'String'));
    case 2
        hand.L_X=1e-9*str2double(get(hand.L2box,'String'));
        hand.L_Y=1e-9*str2double(get(hand.L3box,'String'));
        hand.L_Z=1e-9*str2double(get(hand.L1box,'String'));
    case 3
        hand.L_X=1e-9*str2double(get(hand.L3box,'String'));
        hand.L_Y=1e-9*str2double(get(hand.L1box,'String'));
        hand.L_Z=1e-9*str2double(get(hand.L2box,'String'));
end
if str2double(get(hand.L1box,'String'))~=str2double(get(hand.L2box,'String')) ||...
        str2double(get(hand.L2box,'String'))~=str2double(get(hand.L3box,'String'));
    hand.Aniso=1;
else
    hand.Aniso=0;
end
if get(hand.Check_VaryD,'value')==0
    front=0.5*(hand.L_X*hand.L_Z/hand.L_Y + hand.L_Z*hand.L_Y/hand.L_X + hand.L_X*hand.L_Y/hand.L_Z)^-1;
    hand.c_X=double(front*hand.L_Y*hand.L_Z/hand.L_X);
    hand.c_Y=double(front*hand.L_Z*hand.L_X/hand.L_Y);
    hand.c_Z=double(front*hand.L_X*hand.L_Y/hand.L_Z);
    
    hand.D=1;
end
hand.Qcv=hand.D*(hand.L_Y*b*hand.L_Z*c)/(hand.L_X*a);
[hand]=MBytesRecord(hand,whos,'Preparation1 end'); %Memory

function [hand]=Preparation2(hand)
[a,b,c]=size(hand.Net_Perc);
% Generate maps of nearest neighbours
hand.Map=logical(padarray(hand.Net_Perc, [1,1,1],0));
hand.Net_Perc=1;
% Calculate adjusted map of adjusted nearest neighbours
if get(hand.Check_VaryD,'value')==0
    hand.NN_a=zeros(size(hand.Map),'double');
    hand.NN_a(2:end-1,2:end-1,2:end-1)=...
        hand.c_X*double((hand.Map(1:end-2,2:end-1,2:end-1)+hand.Map(3:end  ,2:end-1,2:end-1)))+...
        hand.c_Y*double((hand.Map(2:end-1,1:end-2,2:end-1)+hand.Map(2:end-1,3:end  ,2:end-1)))+...
        hand.c_Z*double((hand.Map(2:end-1,2:end-1,1:end-2)+hand.Map(2:end-1,2:end-1,3:end  )));
    if hand.Blocked==1
        hand.NN_a(end-1,:,:)=double(hand.Map(end-1,:,:)).*(hand.NN_a(end-1,:,:)+(2*hand.c_X));
    else
        hand.NN_a([2 end-1],:,:)=double(hand.Map([2 end-1],:,:)).*(hand.NN_a([2 end-1],:,:)+(2*hand.c_X));
    end
    if hand.Aniso==0
        hand.NN_a=hand.NN_a/hand.c_X;
    end
    hand.Dmap=1;
    [hand]=MBytesRecord(hand,whos,'Preparation2 NN_A'); %Memory
else
    hand.Dmap(:,[1, end],:)=inf;hand.Dmap(:,:,[1, end])=inf;
    hand.Dmap([1, end],2:end-1,2:end-1)=0;% Easier than changing dx
    hand.NN_a=zeros(size(hand.Dmap));
    [hand]=MBytesRecord(hand,whos,'Dmap'); %Memory
    hand.R.Xm=hand.NN_a;hand.R.Xp=hand.NN_a;
    hand.R.Ym=hand.NN_a;hand.R.Yp=hand.NN_a;
    if c>1
        hand.R.Zm=hand.NN_a;hand.R.Zp=hand.NN_a;
    end
    [hand]=MBytesRecord(hand,whos,'R.XYZ'); %Memory
    hand.NN_a(2:end-1,2:end-1,2:end-1)=2*hand.Map(2:end-1,2:end-1,2:end-1).*(...
        hand.L_Y*hand.L_Z/hand.L_X*(...
        hand.Map(1:end-2,2:end-1,2:end-1)./(hand.Dmap(1:end-2,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1))+...
        hand.Map(3:end  ,2:end-1,2:end-1)./(hand.Dmap(3:end  ,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1)))...
        +hand.L_Z*hand.L_X/hand.L_Y*(...
        hand.Map(2:end-1,1:end-2,2:end-1)./(hand.Dmap(2:end-1,1:end-2,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1))+...
        hand.Map(2:end-1,3:end  ,2:end-1)./(hand.Dmap(2:end-1,3:end  ,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1)))...
        +hand.L_X*hand.L_Y/hand.L_Z*(...
        hand.Map(2:end-1,2:end-1,1:end-2)./(hand.Dmap(2:end-1,2:end-1,1:end-2)+hand.Dmap(2:end-1,2:end-1,2:end-1))+...
        hand.Map(2:end-1,2:end-1,3:end  )./(hand.Dmap(2:end-1,2:end-1,3:end  )+hand.Dmap(2:end-1,2:end-1,2:end-1))));
    hand.NN_a([2 end-1],2:end-1,2:end-1)=2*hand.Map([2 end-1],2:end-1,2:end-1).*(...
        hand.L_Y*hand.L_Z/hand.L_X*(...
        hand.Map([2 end-2],2:end-1,2:end-1)./(hand.Dmap([1 end-2],2:end-1,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1))+...
        hand.Map([3 end-1],2:end-1,2:end-1)./(hand.Dmap([3 end  ],2:end-1,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1)))...
        +hand.L_Z*hand.L_X/hand.L_Y*(...
        hand.Map([2 end-1],1:end-2,2:end-1)./(hand.Dmap([2 end-1],1:end-2,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1))+...
        hand.Map([2 end-1],3:end  ,2:end-1)./(hand.Dmap([2 end-1],3:end  ,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1)))...
        +hand.L_X*hand.L_Y/hand.L_Z*(...
        hand.Map([2 end-1],2:end-1,1:end-2)./(hand.Dmap([2 end-1],2:end-1,1:end-2)+hand.Dmap([2 end-1],2:end-1,2:end-1))+...
        hand.Map([2 end-1],2:end-1,3:end  )./(hand.Dmap([2 end-1],2:end-1,3:end  )+hand.Dmap([2 end-1],2:end-1,2:end-1))));
    
    
    hand.R.Xm(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Y*hand.L_Z/hand.L_X./(hand.Dmap(1:end-2,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    hand.R.Xp(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Y*hand.L_Z/hand.L_X./(hand.Dmap(3:end  ,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    
    hand.R.Ym(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Z*hand.L_X/hand.L_Y./(hand.Dmap(2:end-1,1:end-2,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    hand.R.Yp(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Z*hand.L_X/hand.L_Y./(hand.Dmap(2:end-1,3:end  ,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    if c>1
        hand.R.Zm(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_X*hand.L_Y/hand.L_Z./(hand.Dmap(2:end-1,2:end-1,1:end-2)+hand.Dmap(2:end-1,2:end-1,2:end-1));
        hand.R.Zp(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_X*hand.L_Y/hand.L_Z./(hand.Dmap(2:end-1,2:end-1,3:end  )+hand.Dmap(2:end-1,2:end-1,2:end-1));
    else
        
        hand.R.Zm=1;
        hand.R.Zp=1;
    end
end

hand.Cheq1.P=zeros([a,b,c],'uint8');
% Build checkerboard
hand.Cheq1.P(1:2:end)=1;
if rem(a,2)==0
    hand.Cheq1.P(:,2:2:end,:)=1-hand.Cheq1.P(:,2:2:end,:);
end
if rem(a*b,2)==0
    hand.Cheq1.P(:,:,2:2:end)=1-hand.Cheq1.P(:,:,2:2:end);
end
% Find checkboard neighbours
% N=North, S=South, E=East, W=West, U=Up, D=down
hand.Cheq2.P=1-hand.Cheq1.P;
hand.Cheq1.P=padarray(hand.Cheq1.P>0,[1,1,1]);hand.Cheq1.P(hand.Map==0)=0;
hand.Cheq1.P=uint32(find(hand.Cheq1.P));
hand.Cheq1.P_Xm=hand.Cheq1.P-1;
hand.Cheq1.P_Xp=hand.Cheq1.P+1;
hand.Cheq1.P_Ym=hand.Cheq1.P-(a+2);
hand.Cheq1.P_Yp=hand.Cheq1.P+(a+2);
if c>1
    hand.Cheq1.P_Zm=hand.Cheq1.P-(a+2)*(b+2);
    hand.Cheq1.P_Zp=hand.Cheq1.P+(a+2)*(b+2);
else
    hand.Cheq1.P_Zm=uint32(1);
    hand.Cheq1.P_Zp=uint32(1);
end
[hand]=MBytesRecord(hand,whos,'Prep2 Cheq'); %Memory

hand.Cheq2.P=padarray(hand.Cheq2.P>0,[1,1,1]);hand.Cheq2.P(hand.Map==0)=0;
hand.Cheq2.P=uint32(find(hand.Cheq2.P));
hand.Cheq2.P_Xm=hand.Cheq2.P-1;
hand.Cheq2.P_Xp=hand.Cheq2.P+1;
hand.Cheq2.P_Ym=hand.Cheq2.P-(a+2);
hand.Cheq2.P_Yp=hand.Cheq2.P+(a+2);
if c>1
    hand.Cheq2.P_Zm=hand.Cheq2.P-(a+2)*(b+2);
    hand.Cheq2.P_Zp=hand.Cheq2.P+(a+2)*(b+2);
else
    hand.Cheq2.P_Zm=uint32(1);
    hand.Cheq2.P_Zp=uint32(1);
end

if ~isfield(hand,'MBytesA')
    hand.MBytesA=0;
    hand.MemLoc={''};
end
[hand]=MBytesRecord(hand,whos,'Prep2 Cheq end'); %Memory
%%
% New method
Top=uint32([1:a+2:(a+2)*(b+2)*(c+2)]);
Base=uint32([a+2:a+2:(a+2)*(b+2)*(c+2)]);

[LIA2t]=ismember(hand.Cheq2.P_Xm,Top);
[LIA2b]=ismember(hand.Cheq2.P_Xp,Base);
hand.T2Top=uint32(find(LIA2t));
hand.T2Bot=uint32(find(LIA2b));
[~,hand.Cheq2.P_Xm]=(ismember(hand.Cheq2.P_Xm,hand.Cheq1.P));
hand.Cheq2.P_Xm=uint32(hand.Cheq2.P_Xm);
hand.Cheq2.P_Xm(hand.Cheq2.P_Xm==0)=uint32(length(hand.Cheq1.P)+1);
[~,hand.Cheq2.P_Xp]=(ismember(hand.Cheq2.P_Xp,hand.Cheq1.P));
hand.Cheq2.P_Xp=uint32(hand.Cheq2.P_Xp);
hand.Cheq2.P_Xp(hand.Cheq2.P_Xp==0)=uint32(length(hand.Cheq1.P)+1);
[~,hand.Cheq2.P_Ym]=(ismember(hand.Cheq2.P_Ym,hand.Cheq1.P));
hand.Cheq2.P_Ym=uint32(hand.Cheq2.P_Ym);
hand.Cheq2.P_Ym(hand.Cheq2.P_Ym==0)=uint32(length(hand.Cheq1.P)+1);
[~,hand.Cheq2.P_Yp]=(ismember(hand.Cheq2.P_Yp,hand.Cheq1.P));
hand.Cheq2.P_Yp=uint32(hand.Cheq2.P_Yp);
hand.Cheq2.P_Yp(hand.Cheq2.P_Yp==0)=uint32(length(hand.Cheq1.P)+1);
if c>1
    [~,hand.Cheq2.P_Zm]=ismember(hand.Cheq2.P_Zm,hand.Cheq1.P);
    hand.Cheq2.P_Zm=uint32(hand.Cheq2.P_Zm);
    hand.Cheq2.P_Zm(hand.Cheq2.P_Zm==0)=length(hand.Cheq1.P)+1;
    [~,hand.Cheq2.P_Zp]=ismember(hand.Cheq2.P_Zp,hand.Cheq1.P);
    hand.Cheq2.P_Zp=uint32(hand.Cheq2.P_Zp);
    hand.Cheq2.P_Zp(hand.Cheq2.P_Zp==0)=length(hand.Cheq1.P)+1;
end
hand.Cheq2.P_Xp(LIA2b)=uint32(length(hand.Cheq1.P)+2);

[LIA1t]=ismember(hand.Cheq1.P_Xm,Top);
[LIA1b]=ismember(hand.Cheq1.P_Xp,Base);
hand.T1Top=uint32(find(LIA1t));
hand.T1Bot=uint32(find(LIA1b));
[~,hand.Cheq1.P_Xm]=(ismember(hand.Cheq1.P_Xm,hand.Cheq2.P));
hand.Cheq1.P_Xm=uint32(hand.Cheq1.P_Xm);
hand.Cheq1.P_Xm(hand.Cheq1.P_Xm==0)=uint32(length(hand.Cheq2.P)+1);
[~,hand.Cheq1.P_Xp]=(ismember(hand.Cheq1.P_Xp,hand.Cheq2.P));
hand.Cheq1.P_Xp=uint32(hand.Cheq1.P_Xp);
hand.Cheq1.P_Xp(hand.Cheq1.P_Xp==0)=uint32(length(hand.Cheq2.P)+1);
[~,hand.Cheq1.P_Ym]=(ismember(hand.Cheq1.P_Ym,hand.Cheq2.P));
hand.Cheq1.P_Ym=uint32(hand.Cheq1.P_Ym);
hand.Cheq1.P_Ym(hand.Cheq1.P_Ym==0)=uint32(length(hand.Cheq2.P)+1);
[~,hand.Cheq1.P_Yp]=(ismember(hand.Cheq1.P_Yp,hand.Cheq2.P));
hand.Cheq1.P_Yp=uint32(hand.Cheq1.P_Yp);
hand.Cheq1.P_Yp(hand.Cheq1.P_Yp==0)=uint32(length(hand.Cheq2.P)+1);
if c>1
    [~,hand.Cheq1.P_Zm]=ismember(hand.Cheq1.P_Zm,hand.Cheq2.P);
    hand.Cheq1.P_Zm=uint32(hand.Cheq1.P_Zm);
    hand.Cheq1.P_Zm(hand.Cheq1.P_Zm==0)=length(hand.Cheq2.P)+1;
    [~,hand.Cheq1.P_Zp]=ismember(hand.Cheq1.P_Zp,hand.Cheq2.P);
    hand.Cheq1.P_Zp=uint32(hand.Cheq1.P_Zp);
    hand.Cheq1.P_Zp(hand.Cheq1.P_Zp==0)=length(hand.Cheq2.P)+1;
end
hand.Cheq1.P_Xp(LIA1b)=uint32(length(hand.Cheq2.P)+2);

function [hand]=Preparation3(hand)
%% Third preparation step for Tau calculation where the volume is initialised as linear
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
hand.Area_top=sum(sum(hand.Map(2,:,:)));
hand.Area_bot=sum(sum(hand.Map(end-1,:,:)));
% Specify a relaxation factor

hand.w=2-(pi)/(a*1.5);
if get(hand.Check_VaryD,'value')==1
    %     hand.w=1;
end
% hand.w=1.967;
% hand.w=1;
hand.omw=double(1-hand.w);
% Seperate nearest neighbours into checkerboard and precalculate relaxation
% factor and division
hand.NN_aV.w1=double(hand.w./double(hand.NN_a(hand.Cheq1.P)));
hand.NN_aV.w2=double(hand.w./double(hand.NN_a(hand.Cheq2.P)));
[hand]=MBytesRecord(hand,whos,'Prep3'); %Memory
if get(hand.Check_Impedance,'Value')==0
    hand.NN_a=0;hand=rmfield(hand,'NN_a');
end
[hand]=MBytesRecord(hand,whos,'Prep3 remove NN_a'); %Memory
T=single(hand.Map);
hand.TopStim=0;
hand.BotStim=1;
if get(hand.Check_VaryD,'value')==1
    hand.R.Xm1=hand.R.Xm(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Xp1=hand.R.Xp(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Xm2=hand.R.Xm(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R.Xp2=hand.R.Xp(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R=rmfield(hand.R,{'Xm','Xp'});
    
    hand.R.Ym1=hand.R.Ym(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Yp1=hand.R.Yp(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Ym2=hand.R.Ym(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R.Yp2=hand.R.Yp(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R=rmfield(hand.R,{'Ym','Yp'});
    if c>1
        hand.R.Zm1=hand.R.Zm(hand.Cheq1.P).*hand.NN_aV.w1;
        hand.R.Zp1=hand.R.Zp(hand.Cheq1.P).*hand.NN_aV.w1;
        hand.R.Zm2=hand.R.Zm(hand.Cheq2.P).*hand.NN_aV.w2;
        hand.R.Zp2=hand.R.Zp(hand.Cheq2.P).*hand.NN_aV.w2;
    else
        hand.R.Zm1=1;
        hand.R.Zp1=1;
        hand.R.Zm2=1;
        hand.R.Zp2=1;
    end
    hand.R=rmfield(hand.R,{'Zm','Zp'});
    hand.NN_aV=1;
end
% T(1,2:end-1,2:end-1)=hand.Map(2,2:end-1,2:end-1)*2*hand.TopStim;
% T(end,2:end-1,2:end-1)=hand.Map(end-1,2:end-1,2:end-1)*2*hand.BotStim;
if hand.TopStim==hand.BotStim || hand.Blocked==1 || sum(sum(hand.Map(2,2:end-1,2:end-1)))==0
    T=single(hand.Map*hand.BotStim);
else
    step=(hand.BotStim-hand.TopStim)/a;
    V= (hand.TopStim+step/2:step:hand.BotStim);
    T(2:end-1,2:end-1,2:end-1)=single(ndgrid(V',1:b,1:c).*hand.Map(2:end-1,2:end-1,2:end-1));
end
hand.T1=double(T(hand.Cheq1.P));
hand.T2=double(T(hand.Cheq2.P));
if get(hand.Check_VaryD,'value')==0
    hand.T1(end+2)=2*hand.BotStim;
    hand.T2(end+2)=2*hand.BotStim;
else
    hand.T1(end+2)=hand.BotStim;
    hand.T2(end+2)=hand.BotStim;
end
if hand.InLineMode==0
    T=T(:,:,2);
    [hand]=InitiatePlot1(hand);
    [hand]=InitiatePlot2(hand,T,hand.Map);
end
[hand]=MBytesRecord(hand,whos,'Prep3'); %Memory

function [hand]=Preparation3imp(hand)
%% Alternative third preparation step if impedance is called
% hand.TopStim=0;
% hand.BotStim=1;
% T=complex(hand.Tconv);
% hand.Tconv=1;
% T(1,2:end-1,2:end-1)=2*hand.TopStim;
% T(end,2:end-1,2:end-1)=2*hand.BotStim;
hand.Tconv=complex((hand.Tconv));
hand.T1=double(hand.Tconv(hand.Cheq1.P));
hand.T2=double(hand.Tconv(hand.Cheq2.P));
if get(hand.Check_VaryD,'value')==0
    hand.T1(end+2)=2*hand.BotStim;
    hand.T2(end+2)=2*hand.BotStim;
else
    hand.T1(end+2)=hand.BotStim;
    hand.T2(end+2)=hand.BotStim;
end
if hand.InLineMode==0
    hand.Tconv=hand.Tconv(:,:,2);
    [hand]=InitiatePlot2(hand,(hand.Tconv),hand.Map);
end
hand.Area_top=complex(sum(sum(hand.Map(2,:,:))));
hand.Area_bot=complex(sum(sum(hand.Map(end-1,:,:))));
hand.w=complex(2-(pi)/max(size(hand.Map)*1.3));
if hand.PercFlag==1 && hand.Blocked==0
    if hand.y(hand.freqNo)>0.5
        hand.w=hand.w*0.90^(hand.y(hand.freqNo)-0.5);%-3
    end
else
    if hand.y(hand.freqNo)>-1
        hand.w=hand.w*0.93^(hand.y(hand.freqNo)+1);%-3
    end
end
if hand.Aniso==1
    hand.w=0.95*hand.w;
end
hand.omw=complex(1-hand.w);
hand.NN_aV.w1=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
    hand.D+complex(double(hand.NN_a(hand.Cheq1.P)))) ) );
hand.NN_aV.w2=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
    hand.D+complex(double(hand.NN_a(hand.Cheq2.P)))) ) );

[hand]=MBytesRecord(hand,whos,'Prep3Imp'); %Memory

function [hand]=Iterate_mex(hand)
%% Core iteration function for tortuosity factor calculation
hand.iter=0;
Error=1;
checkNo=1;
hand.whileFlag=1;
[~,~,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
[hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
if hand.impCheck==0
    % ensure variable type is appropriate
else
end
[hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
if hand.impCheck==0 %Not impedance
    if get(hand.Check_VaryD,'value')==0
        if hand.Aniso==0; % if the volume has isotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                [hand.T1,hand.T2]=Mex3DTauIso_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                    hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                    hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
                hand.iter=hand.iter+hand.check_f;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hand);
                end
            end
        elseif hand.Aniso==1 % if the volume has anisotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                [hand.T1,hand.T2]=Mex3DTauAni_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                    hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                    hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
                    hand.c_X,hand.c_Y,hand.c_Z);
                hand.iter=hand.iter+hand.check_f;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hand);
                end
            end
        end
    else % Variable D
        while hand.whileFlag>0 && hand.iter<hand.iter_max
            [hand.T1,hand.T2]=Mex3DTauVar_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,...
                hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
                hand.R.Xm1,hand.R.Xp1,hand.R.Ym1,hand.R.Yp1,hand.R.Zm1,hand.R.Zp1,...
                hand.R.Xm2,hand.R.Xp2,hand.R.Ym2,hand.R.Yp2,hand.R.Zm2,hand.R.Zp2);
            hand.iter=hand.iter+hand.check_f;
            if rem(hand.iter,hand.check_f)==0
                [hand]=Checks(hand);
            end
        end
    end
else %Impedance
    while hand.whileFlag>0 && hand.iter<hand.iter_max
        [hand.T1,hand.T2]=Mex3DTauIso_Imp_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,double(hand.NN_aV.w1),double(hand.NN_aV.w2),...
            hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
            hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
        hand.iter=hand.iter+hand.check_f;
        if rem(hand.iter,hand.check_f)==0
            [hand]=Checks(hand);
        end
    end
end

[hand]=MBytesRecord(hand,whos,'Iterate complete'); %Memory
[hand]=AfterIterate(hand);


function [hand]=AfterIterate(hand)
%%% Future potential to remove voxels with one neighbour at steady state,
%%% then replace for visualisation
% T(hand.NN==1)=T(2:end-1,2:end-1,2:end-1)+T(2:end-1,2:end-1,2:end-1)+...
% T(2:end-1,2:end-1,2:end-1)+T(2:end-1,2:end-1,2:end-1)+...
% T(2:end-1,2:end-1,2:end-1)+T(2:end-1,2:end-1,2:end-1);
[hand]=MBytesRecord(hand,whos,'Iterate complete'); %Memory
hand.Tconv=zeros(size(hand.Map),'single');
% Calculate vertical flux
hand.T1(end-1)=nan;
hand.T2(end-1)=nan;
hand.Tconv(hand.Cheq1.P)=single(hand.T2(hand.Cheq1.P_Xp)-hand.T1(1:end-2));
hand.Tconv(hand.Cheq2.P)=single(hand.T1(hand.Cheq2.P_Xp)-hand.T2(1:end-2));
hand.Tconv(~isfinite(hand.Tconv))=0;
if get(hand.Check_VaryD,'value')==1
    hand.XFlux=sum(sum(hand.Tconv(1:end-1,:,:)./(hand.Dmap(1:end-1,:,:)+hand.Dmap(2:end,:,:)),3),2);
else
    hand.XFlux=sum(sum(hand.Tconv(1:end-1,:,:),3),2);
end
hand.XFlux=flipud(hand.XFlux(2:end-1)/mean(hand.XFlux([2 end-1])));
hand.Tconv=zeros(size(hand.Map),'double');
hand.Tconv(hand.Cheq1.P)=(hand.T1(1:end-2));
hand=rmfield(hand,'T1');
hand.Tconv(hand.Cheq2.P)=(hand.T2(1:end-2));
hand=rmfield(hand,'T2');
if hand.InLineMode==0
    if get(hand.Check_FluxMapSave,'Value')==1 && hand.impCheck==0
        TifSave3D(hand.Tconv(2:end-1,2:end-1,2:end-1),hand,'Cmap');
        z=whos('hand');
        if z.bytes/1024^3<3*2
            save([hand.pathname,hand.fil,'_Cmap'],'-struct','hand','Tconv');
        else
            save([hand.pathname,hand.fil,'_Cmap'],'-struct','hand','Tconv', '-v7.3');
        end
    end
end
hand.Tconv=single(hand.Tconv);
if get(hand.Check_VaryD,'value')==0
    hand.Tconv(end,:,:)=single(2*hand.BotStim);
else
    hand.Tconv(end,:,:)=single(hand.BotStim);
end
if get(hand.Check_Impedance,'Value')==0
    hand.Cheq1=0;hand.Cheq2=0;hand=rmfield(hand,{'Cheq1','Cheq2'});
end
hand.Results.SimTime=TimeString(toc/86400);
hand.Results.Iterations=hand.iter;
if hand.impCheck==0
    hand.Results.Tau=mean([hand.TauFacBot(round(hand.iter/hand.check_f)),hand.TauFacTop(round(hand.iter/hand.check_f))]);
end
[hand]=MBytesRecord(hand,whos,'Tconv'); %Memory
hand.Results.MBytes=ceil(max(hand.MBytesA));
if hand.impCheck==0
    if hand.InLineMode==0
        [hand]=InitiatePlot3(hand);
    end
else
    if hand.Check_FreqPlots==1
        [hand]=InitiatePlot3imp(hand);
    end
end
[hand]=MBytesRecord(hand,whos,'Iterate end'); %Memory

function [hand]=Iterate_mat(hand)
%% Core iteration function for tortuosity factor calculation
hand.iter=0;
Error=1;
checkNo=1;
hand.whileFlag=1;
[~,~,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
[hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
if c>1 % if the volume is 3D
    if get(hand.Check_VaryD,'value')==0
        if hand.Aniso==0; % if the volume has isotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp)+...
                    hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp)+...
                    hand.T2(hand.Cheq1.P_Zm)+hand.T2(hand.Cheq1.P_Zp));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp)+...
                    hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp)+...
                    hand.T1(hand.Cheq2.P_Zm)+hand.T1(hand.Cheq2.P_Zp));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hand);
                end
                
            end
        elseif hand.Aniso==1 % if the volume has anisotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.c_X*(hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp))+...
                    hand.c_Y*(hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp))+...
                    hand.c_Z*(hand.T2(hand.Cheq1.P_Zm)+hand.T2(hand.Cheq1.P_Zp)));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.c_X*(hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp))+...
                    hand.c_Y*(hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp))+...
                    hand.c_Z*(hand.T1(hand.Cheq2.P_Zm)+hand.T1(hand.Cheq2.P_Zp)));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hand);
                end
            end
        end
    else % Variable D
        while hand.whileFlag>0 && hand.iter<hand.iter_max
            hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+...
                hand.R.Xm1.*hand.T2(hand.Cheq1.P_Xm)+...
                hand.R.Xp1.*hand.T2(hand.Cheq1.P_Xp)+...
                hand.R.Ym1.*hand.T2(hand.Cheq1.P_Ym)+...
                hand.R.Yp1.*hand.T2(hand.Cheq1.P_Yp)+...
                hand.R.Zm1.*hand.T2(hand.Cheq1.P_Zm)+...
                hand.R.Zp1.*hand.T2(hand.Cheq1.P_Zp);
            hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+...
                hand.R.Xm2.*hand.T1(hand.Cheq2.P_Xm)+...
                hand.R.Xp2.*hand.T1(hand.Cheq2.P_Xp)+...
                hand.R.Ym2.*hand.T1(hand.Cheq2.P_Ym)+...
                hand.R.Yp2.*hand.T1(hand.Cheq2.P_Yp)+...
                hand.R.Zm2.*hand.T1(hand.Cheq2.P_Zm)+...
                hand.R.Zp2.*hand.T1(hand.Cheq2.P_Zp);
            hand.iter=hand.iter+1;
            if rem(hand.iter,hand.check_f)==0
                [hand]=Checks(hand);
            end
        end
    end
else %% 2D version
    if get(hand.Check_VaryD,'value')==0
        if hand.Aniso==0; % if the volume has isotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp)+...
                    hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp)+...
                    hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hand);
                end
            end
        else % if the volume has anisotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.c_X*(hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp))+...
                    hand.c_Y*(hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp)));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.c_X*(hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp))+...
                    hand.c_Y*(hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp)));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hand);
                end
            end %iterations
        end
    else % Variable D
        while hand.whileFlag>0 && hand.iter<hand.iter_max
            hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+...
                hand.R.Xm1.*hand.T2(hand.Cheq1.P_Xm)+...
                hand.R.Xp1.*hand.T2(hand.Cheq1.P_Xp)+...
                hand.R.Ym1.*hand.T2(hand.Cheq1.P_Ym)+...
                hand.R.Yp1.*hand.T2(hand.Cheq1.P_Yp);
            hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+...
                hand.R.Xm2.*hand.T1(hand.Cheq2.P_Xm)+...
                hand.R.Xp2.*hand.T1(hand.Cheq2.P_Xp)+...
                hand.R.Ym2.*hand.T1(hand.Cheq2.P_Ym)+...
                hand.R.Yp2.*hand.T1(hand.Cheq2.P_Yp);
            hand.iter=hand.iter+1;
            if rem(hand.iter,hand.check_f)==0
                [hand]=Checks(hand);
            end
        end
    end
end
[hand]=MBytesRecord(hand,whos,'Iterate complete'); %Memory
[hand]=AfterIterate(hand);

function [hand]=Checks(hand)
hand.conDwell=5;
hand.conTol=0.005;
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
%% Function to check on the convergence status of the volume
checkNo=hand.iter/hand.check_f;
% Calculate the apparent tortuosity factors at the top and bottom faces
if get(hand.Check_VaryD,'value')==0
    sumT_top=sum(hand.T1(hand.T1Top))+sum(hand.T2(hand.T2Top));
    sumT_bot=sum(hand.T1(hand.T1Bot))+sum(hand.T2(hand.T2Bot));
    
    hand.DeffTop(checkNo)=abs(hand.Area_top*hand.TopStim-sumT_top)*2*hand.D*(a/(b*c));
    hand.DeffBot(checkNo)=abs(hand.Area_bot*hand.BotStim-sumT_bot)*2*hand.D*(a/(b*c));
    
    deltaTB=abs((abs(hand.DeffTop(checkNo))-abs(hand.DeffBot(checkNo)))/mean([hand.DeffTop(checkNo),hand.DeffBot(checkNo)]));
    
    hand.TauFacTop(checkNo)=hand.D*hand.VolFrac./hand.DeffTop(checkNo);
    hand.TauFacBot(checkNo)=hand.D*hand.VolFrac./hand.DeffBot(checkNo);
else
    
    F_top=sum([hand.T1(hand.T1Top);hand.T2(hand.T2Top)]./hand.Dmap([hand.Cheq1.P(hand.T1Top);hand.Cheq2.P(hand.T2Top)]))*(hand.L_Y*hand.L_Z/(0.5*hand.L_X));
    F_bot=sum((hand.BotStim-[hand.T1(hand.T1Bot);hand.T2(hand.T2Bot)])./hand.Dmap([hand.Cheq1.P(hand.T1Bot);hand.Cheq2.P(hand.T2Bot)]))*(hand.L_Y*hand.L_Z/(0.5*hand.L_X));
    deltaTB=abs((abs(F_top)-abs(F_bot))/mean([F_top,F_bot]));
    
    hand.DeffTop(checkNo)=abs(F_top/((hand.L_Y*b*hand.L_Z*c)/(hand.L_X*a)));
    hand.DeffBot(checkNo)=abs(F_bot/((hand.L_Y*b*hand.L_Z*c)/(hand.L_X*a)));
    
    hand.TauFacTop(checkNo)=hand.Qcv/abs(F_top);
    hand.TauFacBot(checkNo)=hand.Qcv/abs(F_bot);
end
% Check for convergence of tau at the two faces
if hand.impCheck~=1
    %     hand.SumT(checkNo)=sum(sum(sum(T(2:end-1,:,:))));
    if checkNo<3
        hand.whileFlag=hand.conDwell;
    else
        %         DsumT=abs((hand.SumT(checkNo)-hand.SumT(checkNo-1))/hand.SumT(checkNo-1))
        if      abs((hand.DeffTop(checkNo)-hand.DeffTop(checkNo-2))/hand.DeffTop(checkNo))<hand.conTol
            if      abs(hand.DeffTop(checkNo)-hand.DeffTop(checkNo-1))/hand.DeffTop(checkNo)<hand.conTol &&...
                    abs(hand.DeffTop(checkNo)-hand.DeffBot(checkNo))<hand.conTol &&...
                    deltaTB<0.04 &&...
                    (mean([hand.TauFacTop(checkNo),hand.TauFacBot(checkNo)])-mean([hand.TauFacTop(checkNo-1),hand.TauFacBot(checkNo-1)]))/hand.TauFacTop(checkNo)<hand.conTol
                hand.whileFlag=hand.whileFlag-1;
                [hand]=MBytesRecord(hand,whos,'Checks start'); %Memory
            else
                if hand.whileFlag~=hand.conDwell
                    hand.whileFlag=hand.conDwell;
                end
            end
        end
        %%%%%%% Increase relaxation factor progressively?
        %         if      abs(hand.TauFacTop(checkNo)-hand.TauFacTop(checkNo-1))<10*hand.conTol &&...
        %                 abs(hand.TauFacTop(checkNo)-hand.TauFacBot(checkNo))<10*hand.conTol
        %             if hand.w==2-(pi)/max(size(T)*1.50);
        %                 %         [hand,~]=Preparation3(hand);
        %                 hand.w=2-(pi)/max(size(T)*1.55);
        %                 hand.omw=double(1-hand.w);
        %                 hand.NN_aV.w1=double(hand.w./double(hand.NN_a(hand.Cheq1.P)));
        %                 hand.NN_aV.w2=double(hand.w./double(hand.NN_a(hand.Cheq2.P)));
        %                 hand.w_increaseFlag=1;
        %             end
        %         end
    end
else % Impedance mode
    TauRat=(hand.L_Y*hand.L_Z)/hand.L_X;
    hand.ImpedanceTop(hand.freqNo,checkNo)=-(0.5/(TauRat*(hand.Area_top*hand.TopStim-sumT_top)/(1/hand.D)));
    hand.ImpedanceBot(hand.freqNo,checkNo)= (0.5/(TauRat*(hand.Area_bot*hand.BotStim-sumT_bot)/(1/hand.D)));
    if checkNo<6
        hand.whileFlag=hand.conDwell-1;
    else
        hand.conTol=0.01;
        %         [mean(abs([hand.T1(:); hand.T2(:)])),...
        %             abs(real(hand.ImpedanceBot(hand.freqNo,checkNo))/real(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1),...
        %             abs(imag(hand.ImpedanceBot(hand.freqNo,checkNo))/imag(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1),...
        %             abs(abs(hand.ImpedanceBot(hand.freqNo,checkNo))/abs(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1),...
        %             abs(angle(hand.ImpedanceBot(hand.freqNo,checkNo))/angle(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1)];
        if      abs(real(hand.ImpedanceBot(hand.freqNo,checkNo))/real(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1)<hand.conTol &&...
                abs(imag(hand.ImpedanceBot(hand.freqNo,checkNo))/imag(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1)<hand.conTol% &&...
            %                abs(real(hand.ImpedanceTop(hand.freqNo,checkNo))/real(hand.ImpedanceTop(hand.freqNo,checkNo-6))-1)<hand.conTol &&...
            %                abs(imag(hand.ImpedanceTop(hand.freqNo,checkNo))/imag(hand.ImpedanceTop(hand.freqNo,checkNo-6))-1)<hand.conTol
            hand.ImpedanceBotConv(hand.freqNo)=hand.ImpedanceBot(hand.freqNo,checkNo);
            hand.ImpedanceTopConv(hand.freqNo)=hand.ImpedanceTop(hand.freqNo,checkNo);
            hand.whileFlag=hand.whileFlag-1;
            
        else
            if hand.whileFlag~=hand.conDwell
                hand.whileFlag=hand.conDwell;
            end
        end
        if hand.iter>=hand.iter_max-hand.check_f
            hand.ImpedanceBotConv(hand.freqNo)=hand.ImpedanceBot(hand.freqNo,checkNo);
            hand.ImpedanceTopConv(hand.freqNo)=hand.ImpedanceTop(hand.freqNo,checkNo);
            Error=1;
            hand.ConvErr=hand.ConvErr+1;
            hand.iter_max=hand.iter_max*1.5;
            display('Simulation not converged')
            return
        end
    end
    %%%%%%     if unstable
    %         hand.w=hand.w*0.99;
    %         hand.omw=complex(1-hand.w);
    %
    %         hand.NN_aV.w1=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
    %             hand.D+complex(double(hand.NN_a(hand.Cheq1.P)))) ) );
    %         hand.NN_aV.w2=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
    %             hand.D+complex(double(hand.NN_a(hand.Cheq2.P)))) ) );
    %     end
end
Error=1;
% hand.TauFacTop(checkNo)=TauRat*hand.VolFrac*hand.Qcv/(2*sum(sum(T(2,:,:)))); %Should be Perhand.VolFrac(counter)
% hand.TauFacBot(checkNo)=TauRat*hand.VolFrac*hand.Qcv/(-2*(sum(sum( T(end-1,:,:)))-hand.Area_bot));

%%%% Damping with Jacobi
% if rem(hand.iter+1,hand.check_f*hand.JacobiRat)==0
%     Tnew=T;
%     Tnew(hand.Cheq1.P)=(hand.NN_aV.w1/hand.w).*(...
%         hand.c_X*(T(hand.Cheq1.P_Xm)+T(hand.Cheq1.P_Xp))+...
%         hand.c_Y*(T(hand.Cheq1.P_Ym)+T(hand.Cheq1.P_Yp))+...
%         hand.c_Z*(T(hand.Cheq1.P_Zm)+T(hand.Cheq1.P_Zp)));
%
%     Tnew(hand.Cheq2.P)=(hand.NN_aV.w2/hand.w).*(...
%         hand.c_X*(T(hand.Cheq2.P_Xm)+T(hand.Cheq2.P_Xp))+...
%         hand.c_Y*(T(hand.Cheq2.P_Ym)+T(hand.Cheq2.P_Yp))+...
%         hand.c_Z*(T(hand.Cheq2.P_Zm)+T(hand.Cheq2.P_Zp)));
%     Error(checkNo)=max(max(max(T-Tnew)));
%     T=Tnew;
% else
%     Error=nan;
% end

% Update figure
if hand.InLineMode==0
    if checkNo>2
        if hand.impCheck==0
            if get(hand.Check_RVA, 'Value')==0
                figure(hand.ResFig);
                hand.axesHandles(6) = subplot(2,4,6);
                ha.Tau=plot(hand.check_f*(1:length(hand.DeffTop)), real(hand.DeffTop),hand.check_f*(1:length(hand.DeffTop)),real(hand.DeffBot));
                set(ha.Tau,'linewidth',1.5)
                leg=legend('Top','Base');
                set(leg,'Location','best','Interpreter','latex')
                title('D$^{\mathrm{eff}}$ Convergence','Interpreter','Latex');%$D^{\mathrm{eff}}$
                xlabel('Iterations','Interpreter','Latex');
                set(hand.axesHandles(6),'TickLabelInterpreter','latex');
                set(hand.axesHandles(6),'position',[0.33 0.16 0.157 0.25]);
                %                 ylabel('D$^{\mathrm{eff}}$','Interpreter','Latex');
                %                     set(hand.axesHandles(6),'YLim',[-0.01 1.01]);
                %                     set(hand.axesHandles(6),'YTick',[0 0.2 0.4 0.6 0.8 1]);
                %                     axis(hand.axesHandles(6),'square')
                drawnow
            end
        else
            %                 if rem(checkNo,5)==0
            %                     plot(hand.check_f*(1:length(hand.ImpedanceBot(1,:))), real(hand.ImpedanceBot(hand.freqNo,:)),hand.check_f*(1:length(hand.ImpedanceTop(1,:))),real(hand.ImpedanceTop(hand.freqNo,:)));
            %                     legend('Base','Top')
            %                     title('Real Imp','Interpreter','Latex');
            %                     xlabel('Iterations','Interpreter','Latex');
            %                     set(hand.axesHandles(6),'TickLabelInterpreter','latex');
            %                     set(hand.axesHandles(6),'position',[0.33 0.16 0.157 0.25])
            %                     %axis(hand.axesHandles(6),'square')
            %                     pause(0.001)
            %                 end
        end
    end
end
%%%% Switch to double precision if not converging
% if checkNo>10 && mean(class(T)=='single')==1 &&...
%         abs(sum(diff(hand.TauFacTop(end-5:end)-hand.TauFacBot(end-5:end))))<0.001
% %     disp('Switching to double')
%     T=double(T);
% end

function [hand]=InitiatePlot1(hand)
%% Open a new figure and plot the phase map
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
[aa,bb,cc]=size(hand.Net_Or);
hand.ResFig=1;
hand.ResFig=figure(...
    'Name',['TF_Results: ','p',hand.Pha(1),'d',num2str(hand.Dir),'_',hand.fil],...
    'Color',[1 1 1],...
    'renderer','painters',...
    'WindowStyle','normal',...
    'PaperPositionMode','auto',...
    'PaperOrientation','landscape',...
    'units','characters',...
    'Position',[30 4 200 48],...
    'PaperType','a4');
if get(hand.Check_RVA, 'Value')==1;
    switch get(hand.Pop_RV, 'Value')
        case 1
            titstr=[' - RVA: ',...
                num2str(round(100*hand.NetVol(end))),' \% (Cubic)'];
        case 2
            titstr=[' - RVA: ',...
                num2str(round(100*hand.NetVol(end))),' \% (L=const.)'];
        case 3
            titstr=[' - RVA: ',...
                num2str(round(100*hand.NetVol(end))),' \% (A=const. from top)'];
        case 4
            titstr=[' - RVA: ',...
                num2str(round(100*hand.NetVol(end))),' \% (A=const. from base)'];
    end
else
    titstr='';
end

annotation(hand.ResFig,'textbox',...
    [0.05 0.90 0.9 0.1],...
    'String',    ['\verb|',hand.filename,'|',titstr],...
    'LineStyle','none',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',16,...
    'FitBoxToText','off');
if get(hand.Check_VaryD,'value')==0
    CentralStr='$D^{\mathrm{eff}}=D\frac{\varepsilon}{\tau}$';
else
    CentralStr='$\displaystyle D^{\mathrm{eff}}=\left(\sum_{i\mathrm{=phase}} D_{i}\varepsilon_{i}\right)\tau^{-1}$';
end
annotation(hand.ResFig,'textbox',...
    [0.4 0.46 0.2 0.1],...
    'String',CentralStr,...
    'LineStyle','none',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',17,...
    'FitBoxToText','off');
if ismac
    set(hand.ResFig,'PaperUnits','normalized',...
        'PaperPosition', [0 0 1 1]);
end
% set(hand.ResFig, 'Visible', 'off');
% set(hand.ResFig,'Position',[20 20 1200 800]);

if a<b
    hand.wh=[0.155 0.235*a/b];
else
    hand.wh=[0.155*b/a 0.235];
end
% hand.axesHandles=1;
hand.axesHandles(1)=subplot(2,4,1);
if get(hand.Check_Reverse, 'Value')==1 && get(hand.Check_Impedance, 'Value')==1
    switch hand.Dir
        case 1
            imagesc(flipud(hand.Net_Or(hand.RVAdims(1,1):hand.RVAdims(1,2),hand.RVAdims(2,1):hand.RVAdims(2,2),hand.RVAdims(3,1))))
        case 2
            if cc~=1
                imagesc(flipud(reshape(hand.Net_Or(hand.RVAdims(1,1),hand.RVAdims(2,1):hand.RVAdims(2,2),hand.RVAdims(3,1):hand.RVAdims(3,2)),a,b)))
            else
                imagesc(rot90(hand.Net_Or(hand.RVAdims(2,1):hand.RVAdims(2,2),hand.RVAdims(1,1):hand.RVAdims(1,2))))
            end
        case 3
            imagesc(rot90(reshape(hand.Net_Or(hand.RVAdims(1,1):hand.RVAdims(1,2),hand.RVAdims(2,2),hand.RVAdims(3,1):hand.RVAdims(3,2)),b,a)))
    end
else
    switch hand.Dir
        case 1
            A=hand.Net_Or;
        case 2
            if cc>1
                A=permute(hand.Net_Or,[2 3 1]);
            else
                A=rot90(hand.Net_Or,3);
            end
        case 3
            A=flip(permute(hand.Net_Or,[3 1 2]),3);
    end
    A=A(hand.RVAdims(1,1):hand.RVAdims(1,2),hand.RVAdims(2,1):hand.RVAdims(2,2),hand.RVAdims(3,1));
    imagesc(A)
end
colormap(hand.myColorMap)
freezeColors;
title('Map','Interpreter','latex');%axis square
% set(hand.axesHandles(1),'position',[0.13 0.635 hand.wh])
axis(hand.axesHandles(1),'image')
set(hand.axesHandles(1), 'XTick', [],'YTick',[]);
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hand]=InitiatePlot2(hand,T,Map)
%% Plot the initial distribution and phase fraction plot
figure(hand.ResFig);
hand.axesHandles(2)=subplot(2,4,2);
hand.h2=imagesc(real(T(2:end-1,2:end-1)));
colormap(jet); freezeColors; colormap(hand.myColorMap);
title('Initialisation (slice)','Interpreter','Latex');
hold on; hand.h2=subimage(double(Map(2:end-1,2:end-1,2)),bone);%axis square;
set(hand.h2, 'AlphaData', 1-Map(2:end-1,2:end-1,2)); hold off
% set(hand.axesHandles(2),'position',[0.33 0.635 hand.wh])
axis(hand.axesHandles(2),'image')
GrayPore=uint8(1-cat(3,single(Map(2:end-1,2:end-1,2)),single(Map(2:end-1,2:end-1,2)),single(Map(2:end-1,2:end-1,2))));
hold on; hand.h2=subimage(128*GrayPore);set(hand.h2, 'AlphaData', 1-single(Map(2:end-1,2:end-1,2)));hold off
set(hand.axesHandles(2), 'XTick', [],'YTick',[]);

if hand.PercFlag==0 || hand.Blocked==1
    title('Real (slice)','Interpreter','Latex');
    hand.axesHandles(3)=subplot(2,4,3);
    hand.h3=imagesc(imag(T(2:end-1,2:end-1)));
    colormap(jet); freezeColors; colormap(hand.myColorMap);
    title('Imaginary (slice)','Interpreter','Latex');
    hold on; hand.h3=subimage(double(Map(2:end-1,2:end-1,2)),bone);%axis square;
    set(hand.h3, 'AlphaData', 1-Map(2:end-1,2:end-1,2)); hold off
    axis(hand.axesHandles(3),'image')
    hold on; hand.h3=subimage(128*GrayPore);
    set(hand.h3, 'AlphaData', 1-single(Map(2:end-1,2:end-1,2)));hold off
    set(hand.axesHandles(3), 'XTick', [],'YTick',[]);
end

[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
if get(hand.Check_VaryD,'value')==0
    hand.PlanePore=single(sum(sum(Map(2:end-1,:,:),2),3)/(b*c));
    hand.axesHandles(5) =subplot(2,4,5);ha.VF=plot(flipud(hand.PlanePore),[0:1:a-1],hand.VolFrac*ones(size([0:1:a])),[0:1:a]);
    legend(vertcat('Local','Perc.'),'Total')
    title('Volume Fraction of Each Plane','Interpreter','Latex');
    xlabel('Phase Volume Fraction','Interpreter','Latex');
    set( hand.axesHandles(5),'XLim',[0 1]);
    if roundsf(hand.VolFrac,2,'round')<0.9 &&...
            roundsf(hand.VolFrac,2,'round')>0.1
        set( hand.axesHandles(5),'XTick',[0 roundsf(hand.VolFrac,2,'round') 1]);
    end
    %     set(hand.axesHandles(5),'position',[0.13 0.16 0.155 0.235])
else
    hand.PlanePore=single(mean(mean(1./hand.Dmap(2:end-1,2:end-1,2:end-1),2),3));
    hand.axesHandles(5) =subplot(2,4,5);ha.VF=plot(flipud(hand.PlanePore),[0:1:a-1],hand.D*ones(size([0:1:a])),[0:1:a]);
    legend(vertcat('Planar','mean D'),vertcat('Global','mean D'))
    title('Diffusivity of Each Plane','Interpreter','Latex');
    xlabel('Mean Diffusivity','Interpreter','Latex');
    %     set( hand.axesHandles(5),'XLim',[0 1]);
end
axis(hand.axesHandles(5),'square');
set(legend,'Interpreter','latex')
set(legend,'Location','North East')
set(ha.VF,'LineWidth',1.5);
ylabel('Voxels from Base','Interpreter','Latex');
set( hand.axesHandles(5),'YLim',[0 inf]);
set(hand.axesHandles(5),'TickLabelInterpreter','latex');
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hand]=InitiatePlot3(hand)
%% Results plot, with tortuoisty iterations and flux map + flux by slice
figure(hand.ResFig);
[a,b,c]=size(hand.Tconv(2:end-1,2:end-1,2:end-1));
itera=[hand.check_f:hand.check_f:hand.iter];
hand.Tconv=single(hand.Tconv);
hand.axesHandles(3) =subplot(2,4,3); imagesc(hand.Tconv(2:a+1,2:b+1,2));
% set(hand.axesHandles(3),'position',[0.53 0.635 hand.wh])
axis(hand.axesHandles(3),'image')
colormap(jet); freezeColors;
title('Steady State (slice)','Interpreter','latex');
hold on; h=subimage(double(hand.Map(2:a+1,2:b+1,2)),bone);%axis square;
set(h, 'AlphaData', 1-single(hand.Map(2:a+1,2:b+1,2)));hold off
GrayPore=uint8(1-cat(3,single(hand.Map(2:end-1,2:end-1,2)),single(hand.Map(2:end-1,2:end-1,2)),single(hand.Map(2:end-1,2:end-1,2))));
hold on; h=subimage(128*GrayPore);
set(h, 'AlphaData', 1-single(hand.Map(2:end-1,2:end-1,2)));hold off
set(hand.axesHandles(3), 'XTick', [],'YTick',[]);

hand.axesHandles(6) =subplot(2,4,6);
ha.Tau=plot(itera, real(hand.DeffTop(1:length(itera))),itera, real(hand.DeffBot(1:length(itera))));
set(ha.Tau,'linewidth',1.5)
legend('Top','Base')
set(legend,'Interpreter','latex');
title('Convergence','Interpreter','latex');
xlabel('Iterations','Interpreter','latex');
ylabel('D$^{\mathrm{eff}}$','Interpreter','Latex');
% TauVal=double(real(hand.TauFacTop(length(itera))));
if get(hand.Check_VaryD,'value')==0
    set(hand.axesHandles(6),'YLim',[0 roundsf(hand.VolFrac,1,'ceil')]);
    %     set(hand.axesHandles(6),'YTick',[0 0.2 0.4 0.6 0.8 1]);
else
    set(hand.axesHandles(6),'YLim',[0 roundsf(hand.D,2,'ceil')]);
end
set(hand.axesHandles(6),'XLim',[0 roundsf(hand.iter,2,'ceil')]);
% set(hand.axesHandles(4),'position',[0.33 0.16 0.155 0.235])
axis(hand.axesHandles(6),'square')
% ax = gca;
% ax.YAxis.TickLabelFormat = '%,.1f';
[hand]=MBytesRecord(hand,whos,'InitPlot3'); %Memory
[Q_thro,Q_plane]=FluxHunter3ani(hand);
% PlaneFlux=sum(sum(Q_thro,2),3)/sum(sum(Q_thro(1,:,:)));%/(2*sum(sum(T(2,:,:))));
[hand]=MBytesRecord(hand,whos,'InitPlot3 afterFlux'); %Memory
if get(hand.Check_FluxMapSave,'Value')==1
    
    FluxFig=figure(...
        'Name',['TF_Flux Projections: ','p',hand.Pha(1),'d',num2str(hand.Dir),'_',hand.fil],...
        'Color',[1 1 1],...
        'PaperOrientation','landscape',...
        'renderer','painters',...
        'units','characters',...
        'position',[170 30 150 22]);
    annotation(FluxFig,'textbox',...
        [0.05 0.90 0.9 0.1],...
        'String',    ['\verb|',hand.filename,'|'],...
        'LineStyle','none',...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'FitBoxToText','off');
    subplot(1,3,1);imagesc(sum((Q_thro+Q_plane),3));
    colormap(hot);axis image;set(gca, 'XTick', [],'YTick',[]);
    title('Total Flux','Interpreter','latex')
    subplot(1,3,2);imagesc(sum((Q_thro),3));
    colormap(hot);axis image;set(gca, 'XTick', [],'YTick',[]);
    title('Through Flux','Interpreter','latex')
    subplot(1,3,3);imagesc(sum((Q_plane),3));
    colormap(hot);axis image;set(gca, 'XTick', [],'YTick',[]);
    title('Plane Flux','Interpreter','latex')
    
    OverNorm=350;
    TifSave3D(uint8((Q_thro+Q_plane)/max(max(max(Q_thro+Q_plane)))*OverNorm),hand,'Flux_Total')
    TifSave3D(uint8(Q_thro/max(Q_thro(:))*OverNorm),hand,'Flux_ThroPln')
    TifSave3D(uint8(Q_plane/max(Q_plane(:))*OverNorm),hand,'Flux_InPln')
    
    print(FluxFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir),'_Flux'],'-dpdf');
    pause(0.5)
    close(FluxFig)
end
figure(hand.ResFig);
hand.axesHandles(4) =subplot(2,4,4);
if abs(hand.Results.Tau-1)>1e-5
    imagesc(sum(Q_thro+Q_plane,3));
else
    imagesc(hand.Tconv(2:end-1,2:end-1,2)>0);
end
axis(hand.axesHandles(4),'image')
colormap(hot); freezeColors;
title('Flux Density (projection)','Interpreter','latex');
if c==1
    GrayPore=uint8(1-cat(3,single(hand.Map(2:a+1,2:b+1,2)),single(hand.Map(2:a+1,2:b+1,2)),single(hand.Map(2:a+1,2:b+1,2))));
    hold on; h=subimage(128*GrayPore);set(h, 'AlphaData', 1-single(hand.Map(2:a+1,2:b+1,round(1+c/2))));hold off
end
figure(hand.ResFig);
hand.axesHandles(7) =subplot(2,4,7);ha.Flux=plot(hand.XFlux,[1:1:a-1],mean(hand.XFlux)*ones(size([1:1:a-1])),[1:1:a-1]);
set(ha.Flux,'linewidth',1.5)
axis(hand.axesHandles(7),'square')
title('Vertical Flux','Interpreter','latex');
xlabel('Comparison to Result','Interpreter','latex');
ylabel('Voxels from Base','Interpreter','latex');
set(hand.axesHandles(7),'XLim',[0.9 1.1]);
set(hand.axesHandles(7),'XTick',[0.9 1 1.1]);
set(hand.axesHandles(7),'YLim',[0 inf]);
set(hand.axesHandles(1:4), 'XTick', [],'YTick',[]);
set(hand.axesHandles(5:7),'TickLabelInterpreter','latex');
% for i=1:length(hand.axesHandles)
%     if i==5
%         i=6;
%     end
%     set(hand.axesHandles(i),'TickLabelInterpreter','latex');
% end
figure(hand.ResFig);
[a b c]=size(hand.Net_Or);
hand.TauFacTop=hand.TauFacTop(end);
hand.TauFacBot=hand.TauFacBot(end);
hand.DeffTop=hand.DeffTop(end);
hand.DeffBot=hand.DeffBot(end);
if get(hand.Check_VaryD,'Value')==1
    if length(hand.NetVals)==3
        VaryD =  {['$D_\mathrm{phase}$ (m$^2$s$^{-1}$):'],...
            ['    * Black = ', get(hand.Edit_D_Black,'String')],...
            ['    * Green = ', get(hand.Edit_D_Green,'String')],...
            ['    * White = ', get(hand.Edit_D_White,'String')],...
            ['']};
    else
        VaryD =  {['$D_\mathrm{phase}$ (m$^2$s$^{-1}$):'],...
            ['    * Black = ', get(hand.Edit_D_Black,'String')],...
            ['    * White = ', get(hand.Edit_D_White,'String')],...
            ['']};
    end
    boxheight=0.51;
    DiffusionFraction=['Mean diffusivity (m$^2$s$^{-1}$) =  ',num2str(roundsf(hand.D,3,'round'))];
else
    DiffusionFraction=['Phase volume fraction  =  ',num2str(roundsf(100*hand.VolFrac,3,'round')),'\%'];
    VaryD={''};
    boxheight=0.45;
end

ResultAnnno=annotation(hand.ResFig,'textbox',...
    [0.72 0.05 0.27 boxheight],'String',{...
    ['Voxel volume: ',num2str(a),'$\times$',num2str(b),'$\times$',num2str(c)],...
    ['Voxel size (nm): ',get(hand.L1box,'String'),'$\times$',get(hand.L2box,'String'),'$\times$',get(hand.L3box,'String')],...
    ['Sample size ($\mu$m): ',...
    num2str(roundsf(a*str2double(get(hand.L1box,'String'))/1e3,3,'round')),'$\times$',...
    num2str(roundsf(b*str2double(get(hand.L2box,'String'))/1e3,3,'round')),'$\times$',...
    num2str(roundsf(c*str2double(get(hand.L3box,'String'))/1e3,3,'round'))],...
    ['Direction ',num2str(hand.Dir)],...
    ['Phase ',hand.Pha],...
    VaryD{:},...
    ['$D^\mathrm{eff}$ (m$^2$s$^{-1}$) =  ',num2str(roundsf(mean([hand.DeffBot,hand.DeffTop]),3,'round'))],...
    ['Tortuosity factor  =  ',num2str(roundsf(mean([hand.TauFacBot,hand.TauFacTop]),3,'round'))],'',...
    DiffusionFraction,...
    ['Directional percolation  =  ',num2str(roundsf(100*hand.VolFrac_Perc,3,'round')),'\%'],'',...
    ['Number of iterations  =  ',num2str(hand.iter)],...
    ['Simulation ',lower(hand.Results.SimTime)]},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'FontSize',11,...
    'Interpreter','Latex');
% sound(1,1000);pause(0.02);sound(1,1000);
set(hand.ResFig,'renderer','painters');
% hand.ImpedanceBot=1;
% hand.ImpedanceTop=1;
[hand.logo_A, hand.logo_map, hand.logo_alpha] = imread('TauFactor_icon2smo.png');
hand.axesHandles(8)=subplot(2,4,8);
set(hand.axesHandles(8),'position',[0.88 0.88 0.12 0.12])
uistack(hand.axesHandles(8),'down',6)
imshow(hand.logo_A, hand.logo_map);
% hand.axesHandles=1;
if get(hand.Check_pdfSave,'Value')==1
    if get(hand.Check_RVA, 'Value')==0
        print(hand.ResFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir)],'-dpdf');
    else
        print(hand.ResFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir),'_RVA',num2str(round(100*hand.NetVol(end)))],'-dpdf');
    end
end

function [hand]=InitiatePlot3imp(hand)
%% Alternative results plot for the impedance simulation
figure(hand.ResFig);
set(hand.ResFig,'name',['TF_y = ',num2str(hand.y(hand.freqNo))])
[a,b,c]=size(hand.Tconv(2:end-1,2:end-1,2:end-1));
itera=[hand.check_f:hand.check_f:hand.iter];
hand.axesHandles(3) =subplot(2,4,3); imagesc(real(hand.Tconv(2:a+1,2:b+1,2)));
% set(hand.axesHandles(3),'position',[0.53 0.635 hand.wh])
axis(hand.axesHandles(3),'image')
colormap(jet); freezeColors; title('Real','Interpreter','latex');
hold on; h=subimage(double(hand.Map(2:a+1,2:b+1,2)),bone);%axis square;
set(h, 'AlphaData', 1-single(hand.Map(2:a+1,2:b+1,2)));hold off
GrayPore=uint8(1-cat(3,single(hand.Map(2:end-1,2:end-1,2)),single(hand.Map(2:end-1,2:end-1,2)),single(hand.Map(2:end-1,2:end-1,2))));
hold on; h=subimage(128*GrayPore);
set(h, 'AlphaData', 1-single(hand.Map(2:end-1,2:end-1,2)));hold off

hand.axesHandles(8) =subplot(2,4,4); imagesc(imag(hand.Tconv(2:a+1,2:b+1,2)));
% set(hand.axesHandles(8),'position',[0.73 0.635 hand.wh])
axis(hand.axesHandles(8),'image')
colormap(jet); freezeColors; title('Imaginary','Interpreter','latex');
hold on; h=subimage(double(hand.Map(2:a+1,2:b+1,2)),bone);%axis square;
set(h, 'AlphaData', 1-single(hand.Map(2:a+1,2:b+1,2)));
h=subimage(128*GrayPore);
set(h, 'AlphaData', 1-single(hand.Map(2:end-1,2:end-1,2)));hold off

hand.axesHandles(6) =subplot(2,4,7); imagesc(abs(hand.Tconv(2:a+1,2:b+1,2)));
% set(hand.axesHandles(6),'position',[0.53 0.16 hand.wh])
axis(hand.axesHandles(6),'image')
colormap(jet); freezeColors; title('Magnitude','Interpreter','latex');
hold on; h=subimage(double(hand.Map(2:a+1,2:b+1,2)),bone);%axis square;
set(h, 'AlphaData', 1-single(hand.Map(2:a+1,2:b+1,2)));
h=subimage(128*GrayPore);
set(h, 'AlphaData', 1-single(hand.Map(2:end-1,2:end-1,2)));hold off

hand.axesHandles(7) =subplot(2,4,8); imagesc(angle(hand.Tconv(2:a+1,2:b+1,2)));
% set(hand.axesHandles(7),'position',[0.73 0.16 hand.wh])
axis(hand.axesHandles(7),'image')
colormap(jet); freezeColors; title('Argument','Interpreter','latex');
hold on; h=subimage(double(hand.Map(2:a+1,2:b+1,2)),bone);%axis square;
set(h, 'AlphaData', 1-single(hand.Map(2:a+1,2:b+1,2)));
h=subimage(128*GrayPore);
set(h, 'AlphaData', 1-single(hand.Map(2:end-1,2:end-1,2)));hold off

hand.axesHandles(5) =subplot(2,4,5);
plot(itera,...
    (real(hand.ImpedanceBot(hand.freqNo,1:length(itera)))-...
    real(hand.ImpedanceBot(hand.freqNo,length(itera))))/...
    real(hand.ImpedanceBot(hand.freqNo,length(itera))),...
    itera,...
    (imag(hand.ImpedanceBot(hand.freqNo,1:length(itera)))-...
    imag(hand.ImpedanceBot(hand.freqNo,length(itera))))/...
    imag(hand.ImpedanceBot(hand.freqNo,length(itera)))...
    );
axis(hand.axesHandles(5),'square');
ylim([-0.5 0.5])
title({'Bottom Surface','Relative Impedance'},'Interpreter','latex');
legend('Real','Imag.')
xlabel('Iterations','Interpreter','latex');

hand.axesHandles(4) =subplot(2,4,6);
plot(itera,...
    (real(hand.ImpedanceTop(hand.freqNo,1:length(itera)))-...
    real(hand.ImpedanceTop(hand.freqNo,length(itera))))/...
    real(hand.ImpedanceTop(hand.freqNo,length(itera))),...
    itera,...
    (imag(hand.ImpedanceTop(hand.freqNo,1:length(itera)))-...
    imag(hand.ImpedanceTop(hand.freqNo,length(itera))))/...
    imag(hand.ImpedanceTop(hand.freqNo,length(itera)))...
    );
axis(hand.axesHandles(4),'square');
ylim([-0.5 0.5])
title({'Top Surface','Relative Impedance'},'Interpreter','latex');
legend('Real','Imag.')
xlabel('Iterations','Interpreter','latex');
% ylim([roundsf(min(real(hand.ImpedanceBot(hand.freqNo,length(itera))),real(hand.ImpedanceTop(hand.freqNo,length(itera))))*0.9,2,'floor'),...
%       roundsf(max(real(hand.ImpedanceBot(hand.freqNo,length(itera))),real(hand.ImpedanceTop(hand.freqNo,length(itera))))*1.1,2,'ceil')]);

set(hand.axesHandles([1:3,6:8]), 'XTick', [],'YTick',[]);

for i=1:length(hand.axesHandles)
    set(hand.axesHandles(i),'TickLabelInterpreter','latex');
end
Metric=['Current  =  ',num2str(roundsf(mean([hand.ImpedanceBot(hand.freqNo,end),...
    hand.ImpedanceTop(hand.freqNo,end)]),3,'round'))];

set(hand.ResFig,'renderer','painters');
% hand.axesHandles=1;
if get(hand.Check_pdfSave,'Value')==1
    print(hand.ResFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),...
        '_p',hand.Pha(1),'d',num2str(hand.Dir),'_',num2str(hand.freqNo)],'-dpdf');
end

function [hand]=ImpPlot(hand)
%% Plot the impedance responce data
figure(hand.impFig)
[a b c]=size(hand.Tconv(2:end-1,2:end-1,2:end-1));
if hand.PercFlag==1 && hand.Blocked==0
    normFactor=hand.D*(b*hand.L_Y*c*hand.L_Z)/(a*hand.L_X);
else
    %     normFactor=(mean(mean(mean(hand.Map(end-1,2:end-1,2:end-1)))))^2/(hand.VolFrac);%/max(real(hand.ImpedanceBotConv)); %hand.VolFrac;
    %     hand.MAA=mean(hand.PlanePore(hand.PlanePore~=0));
    normFactor=3*hand.D*(hand.MAA*hand.L_Y*hand.L_Z)/(hand.L_X*hand.Max_Path);
end
if hand.freqNo==1
    disp(['Mean accessible pore area = ',num2str(round(hand.MAA*hand.L_Y*hand.L_Z*1e18)),' nm^2 (',num2str(round(hand.MAA/(b*c)*100)),'% or CV',')']);
    disp(['Maximum path length = ',num2str(hand.Max_Path*hand.L_X*1e9),' nm (',num2str(round((hand.Max_Path)/a*100)),'% of CV)']);
    annotation(hand.impFig,'textbox',...
        [0.05 0.90 0.9 0.1],...
        'String',    ['\verb|',hand.filename,'|'],...
        'LineStyle','none',...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'FitBoxToText','off');
end
plot(real(hand.ImpedanceBotConv*normFactor),-imag(hand.ImpedanceBotConv*normFactor),'-x','linewidth',1.5);
set(gca,'TickLabelInterpreter','latex')
axis square
if hand.PercFlag==1 && hand.Blocked==0
    hold on; plot(hand.Tau/hand.VolFrac,0,'r*','markersize',7);hold off
    lim=roundsf(max([...
        max( real(hand.ImpedanceBotConv))*normFactor,...
        max(-imag(hand.ImpedanceBotConv))*normFactor,...
        hand.Tau...
        ]),1,'ceil');
else
    lim=round(max(real(hand.ImpedanceBotConv*normFactor))*3);
end
% plot(real(hand.ImpedanceBotConv),-imag(hand.ImpedanceBotConv),'-x')
if hand.y(hand.freqNo)>=0
    hold on
    plot(real(hand.ImpedanceBotConv(hand.y==0)*normFactor),...
        -imag(hand.ImpedanceBotConv(hand.y==0)*normFactor),'or','markersize',7)
    hold off
    if hand.PercFlag==1 && hand.Blocked==0
        legend('Impedance','$\tau/\varepsilon$','$D/L_\mathrm{CV}^2$')
    else
        legend('Impedance','$D/L^2$')
    end
else
    legend('Impedance')
end
set(legend,'Interpreter','latex');
hold on; plot([0 max(lim,1)],[0,max(lim,1)],':k','linewidth',1.5);hold off
xlim([0 max(lim,1)]);
ylim([0 max(lim,1)]);
xlabel('\~{Z}$''$','Interpreter','Latex');
ylabel('-\~{Z}$''''$','Interpreter','Latex');
title('Normalised Impedance','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex',...
    'LineWidth',1.2,...
    'FontSize',16)

varname=[hand.fil,'.Imp','_',num2str(hand.Pha(1)),num2str(num2str(hand.Dir))];
if hand.Blocked==1 %&& hand.PercFlag==1
    varname=[varname,'_B'];
end
if get(hand.Check_Reverse, 'Value')==1;
    varname=[varname,'_R'];
end

if get(hand.Check_pdfSave,'Value')==1
    figure(hand.impFig);
    if hand.freqNo==1
        f1 = getframe(hand.impFig);
        hand.cccmap2=coldat2/255;
        [hand.impFigGif1,hand.cccmap2] = rgb2ind(f1.cdata,hand.cccmap2,'dither');
        f2 = getframe(hand.axesHandles(2));
        [hand.impFigGif2,hand.cccmap2] = rgb2ind(f2.cdata,hand.cccmap2,'dither');
        %         hand.im(1,1,1,length(hand.y)) = 0;
    else
        impFigGif1=hand.impFigGif1;
        impFigGif2=hand.impFigGif2;
    end
    f1 = getframe(hand.impFig);
    f2 = getframe(hand.axesHandles(2));
    hand.impFigGif1(:,:,1,hand.freqNo) = rgb2ind(f1.cdata,hand.cccmap2,'nodither');
    hand.impFigGif2(:,:,1,hand.freqNo) = rgb2ind(f2.cdata,hand.cccmap2,'nodither');
end
assignin(hand.WoSpace,'temp',real(hand.freqSet/hand.freqChar));
evalin(hand.WoSpace,[varname,'.Freq = temp'';']);
assignin(hand.WoSpace,'temp',hand.ImpedanceBotConv);
evalin(hand.WoSpace,[varname,'.Z = temp;']);
evalin(hand.WoSpace,'clear temp');
drawnow

function [Q_thro,Q_plane]=FluxHunter3ani(hand)
% hand.Map=single(hand.Map);
VoxDim=[hand.L_X, hand.L_Y, hand.L_Z];
DimFac(1)=VoxDim(2)*VoxDim(3)/VoxDim(1);
DimFac(2)=VoxDim(3)*VoxDim(1)/VoxDim(2);
DimFac(3)=VoxDim(1)*VoxDim(2)/VoxDim(3);
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
if get(hand.Check_VaryD,'value')==0
    if c~=1 && hand.MBytesA(end)*2<hand.RAM_availableMBs
        Q_thro=0.5*hand.Map(2:a+1,2:b+1,2:c+1)*DimFac(1).*...
            (hand.Map(1:a,  2:b+1,2:c+1).*abs(hand.Tconv(1:a,  2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +hand.Map(3:a+2,2:b+1,2:c+1).*abs(hand.Tconv(3:a+2,2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1)));
        Q_thro(1,:,:)=  Q_thro(1,:,:)+  DimFac(1)*hand.Map(2,2:end-1,2:end-1).*abs(hand.Tconv(2,2:end-1,2:end-1)-hand.TopStim);
        Q_thro(end,:,:)=Q_thro(end,:,:)+DimFac(1)*hand.Map(a+1,2:end-1,2:end-1).*abs(hand.BotStim-hand.Tconv(a+1,2:end-1,2:end-1));
        
        Q_plane=0.5*hand.Map(2:a+1,2:b+1,2:c+1).*...
            (DimFac(2)*hand.Map(2:a+1,1:b,  2:c+1).*abs(hand.Tconv(2:a+1,1:b,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +DimFac(2)*hand.Map(2:a+1,3:b+2,2:c+1).*abs(hand.Tconv(2:a+1,3:b+2,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +DimFac(3)*hand.Map(2:a+1,2:b+1,1:c  ).*abs(hand.Tconv(2:a+1,2:b+1,1:c)-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +DimFac(3)*hand.Map(2:a+1,2:b+1,3:c+2).*abs(hand.Tconv(2:a+1,2:b+1,3:c+2)-hand.Tconv(2:a+1,2:b+1,2:c+1)));
    else
        Q_thro=0.5*hand.Map(2:a+1,2:b+1,2).*...
            (DimFac(1)*hand.Map(1:a,2:b+1,2).*abs(  hand.Tconv(1:a,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2))...
            +DimFac(1)*hand.Map(3:a+2,2:b+1,2).*abs(hand.Tconv(3:a+2,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2)));
        Q_thro(1,:)  =Q_thro(1,:)+  DimFac(1)*hand.Map(2,2:end-1,2).*abs(hand.Tconv(2,2:end-1,2)-hand.TopStim);
        Q_thro(end,:)=Q_thro(end,:)+DimFac(1)*hand.Map(a+1,2:end-1,2).*abs(hand.Tconv(a+1,2:end-1,2)-hand.BotStim);
        
        Q_plane=0.5*hand.Map(2:a+1,2:b+1,2).*...
            (DimFac(2)*hand.Map(2:a+1,1:b,2).*abs(  hand.Tconv(2:a+1,1:b,2)-hand.Tconv(2:a+1,2:b+1,2))...
            +DimFac(2)*hand.Map(2:a+1,3:b+2,2).*abs(hand.Tconv(2:a+1,3:b+2,2)-hand.Tconv(2:a+1,2:b+1,2)));
    end
else
    Q_thro=2*hand.Map(2:a+1,2:b+1,2:c+1)*DimFac(1).*(...
        (hand.Map(1:a,  2:b+1,2:c+1).*abs(hand.Tconv(1:a,  2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(1:a,  2:b+1,2:c+1)+ hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (hand.Map(3:a+2,2:b+1,2:c+1).*abs(hand.Tconv(3:a+2,2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(3:a+2,2:b+1,2:c+1)+ hand.Dmap(2:a+1,2:b+1,2:c+1))));
    Q_thro([1 end],:,:)=2*hand.Map([2 end-1],2:b+1,2:c+1)*DimFac(1).*(...
        (hand.Map([2 end-2],2:b+1,2:c+1).*abs(hand.Tconv([1 end-2],  2:b+1,2:c+1)-hand.Tconv([2 end-1],2:b+1,2:c+1))./...
        (hand.Dmap([1 end-2],  2:b+1,2:c+1)+ hand.Dmap([2 end-1],2:b+1,2:c+1)))+...
        (hand.Map([3 end-1],2:b+1,2:c+1).*abs(hand.Tconv([3 end  ],2:b+1,2:c+1)-  hand.Tconv([2 end-1],2:b+1,2:c+1))./...
        (hand.Dmap([3 end  ],2:b+1,2:c+1)+   hand.Dmap([2 end-1],2:b+1,2:c+1))));
    
    Q_plane=2*hand.Map(2:a+1,2:b+1,2:c+1).*(...
        (DimFac(2)*hand.Map(2:a+1,1:b  ,2:c+1).*abs(hand.Tconv(2:a+1,1:b,2:c+1)-  hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,1:b,2:c+1)+   hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (DimFac(2)*hand.Map(2:a+1,3:b+2,2:c+1).*abs(hand.Tconv(2:a+1,3:b+2,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,3:b+2,2:c+1)+ hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (DimFac(3)*hand.Map(2:a+1,2:b+1,1:c  ).*abs(hand.Tconv(2:a+1,2:b+1,1:c)-  hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,2:b+1,1:c)+   hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (DimFac(3)*hand.Map(2:a+1,2:b+1,3:c+2).*abs(hand.Tconv(2:a+1,2:b+1,3:c+2)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,2:b+1,3:c+2)+ hand.Dmap(2:a+1,2:b+1,2:c+1))));
    %     else
    %         Q_thro=0.5*hand.Map(2:a+1,2:b+1,2).*...
    %             (DimFac(1)*hand.Map(1:a,2:b+1,2).*abs(  hand.Tconv(1:a,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2))...
    %             +DimFac(1)*hand.Map(3:a+2,2:b+1,2).*abs(hand.Tconv(3:a+2,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2)));
    %         Q_thro(1,:)  =Q_thro(1,:)+  DimFac(1)*hand.Map(2,2:end-1,2).*abs(hand.Tconv(2,2:end-1,2)-hand.TopStim);
    %         Q_thro(end,:)=Q_thro(end,:)+DimFac(1)*hand.Map(a+1,2:end-1,2).*abs(hand.Tconv(a+1,2:end-1,2)-hand.BotStim);
    %
    %         Q_plane=0.5*hand.Map(2:a+1,2:b+1,2).*...
    %             (DimFac(2)*hand.Map(2:a+1,1:b,2).*abs(  hand.Tconv(2:a+1,1:b,2)-hand.Tconv(2:a+1,2:b+1,2))...
    %             +DimFac(2)*hand.Map(2:a+1,3:b+2,2).*abs(hand.Tconv(2:a+1,3:b+2,2)-hand.Tconv(2:a+1,2:b+1,2)));
    %     end
end

function TifSave3D(Qa,hand,tag)
%% Function for saving tif stacks
[a b c]=size(Qa);
outputFileName = [hand.pathname,hand.filename(1:find(hand.filename=='.')-1),...
    '_p',hand.Pha(1),'d',num2str(hand.Dir),'_',tag,'.tif'];
Qa(Qa==0)=nan;
for K=1:length(Qa(1, 1, :))
    imwrite(Qa(:,:,K), outputFileName, 'WriteMode', 'append',  'Compression','none');
end

function TextB_Direction_ButtonDownFcn(hObject, eventdata, hand)
% %% Hiden function for revealing the impedance and reversal check boxes
% if ~isfield(hand,'Secrets')
%     hand.Secrets=1;
% else
%     hand.Secrets= hand.Secrets+1;
%     if  hand.Secrets>3
%         a=get(hand.Window,'Position');
%         %         set(hand.Window,'Position',[a(1),a(2),84,39]);
%         set(hand.Check_Impedance,'Value',1);
%         hand.impCheck=1;
%         Check_Impedance_Callback(hObject, eventdata, hand)
%         set(hand.Check_Reverse,'Visible','on')
%         set(hand.Check_Blocked,'Visible','on')
%     end
% end
% guidata(hObject, hand);


function Check_RVA_Callback(hObject, eventdata, hand)
%% RVA checkbox callback + time update
if get(hand.Check_RVA, 'Value')==1;
    set(hand.Pop_RV,'Enable','on')
else
    set(hand.Pop_RV,'Enable','off')
end
ExpectedTime(hObject, eventdata, hand);

function [hand]=RVA_Net_Cube(hand)
%% Calculated the idx for cubic RVA
[a,b,c]=size(hand.Net_full);
s=(1-hand.shrink)^(1/3);
d=(1-s)/2;
hand.RVAdims=[...
    floor(d*a+1),ceil(a*(1-d));...
    floor(d*b+1),ceil(b*(1-d));...
    floor(d*c+1),ceil(c*(1-d))];

function [hand]=RVA_Net_Lconst(hand)
%% Calculated the idx for Lconst RVA
[a,b,c]=size(hand.Net_full);
s=(1-hand.shrink)^(1/2);
d=(1-s)/2;
hand.RVAdims=[...
    1,a;...
    floor(d*b+1),ceil(b*(1-d));...
    floor(d*c+1),ceil(c*(1-d))];

function [hand]=RVA_Net_Aconst(hand)
%% Calculated the idx for Aconst RVA
[a,b,c]=size(hand.Net_full);
s=1-hand.shrink;
d=(1-s);
if get(hand.Pop_RV, 'Value')==3
    hand.RVAdims=[...
        1,ceil(a*(1-d));...
        1,b;...
        1,c];
else
    hand.RVAdims=[...
        floor(a*d+1),a;...
        1,b;...
        1,c];
end

function Window_ButtonDownFcn(hObject, eventdata, hand)
%% Generates a new figure window to show higher resolution preview of large dataset
if isfield(hand,'Net_Or')
    [a b c]=size(hand.Net_Or);
    if max(a,b)>50
        hand.zoomer=figure('name','TF_Phase Preview');
        imagesc(hand.Net_Or(:,:,1))
        axis equal tight
        if length(hand.NetVals)==2
            colormap(gray);
        else
            colormap(hand.myColorMap);
        end
    end
end

function uipanel_ButtonDownFcn(hObject, eventdata, hand)
%% Straight after the app is opened, clicking in the panel 4 times will generate a large random 3phase cube
if strcmp(hand.filename,'Dummy')
    hand.RandomSecrets= hand.RandomSecrets+1;
    if  hand.RandomSecrets==5
        set(hand.uipanel,'Interruptible','off');
        set(hand.LoadData,...
            'String','Loading...',...
            'ForegroundColor',[.9 .9 .9],...
            'Enable','off');
        drawnow
        hand.Net_Or=rand(300,300,300,'single');
        hand.Net_Or(hand.Net_Or<1/3)=1;
        hand.Net_Or(hand.Net_Or<2/3)=2;
        hand.Net_Or(hand.Net_Or<1)=0;
        hand.Net_Or=uint8(hand.Net_Or);
        hand.NetVals=[0 1 2];
        hand.filename='Random300.tif';
        hand.pathname  = [cd,'\'];
        hand.batch_length=1;
        hand.VFs(1)=mean(mean(mean(hand.Net_Or==0)));
        hand.VFs(2)=mean(mean(mean(hand.Net_Or==1)));
        hand.VFs(3)=mean(mean(mean(hand.Net_Or==2)));
        hand.MBytesA=0;
        hand.MemLoc={''};
        %         hand.VFanno(1)= annotation(gcf,'textbox',[0.82 0.84 0.2 0.1],'String',[num2str(100*hand.VFs(1),2),'%']);
        %         hand.VFanno(2)= annotation(gcf,'textbox',[0.82 0.81 0.2 0.1],'String',[num2str(100*hand.VFs(2),2),'%']);
        %         hand.VFanno(3)= annotation(gcf,'textbox',[0.82 0.78 0.2 0.1],'String',[num2str(100*hand.VFs(3),2),'%']);
        %         set(hand.VFanno(2),'Color',[0 1 0]);
        %         set(hand.VFanno(3),'Color',[0.5 0.5 0.5]);
        %         set(hand.VFanno,'LineStyle','none');
        %         hand.batch_current=1;
        [hand]=PlotData_Callback(hObject, eventdata, hand);
        [hand]=PrepWindow_Callback(hObject, eventdata, hand);
        set(hand.LoadData,...
            'String','Load Data',...
            'ForegroundColor',[1 1 1],...
            'Enable','on');
        
    end
end
guidata(hObject, hand);

function Pop_RV_CreateFcn(hObject, eventdata, hand)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Check_pdfSave_Callback(hObject, eventdata, hand)

function Check_FluxMapSave_Callback(hObject, eventdata, hand)

function Check_Reverse_Callback(hObject, eventdata, hand)

function axes1_ButtonDownFcn(hObject, eventdata, hand)

function Check_RVA_Lconst_Callback(hObject, eventdata, hand)

function checkbox14_Callback(hObject, eventdata, hand)

function Pop_RV_Callback(hObject, eventdata, hand)

% --- Executes on button press in Check_VaryD.
function Check_VaryD_Callback(hObject, eventdata, hand)
if get(hand.Check_VaryD,'Value')==1;
    set(hand.Edit_D_Black,'visible','on')
    if length(hand.NetVals)>2
        set(hand.Edit_D_Green,'visible','on')
    end
    set(hand.Edit_D_White,'visible','on')
    
    set(hand.Check_G1,'visible','off')
    set(hand.Check_G2,'visible','off')
    set(hand.Check_G3,'visible','off')
    set(hand.Check_W1,'visible','off')
    set(hand.Check_W2,'visible','off')
    set(hand.Check_W3,'visible','off')
    
    if    get(hand.Check_B1,'Value')==0 && get(hand.Check_B2,'Value')==0 && get(hand.Check_B3,'Value')==0
        set(hand.Check_B1,'Value',1)
    end
    
    set(hand.Check_Impedance, 'Value',0);
    set(hand.Check_Impedance, 'enable','off')
    set(hand.Check_Reverse, 'visible','off')
    set(hand.Check_Blocked, 'visible','off')
else
    set(hand.Edit_D_Black,'visible','off')
    set(hand.Edit_D_Green,'visible','off')
    set(hand.Edit_D_White,'visible','off')
    %     set(hand.Edit_D_Black,'string',1)
    %     set(hand.Edit_D_Green,'string',1)
    %     set(hand.Edit_D_White,'string',1)
    
    set(hand.Check_G1,'visible','on')
    set(hand.Check_G2,'visible','on')
    set(hand.Check_G3,'visible','on')
    set(hand.Check_W1,'visible','on')
    set(hand.Check_W2,'visible','on')
    set(hand.Check_W3,'visible','on')
    
    set(hand.Check_Impedance, 'enable','on')
end
[hand]=ExpectedTime(hObject, eventdata, hand);
guidata(hObject, hand);

function Edit_D_Black_Callback(hObject, eventdata, hand)
Str=get(hand.Edit_D_Black,'string');
if BosanquetCheck(Str)
    set(hand.Edit_D_Black,'backgroundcolor',[0 1 0])
elseif strcmp(Str,'Dseg')==1
    set(hand.Edit_D_Black,'backgroundcolor',[1 0 0])
    pause(0.05)
    if length(hand.NetVals>3)
        set(hand.Edit_D_Green,'string','Dseg')
    end
    set(hand.Edit_D_Black,'string',0)
    set(hand.Edit_D_Black,'backgroundcolor',[1 1 1])
elseif strcmp(Str,'Dmap')==1
    set(hand.Edit_D_Black,'backgroundcolor',[1 0 0])
    pause(0.05)
    set(hand.Edit_D_Green,'string','Dmap')
    set(hand.Edit_D_Black,'string',0)
    set(hand.Edit_D_Black,'backgroundcolor',[1 1 1])
elseif str2double(Str)<0 || length(Str)==0 || isnan(str2double(Str))
    set(hand.Edit_D_Black,'backgroundcolor',[1 0 0])
    pause(0.05)
    set(hand.Edit_D_Black,'string',0)
    set(hand.Edit_D_Black,'backgroundcolor',[1 1 1])
else
    set(hand.Edit_D_Black,'backgroundcolor',[1 1 1])
end

function Edit_D_Black_CreateFcn(hObject, eventdata, hand)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_D_Green_Callback(hObject, eventdata, hand)
Str=get(hand.Edit_D_Green,'string');
if BosanquetCheck(Str)
    set(hand.Edit_D_Green,'backgroundcolor',[0 1 0])
elseif strcmp(Str,'Dmap')==1
    set(hand.Edit_D_Green,'backgroundcolor',[0 1 0])
elseif strcmp(Str,'Dseg')==1
    set(hand.Edit_D_Green,'backgroundcolor',[0 1 0])
elseif str2double(Str)<0 || length(Str)==0 || isnan(str2double(Str))
    set(hand.Edit_D_Green,'backgroundcolor',[1 0 0])
    pause(0.05)
    set(hand.Edit_D_Green,'string',0)
    set(hand.Edit_D_Green,'backgroundcolor',[1 1 1])
else
    set(hand.Edit_D_Green,'backgroundcolor',[1 1 1])
end

function Edit_D_Green_CreateFcn(hObject, eventdata, hand)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Edit_D_White_Callback(hObject, eventdata, hand)
Str=get(hand.Edit_D_White,'string');
if BosanquetCheck(Str)
    set(hand.Edit_D_White,'backgroundcolor',[0 1 0])
elseif strcmp(Str,'Dseg')==1
    set(hand.Edit_D_White,'backgroundcolor',[1 0 0])
    pause(0.05)
    if length(hand.NetVals>3)
        set(hand.Edit_D_Green,'string','Dseg')
    end
    set(hand.Edit_D_White,'string',0)
    set(hand.Edit_D_White,'backgroundcolor',[1 1 1])
elseif strcmp(Str,'Dmap')==1
    set(hand.Edit_D_White,'backgroundcolor',[1 0 0])
    pause(0.05)
    set(hand.Edit_D_Green,'string','Dmap')
    set(hand.Edit_D_White,'string',0)
    set(hand.Edit_D_White,'backgroundcolor',[1 1 1])
elseif str2double(Str)<0 || length(Str)==0 || isnan(str2double(Str))
    set(hand.Edit_D_White,'backgroundcolor',[1 0 0])
    pause(0.05)
    set(hand.Edit_D_White,'string',0)
    set(hand.Edit_D_White,'backgroundcolor',[1 1 1])
else
    set(hand.Edit_D_White,'backgroundcolor',[1 1 1])
end

function Edit_D_White_CreateFcn(hObject, eventdata, hand)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [conf]=BosanquetCheck(Str)
conf=0;
if length(Str)>0
    idx=find(Str=='v');
    if length(idx)==1
        D=str2double(Str(1:idx-1));
        v=str2double(Str(idx+1:end));
        if D*v>0
            conf=1;
        end
    end
end

function Save_XL(hand)
[a, b, c]=size(hand.Net_Or);
XLtemp{1,1}='Filename';XLtemp{1,2}=hand.filename;
XLtemp{2,1}='Number of voxels in d1';XLtemp{2,2}=a;
XLtemp{3,1}='Number of voxels in d2';XLtemp{3,2}=b;
XLtemp{4,1}='Number of voxels in d3';XLtemp{4,2}=c;
XLtemp{5,1}='Voxel Size in d1 (nm)';XLtemp{5,2}=get(hand.L1box,'String');
XLtemp{6,1}='Voxel Size in d2 (nm)';XLtemp{6,2}=get(hand.L2box,'String');
XLtemp{7,1}='Voxel Size in d3 (nm)';XLtemp{7,2}=get(hand.L3box,'String');
XLtemp{8,1}='Black VF';XLtemp{8,2}=hand.VFs(1);
XLtemp{9,1}='Green VF';XLtemp{9,2}=hand.VFs(2);
XLtemp{10,1}='White VF';XLtemp{10,2}=hand.VFs(3);
XLtemp{11,1}='B1 Tau';XLtemp{11,2}=hand.TauSet(1,1);
XLtemp{12,1}='B2 Tau';XLtemp{12,2}=hand.TauSet(1,2);
XLtemp{13,1}='B3 Tau';XLtemp{13,2}=hand.TauSet(1,3);
XLtemp{14,1}='G1 Tau';XLtemp{14,2}=hand.TauSet(2,1);
XLtemp{15,1}='G2 Tau';XLtemp{15,2}=hand.TauSet(2,2);
XLtemp{16,1}='G3 Tau';XLtemp{16,2}=hand.TauSet(2,3);
XLtemp{17,1}='W1 Tau';XLtemp{17,2}=hand.TauSet(3,1);
XLtemp{18,1}='W2 Tau';XLtemp{18,2}=hand.TauSet(3,2);
XLtemp{19,1}='W3 Tau';XLtemp{19,2}=hand.TauSet(3,3);
if length(hand.NetVals)<3
    XLtemp{10,2}=hand.VFs(2);
    XLtemp([9,14:16],:)=[];
end
xlswrite([hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_Tau'],XLtemp)

function [hand]=MBytesRecord(hand,CurrentWoSpa,PlaceName)
hand.MBytesA=[hand.MBytesA sum([CurrentWoSpa.bytes])/1024^2]; %Memory
if nargin<3
    PlaceName='';
end
hand.MemLoc{end+1}=PlaceName;


% --- Executes on button press in Check_mex.
function Check_mex_Callback(hObject, eventdata, hand)
if get(hand.Check_mex,'value')==1
    hand.mexCtrl='mex';
else
    hand.mexCtrl='mat';
end
ExpectedTime(hObject, eventdata, hand);
guidata(hObject, hand);

function [tocStr]=TimeString(tocking)
if tocking<3/86400
    tocStr = ['Time = 3 sec'];
elseif tocking<100/86400
    tocStr = ['Time = ',num2str(round(tocking*86400)),' sec'];
elseif tocking<100*60/86400
    tocStr = ['Time = ',num2str(round(tocking*86400/60)),' min'];
elseif tocking<36*3600/86400
    tocStr = ['Time = ',num2str(round(tocking*86400/3600*2)/2),' hrs'];
else
    tocStr = ['Time = ',num2str(ceil(tocking*2)/2),' days'];
end
