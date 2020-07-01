%function for analysis of the degere to which an analysis volume is
%representative of the bulk through uantification of the volume fractions,
%surface areas and TPB lengths contained in subsampled regions.
function [hand]=MetricsRVA(VoxHeight_nm,VoxWidth_nm,VoxDepth_nm,RVA,hand); % Ignor hand if running independently.
MBytes=whos;MBytesA(1)=sum([MBytes.bytes])/1024^2;
if nargin==0
    VoxHeight_nm=2e3;
    VoxWidth_nm=1e3;
    VoxDepth_nm=1e3;
end
VoxDim=[VoxHeight_nm VoxWidth_nm VoxDepth_nm]*1e-9;
CalledFromTauFactor=1;
%% File imported
if nargin<5
    CalledFromTauFactor=0;
    hand.InLineMode=0;
    hand.WoSpace='base';
    hand.Dir=1;
    RVArepeats=0;
    hand.Check_pdfSave=uicontrol('Parent',[],'Value',1);
    if ~exist('RVA','var')
        RVA=0;
        Mode=0;
    else
        Mode=RVA;
    end
    [hand.filename,hand.pathname]=uigetfile('*.*','Select .tif file.');
    if hand.pathname==0
        hand.pathname  = [cd,'\'];
        hand.filename  = 'OutputResults.tif';
        %     Net_Or=zeros([5 10 3]);Net_Or(1:3,1:6,:)=1;Net_Or(1:3,7:end,:)=2;
        Net_Or=rand(300, 300, 300,'single');
        Net_Or(Net_Or<1/3)=1;
        Net_Or(Net_Or<2/3)=2;
        Net_Or(Net_Or<1)=0;
    else
        imageInfo = imfinfo([hand.pathname,hand.filename]);
        a=imageInfo.Height;
        b=imageInfo.Width;
        c=numel(imageInfo);
        Net_Or = zeros(a,b,c);      % Preallocate the cell array
        if c>1;
            for i = 1:c
                Net_Or(:,:,i) = single(imread([hand.pathname,hand.filename],'Index',i,'Info',imageInfo));
            end
        else
            Net_Or=single(imread([hand.pathname,hand.filename]));
        end
        Net_Or=uint8(Net_Or);
    end
    if length(hand.filename)>53
        hand.fil=[hand.filename(1:50)];
    elseif length(hand.filename)>4
        hand.fil=hand.filename(1:find(hand.filename=='.')-1);
    else
        hand.fil=['OutputResults'];
    end
    hand.fil(~isstrprop(hand.fil,'alphanum'))='_';
    if ~isstrprop(hand.fil(1), 'alpha')
        hand.fil=['X',hand.fil];
    end
    NetVals=unique(Net_Or);
    if length(NetVals)==2
        Net_Or(Net_Or==NetVals(1))=0;
        Net_Or(Net_Or==NetVals(2))=1;
        Net_Or=logical(Net_Or);
    elseif length(NetVals)==3
        Net_Or(Net_Or==NetVals(1))=0;
        Net_Or(Net_Or==NetVals(2))=1;
        Net_Or(Net_Or==NetVals(3))=2;
        Net_Or=uint8(Net_Or);
    else
        Net_Or(Net_Or==NetVals(1))=0;
        Net_Or(Net_Or==NetVals(2))=1;
        Net_Or(Net_Or>=NetVals(3))=2;
        Net_Or=uint8(Net_Or);
        disp('Volume has more than 3 phases');
    end
    clear imageInfo
elseif nargin==5
    NetVals=hand.NetVals;
    RVArepeats=hand.RVArepeats;
    [a b c]=size(hand.Net_Or);
    if hand.Dir==2
        if c>1
            Net_Or=permute(hand.Net_Or,[2 3 1]);
        else
            Net_Or=rot90(hand.Net_Or,3);
        end
    elseif hand.Dir==3
        Net_Or=permute(hand.Net_Or,[3 1 2]);
    else
        Net_Or=hand.Net_Or;
    end
end
xxx=tic;
% Net_Or=uint8(Net);
% Net_Or=padarray(Net_Or,[1 1 1],0);
[a,b,c]=size(Net_Or);

% Print expected calculation time to CW.
if a*b*c>1e7
    TimeCoef=0.45e-6;
    if RVA==0
        TimeCoef=TimeCoef/2;
    end
    disp(['Expected time c. ',num2str(ceil(TimeCoef*a*b*c/10)*10),' s.']);
end

NetP.Net=Net_Or;

MBytes=whos;MBytesA(2)=sum([MBytes.bytes])/1024^2;

%% Calculate total interfacial area of each phase
% Counts how many nearest neighbours are not of the same phase to the node
% Each interface is shared by 2 voxels, so the value is divided by 2 to
% avoid double counting.
if c>1 %3D interfacial areas
    NetP.P1=NetP.Net==0;
    Area.P1=0.5*(...
        VoxDim(2)*VoxDim(3)*(...
        (NetP.P1(2:end-1,2:end-1,2:end-1)~=NetP.P1(1:end-2,2:end-1,2:end-1))+...
        (NetP.P1(2:end-1,2:end-1,2:end-1)~=NetP.P1(3:end  ,2:end-1,2:end-1)))+...
        VoxDim(3)*VoxDim(1)*(...
        (NetP.P1(2:end-1,2:end-1,2:end-1)~=NetP.P1(2:end-1,1:end-2,2:end-1))+...
        (NetP.P1(2:end-1,2:end-1,2:end-1)~=NetP.P1(2:end-1,3:end  ,2:end-1)))+...
        VoxDim(1)*VoxDim(2)*(...
        (NetP.P1(2:end-1,2:end-1,2:end-1)~=NetP.P1(2:end-1,2:end-1,1:end-2))+...
        (NetP.P1(2:end-1,2:end-1,2:end-1)~=NetP.P1(2:end-1,2:end-1,3:end  ))));
    NetP.P1=0;
    
    if length(NetVals)>2
        
        NetP.P2=NetP.Net==1;
        Area.P2=0.5*(...
            VoxDim(2)*VoxDim(3)*(...
            (NetP.P2(2:end-1,2:end-1,2:end-1)~=NetP.P2(1:end-2,2:end-1,2:end-1))+...
            (NetP.P2(2:end-1,2:end-1,2:end-1)~=NetP.P2(3:end  ,2:end-1,2:end-1)))+...
            VoxDim(3)*VoxDim(1)*(...
            (NetP.P2(2:end-1,2:end-1,2:end-1)~=NetP.P2(2:end-1,1:end-2,2:end-1))+...
            (NetP.P2(2:end-1,2:end-1,2:end-1)~=NetP.P2(2:end-1,3:end  ,2:end-1)))+...
            VoxDim(1)*VoxDim(2)*(...
            (NetP.P2(2:end-1,2:end-1,2:end-1)~=NetP.P2(2:end-1,2:end-1,1:end-2))+...
            (NetP.P2(2:end-1,2:end-1,2:end-1)~=NetP.P2(2:end-1,2:end-1,3:end  ))));
        NetP.P2=0;
        
        NetP.P3=NetP.Net==2;
        Area.P3=0.5*(...
            VoxDim(2)*VoxDim(3)*(...
            (NetP.P3(2:end-1,2:end-1,2:end-1)~=NetP.P3(1:end-2,2:end-1,2:end-1))+...
            (NetP.P3(2:end-1,2:end-1,2:end-1)~=NetP.P3(3:end  ,2:end-1,2:end-1)))+...
            VoxDim(3)*VoxDim(1)*(...
            (NetP.P3(2:end-1,2:end-1,2:end-1)~=NetP.P3(2:end-1,1:end-2,2:end-1))+...
            (NetP.P3(2:end-1,2:end-1,2:end-1)~=NetP.P3(2:end-1,3:end  ,2:end-1)))+...
            VoxDim(1)*VoxDim(2)*(...
            (NetP.P3(2:end-1,2:end-1,2:end-1)~=NetP.P3(2:end-1,2:end-1,1:end-2))+...
            (NetP.P3(2:end-1,2:end-1,2:end-1)~=NetP.P3(2:end-1,2:end-1,3:end  ))));
        NetP.P3=0;
    end
else %2D interfacial areas
    NetP.P1=NetP.Net==0;
    Area.P1=0.5*(...
        VoxDim(2)*VoxDim(3)*(...
        NetP.P1(2:end-1,2:end-1)~=NetP.P1(1:end-2,2:end-1)+...
        NetP.P1(2:end-1,2:end-1)~=NetP.P1(3:end  ,2:end-1))+...
        VoxDim(3)*VoxDim(1)*(...
        NetP.P1(2:end-1,2:end-1)~=NetP.P1(2:end-1,1:end-2)+...
        NetP.P1(2:end-1,2:end-1)~=NetP.P1(2:end-1,3:end  )));
    NetP.P1=0;
    
    if length(NetVals)>2
        
        NetP.P2=NetP.Net==1;
        Area.P2=0.5*(...
            VoxDim(2)*VoxDim(3)*(...
            NetP.P2(2:end-1,2:end-1)~=NetP.P2(1:end-2,2:end-1)+...
            NetP.P2(2:end-1,2:end-1)~=NetP.P2(3:end  ,2:end-1))+...
            VoxDim(3)*VoxDim(1)*(...
            NetP.P2(2:end-1,2:end-1)~=NetP.P2(2:end-1,1:end-2)+...
            NetP.P2(2:end-1,2:end-1)~=NetP.P2(2:end-1,3:end  )));
        NetP.P2=0;
        
        NetP.P3=NetP.Net==2;
        Area.P3=0.5*(...
            VoxDim(2)*VoxDim(3)*(...
            NetP.P3(2:end-1,2:end-1)~=NetP.P3(1:end-2,2:end-1)+...
            NetP.P3(2:end-1,2:end-1)~=NetP.P3(3:end  ,2:end-1))+...
            VoxDim(3)*VoxDim(1)*(...
            NetP.P3(2:end-1,2:end-1)~=NetP.P3(2:end-1,1:end-2)+...
            NetP.P3(2:end-1,2:end-1)~=NetP.P3(2:end-1,3:end  )));
        NetP.P3=0;
    end
end
MBytes=whos;MBytesA(3)=sum([MBytes.bytes])/1024^2;
NetP=NetP.Net;

%% Calculate TPB density map if 3 phase data
% Number of edges for cubic volume of dimension n is 3*n^3+6*n^2+3*n
% This approach effectively check all 12 edges of each voxel and then, as
% each edge a shared by 4 voxels, this number is divided by 4 to avoid
% quadruple counting
if length(NetVals)>2
    NetP(NetP==1)=5;
    NetP(NetP==0)=3;
    if c>1 %3D networks
        TPBmap.P1=rem(... %Map of TPBs in dir 1
            NetP(2:end-1,1:end-1,1:end-1).*...
            NetP(2:end-1,2:end  ,1:end-1).*...
            NetP(2:end-1,1:end-1,2:end  ).*...
            NetP(2:end-1,2:end  ,2:end  )...
            ,30)==0;
        TPBmap.P2=rem(... %Map of TPBs in dir 2
            NetP(1:end-1,2:end-1,1:end-1).*...
            NetP(2:end  ,2:end-1,1:end-1).*...
            NetP(1:end-1,2:end-1,2:end  ).*...
            NetP(2:end  ,2:end-1,2:end  )...
            ,30)==0;
        TPBmap.P3=rem(... %Map of TPBs in dir 3
            NetP(1:end-1,1:end-1,2:end-1).*...
            NetP(2:end  ,1:end-1,2:end-1).*...
            NetP(1:end-1,2:end  ,2:end-1).*...
            NetP(2:end  ,2:end  ,2:end-1)...
            ,30)==0;
        
        TPBmap=single(0.25*(... %Map of TPBs in all dirs * 0.25
            VoxDim(1)*(...
            TPBmap.P1(:,1:end-1,1:end-1)+...
            TPBmap.P1(:,1:end-1,2:end  )+...
            TPBmap.P1(:,2:end  ,1:end-1)+...
            TPBmap.P1(:,2:end  ,2:end  ))+...
            VoxDim(2)*(...
            TPBmap.P2(1:end-1,:,1:end-1)+...
            TPBmap.P2(1:end-1,:,2:end  )+...
            TPBmap.P2(2:end  ,:,1:end-1)+...
            TPBmap.P2(2:end  ,:,2:end  ))+...
            VoxDim(3)*(...
            TPBmap.P3(1:end-1,1:end-1,:)+...
            TPBmap.P3(1:end-1,2:end  ,:)+...
            TPBmap.P3(2:end  ,1:end-1,:)+...
            TPBmap.P3(2:end  ,2:end  ,:))));
    else %2D networks
        TPBmap.P3=rem(... %Map of TPBs in dir 3
            NetP(1:end-1,1:end-1).*...
            NetP(2:end  ,1:end-1).*...
            NetP(1:end-1,2:end  ).*...
            NetP(2:end  ,2:end  )...
            ,30)==0;
        
        TPBmap=single(0.25*(...
            VoxDim(3)*(...
            TPBmap.P3(1:end-1,1:end-1,:)+...
            TPBmap.P3(1:end-1,2:end  ,:)+...
            TPBmap.P3(2:end  ,1:end-1,:)+...
            TPBmap.P3(2:end  ,2:end  ,:))));
    end
    MBytes=whos;MBytesA(4)=sum([MBytes.bytes])/1024^2;
end
NetP=1;
MBytes=whos;MBytesA(5)=sum([MBytes.bytes])/1024^2;
%% Setup RVA iterations
if nargin==5
    Mode=get(hand.Pop_RV, 'Value');
else
    Mode=RVA;
end
Dims=[a b c];
if RVA==0
    %% Calculate total voxels in mask region
    Mask_VoxSum=(a-2)*(b-2)*max(1,(c-2));
    Mask_RelVol=1;
    Mask_Vol=Mask_VoxSum*VoxDim(1)*VoxDim(2)*VoxDim(3);
    %% Calculate volume fractions in mask region
    VolFrac(1)=sum(sum(sum(Net_Or==0)))/numel(Net_Or);
    VolFrac(2)=sum(sum(sum(Net_Or==1)))/numel(Net_Or);
    if length(NetVals)>2
        VolFrac(3)=sum(sum(sum(Net_Or==2)))/numel(Net_Or);
    end
    %% Calculate surface areas in mask region
    SA(1)=sum(Area.P1(:))/numel(Area.P1)/prod(VoxDim);
    if length(NetVals)>2
        SA(2)=sum(Area.P2(:))/numel(Area.P2)/prod(VoxDim);
        SA(3)=sum(Area.P3(:))/numel(Area.P3)/prod(VoxDim);
    end
    %% Calculate TPB density in mask region
    if length(NetVals)>2
        TPB_RVA=sum(double(TPBmap(:)))/numel(TPBmap)/prod(VoxDim);
    end
else
    CheckNo=20;
    Mask_RelVol=zeros(CheckNo,1);
    Mask_Vol=zeros(CheckNo,1);
    if length(NetVals)>2
        VolFrac=zeros(CheckNo,3);
        SA=zeros(CheckNo,3);
        TPB_RVA=zeros(CheckNo,1);
    else
        VolFrac=zeros(CheckNo,2);
        SA=zeros(CheckNo,1);
    end
    if c~=1
        Net=Net_Or(2:end-1,2:end-1,2:end-1);
    else
        Net=Net_Or(2:end-1,2:end-1,1);
    end
    for i=1:CheckNo
        shrink=(i-2)/(CheckNo-1);
        %% Get 3D matrix indicies for RVA
        if i>1
            switch Mode
                case 1
                    [SD]=       RVA_Net(Dims-2,shrink);
                case 2
                    [SD]=RVA_Net_Lconst(Dims-2,shrink);
                case 3
                    [SD]=RVA_Net_Aconst(Dims-2,shrink,Mode);
                case 4
                    [SD]=RVA_Net_Aconst(Dims-2,shrink,Mode);
            end
        else
            SD=[1,a;1,b;1,c];
        end
        %% Calculate total voxels in mask region
        Mask_VoxSum=...
            (SD(1,2)-SD(1,1)+1)*...
            (SD(2,2)-SD(2,1)+1)*...
            (SD(3,2)-SD(3,1)+1);
        Mask_RelVol(i)=Mask_VoxSum/(a*b*c);
        Mask_Vol(i)=Mask_VoxSum*prod(VoxDim);
        
        %% Calculate volume fractions in mask region
        if i==1
            VolFrac(i,1)=sum(Net_Or(:)==0)/Mask_VoxSum;
            VolFrac(i,2)=sum(Net_Or(:)==1)/Mask_VoxSum;
            if length(NetVals)>2
                VolFrac(i,3)=sum(Net_Or(:)==2)/Mask_VoxSum;
            end
        else
            VolFrac(i,1)=sum(sum(sum(Net(SD(1,1):SD(1,2),SD(2,1):SD(2,2),SD(3,1):SD(3,2))==0)))/Mask_VoxSum;
            VolFrac(i,2)=sum(sum(sum(Net(SD(1,1):SD(1,2),SD(2,1):SD(2,2),SD(3,1):SD(3,2))==1)))/Mask_VoxSum;
            if length(NetVals)>2
                VolFrac(i,3)=sum(sum(sum(Net(SD(1,1):SD(1,2),SD(2,1):SD(2,2),SD(3,1):SD(3,2))==2)))/Mask_VoxSum;
            end
        end
        %% Calculate surface areas in mask region
        if i==1
            SA(i,1)=nan;
            if length(NetVals)>2
                SA(i,2)=nan;
                SA(i,3)=nan;
            end
        else
            SA(i,1)=sum(sum(sum(Area.P1(SD(1,1):SD(1,2),SD(2,1):SD(2,2),SD(3,1):SD(3,2)))))/Mask_Vol(i);
            if length(NetVals)>2
                SA(i,2)=sum(sum(sum(Area.P2(SD(1,1):SD(1,2),SD(2,1):SD(2,2),SD(3,1):SD(3,2)))))/Mask_Vol(i);
                SA(i,3)=sum(sum(sum(Area.P3(SD(1,1):SD(1,2),SD(2,1):SD(2,2),SD(3,1):SD(3,2)))))/Mask_Vol(i);
            end
        end
        
        %% Calculate TPB density in mask region
        if length(NetVals)>2
            if i==1
                TPB_RVA(i)=nan;
            else
                TPB_RVA(i)=sum(sum(sum(double(TPBmap(SD(1,1):SD(1,2),SD(2,1):SD(2,2),SD(3,1):SD(3,2))))))/Mask_Vol(i);
            end
        end
    end
end
MBytes=whos;MBytesA(6)=sum([MBytes.bytes])/1024^2;
%% Send results to workspace
if hand.InLineMode==0
    varname=[hand.fil,'.Metrics'];
    assignin(hand.WoSpace,'temp',Mask_RelVol(1)');
    evalin(hand.WoSpace,[varname,'.RVA_Steps = temp;']);
    assignin(hand.WoSpace,'temp',VolFrac(1,:));
    evalin(hand.WoSpace,[varname,'.VolFracs = temp'';']);
    if RVA==0
        assignin(hand.WoSpace,'temp',SA*1e-6);
    else
        assignin(hand.WoSpace,'temp',SA(2,:)*1e-6);
    end
    evalin(hand.WoSpace,[varname,'.SurfAreaDens_over_um = temp'';']);
    if length(NetVals)>2
        if RVA==0
            assignin(hand.WoSpace,'temp',TPB_RVA'*1e-12);
        else
            assignin(hand.WoSpace,'temp',TPB_RVA(2)'*1e-12);
        end
        evalin(hand.WoSpace,[varname,'.TPBdens_over_um2 = temp;']);
    end
    if RVA~=0
        if Mode>1
            varname=[hand.fil,'.MetricsRVA',num2str(Mode),'_dir',num2str(num2str(hand.Dir))];
        else
            varname=[hand.fil,'.MetricsRVA',num2str(Mode)];
        end
        assignin(hand.WoSpace,'temp',Mask_RelVol');
        evalin(hand.WoSpace,[varname,'.RVA_Steps = temp;']);
        assignin(hand.WoSpace,'temp',VolFrac);
        evalin(hand.WoSpace,[varname,'.VolFracs = temp'';']);
        assignin(hand.WoSpace,'temp',SA*1e-6);
        evalin(hand.WoSpace,[varname,'.SurfAreaDens_over_um = temp'';']);
        if length(NetVals)>2
            assignin(hand.WoSpace,'temp',TPB_RVA'*1e-12);
            evalin(hand.WoSpace,[varname,'.TPBdens_over_um2 = temp;']);
        end
    end
    
    
    evalin(hand.WoSpace,'clear temp');
else
    if RVA~=0
        if Mode>1
            varname=['OutputResults.MetricsRVA',num2str(Mode),'_dir',num2str(num2str(hand.Dir))];
        else
            varname=['OutputResults.MetricsRVA',num2str(Mode)];
        end
        eval(['hand.',varname,'.RVA_Steps = Mask_RelVol;']);
    else
        varname='OutputResults.Metrics';
    end
    eval(['hand.',varname,'.VolFracs = VolFrac;']);
    eval(['hand.',varname,'.SurfAreaDens_over_um = SA*1e-6;']);
    if length(NetVals)>2
        eval(['hand.',varname,'.TPBdens_over_um2 = TPB_RVA''*1e-12;']);
    end
end

%% Plot results
if hand.InLineMode==0
    if RVA~=0
        % Setup figure window
        RVfig(Mode)=figure(...
            'Name',['TF_Metrics: ',hand.fil],...
            'Color',[1 1 1],...
            'renderer','painters',...
            'WindowStyle','normal',...
            'PaperPositionMode','auto',...
            'PaperOrientation','landscape',...
            'units','characters',...
            'position',[14 4 130 53]);
        if ismac
            set(RVfig(Mode),'PaperUnits','normalized',...
                'PaperPosition', [0 0 1 1]);
        end
        % Volume fraction plot
        h(1)=subplot(2,2,1);ha.V=plot(Mask_RelVol,VolFrac');%,'ColorOrder',[0 0 0;0 0 0;0 1 0]
        ylim([0 1])
        xlim([0 1])
        xlabel('Fraction of Total Volume','Interpreter','Latex');
        ylabel('Phase Volume Fraction','Interpreter','Latex');
        if length(NetVals)<3
            legend('$V_\mathrm{Black}$','$V_\mathrm{White}$');
            set(ha.V(2),'Color',[0 0 0])
            set(ha.V(2),'LineStyle',':')
        else
            legend('$V_\mathrm{Black}$','$V_\mathrm{Green}$','$V_\mathrm{White}$');
            set(ha.V(3),'Color',[0 0 0])
            set(ha.V(3),'LineStyle',':')
            set(ha.V(2),'Color',[0 1 0])
        end
        set(ha.V(1),'Color',[0 0 0])
        
        set(ha.V,'linewidth',1.5)
        set(legend,'Interpreter','latex')
        
        % Surface area plot
        SA_norm=SA;
        SA_norm(:,1)=SA(:,1)/SA(2,1);
        if length(NetVals)>2
            SA_norm(:,2)=SA(:,2)/SA(2,2);
            SA_norm(:,3)=SA(:,3)/SA(2,3);
        end
        
        if length(NetVals)>2
            h(2)=subplot(2,2,2);ha.A=plot(Mask_RelVol,SA_norm'-1);%,Mask_RelVol,SA_norm2'-1);
            ylabel('$\Delta$(Surface Area per Vol.)','Interpreter','Latex');
            set(ha.A(2),'Color',[0 1 0])
            set(ha.A(3),'Color',[0 0 0])
            set(ha.A(3),'LineStyle',':')
            ylim([-0.1 0.1]);
            set(h(2), 'yticklabel', {'-10\%','-5\%','0\%','5\%','10\%'});
            legend('$A_\mathrm{Black}$','$A_\mathrm{Green}$','$A_\mathrm{White}$');
            set(legend,'Interpreter','latex')
        else
            h(2)=subplot(2,2,2);ha.A=plot(Mask_RelVol,SA'./1e6);%,Mask_RelVol(1,:),mean(SA2'));
            ylabel('Surface Area per Vol. / $\mu$m$^{-1}$','Interpreter','Latex');
        end
        xlim([0 1])
        set(ha.A(1),'Color',[0 0 0])
        set(ha.A,'linewidth',1.5)
        xlabel('Fraction of Total Volume','Interpreter','Latex');
        
        % TPB plot
        if length(NetVals)>2
            if TPB_RVA(1)~=0
                h(3)=subplot(2,2,3);ha.T=plot(Mask_RelVol(1:end),TPB_RVA(1:end)*1e-12,'k');
                xlabel('Fraction of Total Volume','Interpreter','Latex');
                ylabel('TPB Density / $\mu$m$^{-2}$','Interpreter','Latex')
                ylim([roundsf(0.9*min(0,TPB_RVA(2))*1e-12,2,'floor') roundsf(1.1*max(1,TPB_RVA(2))*1e-12,2,'ceil')])
                xlim([0 1])
                set(ha.T,'linewidth',1.5)
            end
        end
        set(h,'PlotBoxAspectRatio',[1 1 1]);
        
        if CalledFromTauFactor==1
            set(h,'TickLabelInterpreter','Latex')
            [hand.logo_A, hand.logo_map, hand.logo_alpha] = imread('TauFactor_icon2smo.png');
            h(4)=subplot(2,2,4);
            set(h(4),'position',[0.9 0.9 0.1 0.1])
            uistack(h(4),'down',3)
            imshow(hand.logo_A, hand.logo_map);
            hand.axesHandles=1;
        else
            set(h(2), 'yticklabel', {'-10%','-5%','0%','5%','10%'});
        end
        % Figure title
        
        switch Mode
            case 1
                TitStr=['RVA (Cuboid)'];
            case 2
                TitStr=['RVA (L=const.) in dir ',num2str(hand.Dir)];
            case 3
                TitStr=['RVA (A=const.) in dir ',num2str(hand.Dir),' from top'];
            case 4
                TitStr=['RVA (A=const.) in dir ',num2str(hand.Dir),' from base'];
        end
        
        %% Results annotation
        annotation(RVfig(Mode),'textbox',...
            [0 0.9 1 0.1],...
            'FitBoxToText','off',...
            'EdgeColor','none',...
            'string',{['\textbf{',TitStr,':}   '],['\verb|',hand.filename,'|']},...
            'HorizontalAlignment','center',...
            'Interpreter','Latex',...
            'FontSize',13);
        
        if length(NetVals)>2
            ResultAnnno=annotation(RVfig(Mode),'textbox',...
                [0.54 0.05 0.45 0.42],...
                'String',{...
                ['Voxel Size = ',num2str(roundsf(VoxDim(1)*1e9,4,'round')),' $\times$ ',num2str(roundsf(VoxDim(2)*1e9,4,'round')),' $\times$ ',num2str(roundsf(VoxDim(3)*1e9,4,'round')),' nm'],...
                ['Total Voxels = ',num2str(a),' $\times$ ',num2str(b),' $\times$ ',num2str(c)],...%['ROI Volume = ',num2str(roundsf((a-2)*VoxDim(1)*1e6,3,'round')),' $\times$ ',num2str(roundsf((b-2)*VoxDim(2)*1e6,3,'round')),' $\times$ ',num2str(roundsf(max(c-2,1)*VoxDim(3)*1e6,3,'round')),' $\mu$m'],...
                ['Volume Fractions:'],...
                ['    * Black = ',num2str(roundsf(VolFrac(1,1),3,'round'),'%.3f')],...
                ['    * Green = ',num2str(roundsf(VolFrac(1,2),3,'round'),'%.3f')],...
                ['    * White = ',num2str(roundsf(VolFrac(1,3),3,'round'),'%.3f')],...
                ['Surface Area per Volume / $\mu$m$^{-1}$:'],...
                ['    * Black = ',num2str(roundsf(SA(2,1)*1e-6,3,'round'))],...
                ['    * Green = ',num2str(roundsf(SA(2,2)*1e-6,3,'round'))],...
                ['    * White = ',num2str(roundsf(SA(2,3)*1e-6,3,'round'))],...
                ['Interface per Volume / $\mu$m$^{-1}$:'],...
                ['    * B-G = ',num2str(roundsf((SA(2,1)+SA(2,2)-SA(2,3))/2*1e-6,3,'round'))],...
                ['    * G-W = ',num2str(roundsf((SA(2,2)+SA(2,3)-SA(2,1))/2*1e-6,3,'round'))],...
                ['    * W-B = ',num2str(roundsf((SA(2,3)+SA(2,1)-SA(2,2))/2*1e-6,3,'round'))],...
                ['TPB Density = ',num2str(roundsf(TPB_RVA(2)*1e-12,3,'round')),' $\mu$m$^{-2}$']...
                },...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'FontSize',11,...
                'Interpreter','Latex');
        else
            ResultAnnno=annotation(RVfig(Mode),'textbox',...
                [0.25 0.1 0.5 0.38],...
                'String',{...
                [],['Voxel Size = ',num2str(roundsf(VoxDim(1)*1e9,4,'round')),' $\times$ ',num2str(roundsf(VoxDim(2)*1e9,4,'round')),' $\times$ ',num2str(roundsf(VoxDim(3)*1e9,4,'round')),' nm'],...
                [],['Total Voxels = ',num2str(a),' $\times$ ',num2str(b),' $\times$ ',num2str(c)],...%[],['ROI Volume = ',num2str(roundsf((a-2)*VoxDim(1)*1e6,3,'round')),' $\times$ ',num2str(roundsf((b-2)*VoxDim(2)*1e6,3,'round')),' $\times$ ',num2str(roundsf(max(c-2,1)*VoxDim(3)*1e6,3,'round')),' $\mu$m'],...
                [],['Volume Fractions:'],...
                ['    * Black = ',num2str(roundsf(VolFrac(1,1),3,'round'),'%.3f')],...
                ['    * White = ',num2str(roundsf(VolFrac(1,2),3,'round'),'%.3f')],...
                [],['Surface Area per Volume = ',num2str(roundsf(SA(2,1)*1e-6,3,'round')),' $\mu$m$^{-1}$'],...
                },...
                'FitBoxToText','off',...
                'LineStyle','none',...
                'FontSize',11,...
                'Interpreter','Latex');
        end
        
        if RVA~=0
            if get(hand.Check_pdfSave,'Value')==1
                switch Mode
                    case 1
                        FilStr=['_MetricsRVA_cube'];
                    case 2
                        FilStr=['_MetricsRVA_Lcon_dir',num2str(hand.Dir)];
                    case 3
                        FilStr=['_MetricsRVA_AconTop_dir',num2str(hand.Dir)];
                    case 4
                        FilStr=['_MetricsRVA_AconBase_dir',num2str(hand.Dir)];
                end
                print(RVfig(Mode),[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),FilStr],'-dpdf');
            end
        end
        
    end % No RVA, report to CW
    if RVArepeats==0
        if RVA~=0
            zz=2;
        else
            zz=1;
        end
        if length(NetVals)>2
            disp([...
                '%%%%%%%%%%%% RESULTS %%%%%%%%%%%% - ',hand.filename,char(10),...
                'Voxel Size = ',num2str(roundsf(VoxDim(1)*1e9,4,'round')),' x ',num2str(roundsf(VoxDim(2)*1e9,4,'round')),' x ',num2str(roundsf(VoxDim(3)*1e9,4,'round')),' nm',char(10),...
                'Total Voxels = ',num2str(a),' x ',num2str(b),' x ',num2str(c),char(10),...%'ROI Volume = ',num2str(roundsf(SD(1,2)*VoxDim(1)*1e6,3,'round')),' x ',num2str(roundsf(SD(2,2)*VoxDim(2)*1e6,3,'round')),' x ',num2str(roundsf(SD(3,2)*VoxDim(3)*1e6,3,'round')),' um',char(10),...
                'Volume Fractions:',char(10),...
                char(9),'Black = ',num2str(roundsf(VolFrac(1,1),3,'round')),char(10),...
                char(9),'Green = ',num2str(roundsf(VolFrac(1,2),3,'round')),char(10),...
                char(9),'White = ',num2str(roundsf(VolFrac(1,3),3,'round')),char(10),...
                'Surface Area Densities / um-1:',char(10),...
                char(9),'Black = ',num2str(roundsf(SA(zz,1)*1e-6,3,'round')),char(10),...
                char(9),'Green = ',num2str(roundsf(SA(zz,2)*1e-6,3,'round')),char(10),...
                char(9),'White = ',num2str(roundsf(SA(zz,3)*1e-6,3,'round')),char(10),...
                'TPB Density = ',num2str(roundsf(TPB_RVA(zz)*1e-12,3,'round')),' um-2']);
        else
            disp([...
                '%%%%%%%%%%%% RESULTS %%%%%%%%%%%% - ',hand.filename,char(10),...
                'Voxel Size = ',num2str(roundsf(VoxDim(1)*1e9,4,'round')),' x ',num2str(roundsf(VoxDim(2)*1e9,4,'round')),' x ',num2str(roundsf(VoxDim(3)*1e9,4,'round')),' nm',char(10),...
                'Total Voxels = ',num2str(a),' x ',num2str(b),' x ',num2str(c),char(10),...%'ROI Volume = ',num2str(roundsf(SD(1,2)*VoxDim(1)*1e6,3,'round')),' x ',num2str(roundsf(SD(2,2)*VoxDim(2)*1e6,3,'round')),' x ',num2str(roundsf(SD(3,2)*VoxDim(3)*1e6,3,'round')),' um',char(10),...
                'Volume Fractions:',char(10),...
                char(9),'Black = ',num2str(roundsf(VolFrac(1,1),3,'round')),char(10),...
                char(9),'White = ',num2str(roundsf(VolFrac(1,2),3,'round')),char(10),...
                'Surface Area Density = ',num2str(roundsf(SA(zz)*1e-6,3,'round')), ' um-1']);
        end
    end
    if toc(xxx)<120
        timeStr=[num2str(ceil(toc(xxx))),' s'];
    elseif toc(xxx)<2*3600
        timeStr=[num2str(ceil(toc(xxx)/60)),' mins'];
    else
        timeStr=[num2str(ceil(toc(xxx)/3600)),' hours'];
    end
    if max(round(MBytesA))<1024
        disp(['Calculation for ',num2str(a*b*c),' voxels (Memory < ',num2str(roundsf(max(MBytesA),2,'ceil')),' MB, Time < ',timeStr,').',char(10)])
    else
        disp(['Calculation for ',num2str(round(a*b*c/1000^2)),' million voxels (Memory < ',num2str(roundsf(max(MBytesA)/1024,2,'ceil')),' GB, Time < ',timeStr,').',char(10)])
    end
    if get(hand.Check_pdfSave,'Value')==1
        if RVA~=0
            zz=2;
        else
            zz=1;
        end
        [XLtemp]=Make_XLtemp;
        XLtemp{1,2}=hand.filename;
        XLtemp{2,2}=Dims(1);
        XLtemp{3,2}=Dims(2);
        XLtemp{4,2}=Dims(3);
        XLtemp{5,2}=VoxDim(1)*1e9;
        XLtemp{6,2}=VoxDim(2)*1e9;
        XLtemp{7,2}=VoxDim(3)*1e9;
        XLtemp{8,2}=VolFrac(1,1);
        XLtemp{9,2}=VolFrac(1,2);
        if length(NetVals)==3
            XLtemp{10,2}=VolFrac(1,3);
            XLtemp{11,2}=SA(zz,1)*1e-6;
            XLtemp{12,2}=SA(zz,2)*1e-6;
            XLtemp{13,2}=SA(zz,3)*1e-6;
            XLtemp{14,2}=(SA(zz,1)+SA(zz,2)-SA(zz,3))/2*1e-6;
            XLtemp{15,2}=(SA(zz,2)+SA(zz,3)-SA(zz,1))/2*1e-6;
            XLtemp{16,2}=(SA(zz,3)+SA(zz,1)-SA(zz,2))/2*1e-6;
            XLtemp{17,2}=TPB_RVA(zz)*1e-12;
        else
            XLtemp{9,1}=XLtemp{10,1};
            XLtemp{10,1}='Specific surface area (um-1)'; 
            XLtemp{10,2}=SA(zz,1)*1e-6;
            XLtemp(11:end,:)=[];
        end
        xlswrite([hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_Metrics'],XLtemp)
    end
end

    function [XLtemp]=Make_XLtemp
        XLtemp{1,1}='Filename';
        XLtemp{2,1}='Number of voxels in d1';
        XLtemp{3,1}='Number of voxels in d2';
        XLtemp{4,1}='Number of voxels in d3';
        XLtemp{5,1}='Voxel Size in d1 (nm)';
        XLtemp{6,1}='Voxel Size in d2 (nm)';
        XLtemp{7,1}='Voxel Size in d3 (nm)';
        XLtemp{8,1}='Black VF';
        XLtemp{9,1}='Green VF';
        XLtemp{10,1}='White VF';
        XLtemp{11,1}='Black VSSA (um-1)';
        XLtemp{12,1}='Green VSSA (um-1)';
        XLtemp{13,1}='White VSSA (um-1)';
        XLtemp{14,1}='B-G VSIA (um-1)';
        XLtemp{15,1}='G-W VSIA (um-1)';
        XLtemp{16,1}='W-B VSIA (um-1)';
        XLtemp{17,1}='TPB (um-2)';

    end
% if length(NetVals)>2
%     figure;for i =1:c;imagesc(TPBmap(:,:,i));axis equal tight;pause(0.01);end
% end

%% Functions for finding indicies for RVA
% Cuboid
    function [SD]=RVA_Net(Dims,shrink)
        s=(1-shrink)^(1/3);
        d=(1-s)/2;
        
        SD(1,:)=[floor(d*Dims(1)+1) ceil(Dims(1)*(1-d))];
        SD(2,:)=[floor(d*Dims(2)+1) ceil(Dims(2)*(1-d))];
        SD(3,:)=[floor(d*Dims(3)+1) ceil(Dims(3)*(1-d))];
        
        SD(SD<1)=1;
    end
% L=const.
    function [SD]=RVA_Net_Lconst(Dims,shrink)
        s=(1-shrink)^(1/2);
        d=(1-s)/2;
        
        SD(1,:)=[1 Dims(1)];
        SD(2,:)=[floor(d*Dims(2)+1) ceil(Dims(2)*(1-d))];
        SD(3,:)=[floor(d*Dims(3)+1) ceil(Dims(3)*(1-d))];
        
        SD(SD<1)=1;
    end
% A=const.
    function [SD]=RVA_Net_Aconst(Dims,shrink,Mode)
        s=1-shrink;
        d=(1-s);
        if Mode==3
            SD(1,:)=[1 ceil(Dims(1)*(1-d))];
        else
            SD(1,:)=[floor(1+Dims(1)*d) Dims(1)];
        end
        SD(2,:)=[1 Dims(2)];
        SD(3,:)=[1 Dims(3)];
        
        SD(SD<1)=1;
    end

%% Function for rounding
    function [y]=roundsf(number,sfs,method)
        %opt = {'round','floor','ceil','fix'};R
        og = 10.^(floor(log10(abs(number)) - sfs + 1));
        y = feval(method,number./og).*og;
        y(find(number==0)) = 0;
    end
end

% Alternative isosurface implementation:
% IS=isosurface(padarray(hand.Net_Or,[1 1 1],0))
% a = IS.vertices(IS.faces(:, 2), :) - IS.vertices(IS.faces(:, 1), :);
% b = IS.vertices(IS.faces(:, 3), :) - IS.vertices(IS.faces(:, 1), :);
% c = cross(a, b, 2);
% area = 1/2 * sum(sqrt(sum(c.^2, 2)))* ;
% fprintf('\nThe surface area is %f\n\n', area);