function [T1,T2]=Test_Mex3DTauAni

c=3;
Check_f=100;
T1=ones(102,1);
T2=T1;
omw=1;
NN_aVw1=ones(100,1);
NN_aVw2=ones(100,1);
Cheq1P_Xm=uint32([1:100]');
Cheq1P_Xp=Cheq1P_Xm;
Cheq1P_Ym=Cheq1P_Xm;
Cheq1P_Yp=Cheq1P_Xm;
Cheq1P_Zm=Cheq1P_Xm;
Cheq1P_Zp=Cheq1P_Xm;
Cheq2P_Xm=Cheq1P_Xm;
Cheq2P_Xp=Cheq1P_Xm;
Cheq2P_Ym=Cheq1P_Xm;
Cheq2P_Yp=Cheq1P_Xm;
Cheq2P_Zm=Cheq1P_Xm;
Cheq2P_Zp=Cheq1P_Xm;

c_X=1;
c_Y=1;
c_Z=1;

[T1,T2]=Mex3DTauAni(c,Check_f,T1,T2,omw,NN_aVw1,NN_aVw2,...
    Cheq1P_Xm,Cheq1P_Xp,Cheq1P_Ym,Cheq1P_Yp,Cheq1P_Zm,Cheq1P_Zp,...
    Cheq2P_Xm,Cheq2P_Xp,Cheq2P_Ym,Cheq2P_Yp,Cheq2P_Zm,Cheq2P_Zp,...
    c_X,c_Y,c_Z); %#codegen