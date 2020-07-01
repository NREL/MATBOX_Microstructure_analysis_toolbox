function [T1,T2]=Test_Mex3DTauVar

c=3;
Check_f=100;
T1=ones(102,1);
T2=T1;
omw=1;
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
R_Xm1=ones(100,1);
R_Xp1=R_Xm1;
R_Ym1=R_Xm1;
R_Yp1=R_Xm1;
R_Zm1=R_Xm1;
R_Zp1=R_Xm1;
R_Xm2=R_Xm1;
R_Xp2=R_Xm1;
R_Ym2=R_Xm1;
R_Yp2=R_Xm1;
R_Zm2=R_Xm1;
R_Zp2=R_Xm1;

[T1,T2]=Mex3DTauVar(c,Check_f,T1,T2,omw,...
    Cheq1P_Xm,Cheq1P_Xp,Cheq1P_Ym,Cheq1P_Yp,Cheq1P_Zm,Cheq1P_Zp,...
    Cheq2P_Xm,Cheq2P_Xp,Cheq2P_Ym,Cheq2P_Yp,Cheq2P_Zm,Cheq2P_Zp,...
     R_Xm1,R_Xp1,R_Ym1,R_Yp1,R_Zm1,R_Zp1,...
    R_Xm2,R_Xp2,R_Ym2,R_Yp2,R_Zm2,R_Zp2);
