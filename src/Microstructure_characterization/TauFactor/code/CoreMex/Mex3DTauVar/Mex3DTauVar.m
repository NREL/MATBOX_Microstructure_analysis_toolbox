function [T1,T2]=Mex3DTauVar(c,Check_f,T1,T2,omw,...
    Cheq1P_Xm,Cheq1P_Xp,Cheq1P_Ym,Cheq1P_Yp,Cheq1P_Zm,Cheq1P_Zp,...
    Cheq2P_Xm,Cheq2P_Xp,Cheq2P_Ym,Cheq2P_Yp,Cheq2P_Zm,Cheq2P_Zp,...
    R_Xm1,R_Xp1,R_Ym1,R_Yp1,R_Zm1,R_Zp1,...
    R_Xm2,R_Xp2,R_Ym2,R_Yp2,R_Zm2,R_Zp2) %#codegen

if c>1 %3D
    for i=1:Check_f
        T1(1:end-2)=omw*T1(1:end-2)+...
            R_Xm1.*T2(Cheq1P_Xm)+...
            R_Xp1.*T2(Cheq1P_Xp)+...
            R_Ym1.*T2(Cheq1P_Ym)+...
            R_Yp1.*T2(Cheq1P_Yp)+...
            R_Zm1.*T2(Cheq1P_Zm)+...
            R_Zp1.*T2(Cheq1P_Zp);
        T2(1:end-2)=omw*T2(1:end-2)+...
            R_Xm2.*T1(Cheq2P_Xm)+...
            R_Xp2.*T1(Cheq2P_Xp)+...
            R_Ym2.*T1(Cheq2P_Ym)+...
            R_Yp2.*T1(Cheq2P_Yp)+...
            R_Zm2.*T1(Cheq2P_Zm)+...
            R_Zp2.*T1(Cheq2P_Zp);
    end
else %2D
    for i=1:Check_f
        T1(1:end-2)=omw*T1(1:end-2)+...
            R_Xm1.*T2(Cheq1P_Xm)+...
            R_Xp1.*T2(Cheq1P_Xp)+...
            R_Ym1.*T2(Cheq1P_Ym)+...
            R_Yp1.*T2(Cheq1P_Yp);
        T2(1:end-2)=omw*T2(1:end-2)+...
            R_Xm2.*T1(Cheq2P_Xm)+...
            R_Xp2.*T1(Cheq2P_Xp)+...
            R_Ym2.*T1(Cheq2P_Ym)+...
            R_Yp2.*T1(Cheq2P_Yp);
    end
end