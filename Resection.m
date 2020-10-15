%% Information
%Purpose¡GFinds the camera's internal and external parameters.

%% 1.Input
%Input the control points "object space coordinates" and "image space coordinates" of a single image
%The data from left to right are: Number, X, Y, Z, U, V.
clear all;format long;
[A,B,C,D,E,F]=textread('near.txt','%s %f %f %f %f %f','headerlines',1); 
Int_uc=1624/2;     Int_vc=1234/2;     Int_f=3600;    Int_arfa=3;  Int_sita=3;    Int_tou=1.0;     Int_xc=-0.5;     Int_yc=-1.0;     Int_zc=1.5;
IMG_data=[B C D E F];       clear A B C D E F;
%////////////////////////////Camera parameters//////////////////////////////////////%
UTL=0.0071456;     VTL=0.0054296;  %Length and width of film
UTP=1624;       VTP=1234;    %Total pixels of the camera
UPL=0.0000044;    VPL=0.0000044; %Pixel size
TUC=UTP/2;      TVC=VTP/2;   %Image center point
%////////////////////////////Camera parameters//////////////////////////////////////%
Ndata(:,2:6)=IMG_data(:,1:5);   NR=size(Ndata,1);
%% 2. Resection
syms u v uc vc f x y z xc yc zc r11 r12 r13 r21 r22 r23 r31 r32 r33 sita arfa tou;
r11=-cos(arfa)*cos(sita)-sin(arfa)*cos(tou)*sin(sita);  r12=sin(arfa)*cos(sita)-cos(arfa)*cos(tou)*sin(sita);   r13=-sin(tou)*sin(sita);
r21=cos(arfa)*sin(sita)-sin(arfa)*cos(tou)*cos(sita);   r22=-sin(arfa)*sin(sita)-cos(arfa)*cos(tou)*cos(sita);  r23=-sin(tou)*cos(sita);
r31=-sin(arfa)*sin(tou);                                r32=-cos(arfa)*sin(tou);                                r33=cos(tou);

Fu=(-f*((r11*(x-xc)+r12*(y-yc)+r13*(z-zc))/(r31*(x-xc)+r32*(y-yc)+r33*(z-zc))))+uc-u;
Fv=(-f*((r21*(x-xc)+r22*(y-yc)+r23*(z-zc))/(r31*(x-xc)+r32*(y-yc)+r33*(z-zc))))+vc-v;

dFu_uc=diff(Fu,uc);          dFu_vc=diff(Fu,vc);      dFu_f=diff(Fu,f);     dFu_arfa=diff(Fu,arfa);      dFu_sita=diff(Fu,sita);  dFu_tou=diff(Fu,tou);
dFu_xc=diff(Fu,xc);          dFu_yc=diff(Fu,yc);      dFu_zc=diff(Fu,zc);   dFu_x=diff(Fu,x);            dFu_y=diff(Fu,y);        dFu_z=diff(Fu,z);

dFv_uc=diff(Fv,uc);          dFv_vc=diff(Fv,vc);      dFv_f=diff(Fv,f);     dFv_arfa=diff(Fv,arfa);      dFv_sita=diff(Fv,sita);  dFv_tou=diff(Fv,tou);
dFv_xc=diff(Fv,xc);          dFv_yc=diff(Fv,yc);      dFv_zc=diff(Fv,zc);   dFv_x=diff(Fv,x);            dFv_y=diff(Fv,y);        dFv_z=diff(Fv,z);

CDT1=0;    Times1=1;    CDT2=0;    Times2=1;    CDT3=0;    Times3=1;      
uc=Int_uc;     vc=Int_vc;     f=Int_f;    arfa=Int_arfa;  sita=Int_sita;    tou=Int_tou;     xc=Int_xc;     yc=Int_yc;     zc=Int_zc; 

while CDT1==0;    
    for n=1:NR
        x=Ndata(n,2);   y=Ndata(n,3);    z=Ndata(n,4);      u=Ndata(n,5);   v=Ndata(n,6);
        A(n,1:3)=[eval(dFu_arfa),eval(dFu_sita),eval(dFu_tou)]; 
        A(NR+n,1:3)=[eval(dFv_arfa),eval(dFv_sita),eval(dFv_tou)];    
        B(n,1)=[eval(Fu)*-1];
        B(NR+n,1)=[eval(Fv)*-1];
    end

    RA=rank(A);          
    if RA~=3     
        fprintf('The accuracy is not up to the requirement %d\n');               
        NElement(:,1)=E(:,1);        CDT1=1;
    end         
    C(:,1)=inv(A'*A)*(A'*B);                  
    ETC1(:,Times1)=C(:,1);     
    arfa=arfa+C(1,1);   sita=sita+C(2,1);   tou=tou+C(3,1);
    ETE1(:,Times1)=[uc,vc,f,arfa,sita,tou,xc,yc,zc]; 
    AVC=abs(C(1,1))+abs(C(2,1))+abs(C(3,1));    AVC=AVC/3;

    if AVC<=0.00001 & CDT1==0;
        clear C;
        fprintf('Angle accuracy has reached the requirement %d\n');    
        while CDT2==0;            
            for n=1:NR
                x=Ndata(n,2);   y=Ndata(n,3);    z=Ndata(n,4);      u=Ndata(n,5);   v=Ndata(n,6);
                A(n,1:6)=[eval(dFu_arfa),eval(dFu_sita),eval(dFu_tou),eval(dFu_xc),eval(dFu_yc),eval(dFu_zc)]; 
                A(NR+n,1:6)=[eval(dFv_arfa),eval(dFv_sita),eval(dFv_tou),eval(dFv_xc),eval(dFv_yc),eval(dFv_zc)];    
                B(n,1)=[eval(Fu)*-1];
                B(NR+n,1)=[eval(Fv)*-1];
            end
            RA=rank(A);          
            if RA~=6    
                fprintf('The accuracy is not up to the requirement %d\n');                            
                NElement(:,1)=E(:,1);                CDT1=1;                CDT2=1;
            end         
            C(:,1)=inv(A'*A)*(A'*B);                  
            ETC2(:,Times2)=C(:,1);  
            arfa=arfa+C(1,1);   sita=sita+C(2,1);   tou=tou+C(3,1);
            xc=xc+C(4,1);      yc=yc+C(5,1);    zc=zc+C(6,1);    
            ETE2(:,Times2)=[uc,vc,f,arfa,sita,tou,xc,yc,zc];
            AVC=abs(C(1,1))+abs(C(2,1))+abs(C(3,1))+abs(C(4,1))+abs(C(5,1))+abs(C(6,1));    AVC=AVC/6; 
            if AVC<=0.00001 & CDT2==0;
            clear C;
            fprintf('Angle accuracy has reached the requirement %d\n');  
                while CDT3==0;            
                    for n=1:NR
                        x=Ndata(n,2);   y=Ndata(n,3);    z=Ndata(n,4);      u=Ndata(n,5);   v=Ndata(n,6);
                        A(n,1:7)=[eval(dFu_f),eval(dFu_arfa),eval(dFu_sita),eval(dFu_tou),eval(dFu_xc),eval(dFu_yc),eval(dFu_zc)]; 
                        A(NR+n,1:7)=[eval(dFv_f),eval(dFv_arfa),eval(dFv_sita),eval(dFv_tou),eval(dFv_xc),eval(dFv_yc),eval(dFv_zc)];    
                        B(n,1)=[eval(Fu)*-1];
                        B(NR+n,1)=[eval(Fv)*-1];
                    end
                    RA=rank(A);          
                    if RA~=7     
                        fprintf('The accuracy is not up to the requirement %d\n');                             
                        CDT1=1;                CDT2=1;                CDT3=1;
                    end         
                    C(:,1)=inv(A'*A)*(A'*B);                
                    ETC3(:,Times3)=C(:,1);
                    f=f+C(1,1);        arfa=arfa+C(2,1);   sita=sita+C(3,1);   tou=tou+C(4,1);
                    xc=xc+C(5,1);      yc=yc+C(6,1);         zc=zc+C(7,1);    
                    ETE3(:,Times3)=[uc,vc,f,arfa,sita,tou,xc,yc,zc]; 
                    AVC=abs(C(1,1))+abs(C(2,1))+abs(C(3,1))+abs(C(4,1))+abs(C(5,1))+abs(C(6,1))+abs(C(7,1));       AVC=AVC/7; 
                    if AVC<=0.00001 & CDT3==0;
                        fprintf('Angle accuracy has reached the requirement %d\n');                 
                        CLN_E(1:9,1)=[uc,vc,f,arfa,sita,tou,xc,yc,zc];
                        CLN_E(10,1)=f*UPL; 
                        CLN_E(11:13,1)=[ arfa*180/pi , sita*180/pi , tou*180/pi ];
                        CDT1=1;                CDT2=1;                CDT3=1;
                        delta(:,1)=B-A*C(:,1);
                        Tsigma(1,1)=sqrt((delta(:,1)'*delta(:,1))/(2*NR-10));            
                    end
                    if Times3==200 & CDT2==0;
                        fprintf('The accuracy is not up to the requirement %d\n');                       
                        CDT1=1;                CDT2=1;                CDT3=1;
                    end
                    Times3=Times3+1;   
                end  
                CDT2=1;
            end
            if Times2==500 & CDT2==0;
                fprintf('The accuracy is not up to the requirement %d\n');           
                CDT1=1;                CDT2=1;
            end
            Times2=Times2+1;   
        end
        CDT1=1;
    end

    if Times1==500 & CDT1==0;
        fprintf('The accuracy is not up to the requirement %d\n');                
        NElement(:,1)=E(:,1);
        CDT1=1;
    end
    Times1=Times1+1;    
    clear A B ;    
end
save('CLNE.txt','CLN_E','-ascii');
clear CDT Times Ndata RA AVC NR NL TUC UPL UTL UTP TVC VPL VTL VTP;
clear u v uc vc fu fv x y z xc yc zc r11 r12 r13 r21 r22 r23 r31 r32 r33 sita arfa tou;
clear Fu dFu_uc dFu_vc dFu_fu dFu_fv dFu_arfa dFu_sita dFu_tou dFu_xc dFu_yc dFu_zc dFu_x dFu_y dFu_z;  
clear Fv dFv_uc dFv_vc dFv_fu dFv_fv dFv_arfa dFv_sita dFv_tou dFv_xc dFv_yc dFv_zc dFv_x dFv_y dFv_z;
