%% Information
%Purpose¡GImage analysis of river surface velocity

%% 1.Setting parameters
clear all;format long;
[L,R]=textread('Yufeng1081022CLNLR.txt','%f %f','headerlines',1);
PM_CameraLR=[L,R];
PM_VminA=10;            % Enhanced image parameters A
PM_VminB=10;            % Enhanced image parameters B
PM_Area_percent=1.0;    % Determine whether the image matching window is the threshold percentage of the flowing water body
PM_Climit=0.6;          % Correlation coefficient threshold
PM_WindowsL=50;        % Match the side length of the window
PM_WindowsD=50;        % Match the side length of the window
PM_Limit_WFTSN_min=110;  % The minimum number of successful image matching
PM_UN=1624;             % The horizontal length of the original image
PM_VN=1234;             % Original image longitudinal side length
PM_X=0;                 % Solve the initial value of coordinates of X
PM_Y=0;                 % Solve the initial value of coordinates of Y
PM_T=0.04;              % Image time interval
PM_Range=50;           % Analyze the reserved border of the image
PM_Distance1R=30;       % Limitation of the flow range, the search window moves the maximum distance to the right
PM_Distance2L=0;        % Limitation of the flow range, the search window moves the maximum distance to the left
PM_Distance3U=10;       % Limitation of the flow range, the search window moves the maximum distance to the up
PM_Distance4D=10;       % Limitation of the flow range, the search window moves the maximum distance to the down
PM_Z=684.53;           % Water surface elevation
PM_StationH=693-1.41;
PM_CameraLR(9,1)=PM_CameraLR(9,1)+PM_StationH;    PM_CameraLR(9,2)=PM_CameraLR(9,2)+PM_StationH;
PM_Total_Pixel=PM_UN*PM_VN;

%% 2.LSPIV.

for PM_photo_j=1;
    for PM_photo_k=1:3;
        for PM_photo_i=1:2;
            name=PM_photo_j*100+PM_photo_k*10+PM_photo_i+1;
            fname=[num2str(name),'.tif'];
            Photo_RGB(:,:,:,PM_photo_i)=imread([fname]);
            Photo_Gray(:,:,PM_photo_i)=rgb2gray(Photo_RGB(:,:,:,PM_photo_i));
            CLN_E(:,1)=PM_CameraLR(:,1);
        end
        clear TD_i TD_I Time_Date Ax Ah Bx Bh fname name A B C D RiverAreaData i;
                Photo_Gray(:,:,3)=Photo_Gray(:,:,1)-Photo_Gray(:,:,2)>PM_VminA;
        Time_step=PM_T;
        
        %% Image matching - IA
        VSN_check_times=0;  VS_all_dataN=0; WindowsD=PM_WindowsD;  WindowsL=PM_WindowsL; Photo_Number=1;
        WFTSN=1;
        for WF_M= 1 + PM_Range  : WindowsD : PM_UN-WindowsL - PM_Range
            for WF_N= 1 + PM_Range : WindowsD : PM_VN - WindowsL - PM_Range
                Cmax=PM_Climit;    times=0;
                WF=Photo_Gray(WF_N : WF_N+WindowsL-1 , WF_M : WF_M+WindowsL-1 , 3);    WFT=sum(sum(WF));
                if WFT>(PM_Area_percent/100)*PM_WindowsL*PM_WindowsL
                    WF=Photo_Gray(WF_N : WF_N+WindowsL-1 , WF_M : WF_M+WindowsL-1 , 1);    WFT=sum(sum(WF));
                    
                    %% Image matching - SA
                    WS_Mmin=WF_M-PM_Distance2L; WS_Mmax=WF_M+PM_Distance1R;
                    if WS_Mmin<1;               WS_Mmin=1;               end
                    if WS_Mmax>=PM_UN-WindowsL;  WS_Mmax=PM_UN-WindowsL;      end
                    for WS_M=WS_Mmin: 1 :WS_Mmax
                        WS_Nmin=WF_N-PM_Distance3U; WS_Nmax=WF_N+PM_Distance4D;
                        if WS_Nmin<1;               WS_Nmin=1;               end
                        if WS_Nmax>=PM_VN-WindowsL;  WS_Nmax=PM_VN-WindowsL;      end
                        for WS_N=WS_Nmin: 1 :WS_Nmax
                            WS=Photo_Gray(WS_N : WS_N+WindowsL-1 , WS_M : WS_M+WindowsL-1 , 2);    WST=sum(sum(WS));
                            
                            %% Correlation coefficient
                            avgF=WFT/(WindowsL*WindowsL);      avgS=WST/(WindowsL*WindowsL);        DF=WF-avgF;          DS=WS-avgS;
                            
                            DC=abs(1-(sum(sum(DF.*DS))/sqrt(sum(sum(DF.*DF))*sum(sum(DS.*DS)))));
                            if DC<(1-Cmax) & WF_M-WS_M+WF_N-WS_N~=0
                                if abs(WF_N-WS_N)==1;   WS_N=WF_N;  end
                                CLN_UV(WFTSN,1:2,1)=[WF_M,WF_N];            CLN_UV(WFTSN,1:2,2)=[WS_M,WS_N];
                                times=times+1;                              Cmax=1-DC;
                            end
                        end
                    end
                    if Cmax~=PM_Climit;    WFTSN=WFTSN+1;                end
                end
            end
        end
        clear WF WFT WF_M WF_N WS WST WS_M WS_N WS_Mmax WS_Mmin WS_Nmax WS_Nmin;
        
        %% LSPIV - CLN
        uc=CLN_E(1,1);        vc=CLN_E(2,1);
        f=CLN_E(3,1);         arfa=CLN_E(4,1);   sita=CLN_E(5,1)+4.0/180*pi();      tou=CLN_E(6,1);
        xc=CLN_E(7,1);        yc=CLN_E(8,1);     zc=CLN_E(9,1);        z=PM_Z;
        r11=-cos(arfa)*cos(sita)-sin(arfa)*cos(tou)*sin(sita);  r12=sin(arfa)*cos(sita)-cos(arfa)*cos(tou)*sin(sita);   r13=-sin(tou)*sin(sita);
        r21=cos(arfa)*sin(sita)-sin(arfa)*cos(tou)*cos(sita);   r22=-sin(arfa)*sin(sita)-cos(arfa)*cos(tou)*cos(sita);  r23=-sin(tou)*cos(sita);
        r31=-sin(arfa)*sin(tou);                                r32=-cos(arfa)*sin(tou);                                r33=cos(tou);
        
        clear x y u v UX UY VX VY Uu Uuc Uf Vv Vvc Vf CLN_XYZ d VS dx dy VS_all_data VS_all_dataA VS_all_dataN VS_vector_data;
        for s=1:2
            u=CLN_UV(1:WFTSN-1,1,s);      v=CLN_UV(1:WFTSN-1,2,s);
            UX(1:WFTSN-1,1)=r31*u-r31*uc+r11*f;  UY(1:WFTSN-1,1)=r32*u-r32*uc+r12*f;  VX(1:WFTSN-1,1)=r31*v-r31*vc+r21*f;  VY(1:WFTSN-1,1)=r32*v-r32*vc+r22*f;
            Uu(1:WFTSN-1,1)=xc*r31+yc*r32-z*r33+zc*r33;  Uuc(1:WFTSN-1,1)=-xc*r31-yc*r32+z*r33-zc*r33;    Uf(1:WFTSN-1,1)=-xc*r11-yc*r12+z*r13-zc*r13;
            Vv(1:WFTSN-1,1)=xc*r31+yc*r32-z*r33+zc*r33;  Vvc(1:WFTSN-1,1)=-xc*r31-yc*r32+z*r33-zc*r33;    Vf(1:WFTSN-1,1)=-xc*r21-yc*r22+z*r23-zc*r23;
            
            x(1:WFTSN-1,1) = (v.*(Vv)+vc.*(Vvc)-f.*(Vf)-(u.*(Uu)+uc.*(Uuc)-f.*(Uf))./(UY).*(VY)) ./ ((VX)-(UX)./(UY).*(VY));
            y(1:WFTSN-1,1) = (v.*(Vv)+vc.*(Vvc)-f.*(Vf)-(u.*(Uu)+uc.*(Uuc)-f.*(Uf))./(UX).*(VX)) ./ ((VY)-(UY)./(UX).*(VX));
            
            CLN_XYZ(1:WFTSN-1,1:4,s)=[x(1:WFTSN-1,1),y(1:WFTSN-1,1),CLN_UV(1:WFTSN-1,1,s),CLN_UV(1:WFTSN-1,2,s)];
        end
        
        n=WFTSN-1;
        i=n;
        fprintf('Image coordinate calculation completed¡ATotal success %d Number¡C\n',i);
        clear dFu_x dFu_y dFu_z dFv_x dFv_y dFv_z check_Anan check_Bnan check_nan;
        clear u v uc vc fu fv x y z xc yc zc r11 r12 r13 r21 r22 r23 r31 r32 r33 sita arfa tou;
        clear Fu dFu_uc dFu_vc dFu_fu dFu_fv dFu_arfa dFu_sita dFu_tou dFu_xc dFu_yc dFu_zc;
        clear Fv dFv_uc dFv_vc dFv_fu dFv_fv dFv_arfa dFv_sita dFv_tou dFv_xc dFv_yc dFv_zc;
        
        if WFTSN==1 || i==0
            VSavg(PM_photo_k,PM_photo_j)=0;        VSrmse=0;
            fprintf('Failure¡C\r\n');
        else
            n=i;
            dx=CLN_XYZ(:,1,1)-CLN_XYZ(:,1,2);            dy=CLN_XYZ(:,2,1)-CLN_XYZ(:,2,2);        d(:,1)=sqrt(dx.^2+dy.^2);
            VS(:,1)=d(:,1)/Time_step;
            VS_all_data(:,1:10)=[VS(:,1),d(:,1),CLN_XYZ(:,1,1),CLN_XYZ(:,2,1),CLN_XYZ(:,3,1),CLN_XYZ(:,4,1),CLN_XYZ(:,1,2),CLN_XYZ(:,2,2),CLN_XYZ(:,3,2),CLN_XYZ(:,4,2)];
            VS_all_data(:,11)=Time_step;
            VS_all_dataA=VS_all_data;
            VScheck=0;
            while VScheck==0;
                VSavg(PM_photo_k,PM_photo_j)=sum((VS_all_data(:,1)))/n;        VSrms=sqrt(sum((VS_all_data(:,1)-VSavg(PM_photo_k,PM_photo_j)).^2)/(n-1));       j=0;
                for i=1:n
                    if VS_all_data(i,1)<VSavg(PM_photo_k,PM_photo_j)+VSrms*3 & VS_all_data(i,1)>VSavg(PM_photo_k,PM_photo_j)-VSrms*3 & VS_all_data(i,1) ~= 0
                        j=j+1;
                        VS2_all_data(j,:)=VS_all_data(i,:);
                    end
                end
                clear VS_all_data;            VS_all_data=VS2_all_data;            clear VS2_all_data;
                if j==i;                VScheck=1;            end
                n=j;
            end
            VSavg(PM_photo_k,PM_photo_j)=sum((VS_all_data(:,1)))/n;
        end
        
        if n>=1
            Photo_mark(:,:,1)=Photo_Gray(:,:,1);
            Photo_mark(:,:,2)=Photo_Gray(:,:,1);
            Photo_mark(:,:,3)=Photo_Gray(:,:,1);

            for i=1:n
                if  VS_all_data(i,6)>100
                    MVF=VS_all_data(i,6)+WindowsL/2;                  MUF=VS_all_data(i,5)+WindowsL/2;
                    MVS=VS_all_data(i,10)+WindowsL/2;                 MUS=VS_all_data(i,9)+WindowsL/2;
                    Photo_mark(MVF:MVF+5,MUF:MUF+5,2)=255;              Photo_mark(MVS:MVS+5,MUS:MUS+5,2)=255;
                    Photo_mark(MVF:MVF+5,MUF:MUF+5,3)=0;              Photo_mark(MVS:MVS+5,MUS:MUS+5,3)=0;
                    Photo_mark(MVF:MVF+3,MUF:MUF+3,1)=0;            Photo_mark(MVS:MVS+3,MUS:MUS+3,1)=0;
                    VS_vector_data(i,1:4)=[MUF,MVF,VS_all_data(i,1)*cos(atan((VS_all_data(i,10)-VS_all_data(i,6))/(VS_all_data(i,9)-VS_all_data(i,5)))),VS_all_data(i,1)*sin(atan((VS_all_data(i,10)-VS_all_data(i,6))/(VS_all_data(i,9)-VS_all_data(i,5))))];
                end
            end

            for WF_M= 1+PM_Range : WindowsL : PM_UN-PM_Range+1
                Photo_mark( 1+PM_Range:PM_VN-PM_Range-34 , WF_M-1:WF_M+1,1)=255;
                Photo_mark( 1+PM_Range:PM_VN-PM_Range-34 , WF_M-1:WF_M+1,2)=0;
                Photo_mark( 1+PM_Range:PM_VN-PM_Range-34 , WF_M-1:WF_M+1,3)=0;
            end
            for WF_N= 1 +PM_Range: WindowsL : PM_VN-PM_Range+1
                Photo_mark( WF_N-1:WF_N+1 , 1+PM_Range:PM_UN-PM_Range-24,1)=255;
                Photo_mark( WF_N-1:WF_N+1 , 1+PM_Range:PM_UN-PM_Range-24,2)=0;
                Photo_mark( WF_N-1:WF_N+1 , 1+PM_Range:PM_UN-PM_Range-24,3)=0;
            end
            imshow(Photo_mark);hold on;
            quiver(VS_vector_data(:,1),VS_vector_data(:,2),VS_vector_data(:,3),VS_vector_data(:,4),0.5,'g');
            
        end
    end
end
fprintf('Average surface velocity%7.3f (m/s)¡C\r\n',VSavg);



