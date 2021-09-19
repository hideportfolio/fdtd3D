clear
close all
clc

%%
%%%%%%%%         真空なし

%%%%%%%  センサ20mm   66ピクセル(305e-6)   328ピクセル(61e-6)

%%%%%%%%%%%%%%パラメータ設定
velocity_longitudinal_water=1500;       %水中縦波音速
velocity_share_water=0;                 %水中横波音速
frequency=1.0e+5;                       %周波数
rhow=1000;                              %水密度
dimension=3;                            %次元
Kwl=rhow*velocity_longitudinal_water^2; %縦波　定数
Kws=rhow*velocity_share_water^2;        %横波　定数

density_water=1000;


rhob=2000;                              %骨密度
velocity_longitudinal_bone3=4250;       %骨中縦波音速 3方向axial       
velocity_longitudinal_bone1=3460;       %骨中縦波音速 1方向       これを4250にしたらバグる

l=0.5037;                               %ポアソン比0.33のときの縦波音速、横波音速の比
velocity_share_bone4=velocity_longitudinal_bone3*l;               %骨中横波音速
velocity_share_bone6=velocity_longitudinal_bone1*l;               %骨中横波音速axial方向

dx=305e-6;                            %空間分解能dx
dy=dx;                                  %空間分解能dy
dz=dx;                                  %空間分解能dz

courant=1/(velocity_longitudinal_water*sqrt((1/dx^2)+(1/dy^2)));  %courantの安定条件式
dt=dx/velocity_longitudinal_bone3/sqrt(dimension);                %時間分解能

nt=4000;                                %for文回す回数


%%
%入力波形
N=round(1/frequency/courant);
wave(:)=sin(2*pi*frequency*courant*(1:N));
wd(:)=0.5-0.5*cos(2*pi*frequency*courant*(1:N));
wave2=wave.*wd;

% figure;plot(wave2);

%%

    
%%%%%%%     
    Model2=single(zeros(198,297,2470/5));
    Model3d=single(zeros(198+100,297,2470/5+10));

  for i=11:2471+10
     bone=csvread(sprintf('femurmodel_300μm_%04d.csv',i-10)); %読み込み
     model2=bone>70; %ロジック行列化
     Model(:,:,i)=model2;
     Mod(:,:,i)=bone;
  end
  
%%
  
  for i=1:2470/5
      for Nx=1:198
          for Ny=1:297
              Model2(Nx,Ny,i)=sum(Model(Nx,Ny,10+(i-1)*5+1)+Model(Nx,Ny,10+(i-1)*5+2)+Model(Nx,Ny,10+(i-1)*5+3)+...
          Model(Nx,Ny,10+(i-1)*5+4)+Model(Nx,Ny,10+(i-1)*5+5))/5;
          end
      end
  end
  

Model3d(51:198+50,:,1:2470/5) = Model2;
nx=298;
ny=297;
nz=2470/5+10;

%%
%%%%%%%%%%%%%%応力・粒子速度行列の表記
%XYZ平面
Txx(:,:,:)=single(zeros(nx,ny,nz));
Tyy(:,:,:)=single(zeros(nx,ny,nz));
Tzz(:,:,:)=single(zeros(nx,ny,nz));

% Txy(:,:,:)=single(zeros(nx+1,ny+1,nz));
% Tyz(:,:,:)=single(zeros(nx,ny+1,nz+1));
% Tzx(:,:,:)=single(zeros(nx+1,ny,nz+1));

Txy(:,:,:)=single(zeros(nx+1,ny+1,nz+1));
Tyz(:,:,:)=single(zeros(nx+1,ny+1,nz+1));
Tzx(:,:,:)=single(zeros(nx+1,ny+1,nz+1));

Ux(:,:,:)=single(zeros(nx+1,ny,nz));
Uy(:,:,:)=single(zeros(nx,ny+1,nz));
Uz(:,:,:)=single(zeros(nx,ny,nz+1));

%%
%%%%%  密度パラメータ
density_x=(Model3d(2:nx,:,:)+Model3d(1:nx-1,:,:))/2;  %密度パラメータ
density_y=(Model3d(:,2:ny,:)+Model3d(:,1:ny-1,:))/2;
density_z=(Model3d(:,:,2:nz)+Model3d(:,:,1:nz-1))/2;


% 
% %%%%%%%%%%%%     密度情報を入れていく
for n=1:nx-1
    for m=1:ny
        for l=1:nz
            if density_x(n,m,l)==0
                density_x(n,m,l)=rhow;
                     
            else
                density_x(n,m,l)=rhob;
                
            end
        end
    end
end

for n=1:nx
    for m=1:ny-1
        for l=1:nz
            if density_y(n,m,l)==0
                density_y(n,m,l)=rhow;

            else
                density_y(n,m,l)=rhob;
            end
        end
    end
end

for n=1:nx
    for m=1:ny
        for l=1:nz-1
            if density_z(n,m,l)==0
                density_z(n,m,l)=rhow;
               
            else
                density_z(n,m,l)=rhob;
            end
        end
    end
end



%%
%%%%%%%%%%%%%%モデルに水、骨の情報を反映させていく

%%%%%%%%%%%%%%弾性定数 C44,C55,C66 の設定
c44=zeros(nx,ny,nz);
c55=zeros(nx,ny,nz);
c66=zeros(nx,ny,nz);
C44=zeros(nx-1,ny-1,nz-1);
C55=zeros(nx-1,ny-1,nz-1);
C66=zeros(nx-1,ny-1,nz-1);

for n=1:nx
    for m=1:ny
        for l=1:nz
            if Model3d(n,m,l)==0
                c44(n,m,l)=rhow*velocity_share_water^2;
                c55(n,m,l)=rhow*velocity_share_water^2;
                c66(n,m,l)=rhow*velocity_share_water^2;

            else
                c44(n,m,l)=rhob*velocity_share_bone4^2; %モデル×骨中横波音速
                c55(n,m,l)=rhob*velocity_share_bone4^2; %モデル×骨中横波音速
                c66(n,m,l)=rhob*velocity_share_bone6^2; %モデル×骨中横波音速axial方向
                
            end
        end
    end
end

for n=1:nx-1
    for m=1:ny-1
        for l=1:nz-1
            
                C44(n,m,l)=(c44(n,m,l)+c44(n+1,m,l)+c44(n,m+1,l)+c44(n,m,l+1)+...
                            c44(n+1,m+1,l)+c44(n,m+1,l+1)+c44(n+1,m,l+1)+c44(n+1,m+1,l+1))/8;
                C55(n,m,l)=(c55(n,m,l)+c55(n+1,m,l)+c55(n,m+1,l)+c55(n,m,l+1)+...
                            c55(n+1,m+1,l)+c55(n,m+1,l+1)+c55(n+1,m,l+1)+c55(n+1,m+1,l+1))/8;
                C66(n,m,l)=(c66(n,m,l)+c66(n+1,m,l)+c66(n,m+1,l)+c66(n,m,l+1)+...
                            c66(n+1,m+1,l)+c66(n,m+1,l+1)+c66(n+1,m,l+1)+c66(n+1,m+1,l+1))/8;               
                        
        end
    end
end

%%
%%%%%%%%%%%%%%弾性定数 C11,C22,C33,C12 の設定
C11=zeros(nx,ny,nz);
C22=zeros(nx,ny,nz);
C33=zeros(nx,ny,nz);
C12=zeros(nx,ny,nz);

for n=1:nx
    for m=1:ny
        for l=1:nz
            if Model3d(n,m,l)==0
                C11(n,m,l)=Kwl;
                C22(n,m,l)=Kwl;
                C33(n,m,l)=Kwl;
                C12(n,m,l)=Kwl-2*c66(n,m,l);

            else
                C11(n,m,l)=rhob*velocity_longitudinal_bone1^2; %モデル×radial方向骨中音速
                C22(n,m,l)=rhob*velocity_longitudinal_bone1^2; %モデル×radial方向骨中音速
                C33(n,m,l)=rhob*velocity_longitudinal_bone3^2; %モデル×axial方向骨中音速
                C12(n,m,l)=C11(n,m,l)-2*c66(n,m,l);                    %
            end
        end
    end
end

%%%%%%%%%%%%%%弾性定数 C13,C23 の設定
C13=zeros(nx,ny,nz);
ratio_elastic_C12C13=1.0212; %C12の1.0212倍がC13

for n=1:nx
    for m=1:ny
        for l=1:nz
            if Model3d(n,m,l)==0
                C13(n,m,l)=Kwl-2*c66(n,m,l);    

            else
                C13(n,m,l)=C12(n,m,l)*ratio_elastic_C12C13;
            end
        end
    end
end

C23(:,:,:)=C13(:,:,:);

X=dt/dx;
Y=dt/dy;
Z=dt/dz;
D=X;

%%
%%%%%%%%%%%%%%吸収境界条件パラメータ設定

a1=1; %a1=1/cosθ θは入射角で、1次元の場合θはゼロ
b1=(a1*velocity_longitudinal_water*dt-dx)/(a1*velocity_longitudinal_water*dt+dx);
b2=b1;
b3=b1;
b4=b1;
b5=b1;
b6=b1;

d1=0.00001;
d2=d1;

hig1=zeros(6,ny,nz);  %上端 Txx
hig2=zeros(6,ny,nz);  %下端 Txx
hig3=zeros(nx,6,nz);  %左端 Txx
hig4=zeros(nx,6,nz);  %右端 Txx
hig5=zeros(nx,ny,6);  %奥端 Txx
hig6=zeros(nx,ny,6);  %手前端 Txx

hig7=zeros(6,ny,nz);  %上端 Tyy
hig8=zeros(6,ny,nz);  %下端 Tyy
hig9=zeros(nx,6,nz);  %左端 Tyy
hig10=zeros(nx,6,nz); %右端 Tyy
hig11=zeros(nx,ny,6); %奥端 Tyy
hig12=zeros(nx,ny,6); %手前端 Tyy

hig13=zeros(6,ny,nz); %上端 Tzz
hig14=zeros(6,ny,nz); %下端 Tzz
hig15=zeros(nx,6,nz); %左端 Tzz
hig16=zeros(nx,6,nz); %右端 Tzz
hig17=zeros(nx,ny,6); %奥端 Tzz
hig18=zeros(nx,ny,6); %手前端 Tzz

% 切り口
Txx1=zeros(nx,ny);
Txx2=zeros(nx,nz);  
Txx3=zeros(nx,nz);  
Txx4=zeros(nx,nz);

Tyy1=zeros(nx,ny);
Tyy2=zeros(nx,nz);
Tyy3=zeros(nx,nz);
Tyy4=zeros(nx,nz);

Tzz1=zeros(nx,ny);
Tzz2=zeros(nx,nz);
Tzz3=zeros(nx,nz);
Tzz4=zeros(nx,nz);


%%
centerx = 150;
centery = 150;
centerz = 150;

sensaradial= 0.0427;   %実際の長さ　ピクセル数ではない (140*dx)

xzure45=round((sensaradial/sqrt(2)/dx));   %センサ中心　45度回転でのX軸ずれ
yzure45=round((sensaradial/sqrt(2)/dy));   %センサ中心　45度回転でのY軸ずれ
slength45=round(5e-3/2/sqrt(2)/dx);   %45度のときのセンサ5mmの要素数の半分

xzure30=round(sensaradial*sqrt(3)/2);   %センサ中心　30度回転でのX軸ずれ
yzure30=round(sensaradial/2);                   %センサ中心　30度回転でのY軸ずれ
slength30=round(5e-3/2/dx/2);   %30度のときのセンサ5mmの要素数の半分
% % 
xzure15=round(sensaradial*(2+sqrt(3))/(sqrt(6)+sqrt(2))); %センサ中心　15度回転でのx軸ズレ
yzure15=round(sensaradial/(sqrt(6)+sqrt(2)));             %センサ中心　15度回転でのy軸ズレ
slength15=round((10e-3/(sqrt(6)+sqrt(2)))/dx/2);          %15度の時のセンサ5mmの要素数の半分  
% 
% xzure15_2=round(sensaradial*0.99);
% yzure15_2=round(sensaradial*0.13);
% slength15_2=9; 
% 
% 
% 
sensorradius  = round(5e-3/dx);   %センサ直径5mm


% for m=1:ny
%     for n=1:nz
%         x=m-centery;
%         y=n-centerz;
%         radius=sqrt(x^2+y^2);
%         if radius <=sensorradius
%             sensormodel(1,m,n)=1;
%         end
%     end
% end

%%


Center=zeros(1,nt);
CenterT=zeros(1,nt);

boxelX = zeros(1,80);
boxelY = zeros(1,80);
boxelZ = zeros(1,80);
boxelS = zeros(1,80);
 
%% FDTD
for s=1:nt
    if s <= N
        Txx(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
                  centery+140-round(0.5e-3/dx):centery+140+round(0.5e-3/dx),500)=wave2(1,s);%0度
        Tyy(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
                  centery+140-round(0.5e-3/dx):centery+140+round(0.5e-3/dx),500)=wave2(1,s);%0度
        Tzz(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
                  centery+140-round(0.5e-3/dx):centery+140+round(0.5e-3/dx),500)=wave2(1,s);%0度
              
        Txx(centerx-140-round(0.5e-3/dx):centerx-140+round(0.5e-3/dx),...
                 centery-round(sensorradius/2):centery+round(sensorradius/2),500)=wave2(1,s);%90度
        Tyy(centerx-140-round(0.5e-3/dx):centerx-140+round(0.5e-3/dx),...
                 centery-round(sensorradius/2):centery+round(sensorradius/2),500)=wave2(1,s);%90度
        Tzz(centerx-140-round(0.5e-3/dx):centerx-140+round(0.5e-3/dx),...
                 centery-round(sensorradius/2):centery+round(sensorradius/2),500)=wave2(1,s);%90度
    else
        Txx(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
                  centery+140-round(0.5e-3/dx):centery+140+round(0.5e-3/dx),500)=0;%0度
        Tyy(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
                  centery+140-round(0.5e-3/dx):centery+140+round(0.5e-3/dx),500)=0;%0度
        Tzz(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
                  centery+140-round(0.5e-3/dx):centery+140+round(0.5e-3/dx),500)=0;%0度
              
        Txx(centerx-140-round(0.5e-3/dx):centerx-140+round(0.5e-3/dx),...
                 centery-round(sensorradius/2):centery+round(sensorradius/2),500)=0;%90度
        Tyy(centerx-140-round(0.5e-3/dx):centerx-140+round(0.5e-3/dx),...
                 centery-round(sensorradius/2):centery+round(sensorradius/2),500)=0;%90度
        Tzz(centerx-140-round(0.5e-3/dx):centerx-140+round(0.5e-3/dx),...
                 centery-round(sensorradius/2):centery+round(sensorradius/2),500)=0;%90度
    end
    
    if s <= N
        for mm = -slength45:slength45
                    Txx(centerx+xzure45+mm-round(0.5e-3/dx):centerx+xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45-mm-round(0.5e-3/dx):centery+yzure45-mm+round(0.5e-3/dx),500)=wave2(1,s);%315度
                    Tyy(centerx+xzure45+mm-round(0.5e-3/dx):centerx+xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45-mm-round(0.5e-3/dx):centery+yzure45-mm+round(0.5e-3/dx),500)=wave2(1,s);%315度
                    Tzz(centerx+xzure45+mm-round(0.5e-3/dx):centerx+xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45-mm-round(0.5e-3/dx):centery+yzure45-mm+round(0.5e-3/dx),500)=wave2(1,s);%315度
                    
                    Txx(centerx-xzure45+mm-round(0.5e-3/dx):centerx-xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45+mm-round(0.5e-3/dx):centery+yzure45+mm+round(0.5e-3/dx),500)=wave2(1,s);%45度
                    Tyy(centerx-xzure45+mm-round(0.5e-3/dx):centerx-xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45+mm-round(0.5e-3/dx):centery+yzure45+mm+round(0.5e-3/dx),500)=wave2(1,s);%45度
                    Tzz(centerx-xzure45+mm-round(0.5e-3/dx):centerx-xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45+mm-round(0.5e-3/dx):centery+yzure45+mm+round(0.5e-3/dx),500)=wave2(1,s);%45度 
        end
    else
        for mm = -slength45:slength45
                    Txx(centerx+xzure45+mm-round(0.5e-3/dx):centerx+xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45-mm-round(0.5e-3/dx):centery+yzure45-mm+round(0.5e-3/dx),500)=0;%315度
                    Tyy(centerx+xzure45+mm-round(0.5e-3/dx):centerx+xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45-mm-round(0.5e-3/dx):centery+yzure45-mm+round(0.5e-3/dx),500)=0;%315度
                    Tzz(centerx+xzure45+mm-round(0.5e-3/dx):centerx+xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45-mm-round(0.5e-3/dx):centery+yzure45-mm+round(0.5e-3/dx),500)=0;%315度
                    
                    Txx(centerx-xzure45+mm-round(0.5e-3/dx):centerx-xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45+mm-round(0.5e-3/dx):centery+yzure45+mm+round(0.5e-3/dx),500)=0;%45度
                    Tyy(centerx-xzure45+mm-round(0.5e-3/dx):centerx-xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45+mm-round(0.5e-3/dx):centery+yzure45+mm+round(0.5e-3/dx),500)=0;%45度
                    Tzz(centerx-xzure45+mm-round(0.5e-3/dx):centerx-xzure45+mm+round(0.5e-3/dx),...
                        centery+yzure45+mm-round(0.5e-3/dx):centery+yzure45+mm+round(0.5e-3/dx),500)=0;%45度 
        end
    end

%      if s<=N
%          for nn= 1:16
%              Txx(centerx+round((28/16)*(nn-1))-115,centery-round(sensorradius/2)...
%                  :centery+round(sensorradius/2),500-16+nn)=wave2(1,s);%90度
%              Tyy(centerx+round((28/16)*(nn-1))-115,centery-round(sensorradius/2)...
%                  :centery+round(sensorradius/2),500-16+nn)=wave2(1,s);%90度
%              Tzz(centerx+round((28/16)*(nn-1))-115,centery-round(sensorradius/2)...
%                  :centery+round(sensorradius/2),500-16+nn)=wave2(1,s);%90度
%              
% %              Txx(centerx+round((28/16)*(nn-1))+140,centery+round((28/16)*(nn-1))...
% %                  -round(sensorradius/2):centery+round((28/16)*(nn-1))+round(sensorradius/2),500-nn)=wave2(1,s);%270度
% %              Tyy(centerx+round((28/16)*(nn-1))+140,centery+round((28/16)*(nn-1))...
% %                  -round(sensorradius/2):centery+round((28/16)*(nn-1))+round(sensorradius/2),500-nn)=wave2(1,s);%270度
% %              Tzz(centerx+round((28/16)*(nn-1))+140,centery+round((28/16)*(nn-1))...
% %                  -round(sensorradius/2):centery+round((28/16)*(nn-1))+round(sensorradius/2),500-nn)=wave2(1,s);%270度
%              
%              Txx(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
%                  centery+round((28/16)*(nn-1))+115,500-16+nn)=wave2(1,s);%0度
%             Tyy(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
%                  centery+round((28/16)*(nn-1))+115,500-16+nn)=wave2(1,s);%0度
%             Tzz(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
%                  centery+round((28/16)*(nn-1))+115,500-16+nn)=wave2(1,s);%0度
%              
% %              Txx(centerx+round((28/16)*(nn-1))-round(sensorradius/2):centerx...
% %                  +round((28/16)*(nn-1))+round(sensorradius/2),centery+round((28/16)*(nn-1))-140,500-nn)=wave2(1,s);%180度
% %             Tyy(centerx+round((28/16)*(nn-1))-round(sensorradius/2):centerx...
% %                  +round((28/16)*(nn-1))+round(sensorradius/2),centery+round((28/16)*(nn-1))-140,500-nn)=wave2(1,s);%180度
% %             Tzz(centerx+round((28/16)*(nn-1))-round(sensorradius/2):centerx...
% %                  +round((28/16)*(nn-1))+round(sensorradius/2),centery+round((28/16)*(nn-1))-140,500-nn)=wave2(1,s);%180度
%          end
% 
%          
%             for nn= 1:16
%                 for mm = -slength45:slength45
%                     Txx(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45-mm,500-16+nn)=wave2(1,s);%315度
%                     Tyy(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45-mm,500-16+nn)=wave2(1,s);%315度
%                     Tzz(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45-mm,500-16+nn)=wave2(1,s);%315度
%                     
%                     Txx(centerx+round((28/16)*(nn-1))-xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-16+nn)=wave2(1,s);%45度
%                     Tyy(centerx+round((28/16)*(nn-1))-xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-16+nn)=wave2(1,s);%45度
%                     Tzz(centerx+round((28/16)*(nn-1))-xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-16+nn)=wave2(1,s);%45度
%                     
% %                     Txx(centerx+round((28/16)*(nn-1))-xzure45-mm,centery+round((28/16)*(nn-1))-yzure45+mm,500-nn)=wave2(1,s);%135度
% %                     Tyy(centerx+round((28/16)*(nn-1))-xzure45-mm,centery+round((28/16)*(nn-1))-yzure45+mm,500-nn)=wave2(1,s);%135度
% %                     Tzz(centerx+round((28/16)*(nn-1))-xzure45-mm,centery+round((28/16)*(nn-1))-yzure45+mm,500-nn)=wave2(1,s);%135度
%                     
% %                     Txx(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-nn)=wave2(1,s);%225度
% %                     Tyy(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-nn)=wave2(1,s);%225度
% %                     Tzz(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-nn)=wave2(1,s);%225度
%                
% 
%                 end
%             end
%             
%      else
%          for nn= 1:16
%                 for mm = -slength45:slength45
%                     Txx(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45-mm,500-16+nn)=0;%315度
%                     Tyy(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45-mm,500-16+nn)=0;%315度
%                     Tzz(centerx+round((28/16)*(nn-1))+xzure45+mm,centery+round((28/16)*(nn-1))+yzure45-mm,500-16+nn)=0;%315度
%                     
%                     Txx(centerx+round((28/16)*(nn-1))-xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-16+nn)=0;%45度
%                     Tyy(centerx+round((28/16)*(nn-1))-xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-16+nn)=0;%45度
%                     Tzz(centerx+round((28/16)*(nn-1))-xzure45+mm,centery+round((28/16)*(nn-1))+yzure45+mm,500-16+nn)=0;%45度
%                     
% %                     Txx(centerx+round((57/33)*(nn-1))-xzure45-mm,centery+round((57/33)*(nn-1))-yzure45+mm,500-nn)=0;%135度
% %                     Tyy(centerx+round((57/33)*(nn-1))-xzure45-mm,centery+round((57/33)*(nn-1))-yzure45+mm,500-nn)=0;%135度
% %                     Tzz(centerx+round((57/33)*(nn-1))-xzure45-mm,centery+round((57/33)*(nn-1))-yzure45+mm,500-nn)=0;%135度
%                     
% %                     Txx(centerx+round((57/33)*(nn-1))+xzure45+mm,centery+round((57/33)*(nn-1))+yzure45+mm,500-nn)=0;%225度
% %                     Tyy(centerx+round((57/33)*(nn-1))+xzure45+mm,centery+round((57/33)*(nn-1))+yzure45+mm,500-nn)=0;%225度
% %                     Tzz(centerx+round((57/33)*(nn-1))+xzure45+mm,centery+round((57/33)*(nn-1))+yzure45+mm,500-nn)=0;%225度
%                
% 
%                 end
%          end
%          
%          for nn= 1:16
%              Txx(centerx+round((28/16)*(nn-1))-115,centery-round(sensorradius/2)...
%                  :centery+round(sensorradius/2),500-16+nn)=0;%90度
%              Tyy(centerx+round((28/16)*(nn-1))-115,centery-round(sensorradius/2)...
%                  :centery+round(sensorradius/2),500-16+nn)=0;%90度
%              Tzz(centerx+round((28/16)*(nn-1))-115,centery-round(sensorradius/2)...
%                  :centery+round(sensorradius/2),500-16+nn)=0;%90度
%              
% %              Txx(centerx+round((28/16)*(nn-1))+140,centery+round((28/16)*(nn-1))...
% %                  -round(sensorradius/2):centery+round((28/16)*(nn-1))+round(sensorradius/2),500-nn)=0;%270度
% %              Tyy(centerx+round((28/16)*(nn-1))+140,centery+round((28/16)*(nn-1))...
% %                  -round(sensorradius/2):centery+round((28/16)*(nn-1))+round(sensorradius/2),500-nn)=0;%270度
% %              Tzz(centerx+round((28/16)*(nn-1))+140,centery+round((28/16)*(nn-1))...
% %                  -round(sensorradius/2):centery+round((28/16)*(nn-1))+round(sensorradius/2),500-nn)=0;%270度
%              
%              Txx(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
%                  centery+round((28/16)*(nn-1))+115,500-16+nn)=0;%0度
%             Tyy(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
%                  centery+round((28/16)*(nn-1))+115,500-16+nn)=0;%0度
%             Tzz(centerx-round(sensorradius/2):centerx+round(sensorradius/2),...
%                  centery+round((28/16)*(nn-1))+115,500-16+nn)=0;%0度
%              
% %              Txx(centerx+round((28/16)*(nn-1))-round(sensorradius/2):centerx...
% %                  +round((28/16)*(nn-1))+round(sensorradius/2),centery+round((28/16)*(nn-1))-140,500-nn)=0;%180度
% %             Tyy(centerx+round((28/16)*(nn-1))-round(sensorradius/2):centerx...
% %                  +round((28/16)*(nn-1))+round(sensorradius/2),centery+round((28/16)*(nn-1))-140,500-nn)=0;%180度
% %             Tzz(centerx+round((28/16)*(nn-1))-round(sensorradius/2):centerx...
% %                  +round((28/16)*(nn-1))+round(sensorradius/2),centery+round((28/16)*(nn-1))-140,500-nn)=0;%180度
%          end
%          
%      end
         
         


    
%FDTDの式 
Ux(2:nx,:,:)=Ux(2:nx,:,:)+D/density_x.*((Txx(2:nx,:,:)-Txx(1:nx-1,:,:))+...
                                            (Txy(2:nx,2:ny+1,1:nz)-Txy(2:nx,1:ny,1:nz))+...
                                            (Tzx(2:nx,1:ny,2:nz+1)-Tzx(2:nx,1:ny,1:nz)));

Uy(:,2:ny,:)=Uy(:,2:ny,:)+D/density_y.*((Txy(2:nx+1,2:ny,1:nz)-Txy(1:nx,2:ny,1:nz))+...
                                            (Tyy(:,2:ny,:)-Tyy(:,1:ny-1,:))+...
                                            (Tyz(1:nx,2:ny,2:nz+1)-Tyz(1:nx,2:ny,1:nz)));

Uz(:,:,2:nz)=Uz(:,:,2:nz)+D/density_z.*((Tzx(2:nx+1,1:ny,2:nz)-Tzx(1:nx,1:ny,2:nz))+...
                                            (Tyz(1:nx,2:ny+1,2:nz)-Tyz(1:nx,1:ny,2:nz))...
                                            +(Tzz(:,:,2:nz)-Tzz(:,:,1:nz-1)));


Txx(:,:,:)=Txx(:,:,:)+C11*X.*(Ux(2:nx+1,:,:)-Ux(1:nx,:,:))+...
                                    C12*Y.*(Uy(:,2:ny+1,:)-Uy(:,1:ny,:))+...
                                    C13*Z.*(Uz(:,:,2:nz+1)-Uz(:,:,1:nz));

Tyy(:,:,:)=Tyy(:,:,:)+C12*X.*(Ux(2:nx+1,:,:)-Ux(1:nx,:,:))+...
                                    C22*Y.*(Uy(:,2:ny+1,:)-Uy(:,1:ny,:))+...
                                    C23*Z.*(Uz(:,:,2:nz+1)-Uz(:,:,1:nz));

Tzz(:,:,:)=Tzz(:,:,:)+C13*X.*(Ux(2:nx+1,:,:)-Ux(1:nx,:,:))+...
                                    C23*Y.*(Uy(:,2:ny+1,:)-Uy(:,1:ny,:))+...
                                    C33*Z.*(Uz(:,:,2:nz+1)-Uz(:,:,1:nz));

Tyz(2:nx,2:ny,2:nz)=Tyz(2:nx,2:ny,2:nz)+C44*D.*((Uz(2:nx,2:ny,2:nz)-Uz(2:nx,1:ny-1,2:nz))+...
    (Uy(2:nx,2:ny,2:nz)-Uy(2:nx,2:ny,1:nz-1)));

Tzx(2:nx,2:ny,2:nz)=Tzx(2:nx,2:ny,2:nz)+C55*D.*((Uz(2:nx,2:ny,2:nz)-Uz(1:nx-1,2:ny,2:nz))+...
    (Ux(2:nx,2:ny,2:nz)-Ux(2:nx,2:ny,1:nz-1)));

Txy(2:nx,2:ny,2:nz)=Txy(2:nx,2:ny,2:nz)+C66*D.*((Ux(2:nx,2:ny,2:nz)-Ux(2:nx,1:ny-1,2:nz))+...    
    (Uy(2:nx,2:ny,2:nz)-Uy(1:nx-1,2:ny,2:nz)));

% Txy(2:nx,2:ny,1:nz)=Txy(2:nx,2:ny,1:nz)+dt*0*((Ux(2:nx,2:ny,1:nz)-Ux(2:nx,1:ny-1,1:nz))/dy+...
%     (Uy(2:nx,2:ny,1:nz)-Uy(1:nx-1,2:ny,1:nz))/dx);
% 
% Tyz(1:nx,2:ny,2:nz)=Tyz(1:nx,2:ny,2:nz)+dt*0*((Uz(1:nx,2:ny,2:nz)-Uz(1:nx,1:ny-1,2:nz))/dy+...
%     (Uy(1:nx,2:ny,2:nz)-Uy(1:nx,2:ny,1:nz-1))/dz);
% 
% Tzx(2:nx,1:ny,2:nz)=Tzx(2:nx,1:ny,2:nz)+dt*0*((Uz(2:nx,1:ny,2:nz)-Uz(1:nx-1,1:ny,2:nz))/dx+...
%     (Ux(2:nx,1:ny,2:nz)-Ux(2:nx,1:ny,1:nz-1))/dz);
        
%吸収境界条件式
     Txx(:,ny,:)=(b3+b4)*(Txx(:,ny-1,:)-hig4(:,3,:))-b3*b4*(Txx(:,ny-2,:)-2*hig4(:,2,:)+hig4(:,6,:))...
         -(b3*(1-d2)+b4*(1-d1))*(hig4(:,1,:)-hig4(:,5,:))+((1-d1)+(1-d2))*hig4(:,2,:)...
         -(1-d1)*(1-d2)*hig4(:,4,:); %右端
     Tyy(:,ny,:)=(b1+b2)*(Tyy(:,ny-1,:)-hig10(:,3,:))-b1*b2*(Tyy(:,ny-2,:)-2*hig10(:,2,:)+hig10(:,6,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig10(:,1,:)-hig10(:,5,:))+((1-d1)+(1-d2))*hig10(:,2,:)...
         -(1-d1)*(1-d2)*hig10(:,4,:); %右端
     Tzz(:,ny,:)=(b1+b2)*(Tzz(:,ny-1,:)-hig16(:,3,:))-b1*b2*(Tzz(:,ny-2,:)-2*hig16(:,2,:)+hig16(:,6,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig16(:,1,:)-hig16(:,5,:))+((1-d1)+(1-d2))*hig16(:,2,:)...
         -(1-d1)*(1-d2)*hig16(:,4,:); %右端
    
     Txx(:,1,:)=(b3+b4)*(Txx(:,2,:)-hig3(:,1,:))-b3*b4*(Txx(:,3,:)-2*hig3(:,2,:)+hig3(:,4,:))...
         -(b3*(1-d2)+b4*(1-d1))*(hig3(:,3,:)-hig3(:,5,:))+((1-d1)+(1-d2))*hig3(:,2,:)...
         -(1-d1)*(1-d2)*hig3(:,6,:); %左端
     Tyy(:,1,:)=(b1+b2)*(Tyy(:,2,:)-hig9(:,1,:))-b1*b2*(Tyy(:,3,:)-2*hig9(:,2,:)+hig9(:,4,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig9(:,3,:)-hig9(:,5,:))+((1-d1)+(1-d2))*hig9(:,2,:)...
         -(1-d1)*(1-d2)*hig9(:,6,:); %左端
     Tzz(:,1,:)=(b1+b2)*(Tzz(:,2,:)-hig15(:,1,:))-b1*b2*(Tzz(:,3,:)-2*hig15(:,2,:)+hig15(:,4,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig15(:,3,:)-hig15(:,5,:))+((1-d1)+(1-d2))*hig15(:,2,:)...
         -(1-d1)*(1-d2)*hig15(:,6,:); %左端
     
     Txx(1,:,:)=(b1+b2)*(Txx(2,:,:)-hig1(1,:,:))-b1*b2*(Txx(3,:,:)-2*hig1(2,:,:)+hig1(4,:,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig1(3,:,:)-hig1(5,:,:))+((1-d1)+(1-d2))*hig1(2,:,:)...
         -(1-d1)*(1-d2)*hig1(6,:,:); %上端
     Tyy(1,:,:)=(b1+b2)*(Tyy(2,:,:)-hig7(1,:,:))-b1*b2*(Tyy(3,:,:)-2*hig7(2,:,:)+hig7(4,:,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig7(3,:,:)-hig7(5,:,:))+((1-d1)+(1-d2))*hig7(2,:,:)...
         -(1-d1)*(1-d2)*hig7(6,:,:); %上端
     Tzz(1,:,:)=(b1+b2)*(Tzz(2,:,:)-hig13(1,:,:))-b1*b2*(Tzz(3,:,:)-2*hig13(2,:,:)+hig13(4,:,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig13(3,:,:)-hig13(5,:,:))+((1-d1)+(1-d2))*hig13(2,:,:)...
         -(1-d1)*(1-d2)*hig13(6,:,:); %上端
    
     Txx(nx,:,:)=(b1+b2)*(Txx(nx-1,:,:)-hig2(3,:,:))-b1*b2*(Txx(nx-2,:,:)-2*hig2(2,:,:)+hig2(6,:,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig2(1,:,:)-hig2(5,:,:))+((1-d1)+(1-d2))*hig2(2,:,:)...
         -(1-d1)*(1-d2)*hig2(4,:,:); %下端
     Tyy(nx,:,:)=(b1+b2)*(Tyy(nx-1,:,:)-hig8(3,:,:))-b1*b2*(Tyy(nx-2,:,:)-2*hig8(2,:,:)+hig8(6,:,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig8(1,:,:)-hig8(5,:,:))+((1-d1)+(1-d2))*hig8(2,:,:)...
         -(1-d1)*(1-d2)*hig8(4,:,:); %下端
     Tzz(nx,:,:)=(b1+b2)*(Tzz(nx-1,:,:)-hig14(3,:,:))-b1*b2*(Tzz(nx-2,:,:)-2*hig14(2,:,:)+hig14(6,:,:))...
         -(b1*(1-d2)+b2*(1-d1))*(hig14(1,:,:)-hig14(5,:,:))+((1-d1)+(1-d2))*hig14(2,:,:)...
         -(1-d1)*(1-d2)*hig14(4,:,:); %下端

     Txx(:,:,1)=(b1+b2)*(Txx(:,:,2)-hig6(:,:,1))-b1*b2*(Txx(:,:,3)-2*hig6(:,:,2)+hig6(:,:,4))...
         -(b1*(1-d2)+b2*(1-d1))*(hig6(:,:,3)-hig6(:,:,5))+((1-d1)+(1-d2))*hig6(:,:,2)...
         -(1-d1)*(1-d2)*hig6(:,:,6); %手前
     Tyy(1:nx,1:ny,1)=(b1+b2)*(Tyy(1:nx,1:ny,2)-hig12(:,:,1))-b1*b2*(Tyy(1:nx,1:ny,3)-2*hig12(:,:,2)+hig12(:,:,4))...
         -(b1*(1-d2)+b2*(1-d1))*(hig12(:,:,3)-hig12(:,:,5))+((1-d1)+(1-d2))*hig12(:,:,2)...
         -(1-d1)*(1-d2)*hig12(:,:,6); %手前
     Tzz(1:nx,1:ny,1)=(b1+b2)*(Tzz(1:nx,1:ny,2)-hig18(:,:,1))-b1*b2*(Tzz(1:nx,1:ny,3)-2*hig18(:,:,2)+hig18(:,:,4))...
         -(b1*(1-d2)+b2*(1-d1))*(hig18(:,:,3)-hig18(:,:,5))+((1-d1)+(1-d2))*hig18(:,:,2)...
         -(1-d1)*(1-d2)*hig18(:,:,6); %手前
    
     Txx(:,:,nz)=(b1+b2)*(Txx(:,:,nz-1)-hig5(:,:,3))-b1*b2*(Txx(:,:,nz-2)-2*hig5(:,:,2)+hig5(:,:,6))...
         -(b1*(1-d2)+b2*(1-d1))*(hig5(:,:,1)-hig5(:,:,5))+((1-d1)+(1-d2))*hig5(:,:,2)...
         -(1-d1)*(1-d2)*hig5(:,:,4); %奥
     Tyy(:,:,nz)=(b1+b2)*(Tyy(:,:,nz-1)-hig11(:,:,3))-b1*b2*(Tyy(:,:,nz-2)-2*hig11(:,:,2)+hig11(:,:,6))...
         -(b1*(1-d2)+b2*(1-d1))*(hig11(:,:,1)-hig11(:,:,5))+((1-d1)+(1-d2))*hig11(:,:,2)...
         -(1-d1)*(1-d2)*hig11(:,:,4); %奥
     Tzz(:,:,nz)=(b1+b2)*(Tzz(:,:,nz-1)-hig17(:,:,3))-b1*b2*(Tzz(:,:,nz-2)-2*hig17(:,:,2)+hig17(:,:,6))...
         -(b1*(1-d2)+b2*(1-d1))*(hig17(:,:,1)-hig17(:,:,5))+((1-d1)+(1-d2))*hig17(:,:,2)...
         -(1-d1)*(1-d2)*hig17(:,:,4); %奥
     

     
     
     hig1(4:6,:,:)=hig1(1:3,:,:); %2秒前上端
     hig1(1:3,:,:)=Txx(1:3,:,:);  %1秒前上端
     
     hig2(4:6,:,:)=hig2(1:3,:,:); %2秒前下端
     hig2(1:3,:,:)=Txx(nx-2:nx,:,:);  %1秒前下端
     
     hig3(:,4:6,:)=hig3(:,1:3,:); %2秒前左端
     hig3(:,1:3,:)=Txx(:,1:3,:);  %1秒前左端
     
     hig4(:,4:6,:)=hig4(:,1:3,:); %2秒前右端
     hig4(:,1:3,:)=Txx(:,ny-2:ny,:);  %1秒前右端
     
     hig6(:,:,4:6)=hig6(:,:,1:3); %2秒前手前
     hig6(:,:,1:3)=Txx(:,:,1:3);  %1秒前手前
     
     hig5(:,:,4:6)=hig5(:,:,1:3); %2秒前奥
     hig5(:,:,1:3)=Txx(:,:,nz-2:nz);  %1秒前奥

     hig7(4:6,:,:)=hig7(1:3,:,:); %2秒前上端
     hig7(1:3,:,:)=Tyy(1:3,:,:);  %1秒前上端
     
     hig8(4:6,:,:)=hig8(1:3,:,:); %2秒前下端
     hig8(1:3,:,:)=Tyy(nx-2:nx,:,:);  %1秒前下端
     
     hig9(:,4:6,:)=hig9(:,1:3,:); %2秒前左端
     hig9(:,1:3,:)=Tyy(:,1:3,:);  %1秒前左端
     
     hig10(:,4:6,:)=hig10(:,1:3,:); %2秒前右端
     hig10(:,1:3,:)=Tyy(:,ny-2:ny,:);  %1秒前右端
     
     hig12(:,:,4:6)=hig12(:,:,1:3); %2秒前手前
     hig12(:,:,1:3)=Tyy(:,:,1:3);  %1秒前手前
     
     hig11(:,:,4:6)=hig11(:,:,1:3); %2秒前奥
     hig11(:,:,1:3)=Tyy(:,:,nz-2:nz);  %1秒前奥

     hig13(4:6,:,:)=hig13(1:3,:,:); %2秒前上端
     hig13(1:3,:,:)=Tzz(1:3,:,:);  %1秒前上端
     
     hig14(4:6,:,:)=hig14(1:3,:,:); %2秒前下端
     hig14(1:3,:,:)=Tzz(nx-2:nx,:,:);  %1秒前下端
     
     hig15(:,4:6,:)=hig15(:,1:3,:); %2秒前左端
     hig15(:,1:3,:)=Tzz(:,1:3,:);  %1秒前左端
     
     hig16(:,4:6,:)=hig16(:,1:3,:); %2秒前右端
     hig16(:,1:3,:)=Tzz(:,ny-2:ny,:);  %1秒前右端
     
     hig18(:,:,4:6)=hig18(:,:,1:3); %2秒前手前
     hig18(:,:,1:3)=Tzz(:,:,1:3);  %1秒前手前
     
     hig17(:,:,4:6)=hig17(:,:,1:3); %2秒前奥
     hig17(:,:,1:3)=Tzz(:,:,nz-2:nz);  %1秒前奥
     
     Txx1(:,:)=squeeze(Txx(:,:,centerx)); %XY切り口
     Txx2(:,:)=squeeze(Txx(:,160-40,:)); %XZ切り口
     Txx3(:,:)=squeeze(Txx(:,160,:)); %XZ切り口
     Txx4(:,:)=squeeze(Txx(:,160+40,:)); %XZ切り口
     Txx5(:,:)=squeeze(Txx(:,:,100)); %XY切り口
     Txx6(:,:)=squeeze(Txx(:,:,500)); %XY切り口
     
     Tyy1(:,:)=squeeze(Tyy(:,:,centerx)); %XY切り口
     Tyy2(:,:)=squeeze(Tyy(:,160-40,:)); %XZ切り口
     Tyy3(:,:)=squeeze(Tyy(:,160,:)); %XZ切り口
     Tyy4(:,:)=squeeze(Tyy(:,160+40,:)); %XZ切り口
     Tyy5(:,:)=squeeze(Tyy(:,:,100)); %XY切り口
     Tyy6(:,:)=squeeze(Tyy(:,:,500)); %XY切り口
     
     Tzz1(:,:)=squeeze(Tzz(:,:,centerx)); %XY切り口
     Tzz2(:,:)=squeeze(Tzz(:,160-40,:)); %XZ切り口
     Tzz3(:,:)=squeeze(Tzz(:,160,:)); %XZ切り口
     Tzz4(:,:)=squeeze(Tzz(:,160+40,:)); %XZ切り口
     Tzz5(:,:)=squeeze(Tzz(:,:,100)); %XY切り口
     Tzz6(:,:)=squeeze(Tzz(:,:,500)); %XY切り口
    
     
     
%     音場保存
    T1=sqrt(Txx1.^2+Tyy1.^2+Tzz1.^2); 
    T2=sqrt(Txx2.^2+Tyy2.^2+Tzz2.^2);
    T3=sqrt(Txx3.^2+Tyy3.^2+Tzz3.^2);
    T4=sqrt(Txx4.^2+Tyy4.^2+Tzz4.^2);
    T5=sqrt(Txx5.^2+Tyy5.^2+Tzz5.^2);
    T6=sqrt(Txx6.^2+Tyy6.^2+Tzz6.^2);
    
    Center(s)=Txx(centerx,centery,centerz);
    CenterT(s)=sqrt(Txx(centerx,centery,centerz)^2+Tyy(centerx,centery,centerz)^2+Tzz(centerx,centery,centerz)^2);

    M=getframe; 
%     imagesc(T1)

    if s<=30
        subplot(2,2,1);imagesc(T6);
        axis equal
        caxis([0 0.05])
        
        subplot(2,2,2);imagesc(T2);
        axis equal
        caxis([0 0.05])  
        
        subplot(2,2,3);imagesc(T3);
        axis equal
        caxis([0 0.05])  
        
        subplot(2,2,4);imagesc(T4);
        axis equal
        caxis([0 0.05])
        
    end
    
    if mod(s,100)==0
        subplot(2,2,1);imagesc(T1);
        axis equal
        caxis([0 0.05])
        
        subplot(2,2,2);imagesc(T2);
        axis equal
        caxis([0 0.05])  
        
        subplot(2,2,3);imagesc(T3);
        axis equal
        caxis([0 0.05])  
        
        subplot(2,2,4);imagesc(T4);
        axis equal
        caxis([0 0.05])
        
    end

    
    if mod(s,50) ==0
    csvwrite(sprintf('T1_%04d.csv',s),T1)
    csvwrite(sprintf('T2_%04d.csv',s),T2)
    csvwrite(sprintf('T3_%04d.csv',s),T3)
    csvwrite(sprintf('T4_%04d.csv',s),T4)
    csvwrite(sprintf('T5_%04d.csv',s),T5)
    csvwrite(sprintf('T6_%04d.csv',s),T6)
    
    countx=0;
    sumx=0;
    for z = 1:250
        for gyo = 1:nx
            for retu = 1:110
                XX = sqrt(squeeze(Txx(gyo,retu,z)).^2+squeeze(Tyy(gyo,retu,z)).^2+squeeze(Tzz(gyo,retu,z)).^2);
                if XX~=0
                    sumx = sumx + XX;
                    countx = countx + 1;
                end
            end
        end
    end
    boxelX(s/50) = sumx/countx;
    
    county=0;
    sumy=0;
    for z = 1:250
        for gyo = 1:nx
            for retu = 111:170
                YY = sqrt(squeeze(Txx(gyo,retu,z)).^2+squeeze(Tyy(gyo,retu,z)).^2+squeeze(Tzz(gyo,retu,z)).^2);
                if YY~=0
                    sumy = sumy + YY;
                    county = county + 1;
                end
            end
        end
    end
    boxelY(s/50) = sumy/county;
    
    countz=0;
    sumz=0;
    for z = 1:250
        for gyo = 1:nx
            for retu = 171:ny
                ZZ = sqrt(squeeze(Txx(gyo,retu,z)).^2+squeeze(Tyy(gyo,retu,z)).^2+squeeze(Tzz(gyo,retu,z)).^2);
                if ZZ~=0
                    sumz = sumz + ZZ;
                    countz = countz + 1;
                end
            end
        end
    end
    boxelZ(s/50) = sumz/countz;
    
    counts=0;
    sums=0;
    for z = 251:400
        for gyo = 1:nx
            for retu = 1:110
                SS = sqrt(squeeze(Txx(gyo,retu,z)).^2+squeeze(Tyy(gyo,retu,z)).^2+squeeze(Tzz(gyo,retu,z)).^2);
                if SS~=0
                    sums = sums + SS;
                    counts = counts + 1;
                end
            end
        end
    end
    boxelS(s/50) = sums/counts;
    
    end
    
    clc
    s
end

  csvwrite(('boxelX.csv'),boxelX(:))
  csvwrite(('boxelY.csv'),boxelY(:))
  csvwrite(('boxelZ.csv'),boxelZ(:))
  csvwrite(('boxelS.csv'),boxelS(:))
  
