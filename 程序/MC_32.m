%%% 3(2) %%%
clear,close,clc;
format compact;
warning off;
w=1.53; %导弹偏航和俯仰角速度 rad/s
dC=9.5e3; %末段飞行距离
N=3000;
n=50;
num=50;    %精确发射角区间内的仿真次数
RE=6300; %地球半径 km
Pa=zeros(3,N); %航母位置
va=32*1.852/3.6; %航母航行速度 m/s
X=1000*distance(120.5,27.5,123.75,25.65);  %初始距离 m
k1=1000*RE*abs(pi/180*(27.5-25.65));    %初始位置南北方向距离（沿经线）m
k2=1000*RE*cos(pi/180*25.65)*(pi/180*(123.75-120.5));    %初始位置东西方向距离（沿纬线）m
Xa=sqrt(k1^2+k2^2); %初始距离的修正值 m
Pa(:,1)=[k1;k2;0];
lam=1.8;    %推重比（以RD-0410核火箭发动机为例）
mju=0.86;   %助推器占总重比（以C-301超音速反舰导弹为例）
mA0=1.5e3;  %发射段导弹初始总重（含助推器）
mB=mA0*(1-mju); %中段导弹总重
mC=mB;  %末段导弹总重
g=9.8;  %重力加速度
v0=500; %发射速度
C0=0.1747;  %零升阻力系数
Ci=2.52e-3; %诱导阻力系数
r=1.293;	%干燥空气密度
Dm=0.4; %导弹直径（以黄蜂Ⅲ反舰导弹为例）
S=pi/4*Dm^2;    %导弹迎风面积
alpha=0.5*(C0+Ci)*r*S;  %空气阻力与速度平方的比例系数
Lac=335;    %航母长度
Wac=77;    %航母宽度（以美国福特级航空母舰为例）
dt=1;   %信号延时的均值假设为1s
ss=0;
for i=2:N
    Pa(1,i)=double(k1)+double(va*(i-1));
    Pa(2,i)=k2;
    Pa(3,i)=0;
end
ph=zeros(1,n);  %杀伤率
pk=zeros(1,n);  %命中率
for u=linspace(5.7,6.7,n)
    ss=ss+1;
    nh=0;
    for theta=linspace(0.72,0.81,num)
    Pm=zeros(3,N); %导弹位置
    vm=zeros(3,N); %导弹沿x,y,z轴分速度
    vh=zeros(3,N); %导弹水平分速度
    v=zeros(1,N);   %导弹合速度
    s=zeros(1,N); %导弹水平航行路程
    L=zeros(1,N); %导弹与航母的动态距离
    sig=zeros(1,N); %导弹速度在海平面的投影与经线夹角
    sig(1)=atan(k2/k1);
    Pm(:,1)=[0;0;0];
    vm(:,1)=[v0*cos(theta)*cos(sig(1));v0*cos(theta)*sin(sig(1));v0*sin(theta)];
    vh(1)=v0*sin(theta);
    v(1)=v0;
    s(1)=0;
    L(1)=Xa;
    tA=0;tBp=0;tB=0;tC=0;hmax=0;hhalf=0;hmin=0;

    for t=2:N
        sig(t)=atan(k2/k1);
        vm(1,t)=cos(sig(t))*(vh(t-1)+1*((lam*mA0*(1-mju/X*s(t-1))*g-alpha*v(t-1)^2))*(vh(t-1)/v(t-1))/(mA0*(1-mju/X*s(t-1))));
        vm(2,t)=sin(sig(t))*(vh(t-1)+1*((lam*mA0*(1-mju/X*s(t-1))*g-alpha*v(t-1)^2))*(vh(t-1)/v(t-1))/(mA0*(1-mju/X*s(t-1))));
        vm(3,t)=vm(3,t-1)+1*((lam*mA0*(1-mju/X*s(t-1))*g-alpha*v(t-1)^2)*vm(3,t-1)/v(t-1)/(mA0*(1-mju/X*s(t-1)))-g);
        vh(t)=sqrt(vm(1,t)^2+vm(2,t)^2);
        v(t)=sqrt(vm(1,t)^2+vm(2,t)^2+vm(3,t)^2);
        s(t)=s(t-1)+1*vh(t);
        Pm(1,t)=Pm(1,t-1)+1*vm(1,t);
        Pm(2,t)=Pm(2,t-1)+1*vm(2,t);
        Pm(3,t)=Pm(3,t-1)+1*vm(3,t);
        L(t)=sqrt((Pm(1,t)-Pa(1,t))^2+(Pm(2,t)-Pa(2,t))^2+(Pm(3,t)-Pa(3,t))^2);
        if t~=1
            if Pm(3,t)<=Pm(3,t-1)
                tA=t-1;
                hmax=Pm(3,tA);
                hhalf=hmax/2;
                hmin=hmax/3;
                break;
            end
        end
    end
    temp=0;
    for t=tA+1:N  
        if Pm(3,t-1)>=hhalf
            vnx=vm(1,t-1)-1*alpha*v(t-1)*vm(1,t-1)/mB;
            vny=vm(2,t-1)-1*alpha*v(t-1)*vm(2,t-1)/mB;
            vnz=vm(3,t-1)+1*(alpha*v(t-1)*vm(3,t-1)/mB-g);
            vpxp=Pa(1,t-1-dt)-Pm(1,t-1);
            vpyp=Pa(2,t-1-dt)-Pm(2,t-1);
            vpzp=Pa(3,t-1-dt)-Pm(3,t-1);
            vpx=vpxp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpy=vpyp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpz=vpzp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vn=[vnx,vny,vnz];   %转向前速度矢量
            vp=[vpx,vpy,vpz];   %转向目标速度矢量
            cr=cross(vn,vp);    %Vn矢量与Vp矢量的叉积
            cr=cr/sqrt(cr*cr'); %单位化
            vr=vn*sin(1*w)+cross(cr,vn)*cos(1*w);   %转向后速度矢量
            vm(1,t)=vr(1);vm(2,t)=vr(2);vm(3,t)=vr(3);
            sig(t)=atan(vm(2,t)/vm(1,t));
            vh(t)=sqrt(vm(1,t)^2+vm(2,t)^2);
            v(t)=sqrt(vm(1,t)^2+vm(2,t)^2+vm(3,t)^2);
            s(t)=s(t-1)+1*vh(t);
            Pm(1,t)=Pm(1,t-1)+1*vm(1,t);
            Pm(2,t)=Pm(2,t-1)+1*vm(2,t);
            Pm(3,t)=Pm(3,t-1)+1*vm(3,t);
            L(t)=sqrt((Pm(1,t)-Pa(1,t))^2+(Pm(2,t)-Pa(2,t))^2+(Pm(3,t)-Pa(3,t))^2);
        else
            if temp==0
                tBp=t-1;
                vzhalf=vm(3,tBp);
            end
            temp=1;
            vnx=cos(sig(t-1))*(vh(t-1)+1*((lam*mB*g-alpha*v(t-1)^2)*(vh(t-1)/v(t-1))/mB));
            vny=sin(sig(t-1))*(vh(t-1)+1*((lam*mB*g-alpha*v(t-1)^2)*(vh(t-1)/v(t-1))/mB));
            vnz=vzhalf*((Pm(3,t-1)-hmin)/(hhalf-hmin));
            vpxp=Pa(1,t-1-dt)-Pm(1,t-1);
            vpyp=Pa(2,t-1-dt)-Pm(2,t-1);
            vpzp=Pa(3,t-1-dt)-Pm(3,t-1);
            vpx=vpxp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpy=vpyp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpz=vpzp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vn=[vnx,vny,vnz];   %转向前速度矢量
            vp=[vpx,vpy,vpz];   %转向目标速度矢量
            cr=cross(vn,vp);    %Vn矢量与Vp矢量的叉积
            cr=cr/sqrt(cr*cr'); %单位化
            vr=vn*sin(1*w)+cross(cr,vn)*cos(1*w);   %转向后速度矢量
            vm(1,t)=vr(1);vm(2,t)=vr(2);vm(3,t)=vr(3);
            sig(t)=atan(vm(2,t)/vm(1,t));
            vh(t)=sqrt(vm(1,t)^2+vm(2,t)^2);
            v(t)=sqrt(vm(1,t)^2+vm(2,t)^2+vm(3,t)^2);
            s(t)=s(t-1)+1*vh(t);
            Pm(1,t)=Pm(1,t-1)+1*vm(1,t);
            Pm(2,t)=Pm(2,t-1)+1*vm(2,t);
            Pm(3,t)=Pm(3,t-1)+1*vm(3,t);
            L(t)=sqrt((Pm(1,t)-Pa(1,t))^2+(Pm(2,t)-Pa(2,t))^2+(Pm(3,t)-Pa(3,t))^2);
        end
        if L(t)<dC
           tB=t-1;
           break;
        end
    end
    
    for t=tB+1:N
        sig(t)=atan(k2/k1);
        vnx=cos(sig(t))*(vh(t-1)+1*((lam*mC*g-alpha*v(t-1)^2))*(vh(t-1)/v(t-1))/(mC));
        vny=sin(sig(t))*(vh(t-1)+1*((lam*mC*g-alpha*v(t-1)^2))*(vh(t-1)/v(t-1))/(mC));
        vnz=vm(3,t-1)+1*((-lam*mC*g+alpha*v(t-1)^2)*vm(3,t-1)/v(t-1)/mC-g);
        vpxp=Pa(1,t-1-dt)-Pm(1,t-1);
        vpyp=Pa(2,t-1-dt)-Pm(2,t-1);
        vpzp=Pa(3,t-1-dt)-Pm(3,t-1);
        vpx=vpxp/sqrt(vpxp^2+vpyp^2+vpzp^2);
        vpy=vpyp/sqrt(vpxp^2+vpyp^2+vpzp^2);
        vpz=vpzp/sqrt(vpxp^2+vpyp^2+vpzp^2);
        vn=[vnx,vny,vnz];   %转向前速度矢量
        vp=[vpx,vpy,vpz];   %转向目标速度矢量
        cr=cross(vn,vp);    %Vn矢量与Vp矢量的叉积
        cr=cr/sqrt(cr*cr'); %单位化
        vr=vn*sin(1*w)+cross(cr,vn)*cos(1*w);   %转向后速度矢量
        vm(1,t)=vr(1);vm(2,t)=vr(2);vm(3,t)=vr(3);
        sig(t)=atan(vm(2,t)/vm(1,t));
        vh(t)=sqrt(vm(1,t)^2+vm(2,t)^2);
        v(t)=sqrt(vm(1,t)^2+vm(2,t)^2+vm(3,t)^2);
        s(t)=s(t-1)+1*vh(t);
        Pm(1,t)=Pm(1,t-1)+1*vm(1,t);
        Pm(2,t)=Pm(2,t-1)+1*vm(2,t);
        Pm(3,t)=Pm(3,t-1)+1*vm(3,t);
        L(t)=sqrt((Pm(1,t)-Pa(1,t))^2+(Pm(2,t)-Pa(2,t))^2+(Pm(3,t)-Pa(3,t))^2);
        if Pm(3,t)<0
           tC=t-1;
           break;
        end
    end
    if sqrt((Pm(1,tC)-Pa(1,tC))^2+(Pm(2,tC)-Pa(2,tC))^2)<=u*Lac/2
        nh=nh+1;
    end
    end
    ph(ss)=nh/num;
    pk(ss)=(4*Wac)/(pi*Lac*u^2)*ph(ss);
end

figure
u=linspace(5.7,6.7,n);
[hAx,hLine1,hLine2]=plotyy(u,ph,u,pk);
hLine1.LineStyle='--';
hLine2.LineStyle=':';
hLine1.LineWidth=1.4;
hLine2.LineWidth=1.4;
legend('杀伤率','直接命中率');
xlabel('u=导弹杀伤半径/航母半径');
ylabel(hAx(1),'杀伤率');
ylabel(hAx(2),'命中率');
title('u对杀伤率和命中率的影响');
box on
grid on
