%%% 1(2) %%%
clear,close,clc;
format compact;
warning off;
thmin=0;
thmax=pi/2;
N=3000;
n=10000;
lam=1.8;    %推重比（以RD-0410核火箭发动机为例）
X=1000*distance(120.5,27.5,123.75,25.65)  %初始距离
mju=0.86;   %助推器占总重比（以C-301超音速反舰导弹为例）
w=0.34; %导弹偏航和俯仰角速度 rad/s
dC=111e3;   %末段飞行距离(以美国AN/TPS-79防空雷达作用距离为例)
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
for theta=linspace(pi/16,7*pi/16,n)
    tA=0;tBp=0;tB=0;tC=0;hmax=0;hhalf=0;hmin=0;
    x=zeros(1,N);
    y=zeros(1,N);
    vx=zeros(1,N);
    vy=zeros(1,N);
    v=zeros(1,N);
    delphi=zeros(1,N);
    vx(1)=v0*cos(theta);
    vy(1)=v0*sin(theta);
    v(1)=v0;
    x(1)=0;y(1)=0;
    for t=2:N
        vx(t)=vx(t-1)+1*(lam*mA0*(1-mju/X*x(t-1))*g-alpha*(vx(t-1)^2+vy(t-1)^2))*vx(t-1)/sqrt(vx(t-1)^2+vy(t-1)^2)/(mA0*(1-mju/X*x(t-1)));
        vy(t)=vy(t-1)+1*(lam*mA0*(1-mju/X*x(t-1))*g-alpha*(vx(t-1)^2+vy(t-1)^2))*vy(t-1)/sqrt(vx(t-1)^2+vy(t-1)^2)/(mA0*(1-mju/X*x(t-1)))-g;
        v(t)=sqrt(vx(t)^2+vy(t)^2);
        x(t)=x(t-1)+vx(t)*1;
        y(t)=y(t-1)+vy(t)*1;
        if t~=1
            if y(t)<=y(t-1)
                tA=t-1;
                hmax=y(tA);
                hhalf=hmax/2;
                hmin=hmax/3;
                break;
            end
        end
    end
    temp=0;
    for t=tA+1:N
        
        if y(t-1)>=hhalf
            vx(t)=vx(t-1)+(-alpha*vx(t-1)*sqrt(vx(t-1)^2+vy(t-1)^2))/mB;
            vy(t)=vy(t-1)-g+(alpha*vy(t-1)*sqrt(vx(t-1)^2+vy(t-1)^2))/mB;
            v(t)=sqrt(vx(t)^2+vy(t)^2);
            x(t)=x(t-1)+1*vx(t);
            y(t)=y(t-1)+1*vy(t);
        else
            if temp==0
                tBp=t-1;
            end
            temp=1;
            vyhalf=vy(t-1);
            vx(t)=vx(t-1)+(lam*mB*g-alpha*vx(t-1)*sqrt(vx(t-1)^2+vy(t-1)^2))/mB;
            vy(t)=vyhalf*((y(t-1)-hmin)/(hhalf-hmin));
            v(t)=sqrt(vx(t)^2+vy(t)^2);
            x(t)=x(t-1)+1*vx(t);
            y(t)=y(t-1)+1*vy(t);
        end
        if x(t)>=X-dC
           tB=t-1;
           break;
        end
    end
    for t=tB:N
        phi1=atan((y(t-1)-0)/(X-x(t-1)));
        phi2=atan(-vy(t-1)/vx(t-1));
        delphi(t)=180*(phi1-phi2)/pi;
        vx(t)=vx(t-1)*cos(phi2+w*sign(phi1-phi2));
        vy(t)=vy(t-1)*sin(phi2+w*sign(phi1-phi2));
        vx(t)=vx(t)+1/mC*(-alpha*(vx(t)^2+vy(t)^2)+lam*mC*g)*cos(phi2+w*sign(phi1-phi2));    
        vy(t)=vy(t)+1/mC*(alpha*(vx(t)^2+vy(t)^2)-lam*mC*g)*sin(phi2+w*sign(phi1-phi2))-g;
        v(t)=sqrt(vx(t)^2+vy(t)^2);
        x(t)=x(t-1)+1*vx(t);
        y(t)=y(t-1)+1*vy(t);
        if y(t)<0
            tC=t-1;
            break;
        end
    end
    if abs(x(tC)-X)<=Lac
       if thmin==0
          thmin=theta;
       else thmax=theta;
       end
    end
end
display(['最小发射角度为 ',num2str(thmin),' rad = ',num2str(thmin/pi*180),'°']);
display(['最大发射角度为 ',num2str(thmax),' rad = ',num2str(thmax/pi*180),'°']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=3000;
n=20;

figure
hold on;
for theta=linspace(thmin,thmax,n)
    tA=0;tBp=0;tB=0;tC=0;hmax=0;hhalf=0;hmin=0;
    x=zeros(1,N);
    y=zeros(1,N);
    vx=zeros(1,N);
    vy=zeros(1,N);
    v=zeros(1,N);
    delphi=zeros(1,N);
    vx(1)=v0*cos(theta);
    vy(1)=v0*sin(theta);
    v(1)=v0;
    x(1)=0;y(1)=0;
    for t=2:N
        vx(t)=vx(t-1)+1*(lam*mA0*(1-mju/X*x(t-1))*g-alpha*(vx(t-1)^2+vy(t-1)^2))*vx(t-1)/sqrt(vx(t-1)^2+vy(t-1)^2)/(mA0*(1-mju/X*x(t-1)));
        vy(t)=vy(t-1)+1*(lam*mA0*(1-mju/X*x(t-1))*g-alpha*(vx(t-1)^2+vy(t-1)^2))*vy(t-1)/sqrt(vx(t-1)^2+vy(t-1)^2)/(mA0*(1-mju/X*x(t-1)))-g;
        v(t)=sqrt(vx(t)^2+vy(t)^2);
        x(t)=x(t-1)+vx(t)*1;
        y(t)=y(t-1)+vy(t)*1;
        if t~=1
            if y(t)<=y(t-1)
                tA=t-1;
                hmax=y(tA);
                hhalf=hmax/2;
                hmin=hmax/3;
                break;
            end
        end
    end
    temp=0;
    for t=tA+1:N
        
        if y(t-1)>=hhalf
            vx(t)=vx(t-1)+(-alpha*vx(t-1)*sqrt(vx(t-1)^2+vy(t-1)^2))/mB;
            vy(t)=vy(t-1)-g+(alpha*vy(t-1)*sqrt(vx(t-1)^2+vy(t-1)^2))/mB;
            v(t)=sqrt(vx(t)^2+vy(t)^2);
            x(t)=x(t-1)+1*vx(t);
            y(t)=y(t-1)+1*vy(t);
        else
            if temp==0
                tBp=t-1;
            end
            temp=1;
            vyhalf=vy(t-1);
            vx(t)=vx(t-1)+(lam*mB*g-alpha*vx(t-1)*sqrt(vx(t-1)^2+vy(t-1)^2))/mB;
            vy(t)=vyhalf*((y(t-1)-hmin)/(hhalf-hmin));
            v(t)=sqrt(vx(t)^2+vy(t)^2);
            x(t)=x(t-1)+1*vx(t);
            y(t)=y(t-1)+1*vy(t);
        end
        if x(t)>=X-dC
           tB=t-1;
           break;
        end
    end
    for t=tB:N
        phi1=atan((y(t-1)-0)/(X-x(t-1)));
        phi2=atan(-vy(t-1)/vx(t-1));
        delphi(t)=180*(phi1-phi2)/pi;
        vx(t)=vx(t-1)*cos(phi2+w*sign(phi1-phi2));
        vy(t)=vy(t-1)*sin(phi2+w*sign(phi1-phi2));
        vx(t)=vx(t)+1/mC*(-alpha*(vx(t)^2+vy(t)^2)+lam*mC*g)*cos(phi2+w*sign(phi1-phi2));    
        vy(t)=vy(t)+1/mC*(alpha*(vx(t)^2+vy(t)^2)-lam*mC*g)*sin(phi2+w*sign(phi1-phi2))-g;
        v(t)=sqrt(vx(t)^2+vy(t)^2);
        x(t)=x(t-1)+1*vx(t);
        y(t)=y(t-1)+1*vy(t);
        if y(t)<0
            tC=t-1;
            break;
        end
    end
    plot(x(1:tC),y(1:tC),'.:','linewidth',1.2);
    xlabel('飞行距离x(m)');
    ylabel('飞行高度y(m)');
    title('不同发射角度下的反舰导弹飞行曲线');
end
hold off;