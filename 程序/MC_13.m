%%% 1(3) %%%
clear,close,clc;
format compact;
warning off;
N=3000;
figure;
for i=1:4
n=10^i;
lam=1.8;    %���رȣ���RD-0410�˻��������Ϊ����
X=1000*distance(120.5,27.5,123.75,25.65)  %��ʼ����
mju=0.86;   %������ռ���رȣ���C-301�����ٷ�������Ϊ����
w=0.34; %����ƫ���͸������ٶ� rad/s
dC=111e3;   %ĩ�η��о���(����AN/TPS-79�����״����þ���)
mA0=1.5e3;  %����ε�����ʼ���أ�����������
mB=mA0*(1-mju); %�жε�������
mC=mB;  %ĩ�ε�������
g=9.8;  %�������ٶ�
v0=500; %�����ٶ�
C0=0.1747;  %��������ϵ��
Ci=2.52e-3; %�յ�����ϵ��
r=1.293;	%��������ܶ�
Dm=0.4; %����ֱ�����ԻƷ�󷴽�����Ϊ����
S=pi/4*Dm^2;    %����ӭ�����
alpha=0.5*(C0+Ci)*r*S;  %�����������ٶ�ƽ���ı���ϵ��
Lac=335;    %��ĸ����
s=0;
th=zeros(1,n);
ve=zeros(1,n);
T=zeros(1,n);
for theta=linspace(0.6332,1.3744,n)
    s=s+1;
    th(s)=theta;
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
            ve(s)=v(tC);
            break;
        end
    end
    T(s)=tC;
end
subplot(2,2,i);
[hAx,hLine1,hLine2]=plotyy(th,ve,th,T);
xlabel('�����(rad)');
ylabel(hAx(1),'����ٶ�(m/s)');
ylabel(hAx(2),'������ʱ��(s)');
title(['n=',num2str(n)]);
end
figure
[hAx,hLine1,hLine2]=plotyy(th,6.703./(th).^7+73.47,th,T);
hLine1.LineStyle = '-';
hLine2.LineStyle = ':';
legend('��С����ٶ��������','������ʱ������');
xlabel('�����(rad)');
ylabel(hAx(1),'��С����ٶ�(m/s)');
ylabel(hAx(2),'������ʱ��(s)');
title('��С����ٶ���������������ʱ��ıȽ�');