%%% 1(4) %%%
clear,close,clc;
format compact;
warning off;fa=0.825;    %���ŷ����
display(['���ŷ���Ƕ�Ϊ ',num2str(fa),' rad = ',num2str(fa/pi*180),'��']);
N=3000;
n=50;
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
figure
hold on;
    theta=fa
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
plot(x(1:tC),y(1:tC),'r.-','linewidth',1.8);
xlabel('���о���x(m)');
ylabel('���и߶�y(m)');
title('��ʼ״̬�µ����ľ�̬���');
hold off;

display(['�ܷ���ʱ�䣺',num2str(tC/60),' min']);
display(['�ܷ��о��룺',num2str(x(tC)/1000), ' km']);
display(['���Ŀ���ٶȣ�',num2str(3.6*v(tC)),' km/h']);
display(['��Ŀ�����',num2str(abs(X-x(tC))),' m']);
display(['����γ���ʱ�䣺',num2str(tA/60),' min']);
display(['����η��о��룺',num2str(x(tA)/1000), ' km']);
display(['�жγ���ʱ�䣺',num2str((tB-tA)/60),' min']);
display(['�жη��о��룺',num2str((x(tB)-x(tA))/1000), ' km']);
display(['ĩ�γ���ʱ�䣺',num2str((tC-tB)/60),' min']);
display(['ĩ�η��о��룺',num2str((x(tC)-x(tB))/1000), ' km']);

figure
hold on

subplot(2,2,1);
plot(1:tC,x(1:tC));
ylabel('����(m)');
xlabel('ʱ��(s)');
title('����-ʱ��ͼ��');

subplot(2,2,2);
plot(1:tC,v(1:tC));
ylabel('�ٶ�(m/s)');
xlabel('ʱ��(s)');
title('�ٶ�-ʱ��ͼ��');

subplot(2,2,3);
a=v;
for i=1:tC-1
    a(i)=v(i+1)-v(i);
end
plot(1:tC-1,a(1:tC-1));
ylabel('���ٶ�(m/s)');
xlabel('ʱ��(s)');
title('���ٶ�-ʱ��ͼ��');

subplot(2,2,4);
fy=zeros(1,tC);
for i=1:tC
    fy(i)=180./pi.*atan(vy(i)./vx(i));
end
plot(1:tC,fy);
ylabel('������(��)');
xlabel('ʱ��(s)');
title('������-ʱ��ͼ��');

hold off