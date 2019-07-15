%%% 2(2) %%%
clear,close,clc;
format compact;
warning off;
num=50;
figure
hold on
ttt=0;
for w=linspace(1.45,1.57,num) %%%%%%%%%%%%����ƫ���͸������ٶ� rad/s
for dC=linspace(9e3,17e3,num) %%%%%%%%%%%%%ĩ�η��о���(������AN/TPS-79�����״����þ���Ϊ��)
nh=0;   %���и���
rh=0;   %������
N=3000;
n=50;
RE=6300; %����뾶 km
Pa=zeros(3,N); %��ĸλ��
va=32*1.852/3.6; %��ĸ�����ٶ� m/s
X=1000*distance(120.5,27.5,123.75,25.65);  %��ʼ���� m
k1=1000*RE*abs(pi/180*(27.5-25.65));    %��ʼλ���ϱ�������루�ؾ��ߣ�m
k2=1000*RE*cos(pi/180*25.65)*(pi/180*(123.75-120.5));    %��ʼλ�ö���������루��γ�ߣ�m
Xa=sqrt(k1^2+k2^2); %��ʼ���������ֵ m
Pa(:,1)=[k1;k2;0];
lam=1.8;    %���رȣ���RD-0410�˻��������Ϊ����
mju=0.86;   %������ռ���رȣ���C-301�����ٷ�������Ϊ����

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
for i=2:N
    Pa(1,i)=double(k1)+double(va*(i-1));
    Pa(2,i)=k2;
    Pa(3,i)=0;
end

for theta=linspace(0.6332,1.3744,n) %ѡȡ��1�ʽ���������ʽϸߵķ���Ƿ�Χ
    Pm=zeros(3,N); %����λ��
    vm=zeros(3,N); %������x,y,z����ٶ�
    vh=zeros(3,N); %����ˮƽ���ٶ�
    v=zeros(1,N);   %�������ٶ�
    s=zeros(1,N); %����ˮƽ����·��
    L=zeros(1,N); %�����뺽ĸ�Ķ�̬����
    sig=zeros(1,N); %�����ٶ��ں�ƽ���ͶӰ�뾭�߼н�
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
            vpxp=Pa(1,t-1)-Pm(1,t-1);
            vpyp=Pa(2,t-1)-Pm(2,t-1);
            vpzp=Pa(3,t-1)-Pm(3,t-1);
            vpx=vpxp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpy=vpyp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpz=vpzp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vn=[vnx,vny,vnz];   %ת��ǰ�ٶ�ʸ��
            vp=[vpx,vpy,vpz];   %ת��Ŀ���ٶ�ʸ��
            cr=cross(vn,vp);    %Vnʸ����Vpʸ���Ĳ��
            cr=cr/sqrt(cr*cr'); %��λ��
            vr=vn*sin(1*w)+cross(cr,vn)*cos(1*w);   %ת����ٶ�ʸ��
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
            vpxp=Pa(1,t-1)-Pm(1,t-1);
            vpyp=Pa(2,t-1)-Pm(2,t-1);
            vpzp=Pa(3,t-1)-Pm(3,t-1);
            vpx=vpxp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpy=vpyp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vpz=vpzp/sqrt(vpxp^2+vpyp^2+vpzp^2);
            vn=[vnx,vny,vnz];   %ת��ǰ�ٶ�ʸ��
            vp=[vpx,vpy,vpz];   %ת��Ŀ���ٶ�ʸ��
            cr=cross(vn,vp);    %Vnʸ����Vpʸ���Ĳ��
            cr=cr/sqrt(cr*cr'); %��λ��
            vr=vn*sin(1*w)+cross(cr,vn)*cos(1*w);   %ת����ٶ�ʸ��
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
        vpxp=Pa(1,t-1)-Pm(1,t-1);
        vpyp=Pa(2,t-1)-Pm(2,t-1);
        vpzp=Pa(3,t-1)-Pm(3,t-1);
        vpx=vpxp/sqrt(vpxp^2+vpyp^2+vpzp^2);
        vpy=vpyp/sqrt(vpxp^2+vpyp^2+vpzp^2);
        vpz=vpzp/sqrt(vpxp^2+vpyp^2+vpzp^2);
        vn=[vnx,vny,vnz];   %ת��ǰ�ٶ�ʸ��
        vp=[vpx,vpy,vpz];   %ת��Ŀ���ٶ�ʸ��
        cr=cross(vn,vp);    %Vnʸ����Vpʸ���Ĳ��
        cr=cr/sqrt(cr*cr'); %��λ��
        vr=vn*sin(1*w)+cross(cr,vn)*cos(1*w);   %ת����ٶ�ʸ��
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
    if sqrt((Pm(1,tC)-Pa(1,tC))^2+(Pm(2,tC)-Pa(2,tC))^2)<=5*Lac/2
        nh=nh+1;
    end
end
rh=nh/n;
plot3(w,dC,rh,'b.');
end
ttt=ttt+1;
display(ttt);
end

xlabel('w(rad/s)');ylabel('dC(m)');zlabel('������');
title('��ͬƫ�����������ٶȺ͵з������״����þ����µĵ���������');
box on
grid on
axis vis3d
hold off