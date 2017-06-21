%22.03.2016 U.Y.-T.I
% Spiral Algoritmasý parametre gösterimi
% r=0.95, Theta=pi/4
% r=0.90, Theta=pi/4
% r=0.95, Theta=pi/2
% r=0.95, Theta=pi/3

clc
clear all
D=2;
NP=100;
Xmerkez=zeros(1,D);  
x=zeros(NP,D);
x(1,:)=[10 10];
M=eye(D);
I=eye(D);
r=0.95; Theta=pi/4;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 
figure
subplot(221)
plot(x(:,1),x(:,2),'k');
axis([-15 15 -15 15]) 
xlabel('x_1');
ylabel('x_2');

r=0.9; Theta=pi/4;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 

subplot(222)
plot(x(:,1),x(:,2),'k');
axis([-15 15 -15 15]) 
xlabel('x_1');
ylabel('x_2');

r=0.95; Theta=pi/2;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 

subplot(223)
plot(x(:,1),x(:,2),'k');
axis([-15 15 -15 15]) 
xlabel('x_1');
ylabel('x_2');

r=0.95; Theta=pi/3;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 

subplot(224)
plot(x(:,1),x(:,2),'k');
axis([-15 15 -15 15]) 
xlabel('x_1');
ylabel('x_2');

figure
r=0.95; Theta=pi/4;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 
plot(x(:,1),x(:,2),'k');
axis([-15 15 -15 15]) 
xlabel('x_1');
ylabel('x_2');
hold on
r=0.9; Theta=pi/4;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 
plot(x(:,1),x(:,2),'b');
r=0.95; Theta=pi/6;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 
plot(x(:,1),x(:,2),'r');
r=0.95; Theta=pi/3;
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
for i=2:NP,
    x(i,:)=(r*M)*x(i-1,:)'-(r*M-I)*Xmerkez'; %sonuc 
end 
plot(x(:,1),x(:,2),'g');
legend('r=0.95; Theta=pi/4','r=0.9; Theta=pi/4','r=0.95; Theta=pi/6','r=0.95; Theta=pi/3'); 
legend show