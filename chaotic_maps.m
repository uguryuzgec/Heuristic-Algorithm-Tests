% Chaotic maps uygulamalarý
% 24.08.2016, Bilecik

% Logistic map
% X(n+1)=a.X(n).(1-X(n))
% a katsayýsý 3.56995...
close all
NP = 250;
a = 3.99;
X0 = rand();
popnokta (1,1) = a*X0*(1-X0);
for i=2:NP
   popnokta(i,1) = a*popnokta(i-1,1)*(1-popnokta(i-1,1));	
   % popnokta (i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
figure
subplot(121)
plot(popnokta);
axis tight;
subplot(122)
plot(popnokta(1:NP-1),popnokta(2:NP));
title('Logistic map');
% Henon map
% x(n+1) = 1-a*x(n)^2+y(n);
% y(n+1)=b*x(n);

% Initialize
Y0 = 0;
a = 1.4;
b = 0.3;
popnokta(1,1)=1-a*X0^2+Y0;
y(1,1)=b*X0;
for i=2:NP
   popnokta(i,1) = 1-a*popnokta(i-1,1)^2+y(i-1,1);
   y(i,1)=b*popnokta(i-1,1);
   % popnokta (i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
figure
subplot(121)
plot(y);
axis tight;
subplot(122)
plot(popnokta(1:NP),y(1:NP));
title('Henon map');

% Tent map
% x(n+1) = x(n)/0.7, x(n)<0.7
% x(n+1) = 10/3x(n)(1-x(n)), diðer...
% Initialize

if (X0<0.7)
   popnokta(1,1)=X0/0.7;
else
   popnokta(1,1)=(10/3)*X0*(1-X0);
end
for i=2:NP
   if (popnokta(i-1,1)<0.7)
      popnokta(i,1)=popnokta(i-1,1)/0.7;
   else
      popnokta(i,1)=(10/3)*popnokta(i-1,1)*(1-popnokta(i-1,1));
   end
% popnokta (i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
figure
subplot(121)
plot(popnokta);
axis tight;
subplot(122)
plot(popnokta(1:NP-1),popnokta(2:NP));
title('Tent map');

% Sinusoidal iterator
% x(n+1) = a.x(n)^2.sin(pi.x(n))
% Initialize
%  a = 2.3; X0=0.4323312;
popnokta(1,1)=sin(pi*X0);

for i=2:NP
   popnokta(i,1) = sin(pi*popnokta(i-1,1));
   % popnokta (i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
figure
subplot(121)
plot(popnokta);
axis tight;
subplot(122)
plot(popnokta(1:NP-1),popnokta(2:NP));
title('Sinusoidal iterator');

% Circle map
% X(n+1)=X(n)+ b - (a/2*pi)Sin(2piX(n))
% a = 0.5, b = 0.2
a = 0.5;
b = 0.2;
popnokta (1,1) = mod(X0+b-(a/(2*pi))*sin(2*pi*X0),1);
% mod(x(i)+b-(a/(2*pi))*sin(2*pi*x(i)),1);
for i=2:NP
   popnokta(i,1) = mod(popnokta(i-1,1)+b-(a/(2*pi))*sin(2*pi*popnokta(i-1,1)),1);	
   % popnokta (i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
figure
subplot(121)
plot(popnokta);
axis tight;
subplot(122)
plot(popnokta(1:NP-1),popnokta(2:NP));
title('Circle map');