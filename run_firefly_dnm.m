% Sezgisel Optimizasyon Algoritmalari Nisan 2015
% Duygu ADIYAMAN Ates Bocegi Optimizasyon Algoritmasi (Firefly Optimization Algorithm)
% Self Adaptive yapısı ilave edildi.
% ICNES2015 konferansı çalışması...

clear all;
close all;
prevpath = path;
path(path, genpath(fileparts(mfilename('fullpath'))));
%--------------------------------------------------------------------------------------
% Duygu ADIYAMAN Ates Bocegi Optimizasyon Algoritmasi (Firefly Optimization Algorithm) 
% -------------------------------------------------------------------------------------
alpha = 0.2;      % Randomness 0--1 (highly random)
gamma = 0.5;      % Absorption coefficient
delta = 0.97;     % Randomness reduction (similar to an annealing schedule)
% -------------------------------------------------------------------------------------
% Farksal Gelisim Algoritmasi (Differential Evolution Algorithm) 
% -------------------------------------------------------------------------------------
F = 0.8;      % scale factor DE-step size F ex [0, 2]
CR = 0.5;     % crossover probability constant ex [0, 1]
% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin
strategy = 7;
% -------------------------------------------------------------------------------------

% refresh       intermediate output will be produced after "refresh"
%               iterations. No intermediate output will be produced
%               if refresh is < 1
		refresh = 50; 
% VTR		"Value To Reach" (stop when ofunc < VTR)
		VTR = 1.e-6; 
% itermax       maximum number of iterations (generations)
		itermax = 200; 
% NP 			number of population (NP=10*D)		
        NP = 20;

[dims, lb, ub, solution, minimum, fonk]=RUN_ezimage;
        XVmin = lb; 
		XVmax = ub;
        D=dims;
        tekrar=10;
% -------------------------------------------------------------------------------------
[bestmem1,bestval1,nfeval1,G1,Gi1,Gs1]=firefly(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,alpha,gamma,delta); 

figure
subplot(221)
for i=1:tekrar,
    plot(G1(i,1:Gi1(i)),'color',rand(1,3));
    hold on
end
title('Firefly Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
hold off
subplot(222)
for i=1:tekrar,
    plot((G1(i,1:Gi1(i))-minimum).^2,'color',rand(1,3));
    hold on
end
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN1=mean(G1(:,1:min(Gi1)));
    plot(MEAN1(1:min(Gi1)));
else
    MEAN1=G1;
    plot(MEAN1(1:min(Gi1)));
end
title('Firefly Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow
% -------------------------------------------------------------------------------------
[bestmem2,bestval2,nfeval2,G2,Gi2,Gs2]=firefly_self_adapt(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,alpha,gamma,delta); 

figure
subplot(221)
for i=1:tekrar,
    plot(G2(i,1:Gi2(i)),'color',rand(1,3));
    hold on
end
title('Self adaptive Firefly Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
hold off
subplot(222)
for i=1:tekrar,
    plot((G2(i,1:Gi2(i))-minimum).^2,'color',rand(1,3));
    hold on
end
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN2=mean(G2(:,1:min(Gi2)));
    plot(MEAN2(1:min(Gi2)));
else
    MEAN2=G2;
    plot(MEAN2(1:min(Gi2)));
end
title('Self adaptive Firefly Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow
% -------------------------------------------------------------------------------------   
fprintf(1,'\n%s probleminin sonucu : %3.3f\n',fonk,minimum);
for i=1:size(solution,1),
    fprintf(1,'x1 : %3.3f ve x2 : %3.3f \n',solution(i,1),solution(i,2));
end
% -------------------------------------------------------------------------------------
%.....Toplu sonuclar.....
% -------------------------------------------------------------------------------------
figure
plot(MEAN1(1:min(Gi1))','g'); hold on;
plot(MEAN2(1:min(Gi2))','r');

% legend('Firefly','Differential','Bat','Cuckoo','Brainstorm','Magnetic-Inspired','Forest','Crab','Spiral','Flower','Fruitfly','Hurricane');
legend('FOA','SFOA');
legend show
title('Heuristic Optimization Algorithms')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;