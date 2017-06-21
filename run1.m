% Sezgisel Optimizasyon Algoritmalari Nisan 2015
% Ali CAKMAK Guguk Kusu Optimizasyon Algoritmasi (Cuckoo Optimization Algorithm) XXX
% Ali Sahin BAYOGLU Parlak Fikir Optimizasyon Algoritmasi (Brainstorm Optimization Algorithm) XXX
% Duygu ADIYAMAN Ates Bocegi Optimizasyon Algoritmasi (Firefly Optimization Algorithm) XXX
% Gizem ATAC KALE Yarasa Optimizasyon Algoritmasi (Bat Optimization Algorithm) XXX
% Gurkan Mustafa CAKIR Kasirga temelli Optimizasyon Algoritmasi (Hurricane-Based Optimization Algorithm) XXX
% Merve BULDUR Manyetik Ilhamli Optimizasyon Algoritmasi (Magnetic-Inspired Optimization Algorithm) XXX
% Mustafa KOSTEK Orman Optimizasyon Algoritmasi (Forest Optimization Algorithm) XXX
% Nartan Ayberk ASKIN Yengec Ciftlesme Optimizasyon Algoritmasi (Crab Mating Optimization Algorithm) XXX
% Refika ANBAR Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm) XXX
% Sumeyye SAG Yunus Ekosu ile Yer Tespiti Optimizasyon Algoritmasi (Dolphin Echolocation Optimization Algorithm) ???
% Yusuf BORUCU Cicek Tozlasma Optimizasyon Algoritmasi (Flower Pollination Optimization Algorithm) XXX
% Yusuf ONDER Penguen Arama Optimizasyon Algoritmasi (Penguins Search Optimization Algorithm) ??? 
% Zuhal CALAYOGLU Meyve Sinegi Optimizasyon Algoritmasi (Fruit Fly Optimization Algorithm) XXX 
 
clear all;
close all;
prevpath = path;
path(path, genpath(fileparts(mfilename('fullpath'))));
%--------------------------------------------------------------------------------------
% Duygu ADIYAMAN Ates Bocegi Optimizasyon Algoritmasi (Firefly Optimization Algorithm) 
% -------------------------------------------------------------------------------------
alpha = 0.2;      % Randomness 0--1 (highly random)
gamma = 0.5;       % Absorption coefficient
delta = 0.97;     % Randomness reduction (similar to an annealing schedule)
% -------------------------------------------------------------------------------------
% Gizem ATAC KALE Yarasa Optimizasyon Algoritmasi (Bat Optimization Algorithm)
% -------------------------------------------------------------------------------------
A = 0.5; %A=Loudness  (constant or decreasing)
r = 0.5; %r=Pulse rate (constant or decreasing)
Qmin = 0; %min. frekans
Qmax = 2; %mak.frekans
%--------------------------------------------------------------------------------------
% Ali CAKMAK Guguk Kusu Optimizasyon Algoritmasi (Cuckoo Optimization Algorithm)
% -------------------------------------------------------------------------------------
pa=0.25; % Yabancı Yumurtaları Bulma/Keşif Çözüm Oranı
%--------------------------------------------------------------------------------------
% Ali Sahin BAYOGLU Parlak Fikir Optimizasyon Algoritmasi (Brainstorm Optimization Algorithm)
% -------------------------------------------------------------------------------------
prob_one_cluster = 0.8; % probability for select one cluster to form new individual; 
NC=4; % number of clusters
%--------------------------------------------------------------------------------------
% Merve BULDUR Manyetik Ilhamli Optimizasyon Algoritmasi (Magnetic-Inspired Optimization Algorithm)
% -------------------------------------------------------------------------------------
alfa_moa=100; % the exploitation parameter
ro=100;       % the exploitation parameter
%--------------------------------------------------------------------------------------
% Mustafa KOSTEK Orman Optimizasyon Algoritmasi (Forest Optimization Algorithm)
% -------------------------------------------------------------------------------------
drag=0.1;     % adım parametresi
esik=0.05;    % eşik belirtir
yas_snr=10;   % bu yaştan yüksek olanlar silinir 
%--------------------------------------------------------------------------------------
% Nartan Ayberk ASKIN Yengec Ciftlesme Optimizasyon Algoritmasi (Crab Mating Optimization Algorithm) 
% -------------------------------------------------------------------------------------
receptivity=200;  % alım derecesi... 1 den büyük olmak zorunda
esikdegeri=0.75; 
NoCross=0.75;     % çaprazlama değeri...
NoMut=0.15;       % mutasyon değeri... 
%--------------------------------------------------------------------------------------
% Refika ANBAR Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm)
% -------------------------------------------------------------------------------------
Theta=pi/4;
radius=0.95; 
%--------------------------------------------------------------------------------------
% Yusuf BORUCU Cicek Tozlasma Optimizasyon Algoritmasi (Flower Pollination Optimization Algorithm) 
% -------------------------------------------------------------------------------------
prob_switch = 0.8; % probability switch(olasılık anahtarı)
Levy_beta=3/2; % Levy'de kullanılan formül sabiti
%--------------------------------------------------------------------------------------
% Zuhal CALAYOGLU Meyve Sinegi Optimizasyon Algoritmasi (Fruit Fly Optimization Algorithm)  
% -------------------------------------------------------------------------------------
% No parameter
%--------------------------------------------------------------------------------------
% Yusuf ONDER Penguen Arama Optimizasyon Algoritmasi (Penguins Search Optimization Algorithm)
% -------------------------------------------------------------------------------------
% ???
%--------------------------------------------------------------------------------------
% Gurkan Mustafa CAKIR Kasirga temelli Optimizasyon Algoritmasi (Hurricane-Based Optimization Algorithm)
% -------------------------------------------------------------------------------------
Rmax = 0.2;    % maximum radius
w = pi/4;      % Changed to 8 from 10
R0 = 10^-2;    % initial radius
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
		VTR = 1.e-5; 
% itermax       maximum number of iterations (generations)
		itermax = 200; 
% NP 			number of population (NP=10*D)		
        NP = 20;

[dims, lb, ub, solution, minimum, fonk]=RUN_ezimage;
        XVmin = lb; 
		XVmax = ub;
        D=dims;
        tekrar=1;
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
xlabel('iterations');
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
xlabel('iterations');
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
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow
% -------------------------------------------------------------------------------------
[bestmem2,bestval2,nfeval2,G2,Gi2,Gs2]=devec3_orj(fonk,VTR,D,XVmin,XVmax,NP,itermax,F,CR,strategy,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G2(i,1:Gi2(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G2(i,1:Gi2(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
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
title('Differential Evolution Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
% -------------------------------------------------------------------------------------
[bestmem3,bestval3,nfeval3,G3,Gi3,Gs3]=bat(fonk,VTR,D,XVmin,XVmax,NP,itermax,Qmin,Qmax,A,r,refresh,tekrar);
%function [Gbestmem,Gbestval,nfeval1,Gmin,Giter,GSure] = bat(fname,VTR,D,XVmin,XVmax,NP,itermax,Qmin,Qmax,A,r,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G3(i,1:Gi3(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G3(i,1:Gi3(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN3=mean(G3(:,1:min(Gi3)));
    plot(MEAN3(1:min(Gi3)));
else
    MEAN3=G3;
    plot(MEAN3(1:min(Gi3)));
end
title('Bat Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem4,bestval4,nfeval4,G4,Gi4,Gs4]=cuckoo(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,pa);
%function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure]=cuckoo(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,pa)

figure
subplot(221)
for i=1:tekrar,
    plot(G4(i,1:Gi4(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G4(i,1:Gi4(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN4=mean(G4(:,1:min(Gi4)));
    plot(MEAN4(1:min(Gi4)));
else
    MEAN4=G4;
    plot(MEAN4(1:min(Gi4)));
end
title('Cuckoo Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem5,bestval5,nfeval5,G5,Gi5,Gs5]=brainstorm(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,prob_one_cluster,NC);
%function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = brainstorm(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,prob,NC)

figure
subplot(221)
for i=1:tekrar,
    plot(G5(i,1:Gi5(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G5(i,1:Gi5(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN5=mean(G5(:,1:min(Gi5)));
    plot(MEAN5(1:min(Gi5)));
else
    MEAN5=G5;
    plot(MEAN5(1:min(Gi5)));
end
title('Brainstorm Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem6,bestval6,nfeval6,G6,Gi6,Gs6] = moa(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,alfa_moa,ro);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = moa(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,alfa,ro)

figure
subplot(221)
for i=1:tekrar,
    plot(G6(i,1:Gi6(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G6(i,1:Gi6(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN6=mean(G6(:,1:min(Gi6)));
    plot(MEAN6(1:min(Gi6)));
else
    MEAN6=G6;
    plot(MEAN6(1:min(Gi6)));
end
title('Magnetic-Inspired Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem7,bestval7,nfeval7,G7,Gi7,Gs7] = forest(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,drag,esik,yas_snr);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = forest(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,drag,esik,yas_snr)
figure
subplot(221)
for i=1:tekrar,
    plot(G7(i,1:Gi7(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G7(i,1:Gi7(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN7=mean(G7(:,1:min(Gi7)));
    plot(MEAN7(1:min(Gi7)));
else
    MEAN7=G7;
    plot(MEAN7(1:min(Gi7)));
end
title('Forest Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem8,bestval8,nfeval8,G8,Gi8,Gs8] = crab(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,receptivity,esikdegeri,NoCross,NoMut);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = crab(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,receptivity,esikdegeri,NoCross,NoMut)
figure
subplot(221)
for i=1:tekrar,
    plot(G8(i,1:Gi8(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G8(i,1:Gi8(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN8=mean(G8(:,1:min(Gi8)));
    plot(MEAN8(1:min(Gi8)));
else
    MEAN8=G8;
    plot(MEAN8(1:min(Gi8)));
end
title('Crab Mating Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem9,bestval9,nfeval9,G9,Gi9,Gs9] = spiral(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius)
figure
subplot(221)
for i=1:tekrar,
    plot(G9(i,1:Gi9(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G9(i,1:Gi9(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN9=mean(G9(:,1:min(Gi9)));
    plot(MEAN9(1:min(Gi9)));
else
    MEAN9=G9;
    plot(MEAN9(1:min(Gi9)));
end
title('Spiral Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem10,bestval10,nfeval10,G10,Gi10,Gs10] = flower(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,prob_switch,Levy_beta);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = flower(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,prob_switch,Levy_beta)
figure
subplot(221)
for i=1:tekrar,
    plot(G10(i,1:Gi10(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G10(i,1:Gi10(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN10=mean(G10(:,1:min(Gi10)));
    plot(MEAN10(1:min(Gi10)));
else
    MEAN10=G10;
    plot(MEAN10(1:min(Gi10)));
end
title('Flower Pollination Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem11,bestval11,nfeval11,G11,Gi11,Gs11] = fruitfly(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = fruitfly(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G11(i,1:Gi11(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G11(i,1:Gi11(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN11=mean(G11(:,1:min(Gi11)));
    plot(MEAN11(1:min(Gi11)));
else
    MEAN11=G11;
    plot(MEAN11(1:min(Gi11)));
end
title('Fruit Fly Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem12,bestval12,nfeval12,G12,Gi12,Gs12] = hurricane(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Rmax,w,R0);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = hurricane(fname,VTR,n,itMax,size,lb,ub,tekrar,refresh,Rmax,w,R0)
figure
subplot(221)
for i=1:tekrar,
    plot(G12(i,1:Gi12(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G12(i,1:Gi12(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iterations');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
if tekrar ~=1
    MEAN12=mean(G12(:,1:min(Gi12)));
    plot(MEAN12(1:min(Gi12)));
else
    MEAN12=G12;
    plot(MEAN12(1:min(Gi12)));
end
title('Hurricane-Based Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
% -------------------------------------------------------------------------------------   
fprintf(1,'\n%s probleminin sonucu : %3.3f\n',fonk,minimum);
for i=1:size(solution,1),
    fprintf(1,'x1 : %3.3f ve x2 : %3.3f \n',solution(i,1),solution(i,2));
end
% -------------------------------------------------------------------------------------
%.....Toplu sonuclar.....
% -------------------------------------------------------------------------------------
figure
plot(MEAN1(1:min(Gi1))','r');hold on
plot(MEAN2(1:min(Gi2))','b');
plot(MEAN3(1:min(Gi3))','g');
plot(MEAN4(1:min(Gi4))','m');
plot(MEAN5(1:min(Gi5))','k');
plot(MEAN6(1:min(Gi6))','c');
plot(MEAN7(1:min(Gi7))','Color',[1,0.4,0.6]);
plot(MEAN8(1:min(Gi8))','Color',[0.7,0.5,0.7]);
plot(MEAN9(1:min(Gi9))','Color',[0.25,0.5,0.75]);
plot(MEAN10(1:min(Gi10))','Color',[0.85,0.5,0.15]);
plot(MEAN11(1:min(Gi11))','Color',[0.25,0.75,0.65]);
plot(MEAN12(1:min(Gi12))','Color',[0.5,0.5,0.25]);
legend('Firefly','Differential','Bat','Cuckoo','Brainstorm','Magnetic-Inspired','Forest','Crab','Spiral','Flower','Fruitfly','Hurricane');
% legend('Differential','Bat','Cuckoo','Brainstorm','Magnetic-Inspired','Forest','Crab','Spiral','Flower','Fruitfly','Hurricane');
legend show
title('Heuristic Optimization Algorithms')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
axis tight;