% Sezgisel Optimizasyon Algoritmalari Nisan 2015
% Zuhal CALAYOGLU Meyve Sinegi Optimizasyon Algoritmasi (Fruit Fly Optimization Algorithm) 
% Adaptif karþýtlýk temelli yapi ilave edildi.
% EEB2016 konferans calismasi...

clear all;
close all;
prevpath = path;
path(path, genpath(fileparts(mfilename('fullpath'))));

% refresh       intermediate output will be produced after "refresh"
%               iterations. No intermediate output will be produced
%               if refresh is < 1
		refresh = 100; 
% VTR		"Value To Reach" (stop when ofunc < VTR)
		VTR = 1.e-6; 
% itermax       maximum number of iterations (generations)
		itermax = 500; 
% NP 			number of population (NP=10*D)		
        NP = 20;

[dims, lb, ub, solution, minimum, fonk]=RUN_ezimage;
        XVmin = lb; 
		XVmax = ub;
        D=dims;
        tekrar=50;
% -------------------------------------------------------------------------------------
[bestmem1,bestval1,nfeval1,G1,Gi1,Gs1] = fruitfly(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = fruitfly(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G1(i,1:Gi1(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
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
title('Meyve Sineði Algoritmasý')
grid on;
xlabel('iterasyon');
ylabel('ortalama maliyet')
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem2,bestval2,nfeval2,G2,Gi2,Gs2] = fruitfly_opposit(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = fruitfly_opposit(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G2(i,1:Gi2(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
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
title('Karþýtlýk temelli Meyve Sineði Algoritmasý')
grid on;
xlabel('iterasyon');
ylabel('ortalama maliyet')
axis tight;
drawnow;
% -------------------------------------------------------------------------------------
[bestmem3,bestval3,nfeval3,G3,Gi3,Gs3] = fruitfly_adapt_opposit(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = fruitfly_adapt_opposit(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G3(i,1:Gi3(i)),'color',rand(1,3));
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G3(i,1:Gi3(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
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
title('Adaptif karþýtlýk temelli Meyve Sineði Algoritmasý')
grid on;
xlabel('iterasyon');
ylabel('ortalama maliyet')
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
%legend('Firefly','Differential','Bat','Cuckoo','Brainstorm','Magnetic-Inspired','Forest','Crab','Spiral','Flower','Fruitfly','Hurricane');
% legend('Differential','Bat','Cuckoo','Brainstorm','Magnetic-Inspired','Forest','Crab','Spiral','Flower','Fruitfly','Hurricane');
legend('MSA','KMSA','AKMSA')
legend show
%title('Heuristic Optimization Algorithms')
grid on;
xlabel('iterasyon');
ylabel('ortalama maliyet');  
axis tight;