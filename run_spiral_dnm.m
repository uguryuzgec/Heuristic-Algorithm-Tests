% Sezgisel Optimizasyon Algoritmalari Subat 2016
% Tufan-Ugur calismalar 
% Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm)
% DEA eklendi (24.02.2016)
% ASOA (Adaptive Spiral Optimization Algorithm) eklendi (05.03.2016)

clear all;
close all;
soa=1;asoa1=1;asoa2=1;asoa3=1;
prevpath = path;
path(path, genpath(fileparts(mfilename('fullpath'))));
%--------------------------------------------------------------------------------------
% Refika ANBAR Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm)
% -------------------------------------------------------------------------------------
Theta=pi/4;
radius=0.95; 
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
		VTR = 1.e-6; 
% itermax       maximum number of iterations (generations)
		itermax = 200; 
% NP 			number of population (NP=10*D)		
        NP = 20;

[dims, lb, ub, solution, minimum, fonk]=RUN_ezimage;
        XVmin = lb; 
		XVmax = ub;
        D=dims;
        tekrar=50;
% -------------------------------------------------------------------------------------
if (soa==1)
[bestmem9,bestval9,nfeval9,G9,Gi9,Gs9] = spiral(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius);
% % %      function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius)      figure
figure
subplot(221)
for i=1:tekrar,
	plot(G9(i,1:Gi9(i)),'color',rand(1,3));
	vector(i,1) = sum(G9(i,1:min(Gi9)));
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
[best,best_index]=min(vector);
[worst,worst_index]=max(vector);
if tekrar ~=1
	MEAN9=mean(G9(:,1:min(Gi9)));
	plot(MEAN9(1:min(Gi9)),'b');
	hold on
	plot(G9(best_index,1:min(Gi9)),'r');
	plot(G9(worst_index,1:min(Gi9)),'c');
else
	MEAN9=G9;
	plot(MEAN9(1:min(Gi9)));
end
legend('mean','best','worst'); 
legend show
hold off
title('Spiral Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
drawnow;
end
% -------------------------------------------------------------------------------------
if(asoa3==1)
    [bestmem1,bestval1,nfeval1,G1,Gi1,Gs1] = spiral_adaptive(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral_adaptive(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G1(i,1:Gi1(i)),'color',rand(1,3));
	vector(i,1) = sum(G1(i,1:min(Gi1)));
	hold on	
end
grid on;
xlabel('iterations');
ylabel('cost value f(x)');  
axis tight;
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
[best,best_index]=min(vector);
[worst,worst_index]=max(vector);
if tekrar ~=1
    MEAN1=mean(G1(:,1:min(Gi1)));
    plot(MEAN1(1:min(Gi1)));
	hold on
	plot(G1(best_index,1:min(Gi1)),'r');
	plot(G1(worst_index,1:min(Gi1)),'c');
else
    MEAN1=G1;
    plot(MEAN1(1:min(Gi1)));
end
legend('mean','best','worst'); 
legend show
hold off
title('Adaptive Spiral Optimization Algorithm')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
drawnow;
end
% ------------------------------------------------------------------------------------- 
if(asoa1==1)
    [bestmem2,bestval2,nfeval2,G2,Gi2,Gs2] = spiral_adaptive_r(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral_adaptive_r(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G2(i,1:Gi2(i)),'color',rand(1,3));
	vector(i,1) = sum(G2(i,1:min(Gi2)));
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
[best,best_index]=min(vector);
[worst,worst_index]=max(vector);
if tekrar ~=1
    MEAN2=mean(G2(:,1:min(Gi2)));
    plot(MEAN2(1:min(Gi2)));
	hold on
	plot(G2(best_index,1:min(Gi2)),'r');
	plot(G2(worst_index,1:min(Gi2)),'c');
else
    MEAN2=G2;
    plot(MEAN2(1:min(Gi2)));
end
legend('mean','best','worst'); 
legend show
hold off
title('Adaptive Spiral Optimization Algorithm_1')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
drawnow;
end
% -------------------------------------------------------------------------------------
if(asoa2==1)
    [bestmem3,bestval3,nfeval3,G3,Gi3,Gs3] = spiral_adaptive_theta(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral_adaptive_theta(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
figure
subplot(221)
for i=1:tekrar,
    plot(G3(i,1:Gi3(i)),'color',rand(1,3));
	vector(i,1) = sum(G3(i,1:min(Gi3)));
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
[best,best_index]=min(vector);
[worst,worst_index]=max(vector);
if tekrar ~=1
    MEAN3=mean(G3(:,1:min(Gi3)));
    plot(MEAN3(1:min(Gi3)));
	hold on
	plot(G3(best_index,1:min(Gi3)),'r');
	plot(G3(worst_index,1:min(Gi3)),'c');
else
    MEAN3=G3;
    plot(MEAN3(1:min(Gi3)));
end
legend('mean','best','worst'); 
legend show
hold off
title('Adaptive Spiral Optimization Algorithm_2')
grid on;
xlabel('iterations');
ylabel('mean of cost value f(x)');  
drawnow;
end
% ------------------------------------------------------------------------------------- 
fprintf(1,'\n%s probleminin sonucu : %3.3f\n',fonk,minimum);
for i=1:size(solution,1),
    fprintf(1,'x1 : %3.3f ve x2 : %3.3f \n',solution(i,1),solution(i,2));
end
% -------------------------------------------------------------------------------------
%.....Toplu sonuclar.....
% -------------------------------------------------------------------------------------
figure
if (soa==1) plot(MEAN9(1:min(Gi9))','g'); hold on;
end
if (asoa1==1) plot(MEAN2(1:min(Gi2))','b'); end
if (asoa2==1) plot(MEAN3(1:min(Gi3))','r'); end
if(asoa3==1) plot(MEAN1(1:min(Gi1))','m'); end
% legend('Firefly','Differential','Bat','Cuckoo','Brainstorm','Magnetic-Inspired','Forest','Crab','Spiral','Flower','Fruitfly','Hurricane');
legend('SOA','ASOA_1','ASOA_2','ASOA_3'); %'HOA','DEA','SOA'
legend show
title('Adaptive Spiral Optimization Algorithms')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;