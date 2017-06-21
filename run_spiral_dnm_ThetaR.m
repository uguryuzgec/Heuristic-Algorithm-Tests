% Sezgisel Optimizasyon Algoritmalari Subat 2016
% Tufan-Ugur calismalar 
% Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm)
% DEA eklendi (24.02.2016)
% ASOA (Adaptive Spiral Optimization Algorithm) eklendi (05.03.2016)

clear all;
close all;
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
        tekrar=1;
% -------------------------------------------------------------------------------------
figure
for radius=0.9:0.01:1,
% for Theta=0:(pi/6):pi*2,
rand('seed',2.5);
    [bestmem9,bestval9,nfeval9,G9,Gi9,Gs9] = spiral(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius);
% % %      function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius)      figure

% rand('seed',Theta);
rand('seed',radius*200);
MEAN9=G9;
	plot(MEAN9(1:min(Gi9)),'color',rand(1,3));
hold on
end
legend('0.9','0.91','0.92','0.93','0.94','0.95','0.96','0.97','0.98','0.99','1');
% legend('0','30','60','90','120','150','180','210','240','270','300','330','360');
legend show;
% % %     
% % % subplot(221)
% % % for i=1:tekrar,
% % % 	plot(G9(i,1:Gi9(i)),'color',rand(1,3));
% % % 	vector(i,1) = sum(G9(i,1:min(Gi9)));
% % % 	hold on			   
% % % end
% % % grid on;
% % % xlabel('iterations');
% % % ylabel('cost value f(x)');  
% % % axis tight;
% % % subplot(222)
% % % for i=1:tekrar,
% % % 	plot((G9(i,1:Gi9(i))-minimum).^2,'color',rand(1,3));
% % % 	hold on
% % % end   
% % % hold off
% % % grid on;
% % % xlabel('iterations');
% % % ylabel('squared error');   
% % % axis tight;
% % % subplot(2,2,[3 4])
% % % [best,best_index]=min(vector);
% % % [worst,worst_index]=max(vector);
% % % if tekrar ~=1
% % % 	MEAN9=mean(G9(:,1:min(Gi9)));
% % % 	plot(MEAN9(1:min(Gi9)),'b');
% % % 	hold on
% % % 	plot(G9(best_index,1:min(Gi9)),'r');
% % % 	plot(G9(worst_index,1:min(Gi9)),'c');
% % % else
% % % 	MEAN9=G9;
% % % 	plot(MEAN9(1:min(Gi9)));
% % % end
% % % legend('mean','best','worst'); 
% % % legend show
% % % hold off
% % % title('Spiral Optimization Algorithm')
% % % grid on;
% % % xlabel('iterations');
% % % ylabel('mean of cost value f(x)');  
% % % drawnow;
% % % % -------------------------------------------------------------------------------------
% % % 	%	[bestmem12,bestval12,nfeval12,G12,Gi12,Gs12] = hurricane(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Rmax,w,R0);
% % % 	%	% function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = hurricane(fname,VTR,n,itMax,size,lb,ub,tekrar,refresh,Rmax,w,R0)
% % % 	%	figure
% % % 	%	subplot(221)
% % % 	%	for i=1:tekrar,
% % % 	%		plot(G12(i,1:Gi12(i)),'color',rand(1,3));
% % % 	%		hold on
% % % 	%	end
% % % 	%	grid on;
% % % 	%	xlabel('iterations');
% % % 	%	ylabel('cost value f(x)');  
% % % 	%	axis tight;
% % % 	%	subplot(222)
% % % 	%	for i=1:tekrar,
% % % 	%		plot((G12(i,1:Gi12(i))-minimum).^2,'color',rand(1,3));
% % % 	%		hold on
% % % 	%	end   
% % % 	%	hold off
% % % 	%	grid on;
% % % 	%	xlabel('iterations');
% % % 	%	ylabel('squared error');   
% % % 	%	axis tight;
% % % 	%	subplot(2,2,[3 4])
% % % 	%	if tekrar ~=1
% % % 	%		MEAN12=mean(G12(:,1:min(Gi12)));
% % % 	%		plot(MEAN12(1:min(Gi12)));
% % % 	%	else
% % % 	%		MEAN12=G12;
% % % 	%		plot(MEAN12(1:min(Gi12)));
% % % 	%	end
% % % 	%	title('Hurricane-Based Optimization Algorithm')
% % % 	%	grid on;
% % % 	%	xlabel('iterations');
% % % 	%	ylabel('mean of cost value f(x)');  
% % % 	%	axis tight;
% % % 	%	drawnow;
% % % % -------------------------------------------------------------------------------------
% % % 	%	[bestmem2,bestval2,nfeval2,G2,Gi2,Gs2]=devec3_orj(fonk,VTR,D,XVmin,XVmax,NP,itermax,F,CR,strategy,refresh,tekrar);
% % % 	%	figure
% % % 	%	subplot(221)
% % % 	%	for i=1:tekrar,
% % % 	%		plot(G2(i,1:Gi2(i)),'color',rand(1,3));
% % % 	%		hold on
% % % 	%	end
% % % 	%	grid on;
% % % 	%	xlabel('iterations');
% % % 	%	ylabel('cost value f(x)');  
% % % 	%	axis tight;
% % % 	%	subplot(222)
% % % 	%	for i=1:tekrar,
% % % 	%		plot((G2(i,1:Gi2(i))-minimum).^2,'color',rand(1,3));
% % % 	%		hold on
% % % 	%	end   
% % % 	%	hold off
% % % 	%	grid on;
% % % 	%	xlabel('iterations');
% % % 	%	ylabel('squared error');   
% % % 	%	axis tight;
% % % 	%	subplot(2,2,[3 4])
% % % 	%	if tekrar ~=1
% % % 	%		MEAN2=mean(G2(:,1:min(Gi2)));
% % % 	%		plot(MEAN2(1:min(Gi2)));
% % % 	%	else
% % % 	%		MEAN2=G2;
% % % 	%		plot(MEAN2(1:min(Gi2)));
% % % 	%	end
% % % 	%	title('Differential Evolution Algorithm')
% % % 	%	grid on;
% % % 	%	xlabel('iterations');
% % % 	%	ylabel('mean of cost value f(x)');  
% % % 	%	axis tight;
% % % 	%	drawnow; 
% % % % -------------------------------------------------------------------------------------
% % % [bestmem1,bestval1,nfeval1,G1,Gi1,Gs1] = spiral_adaptive(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
% % % % function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral_adaptive(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
% % % figure
% % % subplot(221)
% % % for i=1:tekrar,
% % %     plot(G1(i,1:Gi1(i)),'color',rand(1,3));
% % % 	vector(i,1) = sum(G1(i,1:min(Gi1)));
% % % 	hold on	
% % % end
% % % grid on;
% % % xlabel('iterations');
% % % ylabel('cost value f(x)');  
% % % axis tight;
% % % subplot(222)
% % % for i=1:tekrar,
% % %     plot((G1(i,1:Gi1(i))-minimum).^2,'color',rand(1,3));
% % %     hold on
% % % end   
% % % hold off
% % % grid on;
% % % xlabel('iterations');
% % % ylabel('squared error');   
% % % axis tight;
% % % subplot(2,2,[3 4])
% % % [best,best_index]=min(vector);
% % % [worst,worst_index]=max(vector);
% % % if tekrar ~=1
% % %     MEAN1=mean(G1(:,1:min(Gi1)));
% % %     plot(MEAN1(1:min(Gi1)));
% % % 	hold on
% % % 	plot(G1(best_index,1:min(Gi1)),'r');
% % % 	plot(G1(worst_index,1:min(Gi1)),'c');
% % % else
% % %     MEAN1=G1;
% % %     plot(MEAN1(1:min(Gi1)));
% % % end
% % % legend('mean','best','worst'); 
% % % legend show
% % % hold off
% % % title('Adaptive Spiral Optimization Algorithm')
% % % grid on;
% % % xlabel('iterations');
% % % ylabel('mean of cost value f(x)');  
% % % drawnow;
% % % % -------------------------------------------------------------------------------------  
% % % fprintf(1,'\n%s probleminin sonucu : %3.3f\n',fonk,minimum);
% % % for i=1:size(solution,1),
% % %     fprintf(1,'x1 : %3.3f ve x2 : %3.3f \n',solution(i,1),solution(i,2));
% % % end
% % % % -------------------------------------------------------------------------------------
% % % %.....Toplu sonuclar.....
% % % % -------------------------------------------------------------------------------------
% % % figure
% % % plot(MEAN9(1:min(Gi9))','g'); hold on;
% % % %plot(MEAN12(1:min(Gi12))','r'); 
% % % %plot(MEAN2(1:min(Gi2))','b');
% % % plot(MEAN1(1:min(Gi1))','m'); 
% % % 
% % % % legend('Firefly','Differential','Bat','Cuckoo','Brainstorm','Magnetic-Inspired','Forest','Crab','Spiral','Flower','Fruitfly','Hurricane');
% % % legend('SOA','ASOA'); %'HOA','DEA','SOA'
% % % legend show
% % % title('Adaptive Spiral Optimization Algorithms')
% % % grid on;
% % % xlabel('iteration');
% % % ylabel('mean of cost value f(x)');  
% % % axis tight;