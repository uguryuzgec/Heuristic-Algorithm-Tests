% Sezgisel Algoritmalar Benchmark Testleri
% PSO, ABC, SA, DE ve TACO Algoritmalarý ile ödev olarak geliþtirilen 
% sezgisel algoritmanýn karþýlaþtýrýlmasý bu kodda yapýlmaktadýr.  
% Yazan : Uður Yüzgeç
% Tarih : Ocak 2017
clear all;
close all;
%%% legend('PSO','ABC','SA','DE','TACO','SOA','ASOA1','ASOA2','ASOA3','OSOA','CSOA');
% 
de=1;   		% Differential Evolution (DE) Algorithm
pso=0;  		% Particle Swarm Optimization (PSO) Algorithm
abc=0;  		% Artificial Bee Colony (ABC) Algorithm
sa=0;   		% Simulated Annealing (SA) Algorithm
taco=0; 		% Touring Ant Colony Optimization (TACO) algorithm 
soa=0;  		% Spiral Optimization Algorithm 
asoa_r=0; 		% Adaptive Spiral Radius based Spiral Optimization Algorithm 
asoa_theta=0; 	% Adaptive Spiral Angle based Spiral Optimization Algorithm 
asoa=0; 		% Full Adaptive Spiral Optimization Algorithm 
osoa=0;			% Opposition based Spiral Optimization Algorithm 
csoa=0;			% Chaotic maps based Spiral Optimization Algorithm 
xxx=0;
%--------------------------------------------------------------------------------------------------------------------------------------------------
% % % Single-objective functions:
% % % [1] ackley : multi dimension, many local minima
% % % [2] beale : two dimension, other
% % % [3] bird :  two dimension, other
% % % [4] booth : two dimension, plate-shaped 
% % % [5] bukin2 : Kullanmayýn...
% % % [6] bukin4 : Kullanmayýn...
% % % [7] bukin6 : two dimension, many local minima
% % % [8] carromtable : two dimension, many local minima
% % % [9] chichinadze : two dimension, other
% % % [10] crossfunc : two dimension, other
% % % [11] crossintray : two dimension, many local minima
% % % [12] crosslegtable : two dimension, other
% % % [13] crownedcross : two dimension, many local minima
% % % [14] cube : two dimension, plate-shaped
% % % [15] easom : Kullanmayýn... (two dimension, steep ridges/drops)
% % % [16] eggholder : multi dimension, many local minima
% % % [17] giunta : two dimension, other
% % % [18] goldsteinprice : two dimension, other
% % % [19] griewank : multi dimension, many local minima
% % % [20] helicalvalley : Kullanmayýn... three dimension, other
% % % [21] himmelblau : two dimension, plate-shaped 
% % % [22] holdertable : two dimension, many local minima
% % % [23] leon : two dimension, valley-shaped
% % % [24] levy : multi dimension, many local minima
% % % [25] matyas : two dimension, plate-shaped
% % % [26] mccormick : two dimension, plate-shaped 
% % % [27] modschaffer1 : two dimension, many local minima
% % % [28] modschaffer2 : two dimension, many local minima
% % % [29] modschaffer3 : two dimension, many local minima
% % % [30] modschaffer4 : two dimension, many local minima
% % % [31] penholder : two dimension, many local minima
% % % [32] powell : Kullanmayýn...
% % % [33] rastrigin : multi dimension, many local minima
% % % [34] rosenbrock : multi dimension, valley-shaped
% % % [35] schweffel : multi dimension, many local minima
% % % [36] sinenvsin : multi dimension, other
% % % [37] sixhumpcamel : two dimension, valley-shaped
% % % [38] sphere : multi dimension, valley-shaped
% % % [39] styblinskitang : multi dimension, other
% % % [40] sum2 : multi dimension, valley-shaped
% % % [41] testtubeholder : two dimension, many local minima
% % % [42] threehumpcamel : two dimension, valley-shaped
% % % [43] trigonometric : Kullanmayýn...
% % % [44] wood : Kullanmayýn...
% % % [45] zakh : multi dimension, valley-shaped
% % % [46] zettl : two dimension, plate-shaped 
%--------------------------------------------------------------------------------------------------------------------------------------------------

prevpath = path;
path(path, genpath(fileparts(mfilename('fullpath'))));
% refresh       intermediate output will be produced after "refresh"
%               iteration. No intermediate output will be produced
%               if refresh is < 1
		refresh = 100; 
% VTR		"Value To Reach" (stop when ofunc < VTR)
		VTR = 1.e-6; 
% itermax       maximum number of iteration (generations)
		itermax = 1000; 

[dims, lb, ub, solution, minimum,fonk, LB, UB, num]=RUN_ezimage;		
       D=4;
% NP number of population (NP=10*D)		
        NP = 10*D;
	   if D == 2
			XVmin = lb; 
			XVmax = ub;
       else
			XVmin(1:D) = lb(1); 
			XVmax(1:D) = ub(1);
			solution(1:D) = solution(1);
	   end
        F=0.8;
        CR=0.5; 
        tekrar=10;
		strategy=7; %7 digeri...
% strategy   1 --> DE/best/1/exp              6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp              7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp              9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp              else  DE/rand/2/bin
%-----------------------------------------------------------------------------------------------------------------
% PSO parametreleri
    c1=2.05;c2=2.05; 
%----------------------------------------------------------------------------------------------------------------- 
% Refika ANBAR Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm)
%-----------------------------------------------------------------------------------------------------------------
Theta=pi/4;
radius=0.95; 
%-----------------------------------------------------------------------------------------------------------------   
if de==1
[bestmem_DE,bestval_DE,nfeval_DE,G_DE,Gi_DE,Gs_DE]=devec3_orj(fonk,VTR,D,XVmin,XVmax,NP,itermax,F,CR,strategy,refresh,tekrar,solution);
figure
subplot(221)
for i=1:tekrar,
    plot(G_DE(i,1:Gi_DE(i)),'color',rand(1,3));
	vector(i,1) = sum(G_DE(i,1:min(Gi_DE))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_DE(i,1:Gi_DE(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_DE=mean(G_DE(:,1:min(Gi_DE)));
    plot(MEAN_DE(1:min(Gi_DE)),'b');
	hold on
	plot(G_DE(best_index,1:min(Gi_DE)),'r'); % eklendi...
	plot(G_DE(worst_index,1:min(Gi_DE)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_DE=G_DE;
    plot(MEAN_DE(1:min(Gi_DE)));
end
title('Differential Evolution Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
if pso==1 
[bestmem_PSO,bestval_PSO,nfeval_PSO,G_PSO,Gi_PSO,Gs_PSO]=PSO(fonk,VTR,D,XVmin,XVmax,NP,itermax,c1,c2,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_PSO(i,1:Gi_PSO(i)),'color',rand(1,3));
	vector(i,1) = sum(G_PSO(i,1:min(Gi_PSO))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_PSO(i,1:Gi_PSO(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_PSO=mean(G_PSO(:,1:min(Gi_PSO)));
    plot(MEAN_PSO(1:min(Gi_PSO)),'b');
	hold on
	plot(G_PSO(best_index,1:min(Gi_PSO)),'r'); % eklendi...
	plot(G_PSO(worst_index,1:min(Gi_PSO)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_PSO=G_PSO;
    plot(MEAN_PSO(1:min(Gi_PSO)));
end
title('Particle Swarm Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow;
end 
%-----------------------------------------------------------------------------------------------------------------
if abc==1 
[bestmem_ABC,bestval_ABC,nfeval_ABC,G_ABC,Gi_ABC,Gs_ABC]=ABC(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_ABC(i,1:Gi_ABC(i)),'color',rand(1,3));
	vector(i,1) = sum(G_ABC(i,1:min(Gi_ABC))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_ABC(i,1:Gi_ABC(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_ABC=mean(G_ABC(:,1:min(Gi_ABC)));
    plot(MEAN_ABC(1:min(Gi_ABC)),'b');
	hold on
	plot(G_ABC(best_index,1:min(Gi_ABC)),'r'); % eklendi...
	plot(G_ABC(worst_index,1:min(Gi_ABC)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_ABC=G_ABC;
    plot(MEAN_ABC(1:min(Gi_ABC)));
end
title('Artificial Bee Colony Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
if sa==1 
[bestmem_SA,bestval_SA,nfeval_SA,G_SA,Gi_SA,Gs_SA]=SA(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_SA(i,1:Gi_SA(i)),'color',rand(1,3));
	vector(i,1) = sum(G_SA(i,1:min(Gi_SA))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_SA(i,1:Gi_SA(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_SA=mean(G_SA(:,1:min(Gi_SA)));
    plot(MEAN_SA(1:min(Gi_SA)),'b');
	hold on
	plot(G_SA(best_index,1:min(Gi_SA)),'r'); % eklendi...
	plot(G_SA(worst_index,1:min(Gi_SA)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_SA=G_SA;
    plot(MEAN_SA(1:min(Gi_SA)));
end
title('Simulated Annealing Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
if taco==1 
[bestmem_TACO,bestval_TACO,nfeval_TACO,G_TACO,Gi_TACO,Gs_TACO]=TACO(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_TACO(i,1:Gi_TACO(i)),'color',rand(1,3));
	vector(i,1) = sum(G_TACO(i,1:min(Gi_TACO))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_TACO(i,1:Gi_TACO(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_TACO=mean(G_TACO(:,1:min(Gi_TACO)));
    plot(MEAN_TACO(1:min(Gi_TACO)),'b');
	hold on
	plot(G_TACO(best_index,1:min(Gi_TACO)),'r'); % eklendi...
	plot(G_TACO(worst_index,1:min(Gi_TACO)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_TACO=G_TACO;
    plot(MEAN_TACO(1:min(Gi_TACO)));
end
title('Touring Ant Colony Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
% Spiral Optimization Algoritmasý - SOA
%-----------------------------------------------------------------------------------------------------------------
if soa==1 
[bestmem_SOA,bestval_SOA,nfeval_SOA,G_SOA,Gi_SOA,Gs_SOA]=spiral(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius);
figure
subplot(221)
for i=1:tekrar,
    plot(G_SOA(i,1:Gi_SOA(i)),'color',rand(1,3));
	vector(i,1) = sum(G_SOA(i,1:min(Gi_SOA))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_SOA(i,1:Gi_SOA(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_SOA=mean(G_SOA(:,1:min(Gi_SOA)));
    plot(MEAN_SOA(1:min(Gi_SOA)),'b');
	hold on
	plot(G_SOA(best_index,1:min(Gi_SOA)),'r'); % eklendi...
	plot(G_SOA(worst_index,1:min(Gi_SOA)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_SOA=G_SOA;
    plot(MEAN_SOA(1:min(Gi_SOA)));
end
title('Spiral Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
% Adaptive Spiral Radius based Spiral Optimization Algorithm - ASOA1
%-----------------------------------------------------------------------------------------------------------------
if asoa_r==1 
[bestmem_ASOA1,bestval_ASOA1,nfeval_ASOA1,G_ASOA1,Gi_ASOA1,Gs_ASOA1]=spiral_adaptive_r(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_ASOA1(i,1:Gi_ASOA1(i)),'color',rand(1,3));
	vector(i,1) = sum(G_ASOA1(i,1:min(Gi_ASOA1))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_ASOA1(i,1:Gi_ASOA1(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_ASOA1=mean(G_ASOA1(:,1:min(Gi_ASOA1)));
    plot(MEAN_ASOA1(1:min(Gi_ASOA1)),'b');
	hold on
	plot(G_ASOA1(best_index,1:min(Gi_ASOA1)),'r'); % eklendi...
	plot(G_ASOA1(worst_index,1:min(Gi_ASOA1)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_ASOA1=G_ASOA1;
    plot(MEAN_ASOA1(1:min(Gi_ASOA1)));
end
title('Adaptive Spiral Radius based Spiral Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
% Adaptive Spiral Angle based Spiral Optimization Algorithm - ASOA2
%-----------------------------------------------------------------------------------------------------------------
if asoa_theta==1 
[bestmem_ASOA2,bestval_ASOA2,nfeval_ASOA2,G_ASOA2,Gi_ASOA2,Gs_ASOA2]=spiral_adaptive_theta(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_ASOA2(i,1:Gi_ASOA2(i)),'color',rand(1,3));
	vector(i,1) = sum(G_ASOA2(i,1:min(Gi_ASOA2))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_ASOA2(i,1:Gi_ASOA2(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_ASOA2=mean(G_ASOA2(:,1:min(Gi_ASOA2)));
    plot(MEAN_ASOA2(1:min(Gi_ASOA2)),'b');
	hold on
	plot(G_ASOA2(best_index,1:min(Gi_ASOA2)),'r'); % eklendi...
	plot(G_ASOA2(worst_index,1:min(Gi_ASOA2)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_ASOA2=G_ASOA2;
    plot(MEAN_ASOA2(1:min(Gi_ASOA2)));
end
title('Adaptive Spiral Angle based Spiral Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
% Full Adaptive Spiral Optimization Algorithm - ASOA3
%-----------------------------------------------------------------------------------------------------------------
if asoa==1 
[bestmem_ASOA3,bestval_ASOA3,nfeval_ASOA3,G_ASOA3,Gi_ASOA3,Gs_ASOA3]=spiral_adaptive(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_ASOA3(i,1:Gi_ASOA3(i)),'color',rand(1,3));
	vector(i,1) = sum(G_ASOA3(i,1:min(Gi_ASOA3))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_ASOA3(i,1:Gi_ASOA3(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_ASOA3=mean(G_ASOA3(:,1:min(Gi_ASOA3)));
    plot(MEAN_ASOA3(1:min(Gi_ASOA3)),'b');
	hold on
	plot(G_ASOA3(best_index,1:min(Gi_ASOA3)),'r'); % eklendi...
	plot(G_ASOA3(worst_index,1:min(Gi_ASOA3)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_ASOA3=G_ASOA3;
    plot(MEAN_ASOA3(1:min(Gi_ASOA3)));
end
title('Full Adaptive Spiral Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
% Opposition based Spiral Optimization Algorithm - OSOA
%-----------------------------------------------------------------------------------------------------------------
if osoa==1 
[bestmem_OSOA,bestval_OSOA,nfeval_OSOA,G_OSOA,Gi_OSOA,Gs_OSOA]=spiral_opposition(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius);
figure
subplot(221)
for i=1:tekrar,
    plot(G_OSOA(i,1:Gi_OSOA(i)),'color',rand(1,3));
	vector(i,1) = sum(G_OSOA(i,1:min(Gi_OSOA))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_OSOA(i,1:Gi_OSOA(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_OSOA=mean(G_OSOA(:,1:min(Gi_OSOA)));
    plot(MEAN_OSOA(1:min(Gi_OSOA)),'b');
	hold on
	plot(G_OSOA(best_index,1:min(Gi_OSOA)),'r'); % eklendi...
	plot(G_OSOA(worst_index,1:min(Gi_OSOA)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_OSOA=G_OSOA;
    plot(MEAN_OSOA(1:min(Gi_OSOA)));
end
title('Opposition based Spiral Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
% Chaotic maps based Spiral Optimization Algorithm - CSOA
%-----------------------------------------------------------------------------------------------------------------
if csoa==1 
%-----------------------------------------------------------------------------------------------------------------
%--chaotic_option is used for selecting the chaotic function---%
%  1 : Logistic map, 2: Henon map, 3: Tent map, 4: Sinusoidal iterator, 5: Circle map %
	option = 2;
%-----------------------------------------------------------------------------------------------------------------
[bestmem_CSOA,bestval_CSOA,nfeval_CSOA,G_CSOA,Gi_CSOA,Gs_CSOA]=spiral_chaotic_maps(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,radius,option,solution);
figure
subplot(221)
for i=1:tekrar,
    plot(G_CSOA(i,1:Gi_CSOA(i)),'color',rand(1,3));
	vector(i,1) = sum(G_CSOA(i,1:min(Gi_CSOA))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_CSOA(i,1:Gi_CSOA(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_CSOA=mean(G_CSOA(:,1:min(Gi_CSOA)));
    plot(MEAN_CSOA(1:min(Gi_CSOA)),'b');
	hold on
	plot(G_CSOA(best_index,1:min(Gi_CSOA)),'r'); % eklendi...
	plot(G_CSOA(worst_index,1:min(Gi_CSOA)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_CSOA=G_CSOA;
    plot(MEAN_CSOA(1:min(Gi_CSOA)));
end
title('Chaotic maps based Spiral Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
%-----------------------------------------------------------------------------------------------------------------
% XXX Algoritmasý - Sizin kodunuzu yukarýdaki gibi buraya yazacaksýnýz... XXX yerine
% kodunuzun kýsaltmasý gelecek.
%-----------------------------------------------------------------------------------------------------------------
if xxx==1 
[bestmem_XXX,bestval_XXX,nfeval_XXX,G_XXX,Gi_XXX,Gs_XXX]=XXX(fonk,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
figure
subplot(221)
for i=1:tekrar,
    plot(G_XXX(i,1:Gi_XXX(i)),'color',rand(1,3));
	vector(i,1) = sum(G_XXX(i,1:min(Gi_XXX))); % eklendi...
    hold on
end
grid on;
xlabel('iteration');
ylabel('cost value f(x)');  
axis tight;
subplot(222)
for i=1:tekrar,
    plot((G_XXX(i,1:Gi_XXX(i))-minimum).^2,'color',rand(1,3));
    hold on
end   
hold off
grid on;
xlabel('iteration');
ylabel('squared error');   
axis tight;
subplot(2,2,[3 4])
[best,best_index]=min(vector); % eklendi...
[worst,worst_index]=max(vector); % eklendi...
if tekrar ~=1
    MEAN_XXX=mean(G_XXX(:,1:min(Gi_XXX)));
    plot(MEAN_XXX(1:min(Gi_XXX)),'b');
	hold on
	plot(G_XXX(best_index,1:min(Gi_XXX)),'r'); % eklendi...
	plot(G_XXX(worst_index,1:min(Gi_XXX)),'c'); % eklendi...
	legend('mean','best','worst');
else
    MEAN_XXX=G_XXX;
    plot(MEAN_XXX(1:min(Gi_XXX)));
end
title('XXX Optimization Algorithm')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
drawnow; 
end
% -------------------------------------------------------------------------------------
fprintf(1,'\n%s probleminin sonucu : %3.3f\n',fonk,minimum);
for i=1:size(solution,1),
    fprintf(1,'x1 : %3.3f ve x2 : %3.3f \n',solution(i,1),solution(i,2));
end
% -------------------------------------------------------------------------------------
% results_print_SOAs;