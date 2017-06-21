% Sezgisel Algoritmalar Benchmark Testleri
% PSO, ABC, SA, DE ve TACO Algoritmalarý ile ödev olarak geliþtirilen 
% sezgisel algoritmanýn karþýlaþtýrýlmasý bu kodda yapýlmaktadýr.  
% Yazan : Uður Yüzgeç
% Tarih : Ocak 2017
clear all;
close all;
%%% legend('PSO','ABC','SA','DE','TACO');
% 
de=1;   % Differential Evolution (DE) Algorithm
pso=1;  % Particle Swarm Optimization (PSO) Algorithm
abc=1;  % Artificial Bee Colony (ABC) Algorithm
sa=0;    % Simulated Annealing (SA) Algorithm
taco=0; % Touring Ant Colony Optimization (TACO) algorithm 
xxx=0;  % sizin algoritmanýn çalýþmasý için 1 yapýlmasý yeterlidir.
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
% % % [24] levi13 : two dimension, many local minima
% % % [25] matyas : multi dimension, plate-shaped
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
% % % [38] styblinskitang : multi dimension, other
% % % [39] testtubeholder : two dimension, many local minima
% % % [40] threehumpcamel : two dimension, valley-shaped
% % % [41] trigonometric : Kullanmayýn...
% % % [42] wood : Kullanmayýn...
% % % [43] zettl : two dimension, plate-shaped 
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
       D=2;
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
        tekrar=5;
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
% -------------------------------------------------------------------------------------
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
%.....Toplu sonuclar.....
% -------------------------------------------------------------------------------------
figure
hold on
if de==1 plot(MEAN_DE(1:min(Gi_DE))','r'); end
if pso==1 plot(MEAN_PSO(1:min(Gi_PSO))','b'); end
if abc==1 plot(MEAN_ABC(1:min(Gi_ABC))','g'); end
if sa==1 plot(MEAN_SA(1:min(Gi_SA))','m'); end
if taco==1 plot(MEAN_TACO(1:min(Gi_TACO))','color',rand(1,3)); end
%-----------------------------------------------------------------------------------------
% Sizin sezgisel kodunuz (XXX) için yazýlacak...
if xxx==1 plot(MEAN_XXX(1:min(Gi_XXX))','color',rand(1,3)); end
%-----------------------------------------------------------------------------------------
legend('DE','PSO','ABC','SA','TACO','XXX');
legend show
title('Heuristic Optimization Algorithms')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
% ----------------------------------------------------------------------------------------
%SONUCLAR...
str = sprintf('Benchmark_%dtekrar_%dboyut_results.xlsx',tekrar,D);
data_baslik = {'Function','PSO','ABC','SA','DE','TACO','XXX'};
fonk_name = sprintf('A%d',num+1);
xlswrite(str,data_baslik, 1, 'A1');
xlswrite(str,data_baslik, 2, 'A1');
xlswrite(str,data_baslik, 3, 'A1');
xlswrite(str,data_baslik, 4, 'A1');
xlswrite(str,{fonk}, 1, fonk_name); xlswrite(str,{fonk}, 2, fonk_name);
xlswrite(str,{fonk}, 3, fonk_name); xlswrite(str,{fonk}, 4, fonk_name);
%-----------------------------------------------------------------------------------------------------------------
%  PSO Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
if pso==1
    Opt_PSO = 1-(sqrt((minimum-mean(bestval_PSO))^2))/(sqrt((UB-LB)^2)); % Optimality
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_PSO(i,j) = 1-(sqrt((solution(k,j)-bestmem_PSO(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_PSO(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_PSO(i,j))
						Acc_PSO(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nPSO Sonuçlarý\n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_PSO));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_PSO));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_PSO(:,1))/tekrar,sum(bestmem_PSO(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_PSO));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_PSO));
    fprintf(1,'Degerlendirilen amac fonk. sayýsý= %f\n',nfeval_PSO);
    fprintf(1,'Optimallik= %1.3f\n',Opt_PSO);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_PSO));
    data_PSO1 =sprintf('%3.3e(%3.3e)',mean(bestval_PSO),std(bestval_PSO));
    data_PSO2 =sprintf('%5.2f(%1.3f)',nfeval_PSO,mean(Gs_PSO));
    str2 = sprintf('B%d',num+1);
    xlswrite(str,{data_PSO1}, 1, str2);
    xlswrite(str,{data_PSO2}, 2, str2);
    xlswrite(str,Opt_PSO, 3, str2);
    xlswrite(str,mean(mean(Acc_PSO)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% ABC Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
if abc==1
    Opt_ABC=1-(sqrt((minimum-mean(bestval_ABC))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_ABC(i,j) = 1-(sqrt((solution(k,j)-bestmem_ABC(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_ABC(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_ABC(i,j))
						Acc_ABC(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nABC Sonuçlarý\n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_ABC));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_ABC));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_ABC(:,1))/tekrar,sum(bestmem_ABC(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_ABC));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_ABC));
    fprintf(1,'Degerlendirilen amac fonk. sayýsý= %f\n',nfeval_ABC);
    fprintf(1,'Optimallik= %1.3f\n',Opt_ABC);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_ABC));
    data_ABC1 =sprintf('%3.3e(%3.3e)',mean(bestval_ABC),std(bestval_ABC));
    data_ABC2 =sprintf('%5.2f(%1.3f)',nfeval_ABC,mean(Gs_ABC));
    str2 = sprintf('C%d',num+1);
    xlswrite(str,{data_ABC1}, 1, str2);
    xlswrite(str,{data_ABC2}, 2, str2);
    xlswrite(str,Opt_ABC, 3, str2);
    xlswrite(str,mean(mean(Acc_ABC)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% SA Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
if sa==1
    Opt_SA=1-(sqrt((minimum-mean(bestval_SA))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_SA(i,j) = 1-(sqrt((solution(k,j)-bestmem_SA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_SA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_SA(i,j))
						Acc_SA(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nSA Sonuçlarý\n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_SA));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_SA));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_SA(:,1))/tekrar,sum(bestmem_SA(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_SA));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_SA));
    fprintf(1,'Degerlendirilen amac fonk. sayýsý= %f\n',nfeval_SA);
    fprintf(1,'Optimallik= %1.3f\n',Opt_SA);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_SA));
    data_SA1 =sprintf('%3.3e(%3.3e)',mean(bestval_SA),std(bestval_SA));
    data_SA2 =sprintf('%5.2f(%1.3f)',nfeval_SA,mean(Gs_SA));
    str2 = sprintf('D%d',num+1);
    xlswrite(str,{data_SA1}, 1, str2);
    xlswrite(str,{data_SA2}, 2, str2);
    xlswrite(str,Opt_SA, 3, str2);
    xlswrite(str,mean(mean(Acc_SA)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% DE Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
if de==1
    Opt_DE=1-(sqrt((minimum-mean(bestval_DE))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_DE(i,j) = 1-(sqrt((solution(k,j)-bestmem_DE(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_DE(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_DE(i,j))
						Acc_DE(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nDE Sonuçlarý\n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_DE));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_DE));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_DE(:,1))/tekrar,sum(bestmem_DE(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_DE));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_DE));
    fprintf(1,'Degerlendirilen amac fonk. sayýsý= %f\n',nfeval_DE);
    fprintf(1,'Optimallik= %1.3f\n',Opt_DE);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_DE));
    data_DE1 =sprintf('%3.3e(%3.3e)',mean(bestval_DE),std(bestval_DE));
    data_DE2 =sprintf('%5.2f(%1.3f)',nfeval_DE,mean(Gs_DE));
    str2 = sprintf('E%d',num+1);
    xlswrite(str,{data_DE1}, 1, str2);
    xlswrite(str,{data_DE2}, 2, str2);
    xlswrite(str,Opt_DE, 3, str2);
    xlswrite(str,mean(mean(Acc_DE)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% TACO Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
if taco==1
    Opt_TACO=1-(sqrt((minimum-mean(bestval_TACO))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_TACO(i,j) = 1-(sqrt((solution(k,j)-bestmem_TACO(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_TACO(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_TACO(i,j))
						Acc_TACO(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nTACO Sonuçlarý\n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_TACO));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_TACO));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_TACO(:,1))/tekrar,sum(bestmem_TACO(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_TACO));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_TACO));
    fprintf(1,'Degerlendirilen amac fonk. sayýsý= %f\n',nfeval_TACO);
    fprintf(1,'Optimallik= %1.3f\n',Opt_TACO);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_TACO));
    data_TACO1 =sprintf('%3.3e(%3.3e)',mean(bestval_TACO),std(bestval_TACO));
    data_TACO2 =sprintf('%5.2f(%1.3f)',nfeval_TACO,mean(Gs_TACO));
    str2 = sprintf('F%d',num+1);
    xlswrite(str,{data_TACO1}, 1, str2);
    xlswrite(str,{data_TACO2}, 2, str2);
    xlswrite(str,Opt_TACO, 3, str2);
    xlswrite(str,mean(mean(Acc_TACO)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% XXX Sonuçlarý - Sizin kodunuzun sonuçlarý burada olacak... XXX yerine
% kodunuzun kýsaltmasý
%-----------------------------------------------------------------------------------------------------------------
if xxx==1
    Opt_XXX=1-(sqrt((minimum-mean(bestval_XXX))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_XXX(i,j) = 1-(sqrt((solution(k,j)-bestmem_XXX(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_XXX(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_XXX(i,j))
						Acc_XXX(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'XXX Sonuçlarý\n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_XXX));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_XXX));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_XXX(:,1))/tekrar,sum(bestmem_XXX(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_XXX));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_XXX));
    fprintf(1,'Degerlendirilen amac fonk. sayýsý= %f\n',nfeval_XXX);
    fprintf(1,'Optimallik= %1.3f\n',Opt_XXX);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_XXX));
    data_XXX1 =sprintf('%3.3e(%3.3e)',mean(bestval_XXX),std(bestval_XXX));
    data_XXX2 =sprintf('%5.2f(%1.3f)',nfeval_XXX,mean(Gs_XXX));
    str2 = sprintf('G%d',num+1);
    xlswrite(str,{data_XXX1}, 1, str2);
    xlswrite(str,{data_XXX2}, 2, str2);
    xlswrite(str,Opt_XXX, 3, str2);
    xlswrite(str,mean(mean(Acc_XXX)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------