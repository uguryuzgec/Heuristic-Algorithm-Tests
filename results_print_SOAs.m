% -------------------------------------------------------------------------------------
%.....Toplu sonuclar.....
% -------------------------------------------------------------------------------------
figure
hold on
if de==1 plot(MEAN_DE(1:min(Gi_DE))','r'); end
if pso==1 plot(MEAN_PSO(1:min(Gi_PSO))','b'); end
if abc==1 plot(MEAN_ABC(1:min(Gi_ABC))','g'); end
if sa==1 plot(MEAN_SA(1:min(Gi_SA))','m'); end
if taco==1;plot(MEAN_TACO(1:min(Gi_TACO))','color',[1,0.5,0.5]); end
%-----------------------------------------------------------------------------------------
% Sizin sezgisel kodunuz (XXX) için yazýlacak...
if soa==1 plot(MEAN_SOA(1:min(Gi_SOA))','color','k'); end
if asoa_r==1 plot(MEAN_ASOA1(1:min(Gi_ASOA1))','color','c'); end
if asoa_theta==1 plot(MEAN_ASOA2(1:min(Gi_ASOA2))','color',rand(1,3)); end
if asoa==1 plot(MEAN_ASOA3(1:min(Gi_ASOA3))','color',rand(1,3)); end
if osoa==1 plot(MEAN_OSOA(1:min(Gi_OSOA))','color',rand(1,3)); end
if csoa==1 plot(MEAN_CSOA(1:min(Gi_CSOA))','color',rand(1,3)); end
%-----------------------------------------------------------------------------------------
% legend('PSO','ABC','SA','DE','TACO','SOA','ASOA1','ASOA2','ASOA3','OSOA','CSOA');
legend('SOA','ASOA1','ASOA2','ASOA3','OSOA','CSOA');
legend show
title('Heuristic Optimization Algorithms')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
axis tight;
magnifyOnFigure(...
        gcf,...
        'units', 'pixels',...
        'magnifierShape', 'rectangle',...
        'initialPositionSecondaryAxes',[150 200 250 150],...
        'initialPositionMagnifier',[0 0 250 100],...    
        'mode','interactive',...    
        'displayLinkStyle','none',...        
        'edgeWidth',0.5,...
        'edgeColor',[0.5 0.5 0.5],...
        'secondaryAxesFaceColor',[0.91 0.91 0.91],'displayLinkStyle','none'); 
% -----------------------------------------------------------------------------------------------------------------------------------
%Logaritmik plots
figure
if de==1;semilogy(MEAN_DE(1:min(Gi_DE))','r');hold on; end
if pso==1;semilogy(MEAN_PSO(1:min(Gi_PSO))','b');hold on; end
if abc==1;semilogy(MEAN_ABC(1:min(Gi_ABC))','g');hold on; end
if sa==1;semilogy(MEAN_SA(1:min(Gi_SA))','m');hold on; end
if taco==1;semilogy(MEAN_TACO(1:min(Gi_TACO))','color',[1,0.5,0.5]); hold on;end
%-----------------------------------------------------------------------------------------
%Sizin sezgisel kodunuz (XXX) için yazýlacak...
if soa==1 semilogy(MEAN_SOA(1:min(Gi_SOA))','color','k'); hold on; end
if asoa_r==1 semilogy(MEAN_ASOA1(1:min(Gi_ASOA1))','color','c'); hold on; end
if asoa_theta==1 semilogy(MEAN_ASOA2(1:min(Gi_ASOA2))','color',rand(1,3)); hold on; end
if asoa==1 semilogy(MEAN_ASOA3(1:min(Gi_ASOA3))','color',rand(1,3)); hold on; end
if osoa==1 semilogy(MEAN_OSOA(1:min(Gi_OSOA))','color',rand(1,3)); hold on; end
if csoa==1 semilogy(MEAN_CSOA(1:min(Gi_CSOA))','color',rand(1,3)); end
%-----------------------------------------------------------------------------------------
legend('DE','PSO','ABC','SA','TACO','ALO','IALO');
% legend('ALO','IALO_1','IALO_2');
%legend show
%title('Heuristic Optimization Algorithms')
grid on;
xlabel('iteration');
ylabel('mean of cost value f(x)');  
%axis tight;
% -----------------------------------------------------------------------------------------------------------------------------------
%SONUCLAR...
str = sprintf('Benchmark_%dtekrar_%dboyut_results.xlsx',tekrar,D);
data_baslik = {'Function','PSO','ABC','SA','DE','TACO','SOA','ASOA1','ASOA2','ASOA3','OSOA','CSOA'};
fonk_name = sprintf('A%d',num+1);
xlswrite(str,data_baslik, 1, 'A1');
xlswrite(str,data_baslik, 2, 'A1');
xlswrite(str,data_baslik, 3, 'A1');
xlswrite(str,data_baslik, 4, 'A1');
xlswrite(str,{fonk}, 1, fonk_name); xlswrite(str,{fonk}, 2, fonk_name);
xlswrite(str,{fonk}, 3, fonk_name); xlswrite(str,{fonk}, 4, fonk_name);
%-----------------------------------------------------------------------------------------------------------------
%  PSO Sonuclari 
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
    fprintf(1,'\nPSO Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_PSO));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_PSO));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_PSO(:,1))/tekrar,sum(bestmem_PSO(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_PSO));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_PSO));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_PSO);
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
% ABC Sonuclari 
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
    fprintf(1,'\nABC Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_ABC));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_ABC));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_ABC(:,1))/tekrar,sum(bestmem_ABC(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_ABC));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_ABC));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_ABC);
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
% SA Sonuclari 
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
    fprintf(1,'\nSA Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_SA));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_SA));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_SA(:,1))/tekrar,sum(bestmem_SA(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_SA));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_SA));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_SA);
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
% DE Sonuclari 
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
    fprintf(1,'\nDE Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_DE));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_DE));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_DE(:,1))/tekrar,sum(bestmem_DE(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_DE));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_DE));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_DE);
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
% TACO Sonuclari 
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
    fprintf(1,'\nTACO Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_TACO));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_TACO));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_TACO(:,1))/tekrar,sum(bestmem_TACO(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_TACO));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_TACO));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_TACO);
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
% SOA Sonuclari 
%-----------------------------------------------------------------------------------------------------------------
if soa==1
    Opt_SOA=1-(sqrt((minimum-mean(bestval_SOA))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_SOA(i,j) = 1-(sqrt((solution(k,j)-bestmem_SOA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_SOA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_SOA(i,j))
						Acc_SOA(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nSOA Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_SOA));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_SOA));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_SOA(:,1))/tekrar,sum(bestmem_SOA(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_SOA));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_SOA));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_SOA);
    fprintf(1,'Optimallik= %1.3f\n',Opt_SOA);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_SOA));
    data_SOA1 =sprintf('%3.3e(%3.3e)',mean(bestval_SOA),std(bestval_SOA));
    data_SOA2 =sprintf('%5.2f(%1.3f)',nfeval_SOA,mean(Gs_SOA));
    str2 = sprintf('G%d',num+1);
    xlswrite(str,{data_SOA1}, 1, str2);
    xlswrite(str,{data_SOA2}, 2, str2);
    xlswrite(str,Opt_SOA, 3, str2);
    xlswrite(str,mean(mean(Acc_SOA)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% ASOA1 Sonuclari 
%-----------------------------------------------------------------------------------------------------------------
if asoa_r==1
    Opt_ASOA1=1-(sqrt((minimum-mean(bestval_ASOA1))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_ASOA1(i,j) = 1-(sqrt((solution(k,j)-bestmem_ASOA1(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_ASOA1(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_ASOA1(i,j))
						Acc_ASOA1(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nASOA1 Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_ASOA1));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_ASOA1));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_ASOA1(:,1))/tekrar,sum(bestmem_ASOA1(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_ASOA1));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_ASOA1));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_ASOA1);
    fprintf(1,'Optimallik= %1.3f\n',Opt_ASOA1);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_ASOA1));
    data_ASOA11 =sprintf('%3.3e(%3.3e)',mean(bestval_ASOA1),std(bestval_ASOA1));
    data_ASOA12 =sprintf('%5.2f(%1.3f)',nfeval_ASOA1,mean(Gs_ASOA1));
    str2 = sprintf('H%d',num+1);
    xlswrite(str,{data_ASOA11}, 1, str2);
    xlswrite(str,{data_ASOA12}, 2, str2);
    xlswrite(str,Opt_ASOA1, 3, str2);
    xlswrite(str,mean(mean(Acc_ASOA1)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% ASOA2 Sonuclari 
%-----------------------------------------------------------------------------------------------------------------
if asoa_theta==1
    Opt_ASOA2=1-(sqrt((minimum-mean(bestval_ASOA2))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_ASOA2(i,j) = 1-(sqrt((solution(k,j)-bestmem_ASOA2(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_ASOA2(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_ASOA2(i,j))
						Acc_ASOA2(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nASOA2 Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_ASOA2));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_ASOA2));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_ASOA2(:,1))/tekrar,sum(bestmem_ASOA2(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_ASOA2));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_ASOA2));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_ASOA2);
    fprintf(1,'Optimallik= %1.3f\n',Opt_ASOA2);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_ASOA2));
    data_ASOA21 =sprintf('%3.3e(%3.3e)',mean(bestval_ASOA2),std(bestval_ASOA2));
    data_ASOA22 =sprintf('%5.2f(%1.3f)',nfeval_ASOA2,mean(Gs_ASOA2));
    str2 = sprintf('I%d',num+1);
    xlswrite(str,{data_ASOA21}, 1, str2);
    xlswrite(str,{data_ASOA22}, 2, str2);
    xlswrite(str,Opt_ASOA2, 3, str2);
    xlswrite(str,mean(mean(Acc_ASOA2)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% ASOA3 Sonuclari 
%-----------------------------------------------------------------------------------------------------------------
if asoa==1
    Opt_ASOA3=1-(sqrt((minimum-mean(bestval_ASOA3))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_ASOA3(i,j) = 1-(sqrt((solution(k,j)-bestmem_ASOA3(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_ASOA3(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_ASOA3(i,j))
						Acc_ASOA3(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nASOA3 Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_ASOA3));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_ASOA3));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_ASOA3(:,1))/tekrar,sum(bestmem_ASOA3(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_ASOA3));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_ASOA3));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_ASOA3);
    fprintf(1,'Optimallik= %1.3f\n',Opt_ASOA3);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_ASOA3));
    data_ASOA31 =sprintf('%3.3e(%3.3e)',mean(bestval_ASOA3),std(bestval_ASOA3));
    data_ASOA32 =sprintf('%5.2f(%1.3f)',nfeval_ASOA3,mean(Gs_ASOA3));
    str2 = sprintf('J%d',num+1);
    xlswrite(str,{data_ASOA31}, 1, str2);
    xlswrite(str,{data_ASOA32}, 2, str2);
    xlswrite(str,Opt_ASOA3, 3, str2);
    xlswrite(str,mean(mean(Acc_ASOA3)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% OSOA Sonuclari 
%-----------------------------------------------------------------------------------------------------------------
if osoa==1
    Opt_OSOA=1-(sqrt((minimum-mean(bestval_OSOA))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_OSOA(i,j) = 1-(sqrt((solution(k,j)-bestmem_OSOA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_OSOA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_OSOA(i,j))
						Acc_OSOA(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nOSOA Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_OSOA));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_OSOA));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_OSOA(:,1))/tekrar,sum(bestmem_OSOA(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_OSOA));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_OSOA));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_OSOA);
    fprintf(1,'Optimallik= %1.3f\n',Opt_OSOA);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_OSOA));
    data_OSOA1 =sprintf('%3.3e(%3.3e)',mean(bestval_OSOA),std(bestval_OSOA));
    data_OSOA2 =sprintf('%5.2f(%1.3f)',nfeval_OSOA,mean(Gs_OSOA));
    str2 = sprintf('K%d',num+1);
    xlswrite(str,{data_OSOA1}, 1, str2);
    xlswrite(str,{data_OSOA2}, 2, str2);
    xlswrite(str,Opt_OSOA, 3, str2);
    xlswrite(str,mean(mean(Acc_OSOA)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% CSOA Sonuclari 
%-----------------------------------------------------------------------------------------------------------------
if csoa==1
    Opt_CSOA=1-(sqrt((minimum-mean(bestval_CSOA))^2))/(sqrt((UB-LB)^2));
	for i=1:tekrar,
		for j=1:D,
		 temp = 0;
			for k=1:size(solution,1),
				if (k==1)
					Acc_CSOA(i,j) = 1-(sqrt((solution(k,j)-bestmem_CSOA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
				else 
					temp = 1-(sqrt((solution(k,j)-bestmem_CSOA(i,j))^2))/(sqrt((XVmax(j)-XVmin(j))^2)); % Accuracy
					if(temp > Acc_CSOA(i,j))
						Acc_CSOA(i,j) = temp;
					end
				end
			end
		end
	end
    fprintf(1,'\nCSOA Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_CSOA));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_CSOA));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_CSOA(:,1))/tekrar,sum(bestmem_CSOA(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_CSOA));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_CSOA));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_CSOA);
    fprintf(1,'Optimallik= %1.3f\n',Opt_CSOA);
	fprintf(1,'Dogruluk: ');
	disp(mean(Acc_CSOA));
    data_CSOA1 =sprintf('%3.3e(%3.3e)',mean(bestval_CSOA),std(bestval_CSOA));
    data_CSOA2 =sprintf('%5.2f(%1.3f)',nfeval_CSOA,mean(Gs_CSOA));
    str2 = sprintf('L%d',num+1);
    xlswrite(str,{data_CSOA1}, 1, str2);
    xlswrite(str,{data_CSOA2}, 2, str2);
    xlswrite(str,Opt_CSOA, 3, str2);
    xlswrite(str,mean(mean(Acc_CSOA)), 4, str2);
end
%-----------------------------------------------------------------------------------------------------------------
% XXX Sonuclari  - Sizin kodunuzun Sonuclari  burada olacak... XXX yerine
% kodunuzun kÄ±saltmasÄ±
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
    fprintf(1,'XXX Sonuclari \n');
    fprintf(1,'Ortalama Sonuc= %e\n',mean(bestval_XXX));
    fprintf(1,'Standart Sapma= %e\n',std(bestval_XXX));
    fprintf(1,'Ortalama Cozum= %f,%f\n',sum(bestmem_XXX(:,1))/tekrar,sum(bestmem_XXX(:,2))/tekrar);
    fprintf(1,'Ortalama iterasyon= %f\n',mean(Gi_XXX));
    fprintf(1,'Ortalama Sure= %fsn\n',mean(Gs_XXX));
    fprintf(1,'Degerlendirilen amac fonk.sayisi = %f\n',nfeval_XXX);
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