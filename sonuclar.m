%SONUCLAR...
%%%bestmem6,bestval6,nfeval6,G6,Gi6,Gs6
%%% GBestmem,GBestval,nfeval1,Gmin,Giter,GSure

fprintf(1,'\n%s probleminin sonucu : %3.3f\n',fonk,minimum);
for i=1:size(solution,1),
    fprintf(1,'x1 : %3.3f ve x2 : %3.3f \n',solution(i,1),solution(i,2));
end
%-----------------------------------------------------------------------------------------------------------------
% SOA Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
fprintf(1,'\nSOA Sonuçlarý\n');
fprintf(1,'\nOrtalama Sonuc= %e\n',mean(bestval9));
fprintf(1,'\nStandart Sapma= %e\n',std(bestval9));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(bestmem9(:,1))/tekrar,sum(bestmem9(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(Gi9));
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(Gs9));
fprintf(1,'\nnumber of function ev= %f\n',nfeval9);
%-----------------------------------------------------------------------------------------------------------------
% ASOA_1 Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
fprintf(1,'\nASOA_1 Sonuçlarý\n');
fprintf(1,'\nOrtalama Sonuc= %e\n',mean(bestval2));
fprintf(1,'\nStandart Sapma= %e\n',std(bestval2));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(bestmem2(:,1))/tekrar,sum(bestmem2(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(Gi2));
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(Gs2));
fprintf(1,'\nnumber of function ev= %f\n',nfeval2);
%-----------------------------------------------------------------------------------------------------------------
% ASOA_2 Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
fprintf(1,'\nASOA_2 Sonuçlarý\n');
fprintf(1,'\nOrtalama Sonuc= %e\n',mean(bestval3));
fprintf(1,'\nStandart Sapma= %e\n',std(bestval3));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(bestmem3(:,1))/tekrar,sum(bestmem3(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(Gi3));
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(Gs3));
fprintf(1,'\nnumber of function ev= %f\n',nfeval3);
%-----------------------------------------------------------------------------------------------------------------
% ASOA_3 Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
fprintf(1,'\nASOA_3 Sonuçlarý\n');
fprintf(1,'\nOrtalama Sonuc= %e\n',mean(bestval1));
fprintf(1,'\nStandart Sapma= %e\n',std(bestval1));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(bestmem1(:,1))/tekrar,sum(bestmem1(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(Gi1));
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(Gs1));
fprintf(1,'\nnumber of function ev= %f\n',nfeval1);
%-----------------------------------------------------------------------------------------------------------------
% OSOA Sonuçlarý
%-----------------------------------------------------------------------------------------------------------------
fprintf(1,'\nOSOA Sonuçlarý\n');
fprintf(1,'\nOrtalama Sonuc= %e\n',mean(bestval4));
fprintf(1,'\nStandart Sapma= %e\n',std(bestval4));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(bestmem4(:,1))/tekrar,sum(bestmem4(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(Gi4));
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(Gs4));
fprintf(1,'\nnumber of function ev= %f\n',nfeval4);
%-----------------------------------------------------------------------------------------------------------------
