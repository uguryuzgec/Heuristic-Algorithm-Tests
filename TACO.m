function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = TACO(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
% 06/02/2016 Touring Ant Colony Optimization Algorithm
%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('devec3 1st argument must be function name'); else 
  if exist(fname)<1; err(1,length(err)+1)=1; end; end;
if nargin<2, VTR = 1.e-6; else 
  if length(VTR)~=1; err(1,length(err)+1)=2; end; end;
if nargin<3, D = 2; else
  if length(D)~=1; err(1,length(err)+1)=3; end; end; 
if nargin<4, XVmin = [-2 -2];else
  if length(XVmin)~=D; err(1,length(err)+1)=4; end; end; 
if nargin<5, XVmax = [2 2]; else
  if length(XVmax)~=D; err(1,length(err)+1)=5; end; end; 
if nargin<6, NP = 10*D; else
  if length(NP)~=1; err(1,length(err)+1)=6; end; end; 
if nargin<7, itermax = 200; else
  if length(itermax)~=1; err(1,length(err)+1)=7; end; end; 
if nargin<8, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=8; end; end; 
if nargin<9, tekrar = 1; else
  if length(tekrar)~=1; err(1,length(err)+1)=9; end; end;   
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('devec3 (string,scalar,scalar,vector,vector,integer,integer,integer,integer)');    	
end

if (NP < 5)
   NP=5;
   fprintf(1,' NP increased to minimal value 5\n');
end
if (itermax <= 0)
   itermax = 200;
   fprintf(1,'itermax should be > 0; set to default value 200\n');
end
refresh = floor(refresh);

%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure
pop = zeros(NP,D); 
nfeval    = 0;                    % number of function evaluations
% BASLAMA %
for r=1:tekrar
fprintf(1,'\n\n%d. tekrar',r);
%------------------------------TACO Baslama-------------------------------
xMax=XVmax;
xMin=XVmin;
KisaYol=0;
os=1;
KarincaSayisi=NP;                                            % Toplam Karýnca Sayýsý (Populasyon)
BitSayisi=18;                                                % Çözüm için Kullanýlan BitSayýsýnýn Çözünürlüðü
Buharlasma=0.1;
Elit=1;
Yollar=(logical(round(rand(KarincaSayisi,BitSayisi*D))));     % Rasgele Üretilen Tüm Karýncalarýn Gittiði Ýlk Yollar
iYollar=zeros(KarincaSayisi,BitSayisi*D);                     % Rasgele Üretilen Tüm Karýncalarýn Gittiði Ýlk Yollar

FeromenMiktari=zeros(BitSayisi*D,2);                          % k. Karýncanýn Giddiði Yolun Toplam Feromen Miktarý
SecilenYol=zeros(KarincaSayisi,BitSayisi*D);
XX=zeros(KarincaSayisi,D);                                  % Herbir Yolun Sayýsal Deðeri
GidilenYol=zeros(KarincaSayisi,BitSayisi*D);                % k. Karýncanýn Giddiði Yolun Bit Karþýlýðý
%-------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                   Baþlangýç Deðerlerinin Oluþturulmasý
%__________________________________________________________________________
iYollar=Yollar;
 for j=1:D
	for ks=1:KarincaSayisi
        for i=0:BitSayisi-1
        XX(ks,j)=XX(ks,j)+ Yollar(ks,BitSayisi*j-i)*2^i;      % Herbir Yolun Sayýsal Deðeri
        end
    end
 end
 
for j=1:D    
  pop(:,j)=XX(:,j)*((xMax(1,j)-xMin(1,j))/((2^BitSayisi)-1))+xMin(1,j);                  % Herbir Yolun Ölçekli Sayýsal Deðeri
end

 for j=1:D
    for ks=1:KarincaSayisi
        for bs=1:BitSayisi
            BitDeger=Yollar(ks,(BitSayisi*(j-1))+bs);
               if (BitDeger==0) BitDeger=2;
               end
            GidilenYol(ks,(BitSayisi*(j-1))+bs)=BitDeger;
        end

        for bs=1:BitSayisi
            FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(ks,(BitSayisi*(j-1))+bs))=FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(ks,(BitSayisi*(j-1))+bs))+2;
        end

    end
 end
 
val       = zeros(1,NP);          % create and reset the "cost array"
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration


%------Evaluate the best member after initialization----------------------

ibest   = 1;                      % start with first population member
val(1)  = feval(fname,pop(ibest,:));
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members
  val(i) = feval(fname,pop(i,:));
  nfeval  = nfeval + 1;
  if (val(i) < bestval)           % if member is better
     ibest   = i;                 % save its location
     bestval = val(i);
  end   
end
bestmemit = pop(ibest,:);         % best member of current iteration
bestmem = bestmemit;              % best member ever

%------TACO-Minimization---------------------------------------------

[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayýcýyý baþlat...
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
 XX(:,:)=0;    
 for j=1:D
	for ks=1:KarincaSayisi
        for bs=1:BitSayisi
         %___________RULET______________________________________________
        
            a=FeromenMiktari((BitSayisi*(j-1))+bs,1)/(FeromenMiktari((BitSayisi*(j-1))+bs,1)+FeromenMiktari((BitSayisi*(j-1))+bs,2));
            b=FeromenMiktari((BitSayisi*(j-1))+bs,2)/(FeromenMiktari((BitSayisi*(j-1))+bs,1)+FeromenMiktari((BitSayisi*(j-1))+bs,2));
            
            ss=rand(1);
            if (ss>0.5) SecilenBit=1; else SecilenBit=0; end

            if (a>b) 
                if (a>rand(1)) SecilenBit=1; end
            end
            if (b>a) 
                if (b>rand(1)) SecilenBit=0; end
            end

            Yollar(ks,(BitSayisi*(j-1))+bs)=SecilenBit;
        end
    end
 end
    
    %--------------------------------------------------------------------------
    %           Seçilen Yola Göre Amaç Fonksiyonunu Tekrar Hesapla
    %--------------------------------------------------------------------------

 for j=1:D
	for ks=1:KarincaSayisi
        for i=0:BitSayisi-1
        XX(ks,j)=XX(ks,j)+ Yollar(ks,BitSayisi*j-i)*2^i; % Herbir Yolun Sayýsal Deðeri
        end
    end
 end
	
for j=1:D    
  pop(:,j)=XX(:,j)*((xMax(1,j)-xMin(1,j))/((2^BitSayisi)-1))+xMin(1,j); % Herbir Yolun Ölçekli Sayýsal Deðeri
end	

for i=1:KarincaSayisi
    % Her bir karýncanýn uygunlugu bulunuyor
    val(i) = feval(fname,pop(i,:));
	nfeval  = nfeval + 1;
end % for i=1:NP...

 for j=1:D
    for ks=1:KarincaSayisi
        for bs=1:BitSayisi
            BitDeger=Yollar(ks,(BitSayisi*(j-1))+bs);
               if (BitDeger==0) BitDeger=2; 
			   end               
            GidilenYol(ks,(BitSayisi*(j-1))+bs)=BitDeger;
        end

        for bs=1:BitSayisi
            FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(ks,(BitSayisi*(j-1))+bs))=Buharlasma*FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(1,(BitSayisi*(j-1))+bs));
            FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(ks,(BitSayisi*(j-1))+bs))= Buharlasma*FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(ks,(BitSayisi*(j-1))+bs))+1/val(ks);
        end
    end  
 end
	
    %___________________________ ELITIZM __________________________________
    
    [ky,s]=sort(val);
    if (ky(1)< bestval)
            bestval = ky(1);
			bestmem = pop(s(1),:);
		for j=1:D
            for bs=1:BitSayisi
            FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(os,(BitSayisi*(j-1))+bs))=Buharlasma*FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(os,(BitSayisi*(j-1))+bs));
            end 
            os=s(1);
            for bs=1:BitSayisi
            FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(s(1),(BitSayisi*(j-1))+bs))= FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(s(1),(BitSayisi*(j-1))+bs))+Elit/val(s(1));
            end
            for bs=1:BitSayisi
            FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(s(2),(BitSayisi*(j-1))+bs))= FeromenMiktari((BitSayisi*(j-1))+bs,GidilenYol(s(2),(BitSayisi*(j-1))+bs))+Elit/val(s(2));
            end
		end
     end
  
%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'\nIteration: %d,  Best: %f, NP: %d\n',iter,bestval,NP);
       for n=1:D
         fprintf(1,'best(%d) = %2.2f, ',n,bestmem(n));
       end
    end
  end

GlobalMins(r,iter)=bestval;
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = iter + 1;
 
end %---end while ((iter < itermax) ...

sonsaat=toc(saat);
Globalbest(r,:)=bestmem;
GlobalBestval(r,1)=bestval;
GlobalIter(r,1)=iter-1;
GlobalSure(r,1)=sonsaat;
end %end for (r=1:tekrar) ...

fprintf(1,'\nOrtalama Sonuc= %e\n',mean(GlobalBestval));
fprintf(1,'\nStandart Sapma= %e\n',std(GlobalBestval));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(Globalbest(:,1))/tekrar,sum(Globalbest(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(GlobalIter));
fprintf(1,'\nOrtalama Sure=%fsn\n',mean(GlobalSure));
nfeval1=nfeval/tekrar;
%save ('ugur.mat','GlobalMins');
 if(tekrar==1) 
%      hold on; 
%      plot(pop(:,1),pop(:,2),'b.','MarkerSize',20); 
	GBestmem = Globalbest;
	GBestval = GlobalBestval;
	Gmin=GlobalMins;
    Giter=GlobalIter;
	GSure=GlobalSure;
 else
	GBestmem = Globalbest;
	GBestval = GlobalBestval;
    Gmin=GlobalMins;
    Giter=GlobalIter;
	GSure=GlobalSure;
     
% 10 tekrar için ortalama ve her tekrara ait maliyet deðiþimleri...
% % % figure %%%HEP AÇIK
% % % subplot(211)
% % % plot(GlobalMins');%%%HEP AÇIK
% % % axis tight;
% % % xlabel('iterasyon sayýsý');
% % % ylabel('maliyet');
% % % title(fname);
% % % legend('TS1','TS2','TS3','TS4','TS5','TS6','TS7','TS8','TS9','TS10');
% % % legend show
% % % MEAN=mean(GlobalMins);
% % % subplot(212)
% % % plot(MEAN(1:floor(GlobalIter/tekrar)));
% % % axis tight;
% % % xlabel('iterasyon sayýsý');
% % % ylabel('maliyet');
% % % title(fname);
 end



 

