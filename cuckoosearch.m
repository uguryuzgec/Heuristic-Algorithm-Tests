function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure]=cuckoosearch(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,pa)

%----------------Cuckoo Optimization Algorithm-----------------%
%------------------------ Ali CAKMAK --------------------------%
%----------------------- 01.05.2015 ---------------------------%

%-----Giriþ Deðerleri Kontrolü---------------------------------------------
err=[ ];
if nargin<1, error('firefly 1st argument must be function name'); else 
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
if ~isempty(err)
  fprintf(stdout,'error in parameter %d\n', err);
  usage('cuckoo (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
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
GlobalMins=zeros(tekrar,itermax-1); % Herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % Herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % Herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % Herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % Herbir tekrar icin bulunan sure
nfeval = 0;
% BASLAMA %
for r=1:tekrar,
fprintf(1,'%d. tekrar\n',r);
nest = zeros(NP,D); %initialize pop to gain speed


ibest  = 1;  
% Rastgele Baþlangýç Çözümleri
for i=1:NP
   nest(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
val = zeros(1,NP);          % create and reset the "cost array"

val(1)  = feval(fname,nest(ibest,:));
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;

for i=2:NP                         % Kalan Üyeleri Kontrol
  val(i) = feval(fname,nest(i,:));
  nfeval  = nfeval + 1;
  if (val(i) < bestval)           % Üyenin Daha Ýyi Olup-Olmadýðý
     bestval = val(i);
	 ibest = i;
  end   
end
bestnest = nest(ibest,:);         % best member of current iteration

[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat... 
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))     %%%%% start iterations

    % Yeni Çözümler Üretmek (En Ýyi Veriyi Güncel Tutmak)
    % Levy Uçuþu
	% Levy Üst Ve Katsayýsý
	beta=3/2;
	sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
		
for j=1:NP,
     
    %% Mantegna Algoritmasý Ýle Levy Uçuþlarý
    u=rand(1,D)*sigma;
    v=rand(1,D);
    step=u./abs(v).^(1/beta);
  
    % Bir Sonraki Denklemde Fark Alýnmasý En Ýyi Deðeri Almak Ýçindir.
    % Eðer Çözüm Yani Sonuç En Ýdeal Ýse  Deðer Deðiþmeden Kalýr.     
    stepsize=0.01*step.*(nest(j,:)-bestnest);
    
    % Þimdi Gerçek Rastgele Yürüyüþ ve Uçuþ Hareketleri
    new_nest(j,:)=nest(j,:)+stepsize.*rand(1,D);
   
end
	% Keþif ve Rastgelelik
    % Keþfedilen yada Durum Vektörü
	K=rand(NP,D)>pa;
 
%% Ýlk Bulunan / Seçici Rastgele Deðerler
	stepsize=rand()*(new_nest(randperm(NP),:)-new_nest(randperm(NP),:));
	new_nest=new_nest+stepsize.*K;
%%%if new_nest() eger sinirlari gecerse ne olacak??? 
  for i=1:NP
      for j=1:D,
         if(new_nest(i,j)<XVmin(j))
            new_nest(i,j)=XVmin(j);
        elseif(new_nest(i,j)>XVmax(j))
            new_nest(i,j)=XVmax(j);
        end
      end
  end
	    
	%%%%%%%%%%%%%%[bestval,bestnest,nest,val,nfeval]=get_best_nest(nest,new_nest,val,fname);
%-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,new_nest(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       nest(i,:) = new_nest(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestnest = new_nest(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'Iteration: %d,  Best: %f,  NP: %d\n',iter,bestval,NP);
       for n=1:D
         fprintf(1,'best(%d) = %f\n',n,bestnest(n));
       end
    end
  end
  
GlobalMins(r,iter)=bestval;  
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);  
iter=iter+1;
end %% Ýterasyonun Bitmesi
sonsaat=toc(saat);
Globalbest(r,:)=bestnest;
GlobalBestval(r,1)=bestval;
GlobalIter(r,1)=iter-1;
GlobalSure(r,1)=sonsaat;

end %end for (rr=1:tekrar) ...

fprintf(1,'\nOrtalama Sonuc= %e\n',mean(GlobalBestval));
fprintf(1,'\nStandart Sapma= %e\n',std(GlobalBestval));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(Globalbest(:,1))/tekrar,sum(Globalbest(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(GlobalIter));
fprintf(1,'\nOrtalama Sure=%fsn\n',mean(GlobalSure));
nfeval1=nfeval/tekrar;

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
 end