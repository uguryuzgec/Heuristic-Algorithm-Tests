function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = forest(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,drag,esik,yas_snr)

%-----------------Forest Optimization Algorithm----------------%
%---------------------- Mustafa KOSTEK ------------------------%
%----------------------- 01.06.2015 ---------------------------%

err=[];
if nargin<1, error('orman 1st argument must be function name'); else 
  if exist(fname)<1; err(1,length(err)+1)=1; end; end;
if nargin<2, VTR = 1.e-6; else 
  if length(VTR)~=1; err(1,length(err)+1)=2; end; end;
if nargin<3, D = 2; else
  if length(D)~=1; err(1,length(err)+1)=3; end; end; 
if nargin<4, XVmin = [-2 -2];else
  if length(XVmin)~=D; err(1,length(err)+1)=4; end; end; 
if nargin<5, XVmax = [2 2]; else
  if length(XVmax)~=D; err(1,length(err)+1)=5; end; end; 
if nargin<7, NP = 10*D; else
  if length(NP)~=1; err(1,length(err)+1)=7; end; end; 
if nargin<8, itermax = 200; else
  if length(itermax)~=1; err(1,length(err)+1)=8; end; end; 
if nargin<12, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=12; end; end; 
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('forest (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
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
%-------------parametreler-----------------------
%%% agacbas=10;
% % % drag=0.1;%adým parametresi
% % % a=0.05;%eþik belirtir
% % % yas_snr=10;%bu yaþtan yüksek olanlar silinir
%% agacson=100;%%düþük belirtilirse programýn yarýda kalmasý muhtemel 
mutasyon=10;%mutasyon kaç döngüde bir olacak
%lmax=2
%lmin=-2
%% v=0.5;
%enbuyuku=10;
%-------------------------------------------------------
%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure
% BASLAMA %
for r=1:tekrar
fprintf(1,'%d. tekrar\n',r);
agac = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

for i=1:NP
   agac(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

agacold   = zeros(size(agac));    % toggle population
val       = zeros(1,NP);          % create and reset the "cost array"
yas       = zeros(1,NP);          % create and reset the "cost array"
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations
val_old   = zeros(1,NP);
%------Evaluate the best member after initialization----------------------

ibest   = 1;                      % start with first population member
val(1)  = feval(fname,agac(ibest,:));
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members
  val(i) = feval(fname,agac(i,:));
  nfeval  = nfeval + 1;
  if (val(i) < bestval)           % if member is better
     ibest   = i;                 % save its location
     bestval = val(i);
  end   
end
bestmemit = agac(ibest,:);        % best member of current iteration
bestvalit = bestval;              % best value of current iteration

bestmem = bestmemit;              % best member ever

[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat... 
   
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))     %%%%% start iterations
agacold = agac;         % save the old population
 if mod(iter,mutasyon)==0 %%mutasyonun tuttuðu deðerde bir olacak 
        sigma=0.8;
     else 
        sigma=0.1;
 end  

for i=1:NP,
  for j=1:D,
    if val_old(i)<val(i) && agacold(i,j)<=XVmax(j) %%bestvalý yerinde v vardý,bestvalý koymak daha etkili  
    %her bir agacýn dalý için geliþme
        agacold(i,j)=agacold(i,j).*(1+drag*rand());
	end		
	if val_old(i)<val(i) && agacold(i,j)>=XVmin(j) 
	%her bir agacýn dalý için solma
        agacold(i,j)=agacold(i,j).*(1-drag*rand());  
    end
  end
end   

ort_val=mean(val);

  for i=1:NP %%aðýrlýklý ortalamasý 0,1 arasýnda rastgele bir deðerle çarpýlarak eþikleniyor(üreme)
     % % % k=size(agacold,1);     
      if (val(i)/ort_val)*rand()> esik %%  &&  agac(i,1)<enbuyuku 
            agacold(i,:)=(1+sigma*(2*rand()-1)).*agacold(i,:);
            yas(i)=0;
      else
            yas(i)=yas(i)+1;
      end
 end 
 for i=1:size(agacold,1),%hudutlarý aþmasýn
  for j=1:D,
    if (agacold(i,j)<XVmin(j))
        agacold(i,j)=XVmin(j);
    elseif (agacold(i,j)>XVmax(j))
        agacold(i,j)=XVmax(j);
    end
  end
end

%   n=size(agacold,1);ömrünü tamamlayan aðaçlar siliniyor
   for i=NP:-1:1
    if yas(i) > yas_snr
        agacold(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
        yas(i)= 0;
    end
   end

 tempval =zeros(1,NP); 
 for i=1:NP,
    tempval(i) = feval(fname,agacold(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
end
[tempval,Index]=sort(tempval);
agac=agacold(Index,:);
val_old = val;
val=tempval(Index);
yas=yas(Index);
if (val(1) < bestval) 
bestval = val(1);      % new best value
bestmem = agacold(1,:);      % new best parameter vector ever
end

bestmemit = bestmem;       % freeze the best member of this iteration for the coming 
                             % iteration. This is needed for some of the strategies.

%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'Iteration: %d,  Best: %f,  NP: %d\n',iter,bestval,NP);
       for n=1:D
         fprintf(1,'best(%d) = %f\n',n,bestmem(n));
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
%  ============== end =====================================