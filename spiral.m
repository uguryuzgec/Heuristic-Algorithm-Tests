function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,r)
 
%----------------Spiral Optimization Algorithm-----------------%
%----------------------- Refika ANBAR -------------------------%
%----------------------- 12.05.2015 ---------------------------%

%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('spiral 1st argument must be function name'); else 
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
  usage('spiral(string,scalar,scalar,vector,vector,scalar,scalar,scalar,scalar,scalar,scalar)');     
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
nfeval    = 0; %sayac                   % number of function evaluations

%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure
% BASLAMA %
for rr=1:tekrar,
fprintf(1,'\n%d. tekrar\n',rr);
popnokta = zeros(NP,D); %initialize pop to gain speed
newpopnokta = zeros(NP,D);
% ------------------------------------------------

% aciklama yazzzzzzzzzz unutma :-)))
%M  matrisi D boyutlu
M=eye(D);
M_guncel=M;
M=eye(D);
M_guncel=M;
for j=2:D,
	i=1;
	while(i<j)
		M(i,i)=cos(Theta);
		M(i,j)=-sin(Theta);
		M(j,i)=sin(Theta);
		M(j,j)=cos(Theta);
		M_guncel=M_guncel*M;
		i=i+1;
        M=eye(D);
	end
end
 I=eye(D);
 M=M_guncel;
% ------------------------------------------------
% generating the initial locations of n spiral(n boyutlu  spiralin ba�lang��  konumunu olu�turur.)
for i=1:NP
   popnokta (i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
 
val       = zeros(1,NP);          % create and reset the "cost array"( her bir bireyin maliyet dizisi olu�tur.)
bestmem   = zeros(1,D);           % best population member ever(�imdiye kadar en  iyi populasyon eleman�)
bestmemit = zeros(1,D);           % best population member in iteration(iterasyondaki  en iyi populasyon elemean�)

 
ibest   = 1;   %en iyi bir birey                   % start with first population member(ilk populasyon eleman� ile birlikte ba�la)
val(1)  = feval(fname,popnokta(ibest,:));%populasyondaki bireylerin ilkini al�yor.

bestval = val(1);                 % best objective function value so far(�imdiye kadar  hedeflenen en iyi  fonksiyon de�eri )
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members(kalan �ye kontrol)
  val(i) = feval(fname,popnokta(i,:));%ikinci bireyin maliyet de�eri 
  nfeval  = nfeval + 1;
  if (val(i) < bestval)           % if member is better(e�er  val(i)<bestval ise)E�er memeber daha iyi ise
     ibest   = i;                 % save its location(konumu kaydet) i indisini g�ncelle
     bestval = val(i);%best val g�ncelle
  end   
end
bestmemit = popnokta(ibest,:);         % best member of current iteration( g�ncel  iterasyonun en iyi eleman�)
bestvalit = bestval;              % best value of current iteration(g�ncel iterasyonun en iyi de�eri)
 
bestmem = bestmemit;
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat... 
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))     %%%%% start iterations
Xmerkez = bestmem; 
 
for i=1:NP,
    newpopnokta(i,:)=(r*M)*popnokta(i,:)'-(r*M-I)*Xmerkez'; %sonuc 
end % end for i
%%%if pop() eger sinirlari gecerse ne olacak??? for i=1:NP,
  for i=1:NP
      for j=1:D,
         if(newpopnokta(i,j)<XVmin(j))
            newpopnokta(i,j)=XVmin(j);
        elseif(newpopnokta(i,j)>XVmax(j))
            newpopnokta(i,j)=XVmax(j);
        end
      end
  end
for i=1:NP,
    tempval = feval(fname,newpopnokta(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    popnokta(i,:) = newpopnokta(i,:);  % replace old vector with new one (for new iteration)
    val(i) = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem =newpopnokta(i,:);      % new best parameter vector ever
       end
end %---end for imember=1:NP   
% plot(popnokta(:,1),popnokta(:,2),'r.','MarkerSize',20);
% drawnow
% pause(); 

bestmemit = bestmem; 
%----Output section----------------------------------------------------------
 
  if (refresh > 0)
    if (rem(iter,refresh) == 0)
        fprintf(1,'Iteration: %d,  Best: %f,  NP: %d\n',iter,bestval,NP);
       for n=1:D
         fprintf(1,'best(%d) = %f, ',n,bestmem(n));
       end
    end
   
  end
 
GlobalMins(rr,iter)=bestval;
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = iter + 1; 	   
end   %---end while ((iter < itermax) ...
sonsaat=toc(saat);
Globalbest(rr,:)=bestmem;
GlobalBestval(rr,1)=bestval;
GlobalIter(rr,1)=iter-1;
GlobalSure(rr,1)=sonsaat;
end %end for (rr=1:tekrar) ...
 
fprintf(1,'\nOrtalama Sonuc= %e\n',mean(GlobalBestval));
fprintf(1,'\nStandart Sapma= %e\n',std(GlobalBestval));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(Globalbest(:,1))/tekrar,sum(Globalbest(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(GlobalIter));
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(GlobalSure));
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