function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = SA(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
% 03/02/2016 Simulated Annealing Algorithm
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
nfeval    = 0;                    % number of function evaluations
% BASLAMA %
for r=1:tekrar
fprintf(1,'\n\n%d. tekrar',r);
pop = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

for i=1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

popold    = zeros(size(pop));     % toggle population
val       = zeros(1,NP);          % create and reset the "cost array"
fitness   = zeros(1,NP); 
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration


%------Evaluate the best member after initialization----------------------

for i=1:NP                        % check the remaining members
  val(i) = feval(fname,pop(i,:));
  nfeval  = nfeval + 1;
 end
[bestval eki]=min(val);
bestmemit = pop(eki,:); % 1. neslin en iyisi.....
bestmem = bestmemit; 

%------SA-Minimization---------------------------------------------
TolFun=1e-4;
eps=0;
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayýcýyý baþlat...
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
	T=iter/itermax; 
    mu=10^(T*itermax);
  % Güncelleme iþlemi...
    for i=1:NP
		y=2*rand(1,D)-1;
		x=(((1+mu).^abs(y)-1)/mu).*sign(y);
	    dpop=x.*(XVmax - XVmin);
		newpop = pop(i,:)+dpop;
		newpop = max(newpop,XVmin);
        newpop = min(newpop,XVmax);
		temp=feval(fname,newpop);
		nfeval  = nfeval + 1;
		df=temp-val(i);
		if (df < 0 || rand < exp(-T*df/(abs(val(i))+eps)/TolFun))==1
            pop(i,:) = newpop;
			val(i)=temp;
        end
	end % for i=1:NP...

	
	[bestval indis]=min(val);
    bestmem = pop(indis,:);
    
  bestmemit = bestmem;   

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



 

