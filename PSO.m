function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = PSO(fname,VTR,D,XVmin,XVmax,NP,itermax,c1,c2,refresh,tekrar)
% 30/01/2016 Particle Swarm Optimization Algorithm
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
if nargin<8, c1 = 2; else
  if length(c1)~=1; err(1,length(err)+1)=8; end; end;
if nargin<9, c2 = 2; else
  if length(c2)~=1; err(1,length(err)+1)=9; end; end; 
if nargin<10, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=10; end; end; 
if nargin<11, tekrar = 1; else
  if length(tekrar)~=1; err(1,length(err)+1)=11; end; end;   
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('devec3 (string,scalar,scalar,vector,vector,integer,integer,scalar,scalar,integer,integer)');    	
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
% PSO parametreleri
fi=c1+c2; %learning factors
ksi=2/abs(2-fi-sqrt(fi^2-4*fi));%constriction factor
v=zeros(NP,D); %ilk parcacik hýzý	
lbest=pop;%ilk nesil local best=sürünün kendisi

%------PSO-Minimization---------------------------------------------

[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayýcýyý baþlat...
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
  popold = pop;                   % save the old population

 % Geçerli sürü icin yerel/global bestler belirle
    for i=1:NP
    % Parcacýk hýzlarý belirleniyor
		a1=rand; a2=rand; %rasgele sayýlar
        v(i,:)=ksi*(v(i)+c1*a1.*(lbest(i,:)-pop(i,:))+c2*a2.*(bestmemit-pop(i,:)));
        % Parcacik hýzlarý sýnýrlanýyor            
        v(i,:)=max(v(i,:),0.1*XVmin);%alt sýnýr/10'a sýnýrlama
        v(i,:)=min(v(i,:),0.1*XVmax);%ust sýnýr/10'a sýnýrlama
        % Parcaciklar guncelleniyor
        pop(i,:)=pop(i,:)+v(i,:);
		% Parcaciklar sýnýrlanýyor            
        pop(i,:)=max(pop(i,:),XVmin);
        pop(i,:)=min(pop(i,:),XVmax);
    end % for i=1:NP...
% Geçerli sürü icin yerel/global bestler belirle
    for i=1:NP
        % Her bir parcacigin uygunlugu bulunuyor
        val_pre(i) = val(i);%onceki uyg deðeri oncekine atanýyor 
        val(i) = feval(fname,pop(i,:));
		nfeval  = nfeval + 1;

        %Yerel (local) ve global en iyi parcacik(lar) belirleniyor 
        if val(i) < val_pre(i) 
		   lbest(i,:) = pop(i,:); 
		end
        if val(i) < bestval
		   bestmem = pop(i,:); 
		   bestval = val(i); 
		end
    end % for i=1:NP...

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



 

