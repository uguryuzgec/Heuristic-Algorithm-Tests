function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = crab(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,receptivity,esikdegeri,NoCross,NoMut)
 
%--------------Crab Mating Optimization Algorithm--------------%
%-------------------- Nartan Ayberk ASKIN----------------------%
%----------------------- 18.05.2015 ---------------------------%

%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('yengec 1st argument must be function name'); else 
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
if nargin<9, tekrar = 10; else
  if length(tekrar)~=1; err(1,length(err)+1)=9; end; end;
if ~isempty(err)% if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('crab (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
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
GBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure

% BASLAMA %
NP=NP/2;  % half of NP for male & female populations...
for r=1:tekrar
fprintf(1,'%d. tekrar\n',r);
malepop = zeros(NP,D); %initialize pop to gain speed
femalepop = zeros(NP,D); %initialize pop to gain speed
newmalepop = zeros(NP,D); %initialize pop to gain speed
newfemalepop = zeros(NP,D); %initialize pop to gain speed

nrmating=0; %çiftlesme sayisi...
% % % receptivity=200; %alým derecesi... 1 den büyük olmak zorunda
% % % esikdegeri=0.75; 
% % % NoCross=0.75; %çaprazlama deðeri...
% % % NoMut=0.15; %mutasyon deðeri...
beta=3; %yakýnsama faktörü...

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

for i=1:NP
   malepop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
   femalepop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

popold    = zeros(size(malepop));     % toggle population
maleval   = zeros(1,NP);          % create and reset the "cost array"
femaleval = zeros(1,NP);          % create and reset the "cost array"
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations

%------Evaluate the best member after initialization----------------------

ibest   = 1;                      % start with first population member
% maleval(1)  = feval(fname,malepop(ibest,:));
% femaleval(1)  = feval(fname,femalepop(ibest,:));
%bestval = val(1);                 % best objective function value so far
%nfeval  = nfeval + 1;
for i=1:NP                        % check the remaining members
  maleval(i) = feval(fname,malepop(i,:));
  nfeval  = nfeval + 1;
  femaleval(i) = feval(fname,femalepop(i,:));
  nfeval  = nfeval + 1;
%    if (val(i) < bestval)           % if member is better
%       ibest   = i;                 % save its location
%       bestval = val(i);
%    end   
end
[maleval,index1]=sort(maleval);
[femaleval,index2]=sort(femaleval);
malepop = malepop(index1,:);
femalepop = femalepop(index2,:);
if maleval(1)<femaleval(1)
   bestmemit = malepop(1,:);         % best member of current iteration
%    bestvalit = maleval(1);              % best value of current iteration
   bestval=maleval(1);
else
   bestmemit = femalepop(1,:);         % best member of current iteration
%    bestvalit = femaleval(1);              % best value of current iteration
   bestval=femaleval(1);
end

bestmem = bestmemit;              % best member ever

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

birles=[maleval femaleval];
[delta_val]=sort(birles);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayýcýyý baþlat...
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
  eggs=[];  
    for i=1:NP
      for j=1:NP

        prob=exp((maleval(i)*nrmating)/receptivity);
    if prob>esikdegeri
        child1=NoCross*malepop(i,:) + (1-NoCross)*femalepop(j,:);
        child2=(1-NoCross)*malepop(i,:) + NoCross*femalepop(j,:);
        nrmating = nrmating + 1;
        neggs = round((beta*NP)/j); %neggs=üretilen yumurtalar
        nfeggs = round((neggs*maleval(i))/(100*nrmating)); % nfeggs=döllenen yumurtalar
        child1M = child1 + (XVmin + rand(1,D).*(XVmax - XVmin))*NoMut;
        child2M = child2 + (XVmin + rand(1,D).*(XVmax - XVmin))*NoMut;
        eggs(size(eggs,1)+1,1:2)=child1M;
        eggs(size(eggs,1)+1,1:2)=child2M;
     end %if
     
    end %for
     nrmating=0;
  end % for

%-----Select which vectors are allowed to enter the new population------------
    for i=1:size(eggs,1),
        tempval(i) = feval(fname,eggs(i,:));  
		nfeval  = nfeval + 1;
    end
    topla_crab_val=[maleval femaleval tempval];
    topla_crab=[malepop;femalepop;eggs];
    [topla_crab_val,ii]=sort(topla_crab_val);
    topla_crab=topla_crab(ii,:);
    malepop=topla_crab(1:NP,:);
    femalepop=topla_crab(NP+1:2*NP,:);
    maleval=topla_crab_val(1:NP);
    femaleval=topla_crab_val(NP+1:NP*2);
    topla_crab_val= [];
    topla_crab =[];
    ii=[];
    tempval=[];

if maleval(1)<femaleval(1)
   bestmem = malepop(1,:);         % best member of current iteration
%    bestvalit = maleval(1);              % best value of current iteration
   bestval=maleval(1);
else
   bestmem = femalepop(1,:);         % best member of current iteration
%    bestvalit = femaleval(1);              % best value of current iteration
   bestval=femaleval(1);
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
 birles=[maleval femaleval];
 [delta_val]=sort(birles);
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
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(GlobalSure));
nfeval1=nfeval/tekrar;
%save ('ugur.mat','GlobalMins');
 if(tekrar==1) 
%      hold on; 
%      plot(malepop(:,1),malepop(:,2),'b.','MarkerSize',20); 
%      plot(femalepop(:,1),femalepop(:,2),'r.','MarkerSize',20); 
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
% % % plot(MEAN(1:floor(GIter/tekrar)));
% % % axis tight;
% % % xlabel('iterasyon sayýsý');
% % % ylabel('maliyet');
% % % title(fname);
 end



 

