function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure]= moa(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,alfa,ro)
 
%-----------Magnetic-Inspired Optimization Algorithm-----------%
%----------------------- Merve BULDUR -------------------------%
%----------------------- 01.05.2015 ---------------------------%

%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('moa 1st argument must be function name'); else 
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
  if length(refresh)~=1; err(1,length(err)+1)=9; end; end; 
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('moa (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
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
%---------------------------------------------------------------------
%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure
% BASLAMA %
for r=1:tekrar
fprintf(1,'%d. tekrar\n',r);
pop = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------


for i=1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

popold    = zeros(size(pop));     % toggle population
newpop    = zeros(size(pop)); 
val       = zeros(1,NP);          % create and reset the "cost array"
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations

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
bestvalit = bestval;              % best value of current iteration

bestmem = bestmemit;              % best member ever

%------MOA-Minimization---------------------------------------------
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat... 
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))     %%%%% start iterations
  popold = pop;
  B=val;                         % fitness manyetik alan icine kaydet...
  Bmin=min(B);
  Bmax=max(B);
  for i=1:NP,
    B(i)=(B(i)-Bmin)/(Bmax-Bmin);
  end
  M=alfa+ro*B;   % kütle hesabý...
  F = zeros(NP,D);
  for i=1:NP,
    for j=1:D,
      if(i==1)
         u=NP;v=2;
      end
      if(i~=NP)
          v=i+1;
      end
      if(i~=1)
          u=i-1;
      end
      if(i==NP)
         v=1;
      end
      %.....................%
      if(j==1)
         y=D;z=2;
      end
      if(j~=D)
          z=j+1;
      end
      if(j~=1)
          y=j-1;
      end
      if(j==D)
         z=1;
      end
  komsular = [popold(u,j);popold(v,j);popold(i,y);popold(i,z)];
  uzaklik = sum(abs(komsular-popold(i,j))/(XVmax(j) - XVmin(j)));
  uzaklik = sum(uzaklik)/D;
  F(i,j) = F(i,j)+sum(komsular-popold(i,j))*B(i)/uzaklik;
  velocity(i,j) = (F(i,j)/M(i))+ XVmin(j) + rand*(XVmax(j) - XVmin(j));
  newpop(i,j) = popold(i,j)+velocity(i,j);
      end %for...
  end %for...
  
%-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,newpop(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       pop(i,:) = newpop(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = newpop(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

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
end   %---end while ((iter < itermax) ...
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



 

