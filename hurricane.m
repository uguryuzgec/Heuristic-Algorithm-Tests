function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = hurricane(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Rmax,w,R0)
																
%-----------Hurricane Based Optimization Algorithm-----------%
%------------------ Gürkan Mustafa ÇAKIR --------------------%
%---------------------- 31.03.2015 --------------------------%

%% catch error
err=[];
if nargin<1, error('firefly 1st argument must be function name'); else 
  if exist(fname)<1; err(1,length(err)+1)=1; end; end;
if nargin<2, VTR = 1.e-6; else 
  if length(VTR)~=1; err(1,length(err)+1)=2; end; end;
if nargin<3, D = 2; else
  if length(D)~=1; err(1,length(err)+1)=3; end; end; 
if nargin<4, itermax = 200; else
  if length(itermax)~=1; err(1,length(err)+1)=4; end; end; 
if nargin<5, NP = 10*D; else
  if length(NP)~=1; err(1,length(err)+1)=5; end; end;
if nargin<6, XVmin = [-2 -2];else
  if length(XVmin)~=D; err(1,length(err)+1)=6; end; end; 
if nargin<7, XVmax = [2 2]; else
  if length(XVmax)~=D; err(1,length(err)+1)=7; end; end;
if nargin<8, tekrar = 1; else
  if length(tekrar)~=1; err(1,length(err)+1)=8; end; end;
if nargin<9, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=9; end; end;

if ~isempty(err)  %if length(err)>0
  fprintf(stdout,'error in parameter %d\D', err);
  usage('hurricane(string,scalar,scalar,vector,vector,scalar,scalar,scalar,scalar)');    	
end

if (NP < 5)
   NP=5;
   fprintf(1,' NP increased to minimal value 5\D');
end

if (itermax <= 0)
   itermax = 200;
   fprintf(1,'itermax should be > 0; set to default value 200\D');
end
refresh = floor(refresh);

%% -----Initialize population and some arrays-------------------------------
GlobalMins = zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest = zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval = zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter = zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure = zeros(tekrar,1); % herbir tekrar icin bulunan sure
nfeval = 0;
% start iteration of tekrar
for rr=1:tekrar
fprintf(1,'%d. tekrar\n',rr);

 
fi = 2*pi*rand(NP,1);
%e  = rand(NP,2); % noktalarin bulundugu koordinatlar
f  = zeros(NP,1);

%% ----------Start normalizasyon------------%

for i=1:NP % normalizasyon yapiliyor
   x(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
 % % %  e(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

%----------End normalizasyon------------%

%% ----------Evaluate the best member after initialization------------%
ibest   = 1;                      % start with first population member
val(1)  = feval(fname,x(ibest,:));
nfeval  = nfeval + 1;
bestval = val(1);                 % best objective function value so far
for i=2:NP                        % check the remaining members
  val(i) = feval(fname,x(i,:));
  nfeval  = nfeval + 1;
  if (val(i) < bestval)           % if member is better
     ibest   = i;                 % save its location
     bestval = val(i);
  end   
end
eye = x(ibest,:);         % best member of current iteration

%------------------------Best Member End------------------------------%
%% stopping condition
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat = tic; % süre basla
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))%for j = 1: itermax
    x_old=x;
	for i=1:NP,
        r(i) = R0*exp(rand()*f(i));
		k = mod(i,D)+1;
        x_old(i,1) = r(i)*cos(fi(i)+f(i))+x_old(i,1);
        x_old(i,2) = r(i)*sin(fi(i)+f(i))+x_old(i,2);
	if r(i) < Rmax 
      f(i) = f(i)+w;
    else
      f(i) = f(i)+ w*((Rmax/r(i))^rand());
    end
	end
%%%if pop() eger sinirlari gecerse ne olacak??? 
  for i=1:NP,
      for j=1:D,
         if(x_old(i,j)<XVmin(j))
            x_old(i,j)=XVmin(j);
			fi(i) = 2*pi*rand(1,1); %fi initial
			f(i)  = 0;
        elseif(x_old(i,j)>XVmax(j))
            x_old(i,j)=XVmax(j);
			fi(i) = 2*pi*rand(1,1); %fi initial
			f(i)  = 0;
        end
      end
  end
  
 for i=1:NP
    tempval = feval(fname,x_old(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       x(i,:) = x_old(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          eye = x_old(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP
  

%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'Iteration: %d,  Best: %f,  NP: %d\n',iter,bestval,NP);
       for n=1:D
         fprintf(1,'best(%d) = %f\n',n,eye(n));
       end
    end
  end

%%     %---------------Output section End--------------------------------------
      
      GlobalMins(rr,iter)= bestval;
      [delta_val]=sort(val);
      durdurma_katsayisi=delta_val(NP)-delta_val(1);
      iter = iter + 1;%iterasyon degiskeni

    end %---end while ((iter < itermax) ... 
    
sonsaat=toc(saat);
Globalbest(rr,:)=eye(1,:);
GlobalBestval(rr,1)= GlobalMins(rr,iter-1);
GlobalIter(rr,1)=iter-1;
GlobalSure(rr,1)=sonsaat;
end %end for (rr=1:tekrar) ...

fprintf(1,'\nOrtalama Sonuc = %e\n',mean(GlobalBestval));
fprintf(1,'\nStandart Sapma = %e\n',std(GlobalBestval));
fprintf(1,'\nOrtalama Cozum = %f,%f\n',sum(Globalbest(:,1))/tekrar,sum(Globalbest(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon = %f\n',mean(GlobalIter));
fprintf(1,'\nOrtalama Sure = %f sn\n',mean(GlobalSure));
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