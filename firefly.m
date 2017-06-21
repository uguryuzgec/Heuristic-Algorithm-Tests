function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = firefly(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,alpha,gamma,delta)

%----------------Firefly Optimization Algorithm----------------%
%---------------------- Duygu ADIYAMAN-------------------------%
%----------------------- 01.05.2015 ---------------------------%

%-----Check input variables---------------------------------------------
err=[];
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
if nargin<10, alpha = 0.5; else
  if length(alpha)~=1; err(1,length(err)+1)=10; end; end; 
if ~isempty(err)  %if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('firefly(string,scalar,scalar,vector,vector,scalar,scalar,scalar,scalar)');    	
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
% % % rand('state',0); %resetleme...
%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure
% BASLAMA %
for rr=1:tekrar,
fprintf(1,'%d. tekrar\n',rr);
pop = zeros(NP,D); %initialize pop to gain speed

% generating the initial locations of n fireflies
for i=1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

val = zeros(1,NP);          % create and reset the "cost array"

nfeval  = 0;                    % number of function evaluations
ibest   = 1;                      % start with first population member
val(1)  = feval(fname,pop(ibest,:));
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members
  val(i) = feval(fname,pop(i,:));
  nfeval  = nfeval + 1;
  if (val(i) < bestval)         % if member is better
     ibest = i;                 % save its location
     bestval = val(i);
  end   
 end
bestmem = pop(ibest,:);         % best member of current iteration

% % % % Iterations or pseudo time marching
% % % %-------------------------------------------------------------------------------------------------------
% % %   x1 = linspace(XVmin(1), XVmax(1), 101);
% % %   x2 = linspace(XVmin(2), XVmax(2), 101);
% % %   x3 = zeros(length(x1), length(x2));
% % %         
% % %         % simply loop through the function (most functions expect 
% % %         % [N x 2] vectors as input, so meshgrid does not work)
% % %         for i = 1:length(x1)
% % %             for j = 1:length(x2)
% % %                 x3(i, j) = feval(fname,[x1(i), x2(j)]);
% % %             end
% % %         end
% % % 		figure(2);contour(x1', x2', x3); hold on; 
% % %         plot(pop(:,1),pop(:,2),'b.','MarkerSize',20);
% % %         drawnow
% % %         hold off
%-------------------------------------------------------------------------------------------------------
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat... 
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))     %%%%% start iterations
% % %  if (refresh > 0)
% % %     if (rem(iter,refresh) == 0)
% % %         figure(2);contour(x1', x2', x3); hold on; 
% % %         plot(pop(:,1),pop(:,2),'b.','MarkerSize',20);
% % %         drawnow
% % %         hold off
% % %     end
% % %  end
% Ranking the fireflies by their light intensity
[Lightn,Index]=sort(val); %sort(val,2,'descend');
pop=pop(Index,:);
Lighto=Lightn;
newpop=pop; %pop kopyalansin...
% Move all fireflies to the better locations
for i=1:NP,
% The attractiveness parameter beta=exp(-gamma*r)
    for j=1:NP,
		r_toplam=0;
		for k=1:D,
			r_toplam=r_toplam+(pop(i,k)-newpop(j,k))^2;
		end	
	r=sqrt(r_toplam);
	if Lightn(i)>Lighto(j), % Brighter and more attractive
		beta0=1;     beta=beta0*exp(-gamma*r.^2);
		pop(i,:) = pop(i,:).*(1-beta)+newpop(j,:).*beta+alpha.*(XVmin + rand(1,D).*(XVmax - XVmin));
		%%%pop(i,:) = pop(i,:).*(1-beta)+newpop(j,:).*beta+alpha.*(rand-0.5); ORJINAL...
		end
    end % end for j
end % end for i
%%%if pop() eger sinirlari gecerse ne olacak??? 
  for i=1:NP
      for j=1:D,
         if(pop(i,j)<XVmin(j))
            pop(i,j)=XVmin(j);
        elseif(pop(i,j)>XVmax(j))
            pop(i,j)=XVmax(j);
        end
      end
  end

for i=1:NP
    val(i) = feval(fname,pop(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
%        %----we update bestval only in case of success to save time-----------
       if (val(i) < bestval)     % if competitor better than the best one ever
          bestval = val(i);      % new best value
          bestmem = pop(i,:);    % new best parameter vector ever
       end
end %---end for imember=1:NP	

% Reduce randomness as iterations proceed
alpha=alpha*delta;

%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'Iteration: %d,  Best: %f,  NP: %d\n',iter,bestval,NP);
       for n=1:D
         fprintf(1,'best(%d) = %f\n',n,bestmem(n));
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
alpha = 0.2; 
end %end for (rr=1:tekrar) ...
% % % figure(2);contour(x1', x2', x3);
% % % hold on
% % % plot(pop(:,1),pop(:,2),'r.','MarkerSize',20);
% % % drawnow
% % % hold off

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