function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = flower(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,prob_switch,Levy_beta)
 
%-----------Flower Pollination Optimization Algorithm-----------%
%----------------------- Yusuf BORUCU --------------------------%
%----------------------- 17.05.2015 ----------------------------%

err=[];
if nargin<1, error('flower 1st argument must be function name'); else 
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
if ~isempty(err)  %if length(err)>0
   fprintf(stdout,'error in parameter %d\n', err);
   usage('flower(string,scalar,scalar,vector,vector,scalar,scalar,scalar,scalar,vector)');    	
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

GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure

for rr=1:tekrar,
fprintf(1,'%d. tekrar\n',rr);
pop = zeros(NP,D);

%Popülasyonu baþlat
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
bestmem = pop(ibest,:);         % Ýlk popülasyondaki en iyi çözümü bul

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
% % %         
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat... 

while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
    oldpop = pop;
% % %     if (refresh > 0)
% % %     if (rem(iter,refresh) == 0)
% % %         figure(2);contour(x1', x2', x3); hold on; 
% % %         plot(pop(:,1),pop(:,2),'b.','MarkerSize',20);
% % %         drawnow
% % %         hold off
% % %     end
% % %     end
for i=1:NP, 
 % Pollens are carried by insects and thus can move in
 % large scale, large distance.
 % This L should replace by Levy flights  
 % Formula: x_i^{t+1}=x_i^t+ L (x_i^t-gbest)
   if rand()>prob_switch,
     sigma=(gamma(1+Levy_beta)*sin(pi*Levy_beta/2)/(gamma((1+Levy_beta)/2)*Levy_beta*2^((Levy_beta-1)/2)))^(1/Levy_beta);
     u=randn(1,D)*sigma;
     v=randn(1,D);
     step=u./abs(v).^(1/Levy_beta);
     L=0.01*step;
     dS=L.*(pop(i,:)-bestmem);
     oldpop(i,:)=pop(i,:)+dS;
          
% If not, then local pollenation of neighbor flowers 
else
    epsilon=rand();
    % Find random flowers in the neighbourhood
    JK=randperm(NP);
    % As they are random, the first two entries also random
    % If the flower are the same or similar species, then
    % they can be pollenated, otherwise, no action.
    % Formula: x_i^{t+1}+epsilon*(x_j^t-x_k^t)
    oldpop(i,:)=oldpop(i,:)+epsilon*(pop(JK(1),:)-pop(JK(2),:));
   end
end % i=1:NP,...
   %%%if pop() eger sinirlari gecerse ne olacak??? 
  for i=1:NP
      for j=1:D,
         if(oldpop(i,j)<XVmin(j))
            oldpop(i,j)=XVmin(j);
        elseif(oldpop(i,j)>XVmax(j))
            oldpop(i,j)=XVmax(j);
        end
      end
  end
  
  % Yeni çözümleri deðerlendir
  %-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,oldpop(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       pop(i,:) = oldpop(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = oldpop(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

  bestmemit = bestmem;       % freeze the best member of this iteration for the coming 
                             % iteration. This is needed for some of the strategies.

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