function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = fruitfly_adapt_opposit(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
 
%---------------Fruit Fly Optimization Algorithm---------------%
%--------------------- Zuhal CALAYOGLU ------------------------%
%----------------------- 18.05.2015 ---------------------------%

%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('fruitfly 1st argument must be function name'); else 
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
  if length(NP)~=1; err(1,length(err)+1)=7; end; end; 
if nargin<7, itermax = 200; else
  if length(itermax)~=1; err(1,length(err)+1)=8; end; end; 
if nargin<8, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=12; end; end; 
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
 usage('fruitfly(string,scalar,scalar,vector,vector,scalar,scalar,scalar,scalar)');    	
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
% Oppositional based parameters...
Jr=0.1;

%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure
% BASLAMA %
for r=1:tekrar
fprintf(1,'%d. tekrar\n',r);

% % % XY_axis1 = XVmin + rand(1,D).*(XVmax - XVmin);
% % % XY_axis2 = XVmin + rand(1,D).*(XVmax - XVmin);
 X_axis = rands(1,D);
 Y_axis = rands(1,D);

X = zeros(NP,D); %initialize pop to gain speed
Y = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

for i=1:NP,
   X(i,:) = X_axis + 2*rand(1,D)-1;
   Y(i,:) = Y_axis + 2*rand(1,D)-1;
  % X(i,:) = XVmax.^-1 + rand(1,D).*(XVmin.^-1 - XVmax.^-1);
  %	Y(i,:) = XVmax.^-1 + rand(1,D).*(XVmin.^-1 - XVmax.^-1);
    for k=1:D,
		Distance(i,k)=sqrt(X(i,k)^2+Y(i,k)^2);
		S(i,k)=1/Distance(i,k);
		if(S(i,k)<XVmin(k))
            S(i,k)=XVmin(k);
        elseif(S(i,k)>XVmax(k))
            S(i,k)=XVmax(k);
        end
    end
end

val = zeros(1,NP);          % create and reset the "cost array"
nfeval  = 0;                    % number of function evaluations
ibest   = 1;                      % start with first population member
val(1)  = feval(fname,S(ibest,:));
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;
for i=1:NP                        % check the remaining members
  val(i) = feval(fname,S(i,:));
  nfeval  = nfeval + 1;
  if (val(i) < bestval)         % if member is better
     ibest = i;                 % save its location
     bestval = val(i);
  end   
end
bestmem = S(ibest,:);         % best member of current iteration

X_axis = X(ibest,:);
Y_axis = Y(ibest,:);

% Iterations or pseudo time marching
%-------------------------------------------------------------------------------------------------------
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
% % %         plot(S(:,1),S(:,2),'b.','MarkerSize',20);
% % %         drawnow
% % %         hold off
%-------------------------------------------------------------------------------------------------------
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat...
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
Sold = S; 
% % % if (refresh > 0)
% % %     if (rem(iter,refresh) == 0)
% % %         figure(2);contour(x1', x2', x3); hold on; 
% % %         plot(S(:,1),S(:,2),'b.','MarkerSize',20);
% % %         drawnow
% % %         hold off
% % %     end
% % %  end
    
%-----Select which vectors are allowed to enter the new population------------

for i=1:NP,
  Xold(i,:) = X_axis +2*rand(1,D)-1;
  Yold(i,:) = Y_axis +2*rand(1,D)-1;
 %  Xold(i,:) = X_axis + rand(1,D).*(XVmin.^-1 - XVmax.^-1);
 %  Yold(i,:) = Y_axis + rand(1,D).*(XVmin.^-1 - XVmax.^-1);
    for k=1:D,
		Distance(i,k) = sqrt(Xold(i,k)^2+Yold(i,k)^2);
        %Sold(i,k)=1/Distance(i,k);
        Sold(i,k)=1/Distance(i,k);
		if(Sold(i,k)<XVmin(k))
            Sold(i,k)=XVmin(k);
        elseif(Sold(i,k)>XVmax(k))
            Sold(i,k)=XVmax(k);
        end
    end
end

%-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,Sold(i,:));   % check cost of smells
    nfeval  = nfeval + 1;
    if (tempval <= val(i))%%true  % if competitor is better than value in "cost array"
       X(i,:) = Xold(i,:);  % replace old vector with new one (for new iteration)
	   Y(i,:) = Yold(i,:);
       val(i) = tempval;  % save value in "cost array"
       S(i,:) = Sold(i,:);
       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = Sold(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

[abc,index]=min(val);
%%% bestmem = S(index,:);

X_axis =  X(index,:);
Y_axis =  Y(index,:);

bestmemit = bestmem;       % freeze the best member of this iteration for the coming 
                             % iteration. This is needed for some of the strategies.
							  % opposition based generation jumping......
  %---------------------------------------------%
  if (rand(1) < Jr),
     %%% fprintf(1,'\n opposition based generation jumping and Jr=%d\n',Jr);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:NP,
      %%%O_pop(i,D) = max(pop(:,D)) + min(pop(:,D)) - pop(i,D);
	  Xold(i,:) = 10-X(i,:);
	  Yold(i,:) = 10-Y(i,:);
    for k=1:D,
		Distance(i,k) = sqrt(Xold(i,k)^2+Yold(i,k)^2);
        %Sold(i,k)=1/Distance(i,k);
        Sold(i,k)=(-1)^randi([1 2])/Distance(i,k);
		if(Sold(i,k)<XVmin(k))
            Sold(i,k)=XVmin(k);
        elseif(Sold(i,k)>XVmax(k))
            Sold(i,k)=XVmax(k);
        end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,Sold(i,:));   % check cost of smells
    nfeval  = nfeval + 1;
    if (tempval <= val(i))%%true  % if competitor is better than value in "cost array"
       X(i,:) = Xold(i,:);  % replace old vector with new one (for new iteration)
	   Y(i,:) = Yold(i,:);
       val(i) = tempval;  % save value in "cost array"
       S(i,:) = Sold(i,:);
       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = Sold(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

[abc,index]=min(val);
%%% bestmem = S(index,:);

X_axis =  X(index,:);
Y_axis =  Y(index,:);
X_best(iter,:) = X_axis; 
Y_best(iter,:) = Y_axis;

bestmemit = bestmem;  
 end %(rand(1) < Jr),
%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'Iteration: %d,  Best: %f,   NP: %d\n',iter,bestval,NP);
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
% % % figure; plot(X_best,Y_best,'b.'); % tekrar = 1 için fruitfly uçuþlarý çizimi...
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