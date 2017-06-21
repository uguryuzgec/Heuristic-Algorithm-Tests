function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = brainstorm(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,prob,NC)

%--------------Brainstorm Optimization Algorithm---------------%
%--------------------- Ali Sahin BAYOGLU ----------------------%
%----------------------- 01.05.2015 ---------------------------%

warning('off', 'all');
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
if ~isempty(err)  %if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('brainstorm(string,scalar,scalar,vector,vector,scalar,scalar,scalar,scalar)');    	
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

prob_one_cluster = prob; % probability for select one cluster to form new individual; 
stepSize = ones(1,D); % effecting the step size of generating new individuals by adding random values
%%% NC=4; % number of clusters
n_iteration=20;

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
centers = zeros(NC,D);
% generating the initial locations
for i=1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
popu_sorted  = pop; % initialize the  population of individuals sorted according to clusters
% initialize cluster probability to be zeros
prob = zeros(NC,1);
best = zeros(NC,1);  % index of best individual in each cluster
for i=1:NC
   centers(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
val = zeros(NP,1);          % create and reset the "cost array"
val_sorted =val;
indi_temp = zeros(1,D);  % store temperary individual

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
% % % %-------------------------------------------------------------------------------------------------------
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
 cluster = kmeans(popu_sorted, NC,'Distance','cityblock','Start',centers,'EmptyAction','singleton'); % k-mean cluster
 % clustering    
 fit_values = 100000000000000000000000000.0*ones(NC,1);  % assign a initial big fitness value  as best fitness for each cluster in minimization problems
 number_in_cluster = zeros(NC,1);  % initialize 0 individual in each cluster

 for idx = 1:NP
        number_in_cluster(cluster(idx,1),1)= number_in_cluster(cluster(idx,1),1) + 1;
               
        % find the best individual in each cluster
        if fit_values(cluster(idx,1),1) > val(idx,1)  % minimization
            fit_values(cluster(idx,1),1) = val(idx,1);
            best(cluster(idx,1),1) = idx;
        end
            
    end
 
% form population sorted according to clusters
    counter_cluster = zeros(NC,1);  % initialize cluster counter to be 0 
    acculate_num_cluster = zeros(NC,1);  % initialize accumulated number of individuals in previous clusters
    
    for idx =2:NC
        acculate_num_cluster(idx,1) = acculate_num_cluster((idx-1),1) + number_in_cluster((idx-1),1);
    end
 %start form sorted population
    for idx = 1:NP
        counter_cluster(cluster(idx,1),1) = counter_cluster(cluster(idx,1),1) + 1 ;
        temIdx = acculate_num_cluster(cluster(idx,1),1) +  counter_cluster(cluster(idx,1),1);
        popu_sorted(temIdx,:) = pop(idx,:);
        val_sorted(temIdx,1) = val(idx,1);
    end
 % record the best individual in each cluster
    for idx = 1:NC
        centers(idx,:) = pop(best(idx,1),:);        
    end
    centers_copy = centers;  % make a copy
    if (rand() < 0.2) %  select one cluster center to be replaced by a randomly generated center
        cenIdx = ceil(rand()*NC);
        centers(cenIdx,:) = XVmin + rand(1,D).*(XVmax - XVmin);
    end 
 % calculate cluster probabilities based on number of individuals in
    % each cluster
    for idx = 1:NC
        prob(idx,1) = number_in_cluster(idx,1)/NP;
        if idx > 1
            prob(idx,1) = prob(idx,1) + prob(idx-1,1);
        end
    end
      % generate n_p new individuals by adding Gaussian random values
                   
    for idx = 1:NP
        r_1 = rand();  % probability for select one cluster to form new individual
        if r_1 < prob_one_cluster % select one cluster
            r = rand();
            for idj = 1:NC
                if r < prob(idj,1)                      
                    if rand() < 0.4  % use the center
                       indi_temp(1,:) = centers(idj,:); 
                    else % use one randomly selected  cluster
                        indi_1 = acculate_num_cluster(idj,1) + ceil(rand() * number_in_cluster(idj,1));
                        indi_temp(1,:) = popu_sorted(indi_1,:);  
                    end
                    break
                end
            end
        else % select two clusters
            % pick two clusters 
            cluster_1 = ceil(rand() * NC);
            indi_1 = acculate_num_cluster(cluster_1,1) + ceil(rand() * number_in_cluster(cluster_1,1));
            
            cluster_2 = ceil(rand() * NC);
            indi_2 = acculate_num_cluster(cluster_2,1) + ceil(rand() * number_in_cluster(cluster_2,1));
            
            tem = rand();
            if rand() < 0.5 %use center
                indi_temp(1,:) = tem * centers(cluster_1,:) + (1-tem) * centers(cluster_2,:); 
            else   % use randomly selected individuals from each cluster            
                indi_temp(1,:) = tem * popu_sorted(indi_1,:) + (1-tem) * popu_sorted(indi_2,:); 
            end
        end         
        
        stepSize = logsig(((0.5*itermax - n_iteration)/20)) * rand(1,D);
        indi_temp(1,:) = indi_temp(1,:)+stepSize .* normrnd(0,1,1,D);
        % if better than the previous one, replace it
        fv = feval(fname,indi_temp);
        if fv < val(idx,1)  % better than the previous one, replace
            val(idx,1) = fv;
            pop(idx,:) = indi_temp(1,:);
        end 
        
    end %%%forrrr
    
      % keep the best for each cluster
    for idx = 1:NC
        pop(best(idx,1),:) = centers_copy(idx,:);  
        val(best(idx,1),1) = fit_values(idx,1);
    end
       
       
    % record the best fitness in each iteration
    [bestval indis] = min(fit_values);
    bestmem = centers_copy(indis,:);
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
fprintf(1,'\nOrtalama Sure%fsn\n',mean(GlobalSure));
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