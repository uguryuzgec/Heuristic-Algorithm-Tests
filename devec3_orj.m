function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = devec3_orj(fname,VTR,D,XVmin,XVmax,NP,itermax,F,CR,strategy,refresh,tekrar,solution)
% minimization of a user-supplied function with respect to x(1:D),
% using the differential evolution (DE) algorithm of Rainer Storn
% (http://www.icsi.berkeley.edu/~storn/code.html)
% 
% Special thanks go to Ken Price (kprice@solano.community.net) and
% Arnold Neumaier (http://solon.cma.univie.ac.at/~neum/) for their
% valuable contributions to improve the code.
% 
% Strategies with exponential crossover, further input variable
% tests, and arbitrary function name implemented by Jim Van Zandt 
% <jrv@vanzandt.mv.com>, 12/97.
%
% Output arguments:
% ----------------
% bestmem        parameter vector with best solution
% bestval        best objective function value
% nfeval         number of function evaluations
%
% Input arguments:  
% ---------------
%
% fname          string naming a function f(x,y) to minimize
% VTR            "Value To Reach". devec3 will stop its minimization
%                if either the maximum number of iterations "itermax"
%                is reached or the best parameter vector "bestmem" 
%                has found a value f(bestmem,y) <= VTR.
% D              number of parameters of the objective function 
% XVmin          vector of lower bounds XVmin(1) ... XVmin(D)
%                of initial population
%                *** note: these are not bound constraints!! ***
% XVmax          vector of upper bounds XVmax(1) ... XVmax(D)
%                of initial population
% y		        problem data vector (must remain fixed during the
%                minimization)
% NP             number of population members
% itermax        maximum number of iterations (generations)
% F              DE-stepsize F from interval [0, 2]
% CR             crossover probability constant from interval [0, 1]
% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin
%                Experiments suggest that /bin likes to have a slightly
%                larger CR than /exp.
% refresh        intermediate output will be produced after "refresh"
%                iterations. No intermediate output will be produced
%                if refresh is < 1
%
%       The first four arguments are essential (though they have
%       default values, too). In particular, the algorithm seems to
%       work well only if [XVmin,XVmax] covers the region where the
%       global minimum is expected. DE is also somewhat sensitive to
%       the choice of the stepsize F. A good initial guess is to
%       choose F from interval [0.5, 1], e.g. 0.8. CR, the crossover
%       probability constant from interval [0, 1] helps to maintain
%       the diversity of the population and is rather uncritical. The
%       number of population members NP is also not very critical. A
%       good initial guess is 10*D. Depending on the difficulty of the
%       problem NP can be lower than 10*D or must be higher than 10*D
%       to achieve convergence.
%       If the parameters are correlated, high values of CR work better.
%       The reverse is true for no correlation.
%
% default values in case of missing input arguments:
% 	VTR = 1.e-6;
% 	D = 2; 
% 	XVmin = [-2 -2]; 
% 	XVmax = [2 2]; 
%	y=[];
% 	NP = 10*D; 
% 	itermax = 200; 
% 	F = 0.8; 
% 	CR = 0.5; 
% 	strategy = 7;
% 	refresh = 10; 
%
% Cost function:  	function result = f(x,y);
%                      	has to be defined by the user and is minimized
%			w.r. to  x(1:D).
%
% Example to find the minimum of the Rosenbrock saddle:
% ----------------------------------------------------
% Define f.m as:
%                    function result = f(x,y);
%                    result = 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%                    end
% Then type:
%
% 	VTR = 1.e-6;
% 	D = 2; 
% 	XVmin = [-2 -2]; 
% 	XVmax = [2 2]; 
% 	[bestmem,bestval,nfeval] = devec3("f",VTR,D,XVmin,XVmax);
%
% The same example with a more complete argument list is handled in 
% run1.m
%
% About devec3.m
% --------------
% Differential Evolution for MATLAB
% Copyright (C) 1996, 1997 R. Storn
% International Computer Science Institute (ICSI)
% 1947 Center Street, Suite 600
% Berkeley, CA 94704
% E-mail: storn@icsi.berkeley.edu
% WWW:    http://http.icsi.berkeley.edu/~storn
%
% devec is a vectorized variant of DE which, however, has a
% propertiy which differs from the original version of DE:
% 1) The random selection of vectors is performed by shuffling the
%    population array. Hence a certain vector can't be chosen twice
%    in the same term of the perturbation expression.
%
% Due to the vectorized expressions devec3 executes fairly fast
% in MATLAB's interpreter environment.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 1, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. A copy of the GNU 
% General Public License can be obtained from the 
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

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
if nargin<8, F = 0.8; else
  if length(F)~=1; err(1,length(err)+1)=8; end; end;
if nargin<9, CR = 0.5; else
  if length(CR)~=1; err(1,length(err)+1)=9; end; end; 
if nargin<10, strategy = 7; else
  if length(strategy)~=1; err(1,length(err)+1)=10; end; end;
if nargin<11, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=11; end; end; 
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('devec3 (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
end

if (NP < 5)
   NP=5;
   fprintf(1,' NP increased to minimal value 5\n');
end
if ((CR < 0) | (CR > 1))
   CR=0.5;
   fprintf(1,'CR should be from interval [0,1]; set to default value 0.5\n');
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
bestvalit = bestval;              % best value of current iteration

bestmem = bestmemit;              % best member ever

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
pm4 = zeros(NP,D);              % initialize population matrix 4
pm5 = zeros(NP,D);              % initialize population matrix 5
bm  = zeros(NP,D);              % initialize bestmember  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rotd= (0:1:D-1);                % rotating index array (size D)
rt  = zeros(NP);                % another rotating index array
rtd = zeros(D);                 % rotating index array for exponential crossover
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array
a4  = zeros(NP);                % index array
a5  = zeros(NP);                % index array
ind = zeros(4);

%-------------------------------------------------------------------------------------------------------
% % % % %   x1 = linspace(XVmin(1), XVmax(1), 101);
% % % % %   x2 = linspace(XVmin(2), XVmax(2), 101);
% % % % %   x3 = zeros(length(x1), length(x2));
% % % % %         
% % % % %         % simply loop through the function (most functions expect 
% % % % %         % [N x 2] vectors as input, so meshgrid does not work)
% % % % %         for i = 1:length(x1)
% % % % %             for j = 1:length(x2)
% % % % %                 x3(i, j) = feval(fname,[x1(i), x2(j)]);
% % % % %             end
% % % % %         end
% % % % % 		figure(2);view(-40, 30);contour(x1', x2', x3);  hold on
% % % % %         plot(solution(:,1), solution(:,2), 'r.', 'MarkerSize', 20);
% % % % %         drawnow; pause(.1);
%-------------------------------------------------------------------------------------------------------

[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayýcýyý baþlat...
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
  popold = pop;                   % save the old population

  ind = randperm(4);              % index pointer array

  a1  = randperm(NP);             % shuffle locations of vectors
  rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
  a2  = a1(rt+1);                 % rotate vector locations
  rt = rem(rot+ind(2),NP);
  a3  = a2(rt+1);                
  rt = rem(rot+ind(3),NP);
  a4  = a3(rt+1);               
  rt = rem(rot+ind(4),NP);
  a5  = a4(rt+1);                

  pm1 = popold(a1,:);             % shuffled population 1
  pm2 = popold(a2,:);             % shuffled population 2
  pm3 = popold(a3,:);             % shuffled population 3
  pm4 = popold(a4,:);             % shuffled population 4
  pm5 = popold(a5,:);             % shuffled population 5

  for i=1:NP                      % population filled with the best member
    bm(i,:) = bestmemit;          % of the last iteration
  end

  mui = rand(NP,D) < CR;          % all random numbers < CR are 1, 0 otherwise

  if (strategy > 5)
    st = strategy-5;		  % binomial crossover
  else
    st = strategy;		  % exponential crossover
    mui=sort(mui');	          % transpose, collect 1's in each column
    for i=1:NP
      n=floor(rand*D);
      if n > 0
         rtd = rem(rotd+n,D);
         mui(:,i) = mui(rtd+1,i); %rotate column i by n
      end
    end
    mui = mui';			  % transpose back
  end % if (strategy > 5)...
  mpo = mui < 0.5;                % inverse mask to mui

  if (st == 1)                      % DE/best/1
    ui = bm + F*(pm1 - pm2);        % differential variation
    ui = popold.*mpo + ui.*mui;     % crossover
  elseif (st == 2)                  % DE/rand/1
    ui = pm3 + F*(pm1 - pm2);       % differential variation
    ui = popold.*mpo + ui.*mui;     % crossover
  elseif (st == 3)                  % DE/rand-to-best/1
    ui = popold + F*(bm-popold) + F*(pm1 - pm2);        
    ui = popold.*mpo + ui.*mui;     % crossover
  elseif (st == 4)                  % DE/best/2
    ui = bm + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
    ui = popold.*mpo + ui.*mui;           % crossover
  elseif (st == 5)                  % DE/rand/2
    ui = pm5 + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
    ui = popold.*mpo + ui.*mui;            % crossover
  end

%-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,ui(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       pop(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = ui(i,:);      % new best parameter vector ever
       end
    end
  end %---end for imember=1:NP

  bestmemit = bestmem;       % freeze the best member of this iteration for the coming 
                             % iteration. This is needed for some of the strategies.

% % % % %    if (refresh > 0)
% % % % %     if (rem(iter,refresh) == 0)
% % % % %         figure(2); contour(x1', x2', x3); 
% % % % % 		hold on; 
% % % % %         plot(solution(:,1), solution(:,2), 'r.', 'MarkerSize', 20);
% % % % %         plot(pop(:,1),pop(:,2),'b*','MarkerSize',8);
% % % % %         drawnow; pause(.1);
% % % % %         hold off
% % % % %     end
% % % % %    end

    %----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'\nIteration: %d,  Best: %f,  F: %f,  CR: %f,  NP: %d\n',iter,bestval,F,CR,NP);
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



 

