function [bestmem,bestval,nfeval,G,Gi] = crab(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar);
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
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('yengec (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
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
GlobalMins=zeros(tekrar,itermax-1);
Globalbest=zeros(tekrar,D);
GBestval=0; 
GIter=0;
GSure=0;
for r=1:tekrar
fprintf(1,'%d. tekrar\n',r);
malepop = zeros(NP,D); %initialize pop to gain speed
femalepop = zeros(NP,D); %initialize pop to gain speed
newmalepop = zeros(NP,D); %initialize pop to gain speed
newfemalepop = zeros(NP,D); %initialize pop to gain speed

nrmating=0; %çiftlesme sayisi...
resepsivite=200; %resepsivite derecesi... 1 den büyük olmak zorunda
esikdegeri=0.75; 
NoCross=0.75; %çaprazlama deðeri...
NoMut=0.15; %mutasyon deðeri...
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
   bestvalit = maleval(1);              % best value of current iteration
   bestval=maleval(1);
else
   bestmemit = femalepop(1,:);         % best member of current iteration
   bestvalit = femaleval(1);              % best value of current iteration
   bestval=femaleval(1);
end

bestmem = bestmemit;              % best member ever

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

iter = 1;
saat=tic; %zamanlayýcýyý baþlat...
while ((iter < itermax) && (abs(bestval) > VTR))
  for i=1:NP
      for j=1:NP

  prob=exp((maleval(i)*nrmating)/resepsivite);
    if prob>esikdegeri
        chield1=NoCross*malepop(i,:) + (1-NoCross)*femalepop(j,:);
        chield2=(1-NoCross)*malepop(i,:) + NoCross*femalepop(j,:);
        nrmating = nrmating + 1;
        neggs = round((beta*NP)/j); %neggs=üretilen yumurtalar
        nfeggs = round((neggs*maleval(i))/(100*nrmating)); % nfeggs=döllenen yumurtalar
        chield1M = chield1 + rand(1,D)*NoMut;
        chield2M = chield2 + rand(1,D)*NoMut;       
        end %if
        newmalepop(i,:) = chield1M;
        newfemalepop(j,:) = chield2M;
      end %for
      nrmating=0;
  end % for

%-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval1 = feval(fname,newmalepop(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    tempval2 = feval(fname,newfemalepop(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval1 <= maleval(i))  % if competitor is better than value in "cost array"
       malepop(i,:) = newmalepop(i,:);  % replace old vector with new one (for new iteration)
       maleval(i)   = tempval1;  % save value in "cost array"
    end   
       if (tempval2 <= femaleval(i))  % if competitor is better than value in "cost array"
       femalepop(i,:) = newfemalepop(i,:);  % replace old vector with new one (for new iteration)
       femaleval(i)   = tempval2;  % save value in "cost array"
       end

       %----we update bestval only in case of success to save time-----------
       if (tempval1 < bestval)     % if competitor better than the best one ever
          bestval = tempval1;      % new best value
          bestmem = newmalepop(i,:);      % new best parameter vector ever
       end
       
        if (tempval2 < bestval)     % if competitor better than the best one ever
          bestval = tempval2;      % new best value
          bestmem = newfemalepop(i,:);      % new best parameter vector ever
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
 iter = iter + 1;
  
end %---end while ((iter < itermax) ...
sonsaat=toc(saat);
Globalbest(r,:)=bestmem;
GIter=GIter+iter-1;
GBestval=GBestval+bestval;
GSure=GSure+sonsaat;
end %end for (r=1:tekrar) ...


fprintf(1,'\nOrtalama Sonuc= %f\n',GBestval/tekrar);
fprintf(1,'\nOrtalama çözüm= %f,%f\n',sum(Globalbest(:,1))/tekrar,sum(Globalbest(:,2))/tekrar);
fprintf(1,'\nOrtalama Ýterasyon= %f\n',GIter/tekrar);
fprintf(1,'\nOrtalama Süre= %fsn\n',GSure/tekrar);
%save ('ugur.mat','GlobalMins');
 if(tekrar==1) 
     hold on; 
     plot(malepop(:,1),malepop(:,2),'b.','MarkerSize',20); 
     plot(femalepop(:,1),femalepop(:,2),'r.','MarkerSize',20); 
     G=GlobalMins;
     Gi=GIter;
 else
     G=GlobalMins;
     Gi=GIter;
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



 

