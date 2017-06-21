function [Gbestmem,Gbestval,nfeval1,Gmin,Giter,GSure] = bat(fname,VTR,D,XVmin,XVmax,NP,itermax,Qmin,Qmax,A,r,refresh,tekrar)
%ps=populasyon size,typically 10 to 40
%N_gen=Number of generations
%A=Loudness  (constant or decreasing)
%r=Pulse rate (constant or decreasing)

%Dýþardan girilecek deðerler.
%Yarasa Algoritmasý (Bat Algorithm)
% ------------------------------------------------------------------------------------
%  A=0.5; %A=Loudness  (constant or decreasing)
%  r=0.5; %r=Pulse rate (constant or decreasing)
%  Qmin=0; %min. frekans
%  Qmax=2; %mak.frekans
% ------------------------------------------------------------------------------------

alfa=0.98;
gama=0.98;


err=[];
if nargin<1, error('Bat Algorithm 1st argument must be function name'); else 
  if exist(fname)<1; err(1,length(err)+1)=1; end; end;
if nargin<2, VTR = 1.e-6; else 
  if length(VTR)~=1; err(1,length(err)+1)=2; end; end;
if nargin<3, D = 2; else
  if length(D)~=1; err(1,length(err)+1)=3; end; end; 
if nargin<4,XVmin = [-2 -2];else
  if length(XVmin)~=D; err(1,length(err)+1)=4; end; end; 
if nargin<5, XVmax = [2 2]; else
  if length(XVmax)~=D; err(1,length(err)+1)=5; end; end; 
if nargin<6 NP = 10*D; else
  if length(NP)~=1; err(1,length(err)+1)=7; end; end; 
if nargin<7, itermax = 200; else
  if length(itermax)~=1; err(1,length(err)+1)=8; end; end; 
% if nargin<8, N_gen = 1000; else
%   if length(N_gen)~=1; err(1,length(err)+1)=8; end; end;
if nargin<8, Qmin = 0; else
  if length(Qmin)~=1; err(1,length(err)+1)=9; end; end;
if nargin<9, Qmax = 2; else
  if length(Qmax)~=1; err(1,length(err)+1)=10; end; end; 
if nargin<10, A = 1; else
  if length(A)~=1; err(1,length(err)+1)=11; end; end;
if nargin<11, r = 1; else
  if length(r)~=1; err(1,length(err)+1)=12; end; end; 
if nargin<12, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=12; end; end; 
if nargin<13, tekrar = 1; else
  if length(tekrar)~=1; err(1,length(err)+1)=12; end; end; 
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('bat algorithm (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer,integer,integer)');    

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
%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar için bulunan tüm en iyi minumumlar
Globalbest=zeros(tekrar,D);         % herbir tekrar için bulunan en iyi çözümler
GBestval=zeros(tekrar,1); %herbir tekrar için bulunan en iyi minimim deðeri
GIter=zeros(tekrar,1); % kaçýncý iterasyonda çýktýðý deðer
GSure=zeros(tekrar,1); % herbir tekrar için bulunann süre
for t=1:tekrar
fprintf(1,'%d. tekrar\n',t);
pop = zeros(NP,D); %initialize pop to gain speed
for i=1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end

val       = zeros(1,NP);          % create and reset the "cost array"
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


bestmem = pop(ibest,:);            % best member ever

%------BAT ALGORITHM-Minimization---------------------------------------------

frequency=zeros(NP,1);   % Frequency
velocities=zeros(NP,D); %Velocities
% Iteration parameters
[p_val]=sort(val);
durdurma_katsayisi=p_val(NP)-p_val(1);
iter=1;
saat=tic; %zamanlayýcýyý baþlat...
% while ((iter < itermax) & (abs(bestval) > VTR))
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))
%     popold = pop;                   % save the old population
     % Start the iterations -- Bat Algorithm (essential part)  %
% 
        S=zeros(NP,D);
        for i=1:NP,
              frequency(i)=Qmin+(Qmin-Qmax)*rand;
              velocities(i,:)=velocities(i,:)+(pop(i,:)-bestmem)*frequency(i);
              S(i,:)=pop(i,:)+velocities(i,:);
              % Apply simple bounds/limits
              S(i,:)=simplebounds(S(i,:),XVmin,XVmax);
              % Pulse rate
              if rand>r
              % The factor 0.001 limits the step sizes of random walks 
                  S(i,:)=bestmem+0.001.*randn(1,D);
              end
            % Evaluate new solutions
               Fnew= feval(fname,S(i,:));
               
             % Update if the solution improves, or not too loud
               if (Fnew<val(i)) && (rand<A) ,
                    pop(i,:)=S(i,:);
                    val(i)=Fnew;
               end

              % Update the current best solution
              if Fnew<bestval,
                    bestmem=S(i,:);
                    bestval=Fnew;
              end
              
        end
%----Output section----------------------------------------------------------

  if (refresh > 0)
    if (rem(iter,refresh) == 0)
       fprintf(1,'Iteration: %d,  Best: %f,  NP: %d\n',iter,bestval,NP);
       for n=1:D
         fprintf(1,'best(%d) = %f\n',n,bestmem(n));
       end
    end
  end

  GlobalMins(t,iter)=bestval;
  [p_val]=sort(val);
  durdurma_katsayisi=p_val(NP)-p_val(1);
   %A ve r güncelle
   A=alfa*A;
   r=r*(1-exp(-gama));
   
   iter = iter + 1;
end %---end while ((iter < itermax) ...
sonsaat=toc(saat);
Globalbest(t,:)=bestmem;
GIter(t,1)=iter-1;
GBestval(t,1)=bestval;
GSure(t,1)=sonsaat;


end%end for (t=1:tekrar) ...

fprintf(1,'\nOrtalama Sonuc= %e\n',mean(GBestval));
fprintf(1,'\nStandart Sapma= %e\n',std(GBestval));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(Globalbest(:,1))/tekrar,sum(Globalbest(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(GIter));
fprintf(1,'\nOrtalama Sure%fsn\n',mean(GSure));

nfeval1=nfeval/tekrar;
%save ('ugur.mat','GlobalMins');
 if(tekrar==1) 
%      hold on; 
%      plot(pop(:,1),pop(:,2),'r.','MarkerSize',20); 
     Gmin=GlobalMins;
     Giter=GIter;
     Gbestmem=Globalbest;
     Gbestval=GBestval;
 else
     Gmin=GlobalMins;
     Giter=GIter;
     Gbestmem=Globalbest;
     Gbestval=GBestval;
 end
 
 end%fonk. sonu
 
% Application of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound vector
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bound vector 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
end



