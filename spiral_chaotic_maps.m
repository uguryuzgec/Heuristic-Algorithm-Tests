function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure] = spiral_chaotic_maps(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar,Theta,r,chaotic_option,solution)
 
%----------------Spiral Optimization Algorithm-----------------%
%----------------------- Refika ANBAR -------------------------%
%----------------------- 12.05.2015 ---------------------------%
%---------Chaotic maps have been implemented into the SOA------%
%----------------------- 18.11.2016 ---------------------------%

%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('spiral 1st argument must be function name'); else 
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
if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('spiral(string,scalar,scalar,vector,vector,scalar,scalar,scalar,scalar,scalar,scalar)');     
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
nfeval    = 0; %sayac                   % number of function evaluations

%-----Initialize population and some arrays-------------------------------
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure
% BASLAMA %
for rr=1:tekrar,
%--chaotic_option is used for selecting the chaotic function---%
%  1 : Logistic map, 2: Henon map, 3: Tent map, 4: Sinusoidal iterator, 5: Circle map %
X0 = rand(); Y0 = 0;
fprintf(1,'%d. tekrar\n',rr);
popnokta = zeros(NP,D); %initialize pop to gain speed
newpopnokta = zeros(NP,D);
% ------------------------------------------------

% aciklama yazzzzzzzzzz unutma :-)))
%M  matrisi D boyutlu
 M=eye(D);
 M(1,1)=cos(Theta);
 M(1,D)=-sin(Theta);
 M(D,1)=sin(Theta);
 M(D,D)=cos(Theta);
 I=eye(D);
% ------------------------------------------------
% generating the initial locations of n spiral(n boyutlu  spiralin baþlangýç  konumunu oluþturur.)

for i=1:NP 
    for j=1:D
        [random(1,j),X0,Y0]=random_chaotic_maps(chaotic_option,X0,Y0);
    end
   popnokta (i,:) = XVmin + random(1,:).*(XVmax - XVmin);
   % popnokta (i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
 
val       = zeros(1,NP);          % create and reset the "cost array"( her bir bireyin maliyet dizisi oluþtur.)
bestmem   = zeros(1,D);           % best population member ever(þimdiye kadar en  iyi populasyon elemaný)
bestmemit = zeros(1,D);           % best population member in iteration(iterasyondaki  en iyi populasyon elemeaný)

 
ibest   = 1;   %en iyi bir birey                   % start with first population member(ilk populasyon elemaný ile birlikte baþla)
val(1)  = feval(fname,popnokta(ibest,:));%populasyondaki bireylerin ilkini alýyor.

bestval = val(1);                 % best objective function value so far(þimdiye kadar  hedeflenen en iyi  fonksiyon deðeri )
nfeval  = nfeval + 1;
for i=2:NP                        % check the remaining members(kalan üye kontrol)
  val(i) = feval(fname,popnokta(i,:));%ikinci bireyin maliyet deðeri 
  nfeval  = nfeval + 1;
  if (val(i) < bestval)           % if member is better(eðer  val(i)<bestval ise)Eðer memeber daha iyi ise
     ibest   = i;                 % save its location(konumu kaydet) i indisini güncelle
     bestval = val(i);%best val güncelle
  end   
end
bestmemit = popnokta(ibest,:);         % best member of current iteration( güncel  iterasyonun en iyi elemaný)
bestvalit = bestval;              % best value of current iteration(güncel iterasyonun en iyi deðeri)
 
bestmem = bestmemit;
[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1;
saat=tic; %zamanlayiciyi baslat... 
% Iterations or pseudo time marching
%-------------------------------------------------------------------------------------------------------
  x1 = linspace(XVmin(1), XVmax(1), 101);
  x2 = linspace(XVmin(2), XVmax(2), 101);
  x3 = zeros(length(x1), length(x2));
        
%         simply loop through the function (most functions expect 
%         [N x 2] vectors as input, so meshgrid does not work)
        for i = 1:length(x1)
            for j = 1:length(x2)
                x3(i, j) = feval(fname,[x1(i), x2(j)]);
            end
        end
		figure;view(-40, 30);contour(x1', x2', x3); hold on; 
        plot(popnokta(:,1),popnokta(:,2),'b.','MarkerSize',20);
        plot(solution(:,1), solution(:,2), 'r.', 'MarkerSize', 20);
        drawnow; pause(.1);
        hold off
        
%-------------------------------------------------------------------------------------------------------
while ((iter < itermax) && (abs(durdurma_katsayisi) > VTR))     %%%%% start iterations
Xmerkez = bestmem; 
oldpopnokta=popnokta;
 
for i=1:NP,
    newpopnokta(i,:)=(r*M)*oldpopnokta(i,:)'-(r*M-I)*Xmerkez'; %sonuc 
end % end for i
%%%if pop() eger sinirlari gecerse ne olacak??? for i=1:NP,
  for i=1:NP
      for j=1:D,
         if(newpopnokta(i,j)<XVmin(j))
            newpopnokta(i,j)=XVmin(j);
        elseif(newpopnokta(i,j)>XVmax(j))
            newpopnokta(i,j)=XVmax(j);
        end
      end
  end
for i=1:NP,
    tempval = feval(fname,newpopnokta(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    popnokta(i,:) = newpopnokta(i,:);  % replace old vector with new one (for new iteration)
    val(i) = tempval;  % save value in "cost array"

% % %        %----we update bestval only in case of success to save time-----------
% % %        if (tempval < bestval)     % if competitor better than the best one ever
% % %           bestval = tempval;      % new best value
% % %           bestmem =newpopnokta(i,:);      % new best parameter vector ever
% % %        end
end %---end for imember=1:NP   
[sort_val,index] = sort(val);
sort_popnokta = popnokta(index,:);
bestmem = sort_popnokta(1,:);
bestval = sort_val(1);

% plot(popnokta(:,1),popnokta(:,2),'r.','MarkerSize',20);
% drawnow
% pause(); 

bestmemit = bestmem; 
%----Output section----------------------------------------------------------
 if (refresh > 0)
     if (rem(iter,refresh) == 0)
        contour(x1', x2', x3); 
        hold on; 
        plot(solution(:,1), solution(:,2), 'r.', 'MarkerSize', 20);
        plot(popnokta(:,1),popnokta(:,2),'b.','MarkerSize',14);
        drawnow; pause(.1);
        hold off
     end
end

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
end
function [chaotic_output,X0,Y0] = random_chaotic_maps(option,X0,Y0)
switch option
	case 1
		% 1 : Logistic map
		a = 3.99;
		chaotic_output = a*X0*(1-X0);
		X0 = chaotic_output;
	case 2
		% 2 : Henon map
		a = 1.4; b = 0.3;
		X1 = 1-a*X0^2+Y0;
		chaotic_output = (X0*b+0.59);
		Y0 =  X0*b;
        X0 = X1;		
	case 3
		% 3 : Tent map
		if (X0<0.7)
			chaotic_output = X0/0.7;
		else
			chaotic_output = (10/3)*X0*(1-X0);
		end
		X0 = chaotic_output;
	case 4
		% 4 : Sinusoidal iterator
		chaotic_output = sin(pi*X0);
		X0 = chaotic_output;
	case 5
		% 5 : Circle map
		a = 0.5; b = 0.2;
		chaotic_output = mod(X0+b-(a/(2*pi))*sin(2*pi*X0),1);
		X0 = chaotic_output;
end
end