function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure]= penguen(fname,VTR,D,XVmin,XVmax,NP,itermax,refresh,tekrar)
 
%------------Penguins Search Optimization Algorithm-------------%
%----------------------- Yusuf ONDER ---------------------------%
%----------------------- 17.05.2015 ----------------------------%

%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('penguen 1st argument must be function name'); else 
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
  usage('penguen (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
end

%%mesela NP 5 den kücük girdi. Hata vermemesi için böyle bi durumda otomatik olarak 5 veriyoruz.
if (NP < 5)
   NP=5;
   fprintf(1,' NP increased to minimal value 5\n');
end

%%ayný	
if (itermax <= 0)
   itermax = 200;
   fprintf(1,'itermax should be > 0; set to default value 200\n');
end
refresh = floor(refresh);

%-----Initialize population and some arrays-------------------------------
%%Herbir tekrar için bulunan tüm en iyi  minimumlar.
GlobalMins=zeros(tekrar,itermax-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GlobalBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GlobalIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GlobalSure=zeros(tekrar,1); % herbir tekrar icin bulunan sure

oxygen = 0.7;
balik=0;

%BASLAMA
for r=1:tekrar
fprintf(1,'%d. tekrar\n',r);
%%pop popülasyon
%%np popülasyon sayýsý
%%Zeros 0 larla bi matris olusturuyo.
pop = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

%% tüm popülasyona bastan sona deger veriyo.
%: üstüstüne nin anlamý sonuna kadar.
% XVmin ve XVmax da oynama.
% rand(1,10) mesela 1*10 luk sýfýrla bir arasýnda float sayýlarla dolduruyo.
% nokta(.) skaler carpýmmýs
% bu for dan sonraNP sayýsý kadar rasgele verilerden oluþmuþ bir matris var 
for i=1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
     % toggle population
val       = zeros(1,NP);          % best population member in iteration
nfeval    = 0;                    % number of function evaluations //kullanýlan fonskiyon sayýsý

%------Evaluate the best member after initialization----------------------
%pop(ibest,:) demek su . birinci satýrdaki tüm üyeler iþte(ibest in 1 olmasý duurumdunda)
ibest   = 1;                      % start with first population member
val(1)  = feval(fname,pop(ibest,:));
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;%feval ý kullandýk o yüzden kullanýlan fonskiyon sayýsýný arttýrdýk.

%
for i=2:NP                        % check the remaining members
  val(i) = feval(fname,pop(i,:));
  nfeval  = nfeval + 1;
  if (val(i) < bestval)           % if member is better / eger üye daha iyisi yani bestValden minimum degeri daha kücükse best vali güncelliyoz
     ibest   = i;                 % save its location// indisi güncelledi en iyi üyenin
     bestval = val(i);			  % en iyi üyeyi güncelliyo.
  end   
end
bestmem = pop(ibest,:);         % best member of current iteration

  x1 = linspace(XVmin(1), XVmax(1), 101);
  x2 = linspace(XVmin(2), XVmax(2), 101);
  x3 = zeros(length(x1), length(x2));
        
        % simply loop through the function (most functions expect 
        % [N x 2] vectors as input, so meshgrid does not work)
        for i = 1:length(x1)
            for j = 1:length(x2)
                x3(i, j) = feval(fname,[x1(i), x2(j)]);
            end
        end
		figure(2);contour(x1', x2', x3); hold on; 
        plot(pop(:,1),pop(:,2),'b.','MarkerSize',20);
        drawnow
        hold off

%% burayý yani zamanlayýcý ve iter i while dan önce baslatcaz 

[delta_val]=sort(val);
durdurma_katsayisi=delta_val(NP)-delta_val(1);
iter = 1; %birinci tekrarýn içine girdik o yüzden 1 yapýyoruz
saat=tic; %zamanlayýcýyý baþlatýyoz

while ((iter < itermax) & (abs(durdurma_katsayisi) > VTR))%%eskiden bestval vardý
  popold = pop;                   % save the old population
  if (refresh > 0)
    if (rem(iter,refresh) == 0)
        figure(2);contour(x1', x2', x3); hold on; 
        plot(pop(:,1),pop(:,2),'b.','MarkerSize',20);
        drawnow
        hold off
    end
  end
  [sirala,indis]=sort(val);
  pop_sorted=pop(indis,:);
for i=1:NP,
    if(rand<oxygen)
    popold(i,:)=pop(i,:)+rand*abs((pop_sorted(1,:)-pop(i,:)));
    %balik(i,:) = rand(1,D);
    
    end
end 

%-----Select which vectors are allowed to enter the new population------------
  for i=1:NP
    tempval = feval(fname,popold(i,:));   % check cost of competitor
    nfeval  = nfeval + 1;
    if (tempval <= val(i))  % if competitor is better than value in "cost array"
       pop(i,:) = popold(i,:);  % replace old vector with new one (for new iteration)
       val(i)   = tempval;  % save value in "cost array"

       %----we update bestval only in case of success to save time-----------
       if (tempval < bestval)     % if competitor better than the best one ever
          bestval = tempval;      % new best value
          bestmem = popold(i,:);      % new best parameter vector ever
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
alpha = 0.2; 
end %end for (rr=1:tekrar) ...
figure(2);contour(x1', x2', x3);
hold on
plot(pop(:,1),pop(:,2),'r.','MarkerSize',20);
drawnow
hold off

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

