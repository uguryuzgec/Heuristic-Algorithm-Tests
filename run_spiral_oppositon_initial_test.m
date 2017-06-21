% Sezgisel Optimizasyon Algoritmalari Haziran 2016
% Tufan-Ugur calismalar 
% Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm)
% Opposition based SOA gelistirildi (29.06.2016)
% Opposition-based and Random Initilization of initial population

clear all;
close all;

prevpath = path;
path(path, genpath(fileparts(mfilename('fullpath'))));
%--------------------------------------------------------------------------------------
% Refika ANBAR Sarmal Optimizasyon Algoritmasi (Spiral Optimization Algorithm)
% -------------------------------------------------------------------------------------
Theta=pi/4;
radius=0.95; 
 
% NP 			number of population (NP=10*D)		
NP = 20;
D = 2;

popnokta = zeros(NP,D);   
O_popnokta = zeros(NP,D); 
val = zeros(1,NP); 
O_val = zeros(1,NP);

[dims, lb, ub, solution, minimum, fonk]=RUN_ezimage;
        XVmin = lb; 
		XVmax = ub;
        D=dims;
        tekrar=1;
		
GlobalMeanFitness=zeros(tekrar,1);
GlobalStdFitness=zeros(tekrar,1);
GlobalMeanFitness_Opp=zeros(tekrar,1);
GlobalStdFitness_Opp=zeros(tekrar,1);
x1 = linspace(lb(1), ub(1), 101);
x2 = linspace(lb(2), ub(2), 101);
x3 = zeros(length(x1), length(x2));
  for i = 1:length(x1)
            for j = 1:length(x2)
                x3(i, j) = feval(fonk,([x1(i), x2(j)]));
            end
  end
  figure
  contour(x1', x2', x3);
  hold on
% -------------------------------------------------------------------------------------
for j=1:tekrar,
for i=1:NP,
   popnokta(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
   O_popnokta(i,:) = XVmin + XVmax - popnokta(i,:);
   val(i) = feval(fonk,popnokta(i,:)); 
   O_val(i) = feval(fonk,O_popnokta(i,:)); 
end
plot(popnokta(:,1), popnokta(:,2), 'r.', 'MarkerSize', 20);
plot(O_popnokta(:,1), O_popnokta(:,2), 'b.', 'MarkerSize', 20);
gec_val=[val O_val];
[gec_val,index]=sort(gec_val);
O_val = gec_val(1:NP);
GlobalMeanFitness(j)=mean(val);
GlobalMeanFitness_Opp(j)=mean(O_val);
GlobalStdFitness(j)=std(val);
GlobalStdFitness_Opp(j)=std(O_val);
end
 fprintf('Random sonuc: Mean-> %f,  Std-> %f\n',mean(GlobalMeanFitness),std(GlobalMeanFitness));
 fprintf('Opposition-based sonuc: Mean-> %f,  Std-> %f\n',mean(GlobalMeanFitness_Opp),std(GlobalMeanFitness_Opp));
