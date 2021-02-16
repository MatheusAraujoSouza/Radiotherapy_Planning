% Cecilia B. Barboza
% Programa preditor-corretor para o metodo de pontos interiores aplicado ao modelo de tratamento de
%radioterapia utilizando a analise absoluta

clear
load dosematrix_15_64_10.mat
%load dosematrix_5_64_32.mat
%load dosematrix_1_64_32.mat
tol= 1e-4;

UniLevel= 2;
p = UniLevel / 100;

erro=sqrt(eps);
erro= 1e-4;

load pre
TumorRowIndex= find(prescription(:,:,2) == 1);
CriticalRowIndex= find(prescription(:,:,2) == 2);
RegularRowIndex= find(prescription(:,:,2) == 3);
%RegularRowIndex= find(prescription(:,:,2) == 3 | prescription(:,:,2) == 0);

At= sparse(A(TumorRowIndex,:));
Ac= sparse(A(CriticalRowIndex,:));
Ag= sparse(A(RegularRowIndex,:));

[mt,n]= size(At);
[mc,n]= size(Ac);
[mg,n]= size(Ag);

CheckColSum = ones(1,mt)*At;
TumorColIndex = find(CheckColSum >= 0.1);
At= sparse(At(:,TumorColIndex));
Ac= sparse(Ac(:,TumorColIndex));
Ag= sparse(Ag(:,TumorColIndex));
AA= A;
A= sparse([At; Ac; Ag]);
[m,n]= size(A);

tau=0.995;
np= m+mt+n+4;
sigma=1/sqrt(np);

tg = prescription(TumorRowIndex);
ut=(1+p)*tg;
lt=(1-p)*tg + tol;
uc= prescription(CriticalRowIndex);
ug= prescription(RegularRowIndex);
w= norm(lt, inf)/ tol;
lambda= min(lt);

%flops(0)

StartTime = clock;

nut=norm(ut);
nlt=norm(lt)+1;
nuc=norm(uc);
nug=norm(ug);

mt1= mt+1;
mtc= mt+mc;
mtc1= mt+mc+1;

%Ponto inicial primal

T= sqrt(6*(mt+3));
su= ut-2*lt-lambda;
M= sparse(A*A' + [speye(mt)/3 zeros(mt,mc+mg); zeros(mc+mg,mt) speye(mc+mg)]);
y= M\([(sum(su)/(mt+3)-su)/6; zeros(mc,1); ones(mg,1)]);
v3= M\[ones(mt,1)/T; ones(mc+mg,1)];
y= y - v3*(sum(y(1:mt))/T+sum(y(mt1:m))/(1+sum(v3(1:mt))/T+sum(v3(mt1:m))));
x= A'*y;
a= y(1:mt);
sc= y(mt1:mtc);
sg= y(mtc1:m);
su= a-su;
su= (su-sum(su)/(mt+3))/3;
sl= (ut+a-su)/2;
a= su+sl-a;
sl= -su;
sigmat= (lambda-sum(su))/2;
T= sigmat + sum(su);
gama= -sum(sc);
beta= -sum(sg);
pp= max(norm([ut; lt; lambda; ones(mg,1)])/100,100);
pp= min(pp,1e4);
pd= max(pp/40,sqrt(w*w+2)*.425);
pd = min(pd,1e3);
pp= max(pp, -min([x; a; sc; gama; sg; beta; sl; su; T; sigmat]));
x= max(pp,x);
a= max(pp,a);
gama= max(pp,gama);
beta= max(pp,beta);
T= max(pp,T);
sc= max(pp,sc);
sg= max(pp,sg);
sl= max(pp,sl);
su= max(pp,ut-a);
sigmat= max(pp,lambda-T);

%Ponto inicial dual

ya=zeros(mt,1);
yu=ones(mt,1);
yl=ones(mt,1);
yc=ones(mc,1);
yg=ones(mg,1);
pit=pd;
zetat=pd;
zetac=pd;
zetag=pd;
zx=ones(n,1);

%Residuos

r1=ut-a-su;
r2=lt-a-T+sl;
r3=a-At*x;
r4=gama-Ac*x-sc;
r5=1-Ag*x+beta-sg;
r6=lambda-T-sigmat;
r7=yu-yl+ya;
r8=w-sum(yl)+pit-zetat;
r9=1-sum(yc)-zetac;
r10=1-sum(yg)-zetag;
r11=-At'*ya+Ac'*yc+Ag'*yg-zx;

gama1=(su'*yu + sl'*yl + sc'*yc + sg'*yg + sigmat*pit +T*zetat + gama*zetac+beta*zetag + x'*zx);
gamarel=gama1/(1+abs(w*T+gama+beta)+abs(-ut'*yu +lt'*yl-sum(yg)-lambda*pit))/np;

merro=max(gamarel, norm(r1/nut));
merro=max(merro, norm(r2/nlt));
merro=max(merro, norm(r3/mt));
merro=max(merro, norm(r4/mc));
merro=max(merro, norm(r5/mg));
merro=max(merro, norm(r6/mt));
merro=max(merro, norm(r7/mt));
merro=max(merro, norm(r8/mt));
merro=max(merro, norm(r9/mc));
merro=max(merro, norm(r10/mg));
merro=max(merro, norm(r11/n));

k=0;
%EE= sparse(zeros(m));
while(merro>erro)
    if k > 50
        disp('Maximo iteracoes')
        return
    end

   %Matrizes Diagonais

D1= yu./su;
D2= yl./sl;
D3= pit/sigmat;
D4= zetat/T;
D5= yc./sc;
D6= gama/zetac;
D7= yg./sg;
D8= beta/zetag;
D9= zx./x;
Di= D1+D2;

Dtil=[Di; D5; D7];
D10= D2./Di;

D2= sl./yl;
D11= 1./(su./yu + D2);
Dch=[-1/(sum(D1.*D10)+D3+D4); -D6; -D8];

%Residuos

r12=-yu.*su;
r13=-yl.*sl;
r14=-yc.*sc;
r15=-yg.*sg;
r16=-pit*sigmat;
r17=-T*zetat;
r18=-gama*zetac;
r19=-beta*zetag;
r20=-x.*zx;
r21= r7+r12./su-r13./sl;
r22= r8+r16/sigmat-r17/T-sum(r13./sl);
r23= r9-sum(r14./sc)-r18/gama;
r24= r10-sum(r15./sg)-r19/beta;
r25= r11+Ac'*(r14./sc)+Ag'*(r15./sg)-r20./x;
r26= r21-yu.*(r1./su);
r27=r22-pit*(r6/sigmat);
r28=r2-r26.*D2;
r29=r27-sum(r26);
r30=r3+r28.*D10;
r31=r29-sum(D11.*r28);
rp=[r30;r4;r5];
rptil=[r31;r23;r24];
Ee(1:mt,1)=D10*Dch(1);
Ee(mt1:mtc)=-Dch(2);
Ee(mtc1:m)=-Dch(3);
Ep(1:mt,1)=Ee(1:mt)*r31;
Ep(mt1:mtc)=Ee(mt1:mtc)*r23;
Ep(mtc1:m)=Ee(mtc1:m)*r24;
rs=rp-Ep;
Etil(1:mt,1)=D10;
Etil(mt1:mtc)=-D5;
Etil(mtc1:m)=-D7;
rr= Etil;
rr(1:mt)= rr(1:mt).*Di;
v3= rr.*Ee;
E3(1,1)= 1/(1 - sum(v3(1:mt)));
E3(2,1)= 1/(1 - sum(v3(mt1:mtc)));
E3(3,1)= 1/(1 - sum(v3(mtc1:m)));

%dx(1:mt,1)= rr(1:mt)*E3(1);
%dx(mt+1:mt+mc)= rr(mt+1:mt+mc)*E3(2);
%dx(mt+mc+1:m)= rr(mt+mc+1:m)*E3(3);
%EE(1:mt,1:mt)= sparse(speye(mt) + Ee(1:mt)*dx(1:mt)');
%EE(mt+1:mt+mc,mt+1:mt+mc)= sparse(speye(mc) + Ee(mt+1:mt+mc)*dx(mt+1:mt+mc)');
%EE(mt+mc+1:m,mt+mc+1:m)= sparse(speye(mg) + Ee(mt+mc+1:m)*dx(mt+mc+1:m)');

v3= rr.*rs;
v3(1:mt)= sum(v3(1:mt))*E3(1);
v3(mt1:mtc)= sum(v3(mt1:mtc))*E3(2);
v3(mtc1:m)= sum(v3(mtc1:m))*E3(3);

r=-A'*(Dtil.*(Ee.*v3+rs)) + r25;

%Direcoes

%dx=-(sparse(A'*sparse(diag(Dtil))*EE*A+sparse(diag(D9)))\r);
%[L,p]=chol((sparse(A'*sparse(diag(Dtil))*EE*A+sparse(diag(D9)))));
[L,p]=chol((sparse(A'*sparse(diag(Dtil))*A+sparse(diag(D9)))));
L=sparse(L);
dx=-(L\(L'\r));

d2= A*dx;
v3= Etil.*d2;
d1(1)= sum(v3(1:mt).*Di);
d1(2,1)= sum(v3(mt1:mtc));
d1(3,1)= sum(v3(mtc1:m));

Ep(1:mt)= Ee(1:mt)*E3(1);
Ep(mt1:mtc)= Ee(mt1:mtc)*E3(2);
Ep(mtc1:m)= Ee(mtc1:m)*E3(3);
Ep= Ep.*Dtil;

v3= Etil.*(A*(L\(L'\(At'*Ep(1:mt)))));
M3(1,1)= 1 + sum(v3(1:mt).*Di);
M3(2,1)= sum(v3(mt1:mtc));
M3(3,1)= sum(v3(mtc1:m));
v3= Etil.*(A*(L\(L'\(Ac'*Ep(mt1:mtc)))));
M3(1,2)= sum(v3(1:mt));
M3(2,2)= 1 + sum(v3(mt1:mtc));
M3(3,2)= sum(v3(mtc1:m));
v3= Etil.*(A*(L\(L'\(Ag'*Ep(mtc1:m)))));
M3(1,3)= sum(v3(1:mt));
M3(2,3)= sum(v3(mt1:mtc));
M3(3,3)= 1 + sum(v3(mtc1:m));

d1= E3.*(M3\d1);

Ep(1:mt)= Ee(1:mt)*d1(1);
Ep(mt1:mtc)= Ee(mt1:mtc)*d1(2);
Ep(mtc1:m)= Ee(mtc1:m)*d1(3);

dx= dx - L\(L'\(A'*(Ep.*Dtil)));

d2= A*dx-rs;
rs= rr.*d2;
rs(1:mt)= sum(rs(1:mt))*E3(1);
rs(mt1:mtc)= sum(rs(mt1:mtc))*E3(2);
rs(mtc1:m)= sum(rs(mtc1:m))*E3(3);
d2= -(Ee.*rs+d2);
dya=d2(1:mt).*Di;
dsc=d2(mt1:mtc);
dsg=d2(mtc1:m);
d1(1)= Etil(1:mt)'*dya;
d1(2)= Etil(mt1:mtc)'*dsc;
d1(3)= Etil(mtc1:m)'*dsg;
d1=(rptil-d1).*Dch;
dT=d1(1);
dgama=d1(2);
dbeta=d1(3);
da=D10.*(r28-dT-(dya.*D2));
dsl=-(r26+D1.*da+dya).*D2;
dsu=r1-da;
dsigmat=r6-dT;
dyu=(r12-yu.*dsu)./su;
dyl=(r13-yl.*dsl)./sl;
dyc=(r14-yc.*dsc)./sc;
dyg=(r15-yg.*dsg)./sg;
dpit=(r16-pit*dsigmat)/sigmat;
dzetat=(r17-zetat*dT)./T;
dzetac=(r18-zetac*dgama)/gama;
dzetag=(r19-zetag*dbeta)/beta;
dzx=(r20-zx.*dx)./x;

%Tamanho do passo primal

  alfap=min(1,-tau/min(-tau,min(dx./x)));
 
  alfap=min(alfap,-tau/min(-tau,min(dsu./su)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsl./sl)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsc./sc)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsg./sg)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsigmat/sigmat)));
    
  alfap=min(alfap,-tau/min(-tau,min(dT/T)));
    
  alfap=min(alfap,-tau/min(-tau,min(dbeta/beta)));
    
  alfap=min(alfap,-tau/min(-tau,min(dgama/gama)));
    
    %Tamanho do passo dual
    
   alfad=min(1,-tau/min(-tau,min(dyu./yu)));
    
   alfad=min(alfad,-tau/min(-tau,min(dyl./yl)));
    
   alfad=min(alfad,-tau/min(-tau,min(dyc./yc)));
    
   alfad=min(alfad,-tau/min(-tau,min(dyg./yg)));
    
   alfad=min(alfad,-tau/min(-tau,min(dpit/pit)));
    
   alfad=min(alfad,-tau/min(-tau,min(dzetat/zetat)));
   
   alfad=min(alfad,-tau/min(-tau,min(dzetac/zetac)));
   
   alfad=min(alfad,-tau/min(-tau,min(dzetag/zetag)));
     
   alfad=min(alfad,-tau/min(-tau,min(dzx./zx)));
    
gamatil=(su + alfap*dsu)'*(yu + alfad*dyu) + (sl + alfap*dsl)'*(yl + alfad*dyl);
gamatil=gamatil + (sc + alfap*dsc)'*(yc + alfad*dyc) + (sg + alfap*dsg)'*(yg + alfad*dyg);
gamatil=gamatil + (sigmat + alfap*dsigmat)'*(pit + alfad*dpit) + (T + alfap*dT)'*(zetat + alfad*dzetat);
gamatil=gamatil + (gama + alfap*dgama)'*(zetac + alfad*dzetac) + (beta + alfap* dbeta)'*(zetag + alfad*dzetag);
 gamatil=gamatil + (x + alfap*dx)'*(zx + alfad*dzx);

 if gama1 > 1
    mi=sigma*gamatil/np*(gamatil/gama1)*(gamatil/gama1);
 else
    %mi=sigma*gamatil/np*(gamatil/gama1);%*(gamatil/gama1);
    mi=sigma*gama1*gama1/np;
 end

   r12=r12+mi;
   r13=r13+mi;
   r14=r14+mi;
   r15=r15+mi;
   r16=r16+mi;
   r17=r17+mi;
   r18=r18+mi;
   r19=r19+mi;
   r20=r20+mi;

r21= r7+r12./su-r13./sl;
r22= r8+r16/sigmat-r17/T-sum(r13./sl);
r23= r9-sum(r14./sc)-r18/gama;
r24= r10-sum(r15./sg)-r19/beta;
r25= r11+Ac'*(r14./sc)+Ag'*(r15./sg)-r20./x;
r26= r21-yu.*(r1./su);
r27=r22-pit*(r6/sigmat);
r28=r2-r26.*D2;
r29=r27-sum(r26);
r30=r3+r28.*D10;
r31=r29-sum(D11.*r28);
rp=[r30;r4;r5];
rptil=[r31;r23;r24];
Ep(1:mt,1)=Ee(1:mt)*r31;
Ep(mt1:mtc)=Ee(mt1:mtc)*r23;
Ep(mtc1:m)=Ee(mtc1:m)*r24;
rs=rp-Ep;

v3= rr.*rs;
v3(1:mt)= sum(v3(1:mt))*E3(1);
v3(mt1:mtc)= sum(v3(mt1:mtc))*E3(2);
v3(mtc1:m)= sum(v3(mtc1:m))*E3(3);

r=-A'*(Dtil.*(Ee.*v3+rs)) + r25;

dx=-(L\(L'\r));

d2= A*dx;
v3= Etil.*d2;
d1(1)= sum(v3(1:mt).*Di);
d1(2)= sum(v3(mt1:mtc));
d1(3)= sum(v3(mtc1:m));

d1= E3.*(M3\d1);

Ep(1:mt)= Ee(1:mt)*d1(1);
Ep(mt1:mtc)= Ee(mt1:mtc)*d1(2);
Ep(mtc1:m)= Ee(mtc1:m)*d1(3);

dx= dx - L\(L'\(A'*(Ep.*Dtil)));

d2= A*dx-rs;
rs= rr.*d2;
rs(1:mt)= sum(rs(1:mt))*E3(1);
rs(mt1:mtc)= sum(rs(mt1:mtc))*E3(2);
rs(mtc1:m)= sum(rs(mtc1:m))*E3(3);
d2= -(Ee.*rs+d2);
dya=d2(1:mt).*Di;
dsc=d2(mt1:mtc);
dsg=d2(mtc1:m);
d1(1)= Etil(1:mt)'*dya;
d1(2)= Etil(mt1:mtc)'*dsc;
d1(3)= Etil(mtc1:m)'*dsg;
d1=(rptil-d1).*Dch;
dT=d1(1);
dgama=d1(2);
dbeta=d1(3);
da=D10.*(r28-dT-(dya.*D2));
dsl=-(r26+D1.*da+dya).*D2;
dsu=r1-da;
dsigmat=r6-dT;
dyu=(r12-yu.*dsu)./su;
dyl=(r13-yl.*dsl)./sl;
dyc=(r14-yc.*dsc)./sc;
dyg=(r15-yg.*dsg)./sg;
dpit=(r16-pit*dsigmat)/sigmat;
dzetat=(r17-zetat*dT)./T;
dzetac=(r18-zetac*dgama)/gama;
dzetag=(r19-zetag*dbeta)/beta;
dzx=(r20-zx.*dx)./x;

%Tamanho do passo primal

  alfap=min(1,-tau/min(-tau,min(dx./x)));
 
  alfap=min(alfap,-tau/min(-tau,min(dsu./su)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsl./sl)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsc./sc)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsg./sg)));
    
  alfap=min(alfap,-tau/min(-tau,min(dsigmat/sigmat)));
    
  alfap=min(alfap,-tau/min(-tau,min(dT/T)));
    
  alfap=min(alfap,-tau/min(-tau,min(dbeta/beta)));
    
  alfap=min(alfap,-tau/min(-tau,min(dgama/gama)));
    
    %Tamanho do passo dual
    
   alfad=min(1,-tau/min(-tau,min(dyu./yu)));
    
   alfad=min(alfad,-tau/min(-tau,min(dyl./yl)));
    
   alfad=min(alfad,-tau/min(-tau,min(dyc./yc)));
    
   alfad=min(alfad,-tau/min(-tau,min(dyg./yg)));
    
   alfad=min(alfad,-tau/min(-tau,min(dpit/pit)));
    
   alfad=min(alfad,-tau/min(-tau,min(dzetat/zetat)));
   
   alfad=min(alfad,-tau/min(-tau,min(dzetac/zetac)));
   
   alfad=min(alfad,-tau/min(-tau,min(dzetag/zetag)));
     
   alfad=min(alfad,-tau/min(-tau,min(dzx./zx)));
    
    %Proximo ponto primal
    
    a=a+alfap*da;
    x=x+alfap*dx;
    su=su+alfap*dsu;
    sl=sl+alfap*dsl;
    sc=sc+alfap*dsc;
    sg=sg+alfap*dsg;
    sigmat=sigmat+alfap*dsigmat;
    T=T+alfap*dT;
    beta=beta+alfap*dbeta;
    gama=gama+alfap*dgama;
    
    %Proximo ponto dual
    
    ya=ya+alfad*dya;
    yu=yu+alfad*dyu;
    yl=yl+alfad*dyl;
    yc=yc+alfad*dyc;
    yg=yg+alfad*dyg;
    pit=pit+alfad*dpit;
    zetat=zetat+alfad*dzetat;
    zetac=zetac+alfad*dzetac;
    zetag=zetag+alfad*dzetag;
    zx=zx+alfad*dzx;

r1=ut-a-su;
r2=lt-a-T+sl;
r3=a-At*x;
r4=gama-Ac*x-sc;
r5=1-Ag*x+beta-sg;
r6=lambda-T-sigmat;
r7=yu-yl+ya;
r8=w-sum(yl)+pit-zetat;
r9=1-sum(yc)-zetac;
r10=1-sum(yg)-zetag;
r11=-At'*ya+Ac'*yc+Ag'*yg-zx;

gama1=(su'*yu + sl'*yl + sc'*yc + sg'*yg + sigmat*pit +T*zetat + gama*zetac+beta*zetag + x'*zx);
gamarel=gama1/(1+abs(w*T+gama+beta)+abs(-ut'*yu +lt'*yl-sum(yg)-lambda*pit))/np;

merro=max(gamarel, norm(r1/nut));
merro=max(merro, norm(r2/nlt));
merro=max(merro, norm(r3/mt));
merro=max(merro, norm(r4/mc));
merro=max(merro, norm(r5/mg));
merro=max(merro, norm(r6/mt));
merro=max(merro, norm(r7/mt));
merro=max(merro, norm(r8/mt));
merro=max(merro, norm(r9/mc));
merro=max(merro, norm(r10/mg));
merro=max(merro, norm(r11/n))

k=k+1;

end
TotalTime = etime(clock,StartTime)
flops
gama1
merro
k
return

[m,n]= size(AA);
ThePlan = zeros(n,1);
ThePlan(TumorColIndex) = x;
RadLevel = AA*ThePlan;
n= 64;
DisplayMatrix = zeros(n,n);

for j=1:n
   for i=1:n
        DisplayMatrix(i,j)=RadLevel(64*(j-1)+i);
   end
end

figure;
ViewPrescription = prescription(:,:,1);
ViewPrescription(RegularRowIndex) = 0;
contour(ViewPrescription)
hold on;
contour(DisplayMatrix)
hold off;

figure;
mesh(DisplayMatrix)

return
