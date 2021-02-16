% Cecilia B. Barboza
% Programa preditor-corretor para o metodo de pontos interiores aplicado ao modelo de tratamento de
%radioterapia da analise media

clear

global CACHE_SIZE;
global LOOP_LEVEL;

CACHE_SIZE = 16;
LOOP_LEVEL = 8;

%load dosematrix_15_64_10.mat
%load dosematrix_5_64_32.mat
load dosematrix_1_64_32.mat
%load /export/home0/aurelio/rad/lex.mat
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

CheckColSum= ones(1,mt)*At;
TumorColIndex= find(CheckColSum >= 0.1);
At= sparse(At(:,TumorColIndex));
Ac= sparse(Ac(:,TumorColIndex));
Ag= sparse(Ag(:,TumorColIndex));
AA= A;
A= sparse([At; Ac; Ag]);
[m,n]= size(A);
%symbfct(sparse(A'*A));

%Numero de variaveis

np=m+6*n;
sigma=1/sqrt(np);
tau=0.99995;

tg = prescription(TumorRowIndex);
ut=(1+p)*tg;
lt=(1-p)*tg + tol;
uc= prescription(CriticalRowIndex);
ug= prescription(RegularRowIndex);
w= norm(lt, inf)/ tol;

flops(0)

StartTime = clock;

nut=norm(ut);
nlt=norm(lt)+1;
nuc=norm(uc);
nug=norm(ug);
nmt=w*mt+1;
nmc=mc+1;
nmg=mg+1;

mt1= mt+1;
mtc= mt+mc;
mtc1= mt+mc+1;

%Ponto inicial primal

y= sparse(A*A' + [3/8*speye(mt) zeros(mt,mc+mg); zeros(mc+mg,mt) 2*speye(mc+mg)])\([(3*ut+lt)/8; zeros(mc+mg,1)]);
x= A'*y;
a= y(1:mt);
sc= y(mt1:mtc);
c= -sc;
sg= y(mtc1:m);
g= - sg;
sl= (5*ut+3*a-lt)/8;
su= (ut-lt-a)/4;
st= (su+lt)/2;
a= sl-su-a;
t= st-su;
pp= max(norm([uc; ug; ut; lt])/100,100);
pp= min(pp,1e4);
pd= max(pp/40,sqrt(w*w*mt+mc+mg)*.425);
pd = min(pd,1e3);
pp= max(pp, -min([x; a; sc; c; sg; g; sl; su; st; t]));
x= max(pp,x);
a= max(pp,a);
c= max(pp,c);
g= max(pp,g);
t= max(pp,t);
sc= max(pp,sc);
sg= max(pp,sg);
sl= max(pp,sl);
st= max(pp,ut-a);
su= max(pp,lt-t);

%Ponto inicial dual

ya=zeros(mt,1);
yu=ones(mt,1);
yl=yu;
yc=ones(mc,1);
yg=ones(mg,1);
yt=yu;
wc=yc;
wg=yg;
zx=ones(n,1);
zt=(w+pd)*ones(mt,1);

%Residuos

r1=ut-a-su;
r2=lt-a-t+sl;
r3=a-At*x;
r4=-Ac*x-sc+c;
r5=ug-Ag*x+g-sg;
r6=lt-t-st;
r7=yu-yl+ya;
r8=-At'*ya+Ac'*yc+Ag'*yg-zx;
r9=w-yl+yt-zt;
r10=1-yc-wc;
r11=1-yg-wg;

gama=(su'*yu + sl'*yl + sc'*yc + sg'*yg + st'*yt + x'*zx + t'*zt + c'*wc + g'*wg);
gamarel=gama/(1+abs(w*sum(t)+sum(c)+sum(g))+abs(-ut'*yu+lt'*yl-ug'*yg-lt'*yt))/np;

merro=max(gamarel, norm(r1/nut));
merro=max(merro, norm(r2/nlt));
merro=max(merro, norm(r3/(norm(a)+1)));
merro=max(merro, norm(r4/(norm(c)+1)));
merro=max(merro, norm(r5/nug));
merro=max(merro, norm(r6/nlt));
merro=max(merro, norm(r7/(norm(ya)+1)));
merro=max(merro, norm(r8/norm(zx)));
merro=max(merro, norm(r9/nmt));
merro=max(merro, norm(r10/nmc));
merro=max(merro, norm(r11/nmg));

k=0;
while(merro>erro)
    if k > 40
        disp('Maximo iteracoes')
        return
    end

   %Matrizes Diagonais

Dc=sc./yc;
Dg=sg./yg;
Du=yu./su;
Dt=yt./st;
Dz=zt./t;
D5= wc./c;
D6= wg./g;
D1= 1+Dc.*D5;
D2= 1+Dg.*D6;
D4= yl./sl;
D3= Du+D4;
D7= zx./x;
D8= D4+Dt+Dz;
D9= D4./D8;
Dtil= [D3-D9.*D4; D5./D1; D6./D2];

%Residuos

r12=-yu.*su;
r13=-yl.*sl;
r14=-yc.*sc;
r15=-yg.*sg;
r16=-yt.*st;
r17=-x.*zx;
r18=-t.*zt;
r19=-c.*wc;
r20=-g.*wg;
r21=r7+r12./su-r13./sl;
r22=r8+Ac'*(r14./sc)+Ag'*(r15./sg)-r17./x;
r23=r9-r13./sl+r16./st-r18./t;
r24=r10-r14./sc-r19./c;
r25=r11-r15./sg-r20./g;
r26=r4+(Dc.*r24);
r27=r5+(Dg.*r25);
r28=r21-(Du.*r1)-(D4.*r2);
r29=r22+Ac'*r24+Ag'*r25;
r30=r23-(D4.*r2)-(Dt.*r6);
r31=r29+At'*r28-At'*(D9.*r30);
rp=[r3;r26;r27];
rptil=A'*(rp.*Dtil)-r31;

%Direcoes

[L,p]=chol(sparse(A'*sparse(diag(Dtil))*A+sparse(diag(D7))));
L=sparse(L);
dx=L\(L'\rptil);
%dx= linsys(sparse(A'*sparse(diag(Dtil))*A+sparse(diag(D7))),rptil);
d=A*dx-rp;
da=d(1:mt);
dc=d(mt1:mtc)./D1;
dg=d(mtc1:m)./D2;
dt=-(r30+D4.*da)./D8;
dya=-D3.*da-D4.*dt-r28;
dsu=r1-da;
dsl=-r2+da+dt;
dsc=-Dc.*(D5.*dc+r24);
dsg=-Dg.*(D6.*dg+r25);
dst=r6-dt;
dyu=r12./su-Du.*dsu;
dyl=r13./sl-D4.*dsl;
dyc=(r14-yc.*dsc)./sc;
dyg=(r15-yg.*dsg)./sg;
dyt=r16./st-Dt.*dst;
dzx=r17./x-D7.*dx;
dzt=r18./t-Dz.*dt;
dwc=r19./c-D5.*dc;
dwg=r20./g-D6.*dg;

%Tamanho do passo primal

  alfap=min(1,-tau/min(-tau,min(dx./x)));
  alfap=min(alfap,-tau/min(-tau,min(dt./t)));
  alfap=min(alfap,-tau/min(-tau,min(dc./c)));
  alfap=min(alfap,-tau/min(-tau,min(dg./g)));
  alfap=min(alfap,-tau/min(-tau,min(dsu./su)));
  alfap=min(alfap,-tau/min(-tau,min(dsl./sl)));
  alfap=min(alfap,-tau/min(-tau,min(dsc./sc)));
  alfap=min(alfap,-tau/min(-tau,min(dsg./sg)));
  alfap=min(alfap,-tau/min(-tau,min(dst./st)));

%Tamanho do passo dual

   alfad=min(1,-tau/min(-tau,min(dyu./yu)));
   alfad=min(alfad,-tau/min(-tau,min(dyl./yl)));
   alfad=min(alfad,-tau/min(-tau,min(dyc./yc)));
   alfad=min(alfad,-tau/min(-tau,min(dyg./yg)));
   alfad=min(alfad,-tau/min(-tau,min(dyt./yt)));
   alfad=min(alfad,-tau/min(-tau,min(dzx./zx)));
   alfad=min(alfad,-tau/min(-tau,min(dzt./zt)));
   alfad=min(alfad,-tau/min(-tau,min(dwc./wc)));
   alfad=min(alfad,-tau/min(-tau,min(dwg./wg)));

   gamatil=(su + alfap*dsu)'*(yu + alfad*dyu) + (sl + alfap*dsl)'*(yl + alfad*dyl);
   gamatil=gamatil + (sc + alfap*dsc)'*(yc + alfad*dyc) + (sg + alfap*dsg)'*(yg + alfad*dyg);
   gamatil=gamatil + (st + alfap*dst)'*(yt + alfad*dyt) + (x + alfap*dx)'*(zx + alfad*dzx);
   gamatil=gamatil + (t + alfap*dt)'*(zt + alfad*dzt) + (c + alfap*dc)'*(wc + alfad*dwc);
   gamatil=gamatil + (g + alfap*dg)'*(wg + alfad*dwg);

   if gama > 1
    mi=sigma*gamatil/np*(gamatil/gama)*(gamatil/gama);
   else
    mi= sigma*gama*gama/np;
    %mi=sigma*gamatil/np*(gamatil/gama)*(gamatil/gama);
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

   r21=r7+r12./su-r13./sl;
   r22=r8+Ac'*(r14./sc)+Ag'*(r15./sg)-r17./x;
   r23=r9-r13./sl+r16./st-r18./t;
   r24=r10-r14./sc-r19./c;
   r25=r11-r15./sg-r20./g;
   r26=r4+(Dc.*r24);
   r27=r5+(Dg.*r25);
   r28=r21-(Du.*r1)-(D4.*r2);
   r29=r22+Ac'*r24+Ag'*r25;
   r30=r23-(D4.*r2)-(Dt.*r6);
   r31=r29+At'*r28-At'*(D9.*r30);
   rp=[r3;r26;r27];
   rptil=A'*(rp.*Dtil)-r31;

   %Direcoes

%dx= linsys([],rptil);
dx=L\(L'\rptil);
d=A*dx-rp;
da=d(1:mt);
dc=d(mt1:mtc)./D1;
dg=d(mtc1:m)./D2;
dt=-(r30+D4.*da)./D8;
dya=-D3.*da-D4.*dt-r28;
dsu=r1-da;
dsl=-r2+da+dt;
dsc=-Dc.*(D5.*dc+r24);
dsg=-Dg.*(D6.*dg+r25);
dst=r6-dt;
dyu=r12./su-Du.*dsu;
dyl=r13./sl-D4.*dsl;
dyc=(r14-yc.*dsc)./sc;
dyg=(r15-yg.*dsg)./sg;
dyt=r16./st-Dt.*dst;
dzx=r17./x-D7.*dx;
dzt=r18./t-Dz.*dt;
dwc=r19./c-D5.*dc;
dwg=r20./g-D6.*dg;

%Tamanho do passo primal

  alfap=min(1,-tau/min(-tau,min(dx./x)));
  alfap=min(alfap,-tau/min(-tau,min(dt./t)));
  alfap=min(alfap,-tau/min(-tau,min(dc./c)));
  alfap=min(alfap,-tau/min(-tau,min(dg./g)));
  alfap=min(alfap,-tau/min(-tau,min(dsu./su)));
  alfap=min(alfap,-tau/min(-tau,min(dsl./sl)));
  alfap=min(alfap,-tau/min(-tau,min(dsc./sc)));
  alfap=min(alfap,-tau/min(-tau,min(dsg./sg)));
  alfap=min(alfap,-tau/min(-tau,min(dst./st)));

    %Tamanho do passo dual

   alfad=min(1,-tau/min(-tau,min(dyu./yu)));
   alfad=min(alfad,-tau/min(-tau,min(dyl./yl)));
   alfad=min(alfad,-tau/min(-tau,min(dyc./yc)));
   alfad=min(alfad,-tau/min(-tau,min(dyg./yg)));
   alfad=min(alfad,-tau/min(-tau,min(dyt./yt)));
   alfad=min(alfad,-tau/min(-tau,min(dzx./zx)));
   alfad=min(alfad,-tau/min(-tau,min(dzt./zt)));
   alfad=min(alfad,-tau/min(-tau,min(dwc./wc)));
   alfad=min(alfad,-tau/min(-tau,min(dwg./wg)));

    %Proximo ponto primal

    a=a+alfap*da;
    x=x+alfap*dx;
    t=t+alfap*dt;
    c=c+alfap*dc;
    g=g+alfap*dg;
    su=su+alfap*dsu;
    sl=sl+alfap*dsl;
    sc=sc+alfap*dsc;
    sg=sg+alfap*dsg;
    st=st+alfap*dst;

    %Proximo ponto dual

    ya=ya+alfad*dya;
    yu=yu+alfad*dyu;
    yl=yl+alfad*dyl;
    yc=yc+alfad*dyc;
    yg=yg+alfad*dyg;
    yt=yt+alfad*dyt;
    zx=zx+alfad*dzx;
    zt=zt+alfad*dzt;
    wc=wc+alfad*dwc;
    wg=wg+alfad*dwg;

r1=ut-a-su;
r2=lt-a-t+sl;
r3=a-At*x;
r4=-Ac*x-sc+c;
r5=ug-Ag*x+g-sg;
r6=lt-t-st;
r7=yu-yl+ya;
r8=-At'*ya+Ac'*yc+Ag'*yg-zx;
r9=w-yl+yt-zt;
r10=1-yc-wc;
r11=1-yg-wg;

gama=(su'*yu + sl'*yl + sc'*yc + sg'*yg + st'*yt + x'*zx + t'*zt + c'*wc + g'*wg);
gamarel=gama/(1+abs(w*sum(t)+sum(c)+sum(g))+abs(-ut'*yu+lt'*yl-ug'*yg-lt'*yt))/np;

merro=max(gamarel, norm(r1/nut));
merro=max(merro, norm(r2/nlt));
merro=max(merro, norm(r3/(norm(a)+1)));
merro=max(merro, norm(r4/(norm(c)+1)));
merro=max(merro, norm(r5/nug));
merro=max(merro, norm(r6/nlt));
merro=max(merro, norm(r7/(norm(ya)+1)));
merro=max(merro, norm(r8/norm(zx)));
merro=max(merro, norm(r9/nmt));
merro=max(merro, norm(r10/nmc));
merro=max(merro, norm(r11/nmg));

k=k+1;

end
TotalTime = etime(clock,StartTime)
flops
gama
merro
k
w*sum(t)+sum(c)+sum(g)-sum(uc)

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
