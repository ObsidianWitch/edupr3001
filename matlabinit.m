function matlabinit

%Etude du pendule pesant a l'aide de Matlab


% 
%- 1 - Calcul de la solution de l'equation differentielle de 
%conditions initiales x0 pour b=1
%

%Cas 1: Oscillations: a=0
figure(1)                               
reptraj([0;1],30,0,2)

%Cas 2: Point en spiral attractif: a=0.2
figure(2)                               
reptraj([0;1],30,0.2,2)     

%Cas 3: Noeud attractif: a=2.2
figure(3)                               
reptraj([0;1],30,2.2,2)

%Cas 4 Spiral repulsif: a=-0.5
figure(4)                               
reptraj([0;5],30,0.2,2)

%
%- 2 - Calcul du portrait de phase (b=1)
%
%
%Methode 1
%
%Cas 1: a=0
figure(5)
[X,Y]=meshgrid(-5:0.01:5,-5:0.01:5);
Z=integprem(X,Y);
contour(X,Y,Z,[-0.9 -0.6 -0.3 0 0.25 0.5 0.75 1 1.5 2 2.5 3 4 5 6 7] )
title('Portrait de phase: a=0  b=1')
%
%Methode 2
%

%Cas 2: a=0.2
%Choix des valeurs initiales dans le carr� de centre O et de c�t� 2
[XX,YY]=meshgrid(-1:1:1,-1:1:1);
Z=[XX(:) YY(:)];
figure(6)
portloc(Z,30,0.2,1)
title('Portrait de phase: a=0.2  b=1')

%Cas 3: a=2.2
figure(7)
%Choix des valeurs initiales dans le cercle de centre O et de rayon r=1
r=1;
t=0:0.2*pi:2*pi;
Z=[r*cos(t') r*sin(t')];
portloc(Z,30,2.2,1)
title('Portrait de phase: a=2.2  b=1')

%
%- 4 - Diagramme de bifurcation
%

figure(8)
a=0:0.03:3;
b=0:0.03:3;
[X,Y]=meshgrid(a,b);
X=X(:);
Y=Y(:);
D=X.^2-4*Y;
S=-X;
P=Y;
plot(X(D<0 & S>0),Y(D<0 & S>0),'.c')
hold on
plot(X(D<0 & S<0),Y(D<0 & S<0),'.m')
plot(X(D>0 & P>0 & S>0),Y(D>0 & P>0 & S>0),'.b')
plot(X(D>0 & P>0 & S<0),Y(D>0 & P>0 & S<0),'.g')
plot(X(D>0 & P<0 ),Y(D>0 & P<0 ),'.r')
hold off
title('Diagramme de bifurcation')
xlabel('a')
ylabel('b')

%
%- 4 - Calcul de la periode (a=0;b=1)
%
X0=0.1:0.1:3.1;
T=zeros(1,length(X0));
for i=1:length(X0)
T(i)=quad(@periode,10^-6,X0(i)-10^-6,[],[],X0(i));
end

figure(9)
plot(X0,T,'b')
xlabel('x0')
ylabel('T')
axis([0 3.2 0 25])
title('Periode/Amplitude')

%
%- 5 - Calcul du nombre de tours avant arr�t en fonction de v0 et a
%

v0=0.1:0.1:10;
a=0.1:0.1:3;
nt=zeros(length(v0),length(a));
for i=1:length(v0)
for j=1:length(a)
   [t,x]=ode45(@pend,[0 30],[0 v0(i)],[],a(j),1);
    if x(length(t),2)>=2 
    nt(i,j)=NaN;
    else 
    nt(i,j)=floor((x(length(t),1)+pi)/(2*pi));
    end
end
end

figure(10)
mesh(a,v0,nt)
xlabel('a')
ylabel('v0')
zlabel('nt')





%-------------------------------------------------------------
%Equation du pendule
%-------------------------------------------------------------
function dxdt=pend(~,x,a,b)
dxdt=[x(2);-a*x(2)-b*sin(x(1))];

%-------------------------------------------------------------
%Representation d'une trajectoire
%--------------------------------------------------------
function reptraj(x0,tf,a,b)
[t,x]=ode45(@pend,[0 tf],x0,[],a,b);

subplot(2,1,1)
plot(t,x)
title(['Parametres: a=' num2str(a) '   b=' num2str(b) '  ;  Conditions initiales: x_0=' num2str(x0(1)) '  v_0=' num2str(x0(2)) ])
xlabel('t')
grid
axis([0 30 -2 2])
subplot(2,1,2)
plot(x(:,1),x(:,2))
xlabel('x(t)')
ylabel('v(t)')


%-------------------------------------------------------------
% Portrait de phase: representation locale au voisinage 
%d'un point d'�quilibre 
%--------------------------------------------------------
function portloc(Z,tf,a,b)
for i=1:size(Z,1)
    [~,x]=ode45(@pend,[0 tf],Z(i,:),[],a,b);
    plot(x(:,1),x(:,2))
    hold on % TODO hold on en dehors de for
    xlabel('position')
    ylabel('vitesse')
    title('Portrait de phase')
end
hold off

%-------------------------------------------------------------
%Integrale premiere
%--------------------------------------------------------
function z=integprem(x,y)
z=0.5*y.^2-cos(x);

%--------------------------------------------
%Fonction a integrer pour la periode 
%--------------------------------------------
function y=periode(x,x0)
y=2*sqrt(2)./sqrt(cos(x)-cos(x0));





