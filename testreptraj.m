function testreptraj()
    close all
    % trajectoire ([x0;w0], temps final, amortissement, mu)
    reptraj([0;2],30,0.5,2)
end

%-------------------------------------------------------------
%Equation du pendule
%-------------------------------------------------------------
function dxdt=pend(~,x,a,mu)
dxdt=[x(2);-a*x(2)-sin(x(1))*(1 - (2/sqrt(mu^2 + 1 -2*mu*cos(x(1)))))];
end

%-------------------------------------------------------------
%Representation d'une trajectoire
%--------------------------------------------------------
function reptraj(x0,tf,a,mu)
[t,x]=ode45(@pend,[0 tf],x0,[],a,mu);

subplot(2,1,1)
% bleu -> x(t) ; vert -> v(t)
plot(t,x)
title(['Parametres: a=' num2str(a) '   mu=' num2str(mu) '  ;  Conditions initiales: x_0=' num2str(x0(1)) '  w_0=' num2str(x0(2)) ])
xlabel('t')
grid
axis([0 30 -10 10])
subplot(2,1,2)
plot(x(:,1),x(:,2))
xlabel('x(t)')
ylabel('w(t)')
end