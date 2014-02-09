function pr3001_B
    close all
    set(0,'DefaultFigureWindowStyle','docked')
    
    % Parameters
    lambda=2;
    lambda_vector=[2 3 4 5];
    mu=[0.5 1.5 2 2.5 4];
    x_min=[pi/2 pi/2 1.3181-(pi/2 - 1.3181)+0.02 -pi/2 -pi/2];
    x_max=[3*pi/2  (1.8235-pi/2)+1.8235 pi/2-0.02 pi/2 pi/2];
    x0=pi/2;
    w0=0;
    
%     % Nombre de points d'équilibre & représentation des points d'équilibre
%     displayNbEqui(lambda_vector, mu);
%     figure
%     plotEqui(lambda);
%     
%     % Portrait de phase & Ep(theta) & période
%     T=zeros(1,length(mu));
%     for i=1:length(mu)
%         figure
%         
%         % portrait de phase
%         subplot(2,1,1)
%         portraitPhase(lambda, mu(i));
%         w0max=vitesseInitMax(lambda, mu(i), x0);
%         line([x0 x0], [-5 5])
%         str=strcat('\leftarrow\omega_{0max} = ', num2str(w0max));
%         text(x0, w0max, str, 'FontSize', 14)
%         
%         % Ep
%         subplot(2,1,2)
%         epPlot(lambda, mu(i));
%         
%         % période pour chaque mu
%         T(i) = quad(@periode, x_min(i), x_max(i),[],[], lambda, mu(i));
%     end
%     
%     % affichage de la période en fonction de mu
%     T = real(T);
%     figure
%     periodePlot(T, mu);
%     
    % diagramme de bifurcation
    mu2=[0 0.5 1.5 2 2.5 3 4];
    figure
    plotBifurcation(mu2);
    
end

%------------------------------------------------------------------------------
% Affiche le nombre de points d'équilibre en fonction de lambda et mu dans
% un tableau.
function displayNbEqui(lambda, mu)
    Z = nbEqui(lambda, mu);
    
    f = figure('name', 'Nombre de points d''équilibre(lambda, mu)', 'Position', [0 0 600 350]);
    t = uitable('Parent', f, 'Position', [50 20 500 150]);
    set(t, 'Data', Z, 'ColumnName', lambda, 'RowName', mu)
end

% Calcul le nombre de points d'équilibre pour chaque lambda et mu
% (vecteurs). Retourne une matrice (ligne: taille de lambda, colonne: taille de 
% mu).
function z=nbEqui(lambda, mu)
    z = zeros(length(mu), length(lambda));
    
    for i=1:length(mu)
        for j=1:length(lambda)
            if (mu(i) < (lambda(j)/(lambda(j)-1)) + 1) && (mu(i) > (lambda(j)/(lambda(j)-1)) - 1)
                z(i,j) = 4;
            else
                z(i,j) = 2;
            end
        end
    end
end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Portrait de phase
function portraitPhase(lambda, mu)
    [X,Y]=meshgrid(-5:0.01:5, -3:0.01:3);
    Z=integPrem(lambda, mu, X, Y);
    contour(X, Y, Z, -5:0.1:5)
    title(['Portrait de phase lambda=', num2str(lambda), ' mu=', num2str(mu)]);
end

% Affiche le terme H(theta) de l'intégrale première (proportionnel à
% l'énergie potentielle)
function epPlot(lambda, mu)
    x=(-5:0.01:5);
    z=H_IntegPrem(lambda, mu, x);
    plot(x,z);
    title(['EP lambda=', num2str(lambda), ' mu=', num2str(mu)]);
end

% Intégrale première
function z=integPrem(lambda, mu, x, y)
    z=0.5*y.^2 + H_IntegPrem(lambda, mu, x);
end

% Terme H(theta) de l'intégrale première
function z=H_IntegPrem(lambda, mu, x)
    z=cos(x) + 0.5*(lambda/mu)*(sqrt(mu.^2 + 1 -2*mu*cos(x)) - 1).^2;
end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function z=vitesseAngulaire(lambda, mu, C, x)
    z=sqrt(2*(C - H_IntegPrem(lambda, mu, x)));
end

% Calcul de la valeur max de la vitesse initiale w0 pour laquelle la
% trajectoire est périodique, sachant que l'on fixe la position initiale
% x0, lambda et mu.
function w0max=vitesseInitMax(lambda, mu, x0)
    H0=H_IntegPrem(lambda, mu, 0);
    Hpi=H_IntegPrem(lambda, mu, pi);
    if H0 > Hpi
        C=H0;
    else
        C=Hpi;
    end
    
    w0max=vitesseAngulaire(lambda, mu, C, x0);
end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Affiche la courbe représentant la période T (vecteur) en fonction de mu (vecteur).
function periodePlot(T, mu)
    plot(mu,T);
    title('T(mu)');
    
    for i=1:length(mu)
        str=strcat('\leftarrow ', num2str(T(i)));
        text(mu(i), T(i), str, 'FontSize', 14)
    end
end

function y=periode(x, lambda, mu)
    C= 0.5 * (lambda/mu) * (sqrt(mu^2 + 1) - 1)^2;
    y=sqrt(2)./sqrt(C - H_IntegPrem(lambda, mu, x));
end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
function plotEqui(lambda)
    x=(0:0.01:5);
    z= real(equiArccos(lambda, x));
    
    hold on
    plot(x, pi,'r');
    plot(x, 0,'r');
    plot(x, -pi,'r');
    plot(x,z,'b');
    plot(x,-z,'b');
    hold off
    title(['Représentation des points d''équilibre en fonction de mu (lambda=', num2str(lambda),')'])
    xlabel('mu')
    ylabel('position')
end

function y=equiArccos(lambda, mu)
    y = acos(((lambda/(lambda - 1))^2 - mu.^2 -1)./(-2*mu));
end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Diagramme de bifurcation
function plotBifurcation(mu)
    lambda = 2;
    alphaFricArray = zeros(1, length(mu));
    
    for i=1:length(mu)
        mu(i)
        ptEqui = [0 pi];

        % points d'équilibres du système
        if (mu(i) < (lambda/(lambda-1)) + 1) && (mu(i) > (lambda/(lambda-1)) - 1) % 4 points
            ptEqui(end + 1) = equiArccos(lambda, mu(i));
        end

        % on choisit un point d'équilibre tel que le système soit stable (det(A) > 0 
        % sachant que tr(A) < 0)
        for j=1:length(ptEqui)
            detPtEqui = det(ptEqui(j), mu(i))
            if detPtEqui >= 0
                ptEquiChoisi = ptEqui(j);
                break;
            end
        end

        alphaFricArray(i) = alphaFriction(ptEquiChoisi, mu(i));
    end
    
    plot(alphaFricArray, mu);
    title('Diagramme de bifurcation')
    xlabel('alpha')
    ylabel('mu')
end

% Calcule le déterminant de la matrice jacobienne pour lambda=2 et en
% fonction de lambda et x1 (theta).
function y=det(x1, mu)
    y = cos(x1)*(1 - (2/sqrt(mu^2 + 1 -2*mu*cos(x1)))) ...
        + 2*mu*((sin(x1)^2)/((mu^2 + 1 -2*mu*cos(x1))^(3/2)));
end

% Calcule alpha en fonction de x1 et mu
function y=alphaFriction(x1, mu)
    y = 2*sqrt(det(x1,mu));
end
%------------------------------------------------------------------------------