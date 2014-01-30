function pr3001_B
    close all
    set(0,'DefaultFigureWindowStyle','docked')
    
    % Parameters
    lambda=2;
    lambda_vector=[1 2 3 4 5];
    mu=[0.5 1.5 2 2.5 4];
    x_min=[pi/2 pi/2 1.3181-(pi/2 - 1.3181)+0.02 -pi/2 -pi/2];
    x_max=[3*pi/2  (1.8235-pi/2)+1.8235 pi/2-0.02 pi/2 pi/2];
    x0=pi/2;
    w0=0;
    
    % Nombre de points d'équilibre
    displayNbEqui(lambda_vector, mu);
    
    % Portrait de phase & Ep(theta) & période
    T=zeros(1,length(mu));
    for i=1:length(mu)
        figure
        
        % portrait de phase
        subplot(2,1,1)
        portraitPhase(lambda, mu(i));
        w0max=vitesseInitMax(lambda, mu(i), x0);
        line([x0 x0], [-5 5])
        str=strcat('\leftarrow\omega_{0max} = ', num2str(w0max));
        text(x0, w0max, str, 'FontSize', 14)
        
        % Ep
        subplot(2,1,2)
        epPlot(lambda, mu(i));
        
        % période pour chaque mu
        T(i) = quad(@periode, x_min(i), x_max(i),[],[], lambda, mu(i));
    end
    
    % affichage de la période en fonction de mu
    T = real(T);
    figure
    periodePlot(T, mu);
    
end

%------------------------------------------------------------------------------
% Affiche le nombre de points d'équilibre en fonction de lambda et mu dans
% un tableau.
function displayNbEqui(lambda, mu)
    Z = nbEqui(lambda, mu);
    
    f = figure('name', 'Nombre de points d''équilibre(lambda, mu)', 'Position', [0 0 600 350]);
    t = uitable('Parent', f, 'Position', [50 700 500 150]);
    set(t, 'Data', Z, 'ColumnName', lambda, 'RowName', mu)
end

% Calcul le nombre de points d'équilibre pour chaque lambda et mu
% (vecteurs). Retourne une matrice (ligne: taille de lambda, colonne: taille de 
% mu).
function z=nbEqui(lambda, mu)
    z = zeros(length(lambda), length(mu));
    
    for i=lambda
        for j=1:length(mu)
            if (mu(j) < (lambda(i)/(lambda(i)-1)) + 1) && (mu(j) > (lambda(i)/(lambda(i)-1)) - 1)
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