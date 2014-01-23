function pr3001_B
    close all
    set(0,'DefaultFigureWindowStyle','docked')
    
    % Parameters
    lambda=2;
    lambda_vector=[1 2 3 4 5];
    mu=[0.5 1.5 2 2.5 4];
    x0=pi/2;
    %w0=0;
    
    % Nombre de points d'équilibre
    displayNbEqui(lambda_vector, mu);
    
    % Portrait de phase & Ep(theta)
    for i=mu
        figure
        
        subplot(2,1,1)
        portraitPhase(lambda, i);
        w0max=vitesseInitMax(lambda, i, x0)
        text(x0, w0max, '\leftarrow\omega_{0max}', 'FontSize', 16)
        
        subplot(2,1,2)
        epPlot(lambda, i);
        
        %plot_h(lambda, i); % permet de déterminer les valeurs max et min de h(x)
    end
    
    %vitesseInitMax(lambda, mu, x0);
end

%------------------------------------------------------------------------------
% Affiche le nombre de points d'équilibres en fonction de lambda et mu dans
% un tableau.
function displayNbEqui(lambda, mu)
    Z = nbEqui(lambda, mu);
    
    f = figure('name', 'Nombre de points d''équilibre(lambda, mu)', 'Position', [0 0 600 350]);
    t = uitable('Parent', f, 'Position', [100 100  700 200]);
    set(t, 'Data', Z, 'ColumnName', lambda, 'RowName', mu)
end

% TODO? autre méthode: résoudre h(theta) = 0 et compter le nombre de
% solutions
%
% Calcul le nombre de points d'équilibres pour chaque lambda et mu
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
% function z=h(x, lambda, mu)
%     z=sin(x).*(-1 + lambda -(lambda./sqrt(mu^2 + 1 - 2*mu*cos(x))));
% end
% 
% function plot_h(lambda, mu)
%     figure
%     x=-10:0.1:10;
%     y=h(x, lambda, mu);
%     plot(x, y);Z
%     title(['h(x) mu=' num2str(mu)]);
% end
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
% Portrait de phase
function portraitPhase(lambda, mu)
    [X,Y]=meshgrid(-5:0.01:5, -5:0.01:5);
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

function w0max=vitesseInitMax(lambda, mu, x0)
    H0=H_IntegPrem(lambda, mu, 0);
    Hpi=H_IntegPrem(lambda, mu, pi);
    if H0 > Hpi
        C=H0;
    else
        C=Hpi;
    end
    
    w0max=vitesseAngulaire(lambda, mu, C, x0)
end


%------------------------------------------------------------------------------