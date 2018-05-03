function [Uw, Dw] = entrainmentWakeModel(x, Ct, E)
%
% [Uw, Dw] = entrainmentWakeModel(x, Ct, E)
%
% where:
% x is the distance downstream of the rotor, normalized by rotor diameter;
% Ct is the thrust coefficient;
% E is the entrainment coefficient (if not supplied, E = 0.16 is assumed);
% Uw is the velocity in the wake, normalized by upstream velocity;
% Dw is the diameter of the wake, normalized by rotor diameter.
%
% To produce sample output and plots, run without any input.
%
% Paolo Luzzatto-Fegiz, May 2018

if ~exist('Ct', 'var')
    Ct = .8;
end
if ~exist('E', 'var')
    E = 0.16;
end
if ~exist('x', 'var')
    x = linspace(0,20,20);
end


X = 6*E*(2/Ct)^0.5 * x + (1-Ct)^(3/4) / (1-(1-Ct)^(1/2))^(3/2);

Uw = X.^(2/3) ./ (X.^(2/3) + 1) ;

Dw = sqrt(Ct/2) * (X.^(2/3) + 1 )./ X.^(1/3);


if nargout == 0
    figure
    subplot(2,1,1)
    plot(x, Uw,'k');hold on
    ylabel('U_w/U_o')
    subplot(2,1,2)
    plot(x, Dw, 'r')
    ylabel('D_w/D')
    xlabel('x/D')
    title('Entrainment model')
    
    % compare momentum conservation to park model
    figure
    k = 0.06; % spreading paramer
    a = 0.5*(1-sqrt(1-Ct));
    Di = sqrt((1 - a)/(1 - 2*a));
    Dpark = Di + 2*k*x;
    Upark = 1 - 2*a./(1+2*k*x./Di).^2;
    thrust = .5*Ct*pi/4; % nondimensionalized by upstream vel., rotor diam.
    momentumDeficitEntrainment = pi/4*Dw.^2 .* Uw.*(1-Uw);
    momentumDeficitPark = pi/4*Dpark.^2 .* Upark.*(1-Upark);
    plot(x,x*0+thrust)
    hold on
    plot(x,momentumDeficitEntrainment,'r.:')
    plot(x,momentumDeficitPark)
    legend('Thrust','Wake momentum deficit, entrainment model',...
        'Wake momentum deficit, park model','location','nw')
    xlabel('x/D')
end
