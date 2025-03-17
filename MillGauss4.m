
% inizializzo campo e parametri
clear; clc; close all;

k = 1; % costante del campo conservativo, arbitraria
soglia_distanza = 10; % soglia per il warning di precisione, a causa dell'elevato utilizzo 
% nel codice dell' integrazione discreta attraverso la funzione trapz reputo
% corretto evidenziarne le criticità. 
% 
% P.s.  i valori ricavati sono formattati nella console 
%  a 6 cifre significative, mentre nel grafico a sole 4, per una
% migliore comprensione del concetto su cui si basa il programma

% disegno le superfici equipotenziali, sulle quali l'energia
% potenziale si mantiene
figure('Name','Campo e percorsi');
axis equal; hold on;
title('Campo conservativo e percorsi');
xlabel('X'); ylabel('Y');

% genero 10 superfici equipotenziali circolari
theta = linspace(0, 2*pi, 100);
for r = linspace(1, 10, 10)
    plot(r*cos(theta), r*sin(theta), 'k:', 'LineWidth', 0.5);
end

% input delle coordinate dei due punti A e B
fprintf('\n=== INSERIMENTO PUNTI ===\n');
ax = input('Coordinata X di A: ');
ay = input('Coordinata Y di A: ');
bx = input('Coordinata X di B: ');
by = input('Coordinata Y di B: ');

% verifico se sono superate le condizioni critiche
distanza = norm([bx-ax, by-ay]);
if any([ax, ay, bx, by] == 0) || distanza > soglia_distanza
    warning(['Attenzione: condizioni critiche rilevate!\n'...
        'Distanza A-B: %.2f (soglia: %.2f)\n'...
        'Coordinate nulle: %s'],...
        distanza, soglia_distanza, mat2str(any([ax,ay,bx,by]==0)));
end

%genero i percorsi
t = linspace(0, 1, 1000).'; % Parametro normalizzato, spezzandolo in 1000 pezzettini uguali, utile per fare successivamente gli integrali

%percorso rettilineo, anche in questo caso spezzo il percorso in 1000
%pezzetti
x_retto = linspace(ax, bx, 1000).'; 
y_retto = linspace(ay, by, 1000).';

% percorso semicircolare 2
vettore_AB = [bx-ax, by-ay];
centro = mean([[ax,ay]; [bx,by]]);
raggio = distanza/2;
angolo_iniziale = atan2(vettore_AB(2), vettore_AB(1)) + pi;

theta_semi = linspace(angolo_iniziale, angolo_iniziale + pi, 1000).';
x_semi = centro(1) + raggio*cos(theta_semi);
y_semi = centro(2) + raggio*sin(theta_semi);

%percorso a spirale 3 (ironicamente questa è quasi la parte più complessa
%di tutto il programma ma volevo un percorso il più irregolare possibile),
%per tale motivo mi permetto di giustificarlo passo passo

% conversione a coordinate polari
[thetaA, rA] = cart2pol(ax, ay); % angolo e raggio di A
[thetaB, rB] = cart2pol(bx, by); %angolo e raggio di B


% regolo l'angolo per garantire un rotazione antioraria
delta_theta = thetaB - thetaA;
if delta_theta < 0
    delta_theta = delta_theta + 2*pi; % correggo il caso di angoli negativi
end
% creiamo una spirale con almeno un giro completo
theta_spirale = thetaA + t*(delta_theta + 2*pi); % forzo +2pi radianti
r_spirale = rA + t*(rB - rA);% interpolo linearmente il raggio

% ritorno a coordinate cartesiane
[x_spirale, y_spirale] = pol2cart(theta_spirale, r_spirale);

%visualizzo i percorsi
hA = plot(ax, ay, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hB = plot(bx, by, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
hRetto = plot(x_retto, y_retto, 'y', 'LineWidth', 2);
hSemi = plot(x_semi, y_semi, 'Color', [1 0.5 0], 'LineWidth', 2);
hSpirale = plot(x_spirale, y_spirale, 'r', 'LineWidth', 2);

legend([hA, hB, hRetto, hSemi, hSpirale],...
    {'Punto A', 'Punto B', 'Retto', 'Semicerchio', 'Spirale'},...
    'Location', 'northwest');

% calcolo i lavori
F = @(x,y) -k*[x(:), y(:)]./(x(:).^2 + y(:).^2).^(3/2);

% faccio l''integrazione integrazione attraverso il metodo dei trapezi
calcola_lavoro = @(x,y) trapz(t, sum(F(x,y).*[gradient(x,t), gradient(y,t)],2)); %calcolo il lavoto scalare, utilizzo la regola del trapezio

% calcolo il lavoro per tutti i percorsi
W_retto = calcola_lavoro(x_retto, y_retto); 
W_semi = calcola_lavoro(x_semi, y_semi);
W_spirale = calcola_lavoro(x_spirale, y_spirale);

% calcolo i lavori cumulativi
Z_retto = cumtrapz(t, sum(F(x_retto,y_retto).*[gradient(x_retto,t), gradient(y_retto,t)],2));%considero quanto lavoro svolgo di istante in istante sempre grazie alla regola dei trapezi e lo sommo di volta in volta 
%in questo modo posso vedere se la somma dei singoli contributi del lavoro
%aumentano o diminusicono, indicandomi quindi un lavoro negativo o positivo, in termini tecnici quello che ottengo è il lavoro cumulativo istantaneo <-- molto
%importante per la piena comprensione di cosa mostri il grafico.

Z_semi = cumtrapz(t, sum(F(x_semi,y_semi).*[gradient(x_semi,t), gradient(y_semi,t)],2));
Z_spirale = cumtrapz(t, sum(F(x_spirale,y_spirale).*[gradient(x_spirale,t), gradient(y_spirale,t)],2));


%% visualizzo i risultati, parte di elaboraazione numerica conlcusa.
% mostro il grafico del lavoro cumulativo
figure('Name','Lavoro cumulativo');
hold on;

% definisco l'area semi-trasparente
fill([t; flipud(t)], [Z_retto; zeros(size(Z_retto))], 'y',...
    'FaceAlpha',0.2, 'EdgeColor','none');
fill([t; flipud(t)], [Z_semi; zeros(size(Z_semi))], [1 0.5 0],...
    'FaceAlpha',0.2, 'EdgeColor','none');
fill([t; flipud(t)], [Z_spirale; zeros(size(Z_spirale))], 'r',...
    'FaceAlpha',0.2, 'EdgeColor','none');

% le linee principali
plot(t, Z_retto, 'y', 'LineWidth', 2);
plot(t, Z_semi, 'Color', [1 0.5 0], 'LineWidth', 2);
plot(t, Z_spirale, 'r', 'LineWidth', 2);
plot([0 1], [0 0], 'k--', 'LineWidth', 1);

title('Andamento del lavoro cumulativo');
xlabel('Parametro percorso normalizzato (t)');
ylabel('Lavoro compiuto (J)');
legend({'Retto','Semicerchio','Spirale','Zero'},...
    'Location','northwest', 'Box','off');
grid on;

% grafico 3D con i percorsi e il lavoro
figure('Name','Visualizzzazione 3D');
hold on; view(3);

% superfici trasparenti
surface([x_retto, x_retto], [y_retto, y_retto], [zeros(size(Z_retto)), Z_retto],...
    'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.3);
surface([x_semi, x_semi], [y_semi, y_semi], [zeros(size(Z_semi)), Z_semi],...
    'FaceColor',[1 0.5 0], 'EdgeColor','none', 'FaceAlpha',0.3);
surface([x_spirale, x_spirale], [y_spirale, y_spirale], [zeros(size(Z_spirale)), Z_spirale],...
    'FaceColor','r', 'EdgeColor','none', 'FaceAlpha',0.3);

% linee 3D
plot3(x_retto, y_retto, Z_retto, 'y', 'LineWidth', 2);
plot3(x_semi, y_semi, Z_semi, 'Color',[1 0.5 0], 'LineWidth', 2);
plot3(x_spirale, y_spirale, Z_spirale, 'r', 'LineWidth', 2);

% formattazione
light('Style','infinite', 'Position',[1 1 1]);
lighting gouraud;
title(sprintf('Lavoro totale\nRetto: %.4f J\nSemicerchio: %.4f J\nSpirale: %.4f J',...
    W_retto, W_semi, W_spirale));
xlabel('X'); ylabel('Y'); zlabel('Lavoro cumulativo');
grid on;
rotate3d on;

%mostro l'output dei risultati
fprintf('\n=== RISULTATI ===\n');
fprintf('Lavoro percorso rettilineo: %.6f J\n', W_retto);
fprintf('Lavoro percorso semicircolare: %.6f J\n', W_semi);
fprintf('Lavoro percorso a spirale: %.6f J\n', W_spirale);