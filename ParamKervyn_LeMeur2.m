%% PARTIE 7.3 - CAS 1
% Vague avec profondeur de l'océan variable.
% question a.
clear all; close all; clc;

% Parametres physiques:
tfin = 101;    % temps physique de la simulation (s)
xL   = 0.0 ; 		     xR   = 5.e3 ; 		     
yL   = 0.0 ;		     yU   = 2.e3;		%dimensions (m)     
pert_amplitude = 0.e0;	 pert_velocity  = 0.e0;	

% Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

% Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 1.e3;      % profondeur maximale du fond marin (m)
h1  = 2.e1;      % profondeur minimale du fond marin (m)
a   = 2000.0; % limite inferieur de l'onde en x (m)
b   = 5000.0; % limite superieur de l'onde en y (m)
Ly  = 2000.0; % longueur characteristique de l'onde en y (m)

% Parametres numeriques :
% nombre de node du maillage en x,y (extreme inclus)
Nx = 128; 			 Ny = 128;

ComputeDt = true;  %ajout PERSO - option pour le calcul du dt
dt  = 0.1;         %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
CFL = 0.95;        %valeur de la condition CFL

type_u2  = 'onde_cas1';	 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;         % if true f est ecrite

% valeur propre en x,y
mode_num_x = 0;			 mode_num_y = 0;			 

% conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'harmonic'; bc_right  = 'neumann';   
bc_lower  = 'neumann';  bc_upper  = 'neumann';   

T = 15; %2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));

impulsion = true;          % si vrais et harmonique une seule onde est emise
type_init = 'homogene';     % type initialisation: harmonic, default: homogene
F0 = 0.e0;			        % amplitude de l'harmonique initiale
A = 1.e0;			        % amplitude condition bord harmonique
omega = 2*pi/T;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrais le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrais la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

% Simulations
filename2  = "Nx_"+ num2str(Nx) + "Ny_" + num2str(Ny);
    Nx_loc = Nx; Ny_loc = Ny; 
    writeConfig;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 

ViewFormat;
mesh = dlmread(filename2+'_mesh.out');
data = load(filename2+'_f.out');
[X,Y]= meshgrid(mesh(1,1:Nx),mesh(2,1:Ny));

% Calcul du profil de la profondeur
hK1 = h0*ones(length(X),length(Y));
for j = 1:length(Y)
    for i=1:length(X)
        x = mesh(1,i);
        if(x>a && x<b) 
            hK1(j,i) = (h0-h1)*x/(a-b)+(a*h1 - b*h0)/(a-b);%h0+(h1-h0)*sin(pi*(x-a)/(b-a));
        end
    end
end

%% Plot du profil de la profondeur
figure('Name','plotProfondeur Cas 1')
for i =1:length(data(:,1))/Ny
    surf(X,Y,-hK1);
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [m]')
%     zlim([-5 5]);
%     view(0,0);
end
SaveIMG("fondBonus");

%% Plot de la simulation au cours du temps
Hmax = 0; delta = 0.140245/2; 
figure('Name','plotsimulation')
for i =1:length(data(:,1))/Ny
%     surf(X,Y,-hK1,'EdgeColor','None');
%     hold on
    surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1),'FaceColor','Interp');
%     hold off
    title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    Max_temp = max(data(Ny*(i-1)+1:Ny*i,2:Nx+1),[],'all'); %returns the largest element of X.
    if(Max_temp > Hmax)
       Hmax = Max_temp; 
    end
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [m]')
    zlim([-2, 2]);  
    view(0,0);
    if((data(Ny*i,1) > 20-delta) && (data(Ny*i,1) < 20+delta))
        SaveIMG("BONUS plot 1");
    end
    if((data(Ny*i,1) > 50-delta) && (data(Ny*i,1) < 50+delta))
        SaveIMG("BONUS2 plot 2");
    end
    pause(1.e-9);
end
disp("Hmax = " + num2str(Hmax));

%% Nouveau maillage en Nx, Amplitude sur une coupe
% Vague avec profondeur de l'océan variable.
% question a.
clear all; clc;

% Parametres physiques:
tfin = 100.0;    % temps physique de la simulation (s)
xL   = 0.0 ; 		     xR   = 5.e3 ; 		     
yL   = 0.0 ;		     yU   = 2.e3;		%dimensions (m)     
pert_amplitude = 0.e0;	 pert_velocity  = 0.e0;	

% Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

% Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 1.e3;      % profondeur maximale du fond marin (m)
h1  = 2.e1;      % profondeur minimale du fond marin (m)
a   = 2000.0; % limite inferieur de l'onde en x (m)
b   = 5000.0; % limite superieur de l'onde en y (m)
Ly  = 2000.0; % longueur characteristique de l'onde en y (m)

% Parametres numeriques :
% nombre de node du maillage en x,y (extreme inclus)
Nx = [100 200 600]; 			 Ny = 4;

ComputeDt = true;  %ajout PERSO - option pour le calcul du dt
dt  = 0.1;         %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
CFL = 0.95;        %valeur de la condition CFL

type_u2  = 'onde_cas1';	 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;         % if true f est ecrite

% valeur propre en x,y
mode_num_x = 0;			 mode_num_y = 0;			 

% conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'harmonic'; bc_right  = 'neumann';   
bc_lower  = 'neumann';  bc_upper  = 'neumann';   

T = 15; %2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));

impulsion = true;           % si vrais et harmonique une seule onde est emise
type_init = 'homogene';     % type initialisation: harmonic, default: homogene
F0 = 0.e0;			        % amplitude de l'harmonique initiale
A = 1.e0;			        % amplitude condition bord harmonique
omega = 2*pi/T;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrais le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrais la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

% Simulations
for n=1:length(Nx)
    filename2  = "Nx_"+ num2str(Nx(n)) + "Ny_" + num2str(Ny);
    Nx_loc = Nx(n); Ny_loc = Ny; 
    writeConfig;
    disp('Ex7 configuration.in');   
    system('Ex7 configuration.in'); 

    ViewFormat;
    mesh = dlmread(filename2+'_mesh.out');
    data = load(filename2+'_f.out');

    % on récupère toutes les tranches à chaque temps pour y fixé (starty)
    starty   = 2; newdata  = data(starty:Ny:end,:);  
    % on trouve le maximum d'amplitude et l'indice en x associé
    AmpliMax = zeros(length(newdata(:,1)),1); IndiceX = zeros(length(newdata(:,1)),1);
    for i=1:length(AmpliMax) 
       [AmpliMax(i) , IndiceX(i)] = max(newdata(i,2:end)); 
    end
    % on fait correspondre l'indice à une valeur de x
    Xreel = zeros(length(IndiceX),1); 
    for i=1:length(IndiceX) 
        Xreel(i,1) = mesh(1,IndiceX(i)); 
    end
    if(n==1)
       A1 = [Xreel, AmpliMax];
    end
    if(n==2)
       A2 = [Xreel, AmpliMax];
    end
    if(n==3)
       A3 = [Xreel, AmpliMax];
    end
end
posy = mesh(2,starty);

U = load(filename2+'_u.out');
newU = U(starty,:);

A0 = (g*h0)^(1/4);

debut = 60; fin = 1234;
figure('Name',"plot Amplitude")
        plot(A1(:,1),A1(:,2),' .','Linewidth',lw);
        hold on
        plot(A2(6:end,1),A2(6:end,2),' .','Linewidth',lw);
        plot(A3(debut:fin,1),A3(debut:fin,2),' .','Linewidth',lw);
        plot(mesh(1,1:Nx(end)),A0./sqrt(newU),'--','Linewidth',lw+1);
        grid minor; set(gca,'fontsize',fs);
        legendStrings = "$n_x = $" + string(Nx');
        legend("$n_x = $" + num2str(Nx(1)),"$n_x = $" + num2str(Nx(2)),"$n_x = $" + num2str(Nx(3)),"$A_0 / \sqrt{u}$",'Location','northwest','NumColumns',1)
        xlabel('$x$ [m]'); ylabel('Relative amplitude [m]');
        ylim([0.8 2.7]); xlim([10,5000]);
SaveIMG("AmplitudeRelative");

%% Plot de la simulation au cours du temps (pour le dernier Nx)
Hmax = 0; delta = 0.140245/2; 

for i =1:length(data(:,1))/Ny
    Max_temp = max(data(Ny*(i-1)+1:Ny*i,2:Nx(n)+1),[],'all'); %returns the largest element of X.
    if(Max_temp > Hmax)
       Hmax = Max_temp; 
    end
    pause(1.e-9);
end
disp("Hmax = " + num2str(Hmax));

%% évolution temporelle en coupe selon x
start = 2;
newdata = data(start:Ny:end,:);
posy = mesh(2,start);

% en Surf
figure('Name',"plot cas 1 en coupe")
        surf(mesh(1,:),newdata(:,1),newdata(:,2:end),'EdgeColor','None');
        grid minor; set(gca,'fontsize',fs);
        xlabel('$x$ [m]'); ylabel('$t$ [s]'); zlabel("$f(x,y=y_{loc},t)$ [m]")
SaveIMG("7_3_CoupeCas1");

% en Pcolor
axes1 = axes('Parent',figure);
for i =1:length(data(:,1))/Ny
    s1 = pcolor(mesh(1,:),newdata(:,1),newdata(:,2:end));
    s1.FaceColor = 'interp';
    s1.EdgeColor = 'None';
    xlabel('$x$ [m]'); ylabel('$t$ [s]'); 
    set(gca,'fontsize',fs)
        h = colorbar(axes1);
        h.Label.String = '$f(x,y=y_{loc},t)$ [m]';
        h.Label.Interpreter = 'Latex';
        h.Label.FontSize = fs;
end
SaveIMG("7_3_PcolorCas1");