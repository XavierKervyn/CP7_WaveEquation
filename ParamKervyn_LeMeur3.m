%% PARTIE 7.3 - CAS 2
% Vague avec profondeur de l'océan variable.
% question a.
clear all; close all; clc;

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
Nx = 256; 			 Ny = 256;

ComputeDt = true;  %ajout PERSO - option pour le calcul du dt
dt  = 0.1;        %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
CFL = 0.95;         % valeur de la condition CFL

type_u2  = 'onde_cas2';	 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;         % if true f est ecrite

% valeur propre en x,y
mode_num_x = 0;			 mode_num_y = 0;			 

%conditions aux bords: dirichlet, neumann, harmonic
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
filename2  = "Nx_"+ num2str(Nx) + "Ny_" + num2str(Ny);
    Nx_loc = Nx; Ny_loc = Ny; 
    writeConfig;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 

ViewFormat;
data = load(filename2+'_f.out');
mesh = load(filename2+'_mesh.out');
[X,Y]= meshgrid(mesh(1,:),mesh(2,:));

% profil de la profondeur
hK2 = h0*ones(length(X),length(Y));
for j = 1:length(Y)
    for i=1:length(X)
        x = mesh(1,i); %disp(x);
        y = mesh(2,j);
        if(x>a && x<b) 
            hK2(j,i) = h0+(h1-h0)*sin(pi*(x-a)/(b-a))*sin(pi*y/(yU-yL));
        end
    end
end

%% Plot du profil de la profondeur
figure('Name','plotProfondeur Cas 2')
    surf(X,Y,-hK2,'EdgeColor','None');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$-h$ [m]')
    xlim([xL,xR]);
    view(-25,30);
SaveIMG("FondCAS2");

%% Plot de la simulation au cours du temps
close all; clc
Hmax = 0; delta = 0.140245/2;
figure('Name','plotsimulation cas 2')
for i =1:length(data(:,1))/Ny
    contourf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1),10);
%     surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1),'EdgeColor','None');
    title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
    grid minor; set(gca,'fontsize',fs);
    h = colorbar;
    h.Label.String = "$f$ [m]";
    h.Label.Interpreter = 'Latex';
    h.Label.FontSize = fs;
%     caxis([0,9])
    Max_temp = max(data(Ny*(i-1)+1:Ny*i,2:Nx+1),[],'all'); %returns the largest element of X.
    if(Max_temp > Hmax)
       Hmax = Max_temp; 
    end
    xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [m]')

    if((data(Ny*i,1) > 25-delta) && (data(Ny*i,1) < 25+delta))
        SaveIMG("7_3_BplotSimulation1");
    end
    if((data(Ny*i,1) > 50-delta) && (data(Ny*i,1) < 50+delta))
        SaveIMG("7_3_BplotSimulation2");
    end
    if((data(Ny*i,1) > 75-delta) && (data(Ny*i,1) < 75+delta))
        SaveIMG("7_3_BplotSimulation3");
    end
    if((data(Ny*i,1) > 90-delta) && (data(Ny*i,1) < 90+delta))
        SaveIMG("7_3_BplotSimulation4");
    end
    pause(1.e-20);
end
disp("Hmax = " + num2str(Hmax));
