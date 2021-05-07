% % Matlab: question d - limite de stabilité
% clear; close all; clc;
% 
% Parametres physiques:
% xL   = 0.0 ; 		     xR   = 10.0 ; 		     
% yL   = 0.0 ;		     yU   = 6.0 ;		%dimensions (m)     
% pert_amplitude = 3.e0;	 pert_velocity  = 0.e0;	
% 
% Parametres physiques milieu uniforme:
% u = 4.e0; % vitesse de propagation (m/s)
% 
% Parametres physiques onde de Belharra:
% g   = 9.81;      % acceleration de gravite
% h0  = 2.e0;      % profondeur maximale du fond marin (m)
% h1  = 1.e0;      % profondeur minimale du fond marin (m)
% a   = 2000.0;    % limite inferieur de l'onde en x (m)
% b   = 5000.0;    % limite superieur de l'onde en y (m)
% Ly  = 2000.0;    % longueur characteristique de l'onde en y (m)
% 
% Parametres numeriques :
% nombre de node du maillage en x,y (extreme inclus)
% Nx = 64; 			 Ny = 64;
% 
% nsimul = 5;
% ComputeDt = false;                   %ajout PERSO - option pour le calcul du dt
% CFL = linspace(0.8,1.2,nsimul);         % valeur de la condition CFL
% hx  = (xR-xL)/(Nx-1); hy  = (yU-yL)/(Ny-1);
% dt  = sqrt(CFL.^2 / (u^2 * (1/(hx)^2 + 1/(hy)^2)));
% 
% type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
% ecrire_f = true;		 % if true f est ecrite
% 
% valeur propre en x,y
% mode_num_x = 3;			 mode_num_y = 5;			 
% 
% conditions aux bords: dirichlet, neumann, harmonic
% bc_left   = 'dirichlet';  bc_right  = 'dirichlet';   
% bc_lower  = 'dirichlet';  bc_upper  = 'dirichlet';   
% 
% impulsion = false;          % si vrais et harmonique une seule onde est emise
% type_init = 'harmonic';     % type initialisation: harmonic, default: homogene
% F0 = 1.e0;			        % amplitude de l'harmonique initiale
% A = 1.e0;			        % amplitude condition bord harmonique
% omega = 10.e0;			    % frequence condition bord harmonique
% write_mesh = true;		    % si vrai le maillage est ecrite dans output_mesh.out
% write_f = true;		 	    % si vrai la solution est ecrite dans output_f.out
% n_stride = 0;			    % nombre de stride
% 
% période d'oscillation th.
% T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
% tfin = 5*T;
% omegaMN = 2*pi/T;
% 
% Simulations
% filename2  = "CFL_"+ string(CFL);
% 
% for i=1:nsimul
%     CFL_loc = CFL(i);
%     dt_loc  = dt(i);
%     writeConfig3;
%     disp('Exercice7_Kervyn_LeMeur configuration.in');   
%     system('Exercice7_Kervyn_LeMeur configuration.in'); 
% end
% 
% ViewFormat; facteurTemps = 8;
% for k=1:nsimul
%     data = load(filename2(k)+'_f.out');
%     mesh = load(filename2(k)+'_mesh.out');
%     figure('Name',"plotsimulation CFL="+num2str(CFL(k)))
%     for i =1:length(data(:,1))/Ny
%         surfc(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1));
%         title("$t =$"+num2str(data(Ny*i,1))+"s, $\beta_{CFL} =$"+num2str(CFL(k)),'Interpreter','latex');
%         grid minor; set(gca,'fontsize',fs);
%         xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
%         zlim([-3 3]);
%         view(0,0);
%         pause((1/facteurTemps)*dt);
%         pause();
%     end
% end

%% limite de stabilité avec epsilon
clear; clc;

%Parametres physiques:
xL   = 0.0 ; 		     xR   = 10.0 ; 		     
yL   = 0.0 ;		     yU   = 6.0 ;		%dimensions (m)     
pert_amplitude = 3.e0;	 pert_velocity  = 0.e0;	

%Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

%Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 2.e0;      % profondeur maximale du fond marin (m)
h1  = 1.e0;      % profondeur minimale du fond marin (m)
a   = 2000.0;    % limite inferieur de l'onde en x (m)
b   = 5000.0;    % limite superieur de l'onde en y (m)
Ly  = 2000.0;    % longueur characteristique de l'onde en y (m)

%Parametres numeriques :
%nombre de node du maillage en x,y (extreme inclus)
Nx = 64; 			 Ny = 64;

% nsimul = 2;
ComputeDt = false;      % ajout PERSO - option pour le calcul du dt
CFL = [1 1.03];         % valeur de la condition CFL
nsimul = length(CFL);
hx  = (xR-xL)/(Nx-1); hy  = (yU-yL)/(Ny-1);
dt  = sqrt(CFL.^2 / (u^2 * (1/(hx)^2 + 1/(hy)^2)));

type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;		 % if true f est ecrite

%valeur propre en x,y
mode_num_x = 3;			 mode_num_y = 5;			 

%conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'dirichlet';  bc_right  = 'dirichlet';   
bc_lower  = 'dirichlet';  bc_upper  = 'dirichlet';   

impulsion = false;          % si vrais et harmonique une seule onde est emise
type_init = 'harmonic';     % type initialisation: harmonic, default: homogene
F0 = 1.e0;			        % amplitude de l'harmonique initiale
A = 1.e0;			        % amplitude condition bord harmonique
omega = 10.e0;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrai le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrai la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

%période d'oscillation th.
T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
tfin = 5*T;
omegaMN = 2*pi/T;

%Simulations
filename2  = "CFL_"+ string(CFL);

for i=1:nsimul
    CFL_loc = CFL(i);
    dt_loc  = dt(i);
    writeConfig3;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 
end

ViewFormat; 
for k=1:nsimul
    delta = dt(k)/2;
    data = load(filename2(k)+'_f.out');
    mesh = load(filename2(k)+'_mesh.out');
    figure('Name',"plotsimulation CFL="+num2str(CFL(k)))
    for i =1:length(data(:,1))/Ny
        F = data(Ny*(i-1)+1:Ny*i,2:Nx+1);
        surfc(mesh(1,:),mesh(2,:),F);
        title("$t =$"+num2str(data(Ny*i,1))+"s, $\beta_{CFL} =$"+num2str(CFL(k)),'Interpreter','latex');
        grid minor; set(gca,'fontsize',fs);
        xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
        zlim([-1.8 1.8]);
%         view(0,0);
        if((data(Ny*i,1) > 1.49-delta) && (data(Ny*i,1) < 1.49+delta))
            SaveIMG("DsurfaceCFL"+num2str(k)+"_1");
        end
        if((data(Ny*i,1) > 1.68-delta) && (data(Ny*i,1) < 1.68+delta))
            SaveIMG("DsurfaceCFL"+num2str(k)+"_2");
        end
        if((data(Ny*i,1) > 1.85-delta) && (data(Ny*i,1) < 1.85+delta))
            SaveIMG("DsurfaceCFL"+num2str(k)+"_3");
        end
        pause(1.e-9);
%         pause();
    end
end

%% limite de stabilité sur une coupe
clear; clc;

%Parametres physiques:
xL   = 0.0 ; 		     xR   = 10.0 ; 		     
yL   = 0.0 ;		     yU   = 6.0 ;		%dimensions (m)     
pert_amplitude = 3.e0;	 pert_velocity  = 0.e0;	

%Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

%Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 2.e0;      % profondeur maximale du fond marin (m)
h1  = 1.e0;      % profondeur minimale du fond marin (m)
a   = 2000.0;    % limite inferieur de l'onde en x (m)
b   = 5000.0;    % limite superieur de l'onde en y (m)
Ly  = 2000.0;    % longueur characteristique de l'onde en y (m)

%Parametres numeriques :
%nombre de node du maillage en x,y (extreme inclus)
Nx = 64; 			 Ny = 64;

% nsimul = 2;
ComputeDt = false;      % ajout PERSO - option pour le calcul du dt
CFL = 1.01;         % valeur de la condition CFL
nsimul = length(CFL);
hx  = (xR-xL)/(Nx-1); hy  = (yU-yL)/(Ny-1);
dt  = sqrt(CFL.^2 / (u^2 * (1/(hx)^2 + 1/(hy)^2)));

type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;		 % if true f est ecrite

%valeur propre en x,y
mode_num_x = 3;			 mode_num_y = 5;			 

%conditions aux bords: dirichlet, neumann, harmonic
bc_left   = 'dirichlet';  bc_right  = 'dirichlet';   
bc_lower  = 'dirichlet';  bc_upper  = 'dirichlet';   

impulsion = false;          % si vrais et harmonique une seule onde est emise
type_init = 'harmonic';     % type initialisation: harmonic, default: homogene
F0 = 1.e0;			        % amplitude de l'harmonique initiale
A = 1.e0;			        % amplitude condition bord harmonique
omega = 10.e0;			    % frequence condition bord harmonique
write_mesh = true;		    % si vrai le maillage est ecrite dans output_mesh.out
write_f = true;		 	    % si vrai la solution est ecrite dans output_f.out
n_stride = 0;			    % nombre de stride

%période d'oscillation th.
T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
tfin = 5.5*T;
omegaMN = 2*pi/T;

%Simulations
filename2  = "CFL_"+ num2str(CFL);
    CFL_loc = CFL;
    dt_loc  = dt;
    writeConfig3;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 

ViewFormat; 
% delta = dt/2;
data = load(filename2+'_f.out');
mesh = load(filename2+'_mesh.out');

start = 2;
newdata = data(start:Ny:end,:);
posy = mesh(2,start);

figure('Name',"plotStabilité CFL="+num2str(CFL))
        surf(mesh(1,:),newdata(:,1),newdata(:,2:end));
        title("$\beta_{CFL} =$"+num2str(CFL),'Interpreter','latex');
        grid minor; set(gca,'fontsize',fs);
        xlabel('$x$ [m]'); ylabel('$t$ [s]'); zlabel("$f(x,y=y_{loc},t)$ [a.u.]")
%         zlim([-1.8 1.8]);
%         view(0,0);
%         pause();
SaveIMG("StabiliteCFL")
