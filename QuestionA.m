%% Question 7.2.a - première partie (bord droit libre)
clear all; close all; clc;

% Parametres physiques:
tfin = 8.0;		         % temps physique de la simulation (s)
xL   = 0.0 ; 		     xR   = 10.0 ; 		     
yL   = 0.0 ;		     yU   = 6.0 ;		     
pert_amplitude = 0.e0;	 pert_velocity  = 0.e0;	 

% Parametres physiques milieu uniforme:
u = 4.e0; % vitesse de propagation (m/s)

% Parametres physiques onde de Belharra:
g   = 9.81;      % acceleration de gravite
h0  = 2.e0;      % profondeur maximale du fond marin (m)
h1  = 1.e0;      % profondeur minimale du fond marin (m)
a   = 2000.0; % limite inferieur de l'onde en x (m)
b   = 5000.0; % limite superieur de l'onde en y (m)
Ly  = 2000.0; % longueur characteristique de l'onde en y (m)

% Parametres numeriques :
Nx = 64; 			 % nombre de node du maillage en x (extreme inclus)
Ny = 64; 			 % nombre de node du maillage en y (extreme inclus)

ComputeDt = false;   %ajout PERSO - option pour le calcul du dt
dt  = 0.02;           %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
CFL = 1.0;           % valeur de la condition CFL

type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
ecrire_f = true;		 % if true f est ecrite
mode_num_x = 2;			 % valeur propre en x
mode_num_y = 2;			 % valeur propre en y
bc_left   = 'harmonic';  % condition au bord gauche: dirichlet, neumann, harmonic
bc_right  = 'neumann';   % condition au bord droite: dirichlet, neumann, harmonic
bc_lower  = 'neumann';   % condition au bord inferieur: dirichlet, neumann, harmonic
bc_upper  = 'neumann';   % condition au bord superieur: dirichlet, neumann, harmonic
impulsion = true;		 % si vrais et harmonique une seule onde est emise
type_init = 'homogene';  % type initialisation: harmonic, default: homogene
F0 = 1.e0;			           % amplitude de l'harmonique initiale
A = 1.e0;			           % amplitude condition bord harmonique
omega = 5.e0;			       % frequence condition bord harmonique
write_mesh = true;		     % si vrais le maillage est ecrite dans output_mesh.out
write_f = true;		 	     % si vrais la solution est ecrite dans output_f.out
n_stride = 0;			     % nombre de stride

%période d'oscillation th.
T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
omegaMN = 2*pi/T;

% Simulations
    filename2  = ['AA_Nx_',num2str(Nx),'Ny_',num2str(Ny)];
    Nx_loc = Nx; Ny_loc = Ny; fname_loc = filename2;
    writeConfig1;
    disp('Exercice7_Kervyn_LeMeur configuration.in');   
    system('Exercice7_Kervyn_LeMeur configuration.in'); 

%% Plot simulation
% ViewFormat;
% data = load([filename2,'_f.out']);
% mesh = load([filename2,'_mesh.out']);
% 
% start = 2;
% newdata = data(start:Ny:end,:);
% posy = mesh(2,start);
% 
% delta = dt/2;
% figure('Name','plotsimulationAA')
% for i =1:length(newdata(:,1))
%     plot(mesh(1,:),newdata(i,2:end),'-','Linewidth',lw);
%     xlabel('$x$ [m]'); ylabel('$f$ [a.u.]'); 
%     title("$t =$ "+num2str(newdata(i,1))+"s",'Interpreter','latex');
%     set(gca,'fontsize',fs); grid on;
%     ylim([-2,2]);
%     if((newdata(i,1) > 2-delta) && (newdata(i,1) < 2+delta))
%         SaveIMG("AAsurface1");
%     end
%     if((newdata(i,1) > 2.7-delta) && (newdata(i,1) < 2.7+delta))
%         SaveIMG("AAsurface2");
%     end
%     if((newdata(i,1) > 4-delta) && (newdata(i,1) < 4+delta))
%         SaveIMG("AAsurface3");
%     end
%     if((newdata(i,1) > 5-delta) && (newdata(i,1) < 5+delta))
%         SaveIMG("AAsurface4");
%     end
%     pause(1.e-9);
% end

%
% delta = dt/2;
% % figure('Name','plotsimulationAA')
% [X, Y]  = meshgrid(mesh(1,:),mesh(2,:));
% axes1 = axes('Parent',figure);
% for i =1:length(data(:,1))/Ny
%     s1 = pcolor(X,Y,data(Ny*(i-1)+1:Ny*i,2:Nx+1));
%     s1.FaceColor = 'interp';
%     xlabel('$x$ [m]'); ylabel('$y$ [m]'); 
%     title("$t =$ "+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
%     set(gca,'fontsize',fs)
%     caxis([-(A+0.2) A+0.2]);
%         h = colorbar(axes1);
%         h.Label.String = '$f$ [a.u.]';
%         h.Label.Interpreter = 'Latex';
%         h.Label.FontSize = fs;
%     if((data(Ny*i,1) > 2-delta) && (data(Ny*i,1) < 2+delta))
%         SaveIMG("AAsurface1");
%     end
%     if((data(Ny*i,1) > 2.7-delta) && (data(Ny*i,1) < 2.7+delta))
%         SaveIMG("AAsurface2");
%     end
%     if((data(Ny*i,1) > 4-delta) && (data(Ny*i,1) < 4+delta))
%         SaveIMG("AAsurface3");
%     end
%     if((data(Ny*i,1) > 5-delta) && (data(Ny*i,1) < 5+delta))
%         SaveIMG("AAsurface4");
%     end
%     pause(1.e-9);
% end

% for i =1:length(data(:,1))/Ny
%     pcolor(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1),'FaceColor','interp');
%     title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
%     grid minor; set(gca,'fontsize',fs);
%     xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
%     zlim([-5 5]);
% %     view(0,0);
%     if((data(Ny*i,1) > 2-delta) && (data(Ny*i,1) < 2+delta))
%         SaveIMG("AAsurface1");
%     end
%     if((data(Ny*i,1) > 2.7-delta) && (data(Ny*i,1) < 2.7+delta))
%         SaveIMG("AAsurface2");
%     end
%     if((data(Ny*i,1) > 4-delta) && (data(Ny*i,1) < 4+delta))
%         SaveIMG("AAsurface3");
%     end
%     if((data(Ny*i,1) > 5-delta) && (data(Ny*i,1) < 5+delta))
%         SaveIMG("AAsurface4");
%     end
%     pause(1.e-9);
% end

%% Plot en coupe en Pcolor
ViewFormat;
data = load([filename2,'_f.out']);
mesh = load([filename2,'_mesh.out']);

start = 2;
newdata = data(start:Ny:end,:);
posy = mesh(2,start)

% axes1 = axes('Parent',figure);
% for i =1:length(data(:,1))/Ny
%     s1 = contourf(mesh(1,:),newdata(:,1),newdata(:,2:end));
% %     s1.FaceColor = 'interp';
% %     s1.EdgeColor = 'None';
%     colormap cool
%     xlabel('$x$ [m]'); ylabel('$t$ [s]'); 
%     set(gca,'fontsize',fs)
%         h = colorbar(axes1);
%         h.Label.String = '$f(x,y=y_{loc},t)$ [m]';
%         h.Label.Interpreter = 'Latex';
%         h.Label.FontSize = fs;
% end
% SaveIMG("CoupeAA");

axes1 = axes('Parent',figure);
for i =1:length(data(:,1))/Ny
    s1 = contourf(mesh(1,:),newdata(:,1),newdata(:,2:end),5);
%     s1.FaceColor = 'interp';
%     s1.EdgeColor = 'None';
%     colormap cool
    xlabel('$x$ [m]'); ylabel('$t$ [s]'); 
    set(gca,'fontsize',fs)
        h = colorbar(axes1);
        h.Label.String = '$f(x,y=y_{loc},t)$ [a.u.]';
        h.Label.Interpreter = 'Latex';
        h.Label.FontSize = fs;
end
SaveIMG("CoupeAA");

%% Matlab - question a (bord droit fixe)
% clear all; close all; clc;
% 
% % Parametres physiques:
% tfin = 8.0;		     % temps physique de la simulation (s)
% xL   = 0.0 ; 		     % limite gauche du fond marin (m) const: 0.
% xR   = 10.0 ; 		     % limite droite du fond marin (m) const: 15.
% yL   = 0.0 ;		     % limite inferieur du fond marin (m) const: 0.
% yU   = 6.0 ;		     % limite superieur du fond marin (m) const: 15.
% pert_amplitude = 0.e0;	 % amplitude de la perturbation du terme de droite (m)
% pert_velocity  = 0.e0;	 % vitesse angulaire de la perturbation du terme de droite (m)
% 
% % Parametres physiques milieu uniforme:
% u = 4.e0; % vitesse de propagation (m/s)
% 
% % Parametres physiques onde de Belharra:
% g   = 9.81;      % acceleration de gravite
% h0  = 2.e0;      % profondeur maximale du fond marin (m)
% h1  = 1.e0;      % profondeur minimale du fond marin (m)
% a   = 2000.0; % limite inferieur de l'onde en x (m)
% b   = 5000.0; % limite superieur de l'onde en y (m)
% Ly  = 2000.0; % longueur characteristique de l'onde en y (m)
% 
% % Parametres numeriques :
% Nx = 64; 			 % nombre de node du maillage en x (extreme inclus)
% Ny = 64; 			 % nombre de node du maillage en y (extreme inclus)
% 
% ComputeDt = false;   %ajout PERSO - option pour le calcul du dt
% dt  = 0.02;           %ajout PERSO - si le dt n'est pas calculé, cette valeur est prise
% CFL = 1.0;           % valeur de la condition CFL
% 
% type_u2 = 'const';		 % type de simulation: "const", "onde_cas1" et "onde_cas2"
% ecrire_f = true;		 % if true f est ecrite
% mode_num_x = 2;			 % valeur propre en x
% mode_num_y = 2;			 % valeur propre en y
% bc_left   = 'harmonic';  % condition au bord gauche: dirichlet, neumann, harmonic
% bc_right  = 'dirichlet';   % condition au bord droite: dirichlet, neumann, harmonic
% bc_lower  = 'neumann';   % condition au bord inferieur: dirichlet, neumann, harmonic
% bc_upper  = 'neumann';   % condition au bord superieur: dirichlet, neumann, harmonic
% impulsion = true;		 % si vrais et harmonique une seule onde est emise
% type_init = 'homogene';  % type initialisation: harmonic, default: homogene
% F0 = 1.e0;			           % amplitude de l'harmonique initiale
% A = 1.e0;			           % amplitude condition bord harmonique
% omega = 5.e0;			       % frequence condition bord harmonique
% write_mesh = true;		     % si vrais le maillage est ecrite dans output_mesh.out
% write_f = true;		 	     % si vrais la solution est ecrite dans output_f.out
% n_stride = 0;			     % nombre de stride
% 
% %période d'oscillation th.
% T = 2/(u*sqrt((mode_num_x/(xR-xL))^2 + (mode_num_y/(yU-yL))^2));
% nbperiodes = 5;
% omegaMN = 2*pi/T;
% 
% % Simulations
%     filename2  = ['AB_Nx_',num2str(Nx),'Ny_',num2str(Ny)];
%     Nx_loc = Nx; Ny_loc = Ny; fname_loc = filename2;
%     writeConfig1;
%     disp('Exercice7_Kervyn_LeMeur configuration.in');   
%     system('Exercice7_Kervyn_LeMeur configuration.in'); 
% 
% % Plot simulation
% ViewFormat; delta = dt/2;
% data = load([filename2,'_f.out']);
% mesh = load([filename2,'_mesh.out']);
% 
% start = 2;
% newdata = data(start:Ny:end,:);
% posy = mesh(2,start);
% 
% delta = dt/2;
% figure('Name','plotsimulationAA')
% for i =1:length(newdata(:,1))
%     plot(mesh(1,:),newdata(i,2:end),'-','Linewidth',lw);
%     xlabel('$x$ [m]'); ylabel('$f$ [a.u.]'); 
%     title("$t =$ "+num2str(newdata(i,1))+"s",'Interpreter','latex');
%     set(gca,'fontsize',fs); grid on;
%     ylim([-2,2]);
%     if((newdata(i,1) > 2.08-delta) && (newdata(i,1) < 2.08+delta))
%         SaveIMG("ABsurface1");
%     end
%     if((newdata(i,1) > 3.14-delta) && (newdata(i,1) < 3.14+delta))
%         SaveIMG("ABsurface2");
%     end
%     if((newdata(i,1) > 4.4-delta) && (newdata(i,1) < 4.4+delta))
%         SaveIMG("ABsurface3");
%     end
%     if((newdata(i,1) > 5.9-delta) && (newdata(i,1) < 5.9+delta))
%         SaveIMG("ABsurface5");
%     end
%     if((newdata(i,1) > 6.9-delta) && (newdata(i,1) < 6.9+delta))
%         SaveIMG("ABsurface6");
%     end
%     pause(1.e-9);
% end

% [X, Y]  = meshgrid(mesh(1,:),mesh(2,:));
% axes1 = axes('Parent',figure);
% for i =1:length(data(:,1))/Ny
%     s1 = pcolor(X,Y,data(Ny*(i-1)+1:Ny*i,2:Nx+1));
%     s1.FaceColor = 'interp';
%     xlabel('$x$ [m]'); ylabel('$y$ [m]'); 
%     title("$t =$ "+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
%     set(gca,'fontsize',fs)
%     caxis([-(A+0.2) A+0.2]);
%         h = colorbar(axes1);
%         h.Label.String = '$f$ [a.u.]';
%         h.Label.Interpreter = 'Latex';
%         h.Label.FontSize = fs;
%     if((data(Ny*i,1) > 2.08-delta) && (data(Ny*i,1) < 2.08+delta))
%         SaveIMG("ABsurface1");
%     end
%     if((data(Ny*i,1) > 3.14-delta) && (data(Ny*i,1) < 3.14+delta))
%         SaveIMG("ABsurface2");
%     end
%     if((data(Ny*i,1) > 4.4-delta) && (data(Ny*i,1) < 4.4+delta))
%         SaveIMG("ABsurface3");
%     end
% %     if((data(Ny*i,1) > 2.08-delta) && (data(Ny*i,1) < 2.08+delta))
% %         SaveIMG("ABsurface4");
% %     end
%     if((data(Ny*i,1) > 5.9-delta) && (data(Ny*i,1) < 5.9+delta))
%         SaveIMG("ABsurface5");
%     end
%     if((data(Ny*i,1) > 6.9-delta) && (data(Ny*i,1) < 6.9+delta))
%         SaveIMG("ABsurface6");
%     end
%     pause(1.e-7);
% %     pause();
% end

%%
% figure('Name','plotsimulation')
% for i =1:length(data(:,1))/Ny
%     surf(mesh(1,:),mesh(2,:),data(Ny*(i-1)+1:Ny*i,2:Nx+1));
%     title("$t =$"+num2str(data(Ny*i,1))+"s",'Interpreter','latex');
%     grid minor; set(gca,'fontsize',fs);
%     xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$f$ [a.u.]')
%     zlim([-5 5]);
% %     view(0,0);
%     pause(1.e-9);
% end