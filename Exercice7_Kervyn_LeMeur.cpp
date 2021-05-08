#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.h"

using namespace std;

double const pi=3.14159265358979323846264338327950288419716939937510582097494459230e0;

typedef enum{not_defined,dirichlet,neumann,harmonic} bound_cond;

// methode pour lire les conditions aux limites
bound_cond read_boundary_condition(string const& bnd_in){
  if(bnd_in == "dirichlet")
    return dirichlet;
  else if(bnd_in == "neumann")
    return neumann;
  else if(bnd_in == "harmonic")
    return harmonic;
  else
  {
    cerr << "Please select "+bnd_in+"=""dirichlet"", ""neumann"", ""harmonic""." << endl;
    return not_defined;
  }
}

//
// Objets pour la fonction u^2(x,y)
//
class U2 {
public:

  U2(ConfigFile const& configFile){
    xL=configFile.get<double>("xL");
    xR=configFile.get<double>("xR");
    yL=configFile.get<double>("yL");
    yU=configFile.get<double>("yU");
  }

  double get_left_extremum()  {return xL;}
  double get_right_extremum() {return xR;}
  double get_upper_extremum() {return yU;}
  double get_lower_extremum() {return yL;}

  // Methodes virtuelles pures => classe abstraite
  virtual double operator()(double const& x, double const& y) const = 0; // Evalue u^2 au point (x,y)
  virtual double max() const = 0; // Renvoie max(u^2(x,y))

protected:
  double xL,xR,yL,yU;

};

class U2_const: public U2 {
public:
  // Pas de constructeur par defaut => on force a specifier une valeur
  U2_const(ConfigFile const& configFile) :
  U2(configFile), u2(pow(configFile.get<double>("u"),2))
  {}

  // Definition des methodes virtuelles pures :
  double operator()(double const& x,double const& y) const
  {
    return u2;
  }

  double max() const
  {
    return u2;
  }

private:
  double u2;
};

class U2_onde: public U2 {
public:

  U2_onde(ConfigFile const& configFile) :
  U2(configFile),
  g(configFile.get<double>("g")),
  h0(configFile.get<double>("h0")),
  h1(configFile.get<double>("h1")),
  a(configFile.get<double>("a")),
  b(configFile.get<double>("b"))
  {}

protected:
  double g, h0, h1, a, b, Lx;
};

class U2_onde_cas1: public U2_onde {
public:
  U2_onde_cas1(ConfigFile const& configFile) :
  U2_onde(configFile) {}

  double h(double const& x,double const& y) const {
    double h_val(h0);
    if(x>a && x<b) {h_val = h0+(h1-h0)*sin(pi*(x-a)/(b-a));}
    return h_val;
  }

  double h_prime(double const& x) const {
    return (pi*(h1-h0)*cos(pi*(x-a)/(b-a)))/(b-a);
  }

  double h_prime2(double const& x) const {
    return -(pi*pi*(h1-h0)*sin(pi*(x-a)/(b-a)))/((b-a)*(b-a));
  }

  double operator()(double const& x,double const& y) const
  {
    return g * h(x,y);
  }

  double max() const
  {
    unsigned int iteration(0),maxit(1000);
    double error(1.e20),tol(1.e-10);
    double x(5.e-1*(a+b));

   // use the newton method for finding the maximum / minimum velocity
   while((error >= tol) && (iteration <= maxit)){
     // compute new solution
     x = x - (h_prime(x)/h_prime2(x));
     // compute new error
     error = abs(h_prime(x));
     // increase iterator
     iteration += 1;
   }

   // send an error if an extremum has been found
   if(iteration > maxit){
     cout << "Newton method: convergence failed!" << endl;
   }
   // check if we are at a local maximum
   if(error>=tol){
     cout << "Newton method: stationary point not found!" << endl;
   } else {
     if(h_prime2(x)<0.e0){
       cout << "Newton method: maximum not found!" << endl;
     }
   }
   // calculer la vitesse maximale
   return g*(std::max(h0,std::max(h(std::max(a,std::min(x,b)),0.),h1)));
  }

};

class U2_onde_cas2: public U2_onde {
protected:
  double Ly;
public:
  U2_onde_cas2(ConfigFile const& configFile) :
  U2_onde(configFile), Ly(configFile.get<double>("Ly")) {}

  double h(double const& x,double const& y) const {
  double h_value(h0);
  if(x>a && x<b ) h_value = h0+(h1-h0)*sin(pi*(x-a)/(b-a))*sin(pi*y/Ly);
    return h_value;
  }

  vector<double> h_prime(double const& x,double const& y) const
  {
    vector<double> grad(2);
    grad[0] = (pi*(h1-h0)*cos(pi*(x-a)/(b-a))*sin(pi*y/Ly))/(b-a);
    grad[1] = (pi*(h1-h0)*sin(pi*(x-a)/(b-a))*cos(pi*y/Ly))/Ly;
    return grad;
  }

  vector<double> h_prime2(double const& x,double const& y) const {
    vector<double> d2h(3,0.0);
    d2h[0] = -(pi*pi*(h1-h0)*sin(pi*(x-a)/(b-a))*sin(pi*y/Ly))/((b-a)*(b-a));
    d2h[1] = (pi*pi*(h1-h0)*cos(pi*(x-a)/(b-a))*cos(pi*y/Ly))/((b-a)*Ly);
    d2h[2] = -(pi*pi*(h1-h0)*sin(pi*(x-a)/(b-a))*sin(pi*y/Ly))/(Ly*Ly);
    return d2h;
  }

  vector<vector<double>> h_prime2_inv(double const& x,double const& y) const
  {
    vector<vector<double>> jac_inv(2,vector<double>(2));
    vector<double> d2h(h_prime2(x,y));
    double den=0.0;

    // calculer l inverse du jacobian
    den = d2h[0]*d2h[2]-d2h[1]*d2h[1];
    jac_inv[0][0] = d2h[2]/den;
    jac_inv[0][1] = -d2h[1]/den;
    jac_inv[1][0] = -d2h[1]/den;
    jac_inv[1][1] = d2h[0]/den;
    return jac_inv;
  }

  double operator()(double const& x,double const& y) const
  {
    return g * h(x,y);
  }

  double max() const
  {
    unsigned int iteration(0),maxit(1000);
    double error(1.e20),tol(1.e-10),x_const(0.0),y_const(0.0);
    vector<double> x(2,0.e0);
    vector<double> grad(2);
    vector<double> d2h(3);
    vector<vector<double>> jac_inv(2,vector<double>(2));

    // initialisation des vitesses aux milieu du maillage
    x[0] = 0.5*(b+a);
    x[1] = 0.5*Ly;

   // use the newton method for finding the maximum / minimum velocity
   while((error >= tol) && (iteration <= maxit)){
     // compute the gradient
     grad = h_prime(x[0],x[1]);
     // compute the inverse of the jacobian
     jac_inv = h_prime2_inv(x[0],x[1]);
     // compute new solution
     x[0] = x[0]-jac_inv[0][0]*grad[0] - jac_inv[0][1]*grad[1];
     x[1] = x[1]-jac_inv[1][0]*grad[0] - jac_inv[1][1]*grad[1];
     // compute new error
     grad = h_prime(x[0],x[1]);
     error = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
     // increase iterator
     iteration += 1;
   }

   // send an error if an extremum has been found
   if(iteration > maxit){
     cout << "Newton method: convergence failed!" << endl;
   }
   // check if we are at a local maximum
   if(error>=tol){
     cout << "Newton method: stationary point not found!" << endl;
   } else {
       d2h = h_prime2(x[0],x[1]);
       if(!(d2h[0]<0 && d2h[0]*d2h[2]-d2h[1]*d2h[1]>0)){
         cout << "Newton method: maximum not found!" << endl;
       }
     }
   // calculer la vitesse maximale
   x_const = std::max(a,std::min(x[0],b));
   y_const = std::max(0.,std::min(x[1],Ly));
   return g * std::max(h0,std::max(h(x_const,y_const),h1));
  };
};

//
//  Calcul de l'energie de l'onde
//
double energy(vector<vector<double>> const& f, double const& dx, double const& dy,\
unsigned int const& start_xiD,unsigned int const& start_yiD,\
unsigned int const& end_xiD,unsigned int const& end_yiD)
{
  double energy_(0.);
  // en utilisant la formule des trapezes
  for(unsigned int i(start_xiD); i < end_xiD; ++i)
  {
      for(unsigned int j(start_yiD); j < end_yiD; ++j)
      {
          energy_ += pow(f[i][j],2) + pow(f[i+1][j],2) + pow(f[i][j+1],2) + pow(f[i+1][j+1],2);
      }
  }
  return 0.25*dx*dy*energy_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau 1D
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());

  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";

  if(len > 0)
    o << v[len-1];

  return o;
}
//
// Surcharge de l'operateur pour ecrire les elements d'un tableau 2D
//
template<typename T> void print_table(ostream& o, double const& t , vector<vector<T>> const& v,\
unsigned int const& start_xiD,unsigned int const& start_yiD,\
unsigned int const& end_xiD,unsigned int const& end_yiD)
{
  unsigned int i,j;
  for(i=start_yiD; i<end_yiD+1; ++i){
      o << t << " ";
      for(j=start_xiD; j<end_xiD; ++j){
        o << v[i][j] << " ";
      }
      o << v[i][end_xiD] << endl;
  }
}

// TODO: calculer la perturbation
double perturbation(double const& t,double const& x,double const& y,\
double const& xL, double const& xR,double const& yL, double const& yU, \
double const& pert_velocity, double const& pert_amplitude, \
int const mode_num_x,int const mode_num_y) {

  return pert_amplitude*sin(mode_num_x*pi*x/(xR-xL))*sin(mode_num_y*pi*y/(yU-yL))*sin(pert_velocity*t);
}

//
// Main
//
int main(int argc, char* argv[])
{
  unsigned int i,j;
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int k=2; k<argc; ++k) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[k]);

  // Parametres de simulation :
  double tfin		        = configFile.get<double>("tfin");
  unsigned int Nx_real	= configFile.get<int>("Nx");
  unsigned int Ny_real	= configFile.get<int>("Ny");
  double CFL		        = configFile.get<double>("CFL");
  //ajout personnel
  bool ComputeDt        = configFile.get<bool>("ComputeDt");
  //
  string type_u2 	      = configFile.get<string>("type_u2").c_str();
  double pert_amplitude = configFile.get<double>("pert_amplitude");
  double pert_velocity  = configFile.get<double>("pert_velocity");
  int mode_num_x	      = configFile.get<int>("mode_num_x");
  int mode_num_y	      = configFile.get<int>("mode_num_y");

  U2* u2;
  if(type_u2 == "const")
    u2 = new U2_const(configFile);
  else if(type_u2 == "onde_cas1")
    u2 = new U2_onde_cas1(configFile);
  else if(type_u2 == "onde_cas2")
    u2 = new U2_onde_cas2(configFile);
  else
  {
    cerr << "Please select type_u2=""const"" or ""onde_cas1"" or ""onde_cas2""." << endl;
    return -1;
  }

  // TODO calculer le maillage et le pas de temps
  double dx = abs(u2->get_right_extremum()-u2->get_left_extremum()) /double(Nx_real-1);
  double dy = abs(u2->get_upper_extremum()-u2->get_lower_extremum())/double(Ny_real-1);
  double dt = configFile.get<double>("dt");
  if(ComputeDt) { //si l'option est valide, on détermine dt à partir de la formule de l'énoncé
    dt  = (CFL*dx*dy)/sqrt((pow(dx,2)+pow(dy,2))*u2->max()); cout << dt << endl;
  }

  bool write_mesh = configFile.get<bool>("write_mesh"); // Exporter x,y
  bool write_f    = configFile.get<bool>("write_f");    // Exporter f(x,y,t) ou non

  //Conditions aux bords (les strings sont converties en valeurs numeriques a l'aide d'un enumerateur) :
  bound_cond bc_left, bc_right, bc_upper, bc_lower;
  unsigned int Nx=Nx_real,Ny=Ny_real,start_xiD=0,end_xiD=Nx_real-1,start_yiD=0,end_yiD=Ny_real-1;
  bc_left = read_boundary_condition(configFile.get<string>("bc_left").c_str());
  if(bc_left==neumann){
    ++start_xiD;
    ++end_xiD;
    ++Nx;
  }
  bc_right = read_boundary_condition(configFile.get<string>("bc_right").c_str());
  if(bc_right==neumann) ++Nx; // TODO
  bc_lower = read_boundary_condition(configFile.get<string>("bc_lower").c_str());
  if(bc_lower==neumann){
    ++start_yiD;
    ++end_yiD;
    ++Ny;
  }
  bc_upper = read_boundary_condition(configFile.get<string>("bc_upper").c_str());
  if(bc_upper==neumann) ++Ny; // TODO

  // lire si un impulsion doit etre execute:
  // c'est une option facultative, qui veut dire qu'on excite une seule periode d'excitation
  bool impulsion(configFile.get<bool>("impulsion"));

  double A, omega; // Parametres d'excitation
  if(bc_left == harmonic || bc_right == harmonic)
  {
    A = configFile.get<double>("A");
    omega = configFile.get<double>("omega");
  }

  // Fichiers de sortie :
  ofstream fichier_mesh(configFile.get<string>("output_mesh").c_str());
  fichier_mesh.precision(15);

  ofstream fichier_f(configFile.get<string>("output_file").c_str());
  fichier_f.precision(15);

  ofstream fichier_E(configFile.get<string>("output_energy").c_str());
  fichier_E.precision(15);

  ofstream fichier_u(configFile.get<string>("output_velocity").c_str());
  fichier_u.precision(15);
  for(double y(u2->get_lower_extremum()); y<=u2->get_upper_extremum()+.5*dy; y+=dy){
    for(double x(u2->get_left_extremum()); x<=u2->get_right_extremum()-.5*dx; x+=dx){
      fichier_u << sqrt((*u2)(x,y)) << " ";
    }
    fichier_u << sqrt((*u2)(u2->get_right_extremum(),y)) << endl;
  }
  fichier_u.close();

  // lire les donnees pour initialiser les tableaux
  string type_init(configFile.get<string>("type_init").c_str());
  double F0(configFile.get<double>("F0"));
  double u2_loc(0);

  // TODO: calcul de l'inverse de la longueur d'onde en x et y: k_x=m*pi/L_x; k_y=n*pi/L_y;
  double k_wave_x(pi*mode_num_x/abs(u2->get_right_extremum()-u2->get_left_extremum()) );
  double k_wave_y(pi*mode_num_y/abs(u2->get_upper_extremum()-u2->get_lower_extremum()));

  // Initialisation des tableaux du schema numerique :
  vector<vector<double>> fpast(Ny,vector<double>(Nx)), fnow(Ny,vector<double>(Nx)),\
  fnext(Ny,vector<double>(Nx)), betax2(Ny,vector<double>(Nx)),betay2(Ny,vector<double>(Nx));

  if(type_init=="harmonic"){
    // On initialise alors un mode propre (m,n);
    // TODO: completer fnow, fpast, beta_x^2, beta_y^2
    // Note : La syntaxe pour evaluer u^2 au point x est (*u2)(x,y)
    for(unsigned int i(start_yiD); i < end_yiD+1; ++i) {
        for(unsigned int j(start_xiD); j < end_xiD+1; ++j) {
          fpast[i][j]   = F0*sin(k_wave_x*(u2->get_left_extremum()+j*dx))*sin(k_wave_y*(u2->get_lower_extremum()+i*dy));
          fnow[i][j]    = fpast[i][j];
          betax2[i][j]  = (*u2)(u2->get_left_extremum()+j*dx,u2->get_lower_extremum()+i*dy) * pow(dt,2)/pow(dx,2);
          betay2[i][j]  = (*u2)(u2->get_left_extremum()+j*dx,u2->get_lower_extremum()+i*dy) * pow(dt,2)/pow(dy,2);
        }
    }

  } else {
    // On initialise alors une perturbation nulle.
    // TODO: completer fnow, fpast, beta_x^2, beta_y^2
    // Note : La syntaxe pour evaluer u^2 au point x est (*u2)(x,y)
    for(unsigned int i(start_yiD); i < end_yiD+1; ++i) {
        for(unsigned int j(start_xiD); j < end_yiD+1; ++j) {
          fpast[i][j]   = 0.;
          fnow[i][j]    = 0.;
          betax2[i][j]  = (*u2)(u2->get_left_extremum()+j*dx,u2->get_lower_extremum()+i*dy) * pow(dt,2)/pow(dx,2);
          betay2[i][j]  = (*u2)(u2->get_left_extremum()+j*dx,u2->get_lower_extremum()+i*dy) * pow(dt,2)/pow(dy,2);
        }
    }
  }

  // Boucle temporelle :
  double t, amplitude_bnd;
  unsigned int stride(0);
  unsigned int n_stride(configFile.get<unsigned int>("n_stride"));
  // put has first line the position vector
  vector<double> x_mesh(Nx_real),y_mesh(Ny_real);
  for(i=0; i<Nx_real; ++i){
    x_mesh[i] = u2->get_left_extremum()+(i)*dx;
  }
  for(i=0; i<Ny_real; ++i){
    y_mesh[i] = u2->get_lower_extremum()+(i)*dy;
  }
  /*for(i=start_xiD; i<end_xiD+1; ++i){
    x_mesh[i] = u2->get_left_extremum()+(i-start_xiD)*dx;
  }
  for(i=start_yiD; i<end_yiD+1; ++i){
    y_mesh[i] = u2->get_lower_extremum()+(i-start_yiD)*dy;
  }*/
  if(write_mesh){
    fichier_mesh << x_mesh << endl;
    fichier_mesh << y_mesh << endl;
    fichier_mesh.close();
  }

  for(t=0.; t<=tfin-dt*9./10; t+=dt)
  {
    // Ecriture :
    if(stride >= n_stride )
    {
      if(write_f) print_table(fichier_f,t,fnow,start_xiD,start_yiD,end_xiD,end_yiD);
      fichier_E << t << " " << energy(fnow,dx,dy,start_xiD,start_yiD,end_xiD,end_yiD) << endl;
      stride = 0;
    }
    ++stride;
    // TODO: calculer fnext selon le schema explicite a 3 niveaux
    for(i=1; i<Ny-1; ++i){
      for(j=1; j<Nx-1; ++j){
        fnext[i][j] =     betay2[i][j]*(fnow[i+1][j] - 2*fnow[i][j] + fnow[i-1][j])
                        + betax2[i][j]*(fnow[i][j+1] - 2*fnow[i][j] + fnow[i][j-1])
                        /*+ 0.25*(betay2[i+1][j]-betay2[i-1][j])*(fnow[i+1][j] - fnow[i-1][j])
                        + 0.25*(betax2[i][j+1]-betax2[i][j-1])*(fnow[i][j+1] - fnow[i][j-1])*/
                        + pow(dt,2) * perturbation(t,x_mesh[j],y_mesh[i],u2->get_left_extremum(),u2->get_right_extremum(),u2->get_lower_extremum(),u2->get_upper_extremum(),pert_velocity,pert_amplitude,mode_num_x,mode_num_y)
                        + 2*fnow[i][j] - fpast[i][j];
      }
    }

    // Conditions aux bords:
    switch(bc_left)
    { // condition au bord "gauche" (x=0)
      case dirichlet: for(unsigned int j(0); j<Ny; ++j) fnext[j][0] = 0;    // ("fixe")
        break;
      case neumann:   for(unsigned int j(0); j<Ny; ++j) fnext[j][0] = fnext[j][2]; // ("libre")
        break;
      case harmonic:  for(unsigned int j(0); j<Ny; ++j) fnext[j][0] = A*sin(omega*(t+dt));
      if(impulsion){
        if(t > 2*pi/omega)
          bc_left = neumann;
      }
        break;
      default:
        throw "Invalid left boundary condition!";
    }

    switch(bc_right) // condition au bord droite (x=L_x)
    {
      case dirichlet: for(unsigned int j(0); j<Ny; ++j) fnext[j][Nx-1] = 0; // ("fixe")
        break;
      case neumann:   for(unsigned int j(0); j<Ny; ++j) fnext[j][Nx-1] = fnext[j][Nx-3];// ("libre")
        break;
      case harmonic:  for(unsigned int j(0); j<Ny; ++j) fnext[j][Nx-1] = A*sin(omega*(t+dt)); // TODO : Completer la condition au bord droit harmonique
      if(impulsion){
        if(t > 2*pi/omega)
          bc_right = neumann;
      }
        break;
      default:
        throw "Invalid right boundary condition!";
    }

    switch(bc_lower) // condition au bord inferieur (y=0)
    {
      case dirichlet: for(unsigned int i(0); i<Nx; ++i) fnext[0][i] = 0; // ("fixe")
        break;
      case neumann:   for(unsigned int i(0); i<Nx; ++i) fnext[0][i] = fnext[2][i]; // ("libre")
        break;
      case harmonic:  for(unsigned int i(0); i<Nx; ++i) fnext[0][i] = A*sin(omega*(t+dt)); // TODO : Completer la condition au bord inferieur harmonique
      if(impulsion){
        if(t > 2*pi/omega)
          bc_lower = neumann;
      }
        break;
      default:
        throw "Invalid lower boundary condition!";
    }
    switch(bc_upper) // condition au bord superieur (y=L_y)
    {
      case dirichlet: for(unsigned int i(0); i<Nx; ++i) fnext[Ny-1][i] = 0;// ("fixe")
        break;
      case neumann:   for(unsigned int i(0); i<Nx; ++i) fnext[Ny-1][i] = fnext[Ny-3][i]; // ("libre")
        break;
      case harmonic:  for(unsigned int i(0); i<Nx; ++i) fnext[Ny-1][i] = A*sin(omega*(t+dt)); // TODO : Completer la condition au bord superieur harmonique
      if(impulsion){
        if(t > 2*pi/omega)
          bc_upper = neumann;
      }
        break;
      default:
        throw "Invalid upper boundary condition!";
    }

    // Mise a jour :
    fpast = fnow;
    fnow  = fnext;
  }

  // Ecrire la solution
  if(write_f) print_table(fichier_f,t,fnow,start_xiD,start_yiD,end_xiD,end_yiD);
  fichier_E << t << " " << energy(fnow,dx,dy,start_xiD,start_yiD,end_xiD,end_yiD) << endl;

  fichier_f.close();
  fichier_E.close();

  return 0;
}
