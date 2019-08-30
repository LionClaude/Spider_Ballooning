#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
using namespace std;
const int N_beads=15;
const int N_lines=3;
const int N_steps=1000000;



class TemporalSeries{
public:
  void initializeValues();
  void genNumber();
  void genPoint();
  void genGaussian();
  void genSeries();
  void printSeries();
  double T_tot;
  int N_events;
  int i;
  int j;
  int seed;
  double lambda;
  double mu;
  double sigma;
  double sigma_x;
  double thr;
  double w;
  double t;
  double dt;
  double f[N_steps];
  double psi(double x);
  double d_psi(double x);
  default_random_engine g;
} T1;

void TemporalSeries::initializeValues(){
  dt=0.00001;
  T_tot=N_steps*dt;
  lambda=0.4;
  sigma=1.0;
  sigma_x=1.0;
  thr=4*sigma_x;
  w=0.8;
  for (i=0; i<N_steps; i++){
    f[i]=0;
  }
  seed=1234;
};

void TemporalSeries::genNumber(){
  //std::default_random_engine m(1021);
  std::poisson_distribution<int> poisson(lambda*T_tot);
  N_events=poisson(g);
  cout << N_events << "\n";
};

void TemporalSeries::genPoint(){
  std::uniform_real_distribution<double> uniform(0.0,T_tot);
  mu=uniform(g);
};

void TemporalSeries::genGaussian(){
  std::normal_distribution<double> gauss(0.0,1.0);
  double a=gauss(g);
  for (i=0; i<N_steps; i++){
    t=i*dt;
    f[i]+=0.1*a*exp(-pow(t-mu,2)/(2*pow(sigma,2)));
  }
};

void TemporalSeries::genSeries(){
  initializeValues();
  genNumber();
  for (j=0; j<N_events; j++){
    genPoint();
    genGaussian();
  }
  printSeries();
};

void TemporalSeries::printSeries(){
  ofstream outtemp;
  outtemp.open("TSeries.txt");
  for (i=0; i<N_steps; i++){
    t=i*dt;
    outtemp << t << " " << f[i] << "\n";
  }
  outtemp.close();
};

double TemporalSeries::d_psi(double x){
  return exp((-pow(x,2))/(2*pow(sigma_x,2)));
};

double TemporalSeries::psi(double x){
  double delta_x=0.01;
  double y=-thr;
  double s=0;
  while (y<x){
    s+=d_psi(y)*delta_x;
    y+=delta_x;
  }
  return s;
};



class Bead{
  public:
    double q_old[3];
    double q[3];
    double q_new[3];
    double vel[3];
    double wind_vel[3];
    double F_el[3];
    double F_KP[3];
    double F_drag[3];
    double F_weight[3];
    double r;                  //distances between beads
    double *r_left;
    double *r_right;
    double *r_fleft;
    double sp;                 //scalar products
    double *sp_left;
    double *sp_right;
    double theta;
    double phi;
    double *theta_neigh;
    double *phi_neigh;
    double k;
    double J;
    double s0;
    double h;
    double b;
    double m;
    double g;
    double dt;
    void getDistances();
    void getAngles();
    void getFel();
    void getScalarProd();
    void getFKP();
    void getVel();
    void getWindVel();
    void getDrag();
    void getWeight();
    void oneStepProp();
    void firstStepProp();
    int i;
    int n;
    int l;
    double *qlleft[3];
    double *qleft[3];
    double *qright[3];
    double *qrright[3];
};

void Bead::getDistances(){
  r=0;
  for (i=0; i<3; i++){
    r+=pow((*qright[i]-q[i]),2);
  }
  r=sqrt(r);
};

void Bead::getAngles(){
  theta=acos((*qright[2]-q[2])/r);
  phi=atan2((*qright[1]-q[1]),(*qright[0]-q[0]));
};

void Bead::getFel(){
  for (i=0; i<3; i++){
    F_el[i]=0.0;
  }
  if (n<N_beads-1){
    F_el[0]=k*(*qright[0]-q[0]-s0*cos(phi)*sin(theta));
    F_el[1]=k*(*qright[1]-q[1]-s0*sin(phi)*sin(theta));
    F_el[2]=k*(*qright[2]-q[2]-s0*cos(theta));
  }
  if (n>0){
    F_el[0]-=k*(q[0]-*qleft[0]-s0*cos(*phi_neigh)*sin(*theta_neigh));
    F_el[1]-=k*(q[1]-*qleft[1]-s0*sin(*phi_neigh)*sin(*theta_neigh));
    F_el[2]-=k*(q[2]-*qleft[2]-s0*cos(*theta_neigh));
  }
};

void Bead::getScalarProd(){
  sp=0;
  if ((n>0)&&(n<N_beads-1)){
    for (i=0; i<3; i++){
      sp+=(*qleft[i]-q[i])*(*qright[i]-q[i]);
    }
  }
};

void Bead::getFKP(){
  for (i=0; i<3; i++){
    F_KP[i]=0;
  }
  for (i=0; i<3; i++){
    if (n>1){
      F_KP[i]=-J*((*qlleft[i]-*qleft[i])/(*r_fleft*(*r_left))-*sp_left*(q[i]-*qleft[i])/(*r_fleft*pow(*r_left,3)));
    }
    if ((n>0) && (n<N_beads-1)){
      F_KP[i]-=J*((2*q[i]-*qleft[i]-*qright[i])/(*r_left*r)+sp*((*qleft[i]-q[i])/(pow(*r_left,3)*r)+(*qright[i]-q[i])/(*r_left*pow(r,3))));
    }
    if (n<N_beads-2){
      F_KP[i]-=J*((*qrright[i]-*qright[i])/(r*(*r_right))-*sp_right*(q[i]-*qright[i])/(pow(r,3)*(*r_right)));
    }
  }
};

void Bead::getVel(){
  for (i=0; i<3; i++){
    vel[i]=(q_new[i]-q_old[i])/(2*dt);
  }
};

void Bead::getWindVel(){
  wind_vel[0]=T1.w*q[2]-h*T1.psi(q[0])*q[2];
  wind_vel[1]=0.0;
  wind_vel[2]=h*T1.d_psi(q[0])*pow(q[2],2)/2;
};

void Bead::getWeight(){
  for (i=0; i<3; i++){
    if ((i==2)&&(n==0)){
      F_weight[i]=-m*g;
    } else{
      F_weight[i]=0;
    }
  }
};

void Bead::getDrag(){
  for (i=0; i<3; i++){
    F_drag[i]=-b*(vel[i]-wind_vel[i]);
  }
};

void Bead::firstStepProp(){
  for (i=0; i<3; i++){
    q[i]=q_old[i]+vel[i]*dt+(pow(dt,2)/(2*m))*(b*(wind_vel[i]-vel[i])+F_el[i]+F_KP[i]+F_weight[i]);
  }
};

void Bead::oneStepProp(){
  for (i=0; i<3; i++){
    q_new[i]=(1/(m+b/2*dt))*(2*m*q[i]+(b/2*dt-m)*q_old[i]+pow(dt,2)*(F_el[i]+b*wind_vel[i]+F_KP[i]+F_weight[i]));
    q_old[i]=q[i];
    q[i]=q_new[i];
  }
};



class System: public Bead{
  public:
    int N_springs;
    int N_steps_therm;
    double dt;
    double g;
    Bead p[N_lines][N_beads];
    double z0;
    int i;
    int j;
    int k;
    int t;
    void getValues();
    void computeForces();
    void initialize();
    void thermalize();
    void evolve();
    void printCoord();
};

void System::getValues(){
  N_springs=N_beads-1;
  z0=2.0;
  T1.genSeries();

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      p[i][j].k=0.5;
      p[i][j].J=0.1;
      p[i][j].s0=0.1;
      p[i][j].dt=T1.dt;
      p[i][j].h=T1.f[0];

      if (j!=0){
        p[i][j].m=0.0005;
      } else{
        p[i][j].m=0.01;
      }

      if ((j==0)&&(i>0)){
        p[i][0].g=0.0;
        p[i][0].b=0.0;
      } else{
        p[i][j].b=0.1;
        p[i][j].g=9.81;
      }
    }
  }
};

void System::initialize(){
  getValues();

  std::default_random_engine generator(1234);
  std::normal_distribution<double> distribution(1.0, 0.1);
  std::normal_distribution<double> distribution2(0.0, 0.1);

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      if ((j==0) && (i>0)){
        for (k=0; k<3; k++){
          p[i][j].q_old[k]=p[0][0].q_old[k];
        }
      } else{
        double a=distribution(generator);
        double b=distribution2(generator);
        double c=distribution2(generator);
        p[i][j].q_old[0]=p[i][j].s0*j+p[i][j].s0*a;                       //initialize positions
        p[i][j].q_old[1]=1.0+p[i][j].s0*b;
        p[i][j].q_old[2]=z0+p[i][j].s0*c;
      }
      for (k=0; k<3; k++){
        p[i][j].q[k]=p[i][j].q_old[k];
        p[i][j].vel[k]=0;
      }
      p[i][j].n=j;                                //pass position of the bead
    }
    p[i][j].l=i;                                  //pass line of the bead
  }

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      for (k=0; k<3; k++){
        p[i][j].qleft[k]=&p[i][(j-1+N_beads)%N_beads].q[k];       //initialize pointers to neighbors with PBC
        p[i][j].qlleft[k]=&p[i][(j-2+N_beads)%N_beads].q[k];
        p[i][j].qright[k]=&p[i][(j+1)%N_beads].q[k];
        p[i][j].qrright[k]=&p[i][(j+2)%N_beads].q[k];
      }
    p[i][j].getDistances();
    p[i][j].getScalarProd();
    }
  }

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      p[i][j].r_left=&p[i][(j-1+N_beads)%N_beads].r;
      p[i][j].r_right=&p[i][(j+1)%N_beads].r;
      p[i][j].r_fleft=&p[i][(j-2)%N_beads].r;
      p[i][j].sp_left=&p[i][(j-1+N_beads)%N_beads].sp;
      p[i][j].sp_right=&p[i][(j+1)%N_beads].sp;
      p[i][j].getAngles();
    }
  }

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      p[i][j].theta_neigh=&p[i][(j-1+N_beads)%N_beads].theta;
      p[i][j].phi_neigh=&p[i][(j-1+N_beads)%N_beads].phi;
    }
  }

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      p[i][j].getWindVel();
      p[i][j].getWeight();
      p[i][j].getFKP();
      p[i][j].getFel();
    }
  }

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      if ((j==0) && (i>0)){
        for (k=0; k<3; k++){
          p[i][0].q_old[k]=p[i-1][0].q[k];
        }
      }
      p[i][j].firstStepProp();
    }
  }

  for (i=0; i<N_lines-1; i++){
    for (k=0; k<3; k++){
      p[i][0].q_old[k]=p[N_lines-1][0].q_old[k];
      p[i][0].q[k]=p[N_lines-1][0].q[k];
    }
  }
};

void System::thermalize(){
  N_steps_therm=100000;

  for (t=0; t<N_steps_therm; t++){
    computeForces();
    for (i=0; i<N_lines; i++){
      for (j=1; j<N_beads; j++){
        p[i][j].oneStepProp();
      }
    }
  }
};

void System::computeForces(){
  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      p[i][j].getDistances();
      p[i][j].getScalarProd();
      p[i][j].getAngles();
      p[i][j].getWindVel();
      p[i][j].getDrag();
      p[i][j].getWeight();
      p[i][j].getFKP();
      p[i][j].getFel();
    }
  }
};

void System::evolve(){
  int flag=0;
  initialize();
  thermalize();

  for (t=0; t<N_steps-N_steps_therm; t++){
    computeForces();
    for (i=0; i<N_lines; i++){
      for (j=0; j<N_beads; j++){
        p[i][j].h=T1.f[t];
        if ((j==0) && (i>0)){
          for (k=0; k<3; k++){
            p[i][0].q_old[k]=p[i-1][0].q_old[k];
            p[i][0].q[k]=p[i-1][0].q[k];
          }
        }
        p[i][j].oneStepProp();
        if (p[i][j].q[2]<0.0){
          flag=1;
        }
      }
    }
    for (i=0; i<N_lines-1; i++){
      for (k=0; k<3; k++){
        p[i][0].q_old[k]=p[N_lines-1][0].q_old[k];
        p[i][0].q[k]=p[N_lines-1][0].q[k];
      }
    }
    if (flag==1){
      cout << t << "\n";
      break;
    }
    for (i=0; i<N_lines; i++){
      for (j=0; j<N_beads; j++){
        p[i][j].getVel();
      }
    }
  }
};

void System::printCoord(){
  ofstream outdata;
  outdata.open("coord.txt");
  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      for (k=0; k<3; k++){
        outdata << p[i][j].q[k] << " ";
        cout << p[i][j].q[k] << " ";
      }
      outdata << "\n";
      cout << "\n";
    }
    cout << "\n";
    outdata << "\n\n";
  }
  outdata.close();
};

int main(){
  System S1;
  S1.evolve();
  //S1.initialize();
  //S1.thermalize();
  S1.printCoord();
  return 0;
};
