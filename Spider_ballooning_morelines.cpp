#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
using namespace std;
const int N_beads=15;
const int N_lines=3;




class Bead{
  public:
    double q_old[3];
    double q[3];
    double q_new[3];
    double vel[3];
    double wind_vel[3];
    double vel_in[3];
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
    double Kin_En;
    double Pot_En;
    double k;
    double J;
    double s0;
    double b;
    double w;
    double thr;
    double m;
    double g;
    double dt;
    double tfunc;
    void getDistances();
    void getAngles();
    void getFel();
    void getScalarProd();
    void getFKP();
    void getVel();
    void getWindVel();
    double psi(double x);
    double d_psi(double x);
    double f(double t_in);
    void getDrag();
    void getWeight();
    void oneStepProp();
    void firstStepProp();
    int i;
    int n;
    int l;
    int time;
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
  tfunc=time*dt;
  wind_vel[0]=w*q[2]-f(tfunc)*psi(q[0])*q[2];
  wind_vel[1]=0.0;
  wind_vel[2]=f(tfunc)*d_psi(q[0])*pow(q[2],2)/2;
};

double Bead::f(double t_in){
  return pow(sin(t_in),2);
};

double Bead::d_psi(double x){
  if ((x<thr) && (x>-thr)){
    return 1;
  } else{
    return 0;
  }
};

double Bead::psi(double x){
  if ((x<thr) && (x>-thr)){
    return x+thr;
  } else{
    if (x<-thr){
      return 0;
    } else{
      return 2*thr;
    }
  }
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
    int N_steps;
    double dt;
    double g;
    Bead p[N_lines][N_beads];
    double z0;
    double K_tot;
    double U_tot;
    int i;
    int j;
    int k;
    int t;
    void getValues();
    void computeForces();
    void initialize();
    void evolve();
    void printCoord();
    void printForces();
};

void System::getValues(){
  N_springs=N_beads-1;
  N_steps=10;
  z0=2.0;

  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      p[i][j].k=0.5;
      p[i][j].J=0.1;
      p[i][j].s0=0.1;
      p[i][j].dt=0.00001;
      p[i][j].time=0;
      p[i][j].w=1.0;

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
  std::normal_distribution<double> distribution2(1.0, 0.1);

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


void System::computeForces(){
  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      p[i][j].time=t+1;
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
}

void System::printForces(){
  double mod_KP;
  double mod_el;
  double mod_w;
  double mod_drag;
  ofstream outdata3;
  //outdata3.open("ModForces.txt");
  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      mod_KP=0;
      mod_el=0;
      mod_drag=0;
      mod_w=p[i][j].F_weight[2];
      for (k=0; k<3; k++){
        mod_el+=pow(p[i][j].F_el[k],2);
        mod_KP+=pow(p[i][j].F_KP[k],2);
        mod_drag+=pow(p[i][j].F_drag[k],2);
      }
      cout << t << " " << sqrt(mod_el) << " " << sqrt(mod_KP) << " " << sqrt(mod_drag) << " " << mod_w << " ";
      for (k=0; k<3; k++){
        cout << p[i][j].q_new[k] << " ";
      }
      cout << "\n";
    }
  }
  //outdata3.close();
};

void System::evolve(){
  ofstream outdata3;
  outdata3.open("ModForces.txt");
  int flag=0;
  initialize();

  for (t=0; t<N_steps; t++){
    computeForces();
    for (i=0; i<N_lines; i++){
      for (j=0; j<N_beads; j++){
        if ((j==0) && (i>0)){
          for (k=0; k<3; k++){
            p[i][0].q_old[k]=p[0][0].q_old[k];
            p[i][0].q[k]=p[0][0].q[k];
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

    //if (t%100000==0){
      //cout << "banane\n";
      printForces();
    //}
    if (flag==1){
      cout << t << "\n";
      break;
    }
    for (i=0; i<N_lines; i++){
      for (j=0; j<N_beads; j++){
        p[i][j].getVel();
        //cout << "dionudo\n";
      }
    }
  }
  outdata3.close();
};

void System::printCoord(){
  ofstream outdata;
  ofstream outdata2;
  outdata.open("coord.dat");
  outdata2.open("forces.txt");
  for (i=0; i<N_lines; i++){
    for (j=0; j<N_beads; j++){
      for (k=0; k<3; k++){
        outdata << p[i][j].q[k] << " ";
        outdata2 << p[i][j].F_el[k] << " ";
        cout << p[i][j].q[k] << " ";
      }
      outdata << "\n";
      outdata2 << "\n";
      cout << "\n";
    }
    cout << "\n";
    outdata << "\n\n";
  }
  outdata.close();
  outdata2.close();
};

int main(){
  System S1;
  //S1.initialize();
  S1.evolve();
  S1.printCoord();
  return 0;
};
