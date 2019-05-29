#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
using namespace std;
const int N_beads=21;

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
    double *r_fright;
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
    void firstStepProp();#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
using namespace std;
const int N_beads=21;

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
    double *r_fright;
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
    double m;
    double g;
    double dt;
    void getDistances();
    void getAngles();
    void getFel();
    void getScalarProd();
    void getFKP();
    void getVel();
    void getKinEn();
    void getPotEn();
    void getWindVel();
    void getDrag();
    void getWeight();
    void oneStepProp();
    void firstStepProp();
    int i;
    int n;
    double *qlleft[3];
    double *qleft[3];
    double *qright[3];
    double *qrright[3];
};

void Bead::getDistances(){
  r=0;
  for (i=0; i<3; i++){
    r+=pow((q[i]-*qleft[i]),2);
  }
  r=sqrt(r);
  //cout << r << "\n";
};

void Bead::getAngles(){
  getDistances();
  if (abs(*qright[2]-q[2])<*r_right){
    theta=acos((*qright[2]-q[2])/(*r_right));
  } else{
    theta=0.0;
  //theta=M_PI/2;
  }
  //cout << theta << "\n";
  //cout << *qright[2]-q[2] << " " << *r_right << "\n";
  //phi=atan2((*qright[1]-q[1]),(*qright[0]-q[0]));
  //phi=0.0;
};

void Bead::getFel(){
  getAngles();
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

  //cout << *sp_right << " " << (n+1)%N_beads << "\n";
};

void Bead::getFKP(){
  getDistances();
  getScalarProd();
  for (i=0; i<3; i++){
    F_KP[i]=0;
  }
  //cout << r << "\n";
  for (i=0; i<3; i++){
    if (n>1){
      F_KP[i]=-J*((*qlleft[i]-*qleft[i])/(*r_left*r)-*sp_left*(q[i]-*qleft[i])/(*r_left*pow(r,3)));
    }
    if ((n>0) && (n<N_beads-1)){
      F_KP[i]-=J*((2*q[i]-*qleft[i]-*qright[i])/(r*(*r_right))+sp*((*qleft[i]-q[i])/(pow(r,3)*(*r_right))+(*qright[i]-q[i])/(r*pow(*r_right,3))));
    }
    if (n<N_beads-2){
      F_KP[i]-=J*((*qrright[i]-*qright[i])/(*r_right*(*r_fright))-*sp_right*(q[i]-*qright[i])/(pow(*r_right,3)*(*r_fright)));
    }
  }

  //for (i=0; i<3; i++){
    //cout << theta << " " << n << "\n";
  //}
};

void Bead::getVel(){
  for (i=0; i<3; i++){
    vel[i]=(q_new[i]-q_old[i])/(2*dt);
  }
};

void Bead::getKinEn(){
  getVel();
  Kin_En=0;
  for (i=0; i<3; i++){
    Kin_En+=0.5*m*pow(vel[i],2);
  }
};

void Bead::getPotEn(){
  Pot_En=0.5*k*pow(r-s0,2);
};

void Bead::getWindVel(){
  wind_vel[0]=1.0;
  wind_vel[1]=0.0;
  wind_vel[2]=0.0;
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

void Bead::firstStepProp(){
  getWindVel();
  getWeight();
  getFKP();
  getFel();
  for (i=0; i<3; i++){
    q[i]=q_old[i]+vel_in[i]*dt+(pow(dt,2)/(2*m))*(b*(wind_vel[i]-vel_in[i])+F_el[i]+F_KP[i]+F_weight[i]);
  }
};

void Bead::oneStepProp(){
  getWindVel();
  getFKP();
  getFel();
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
    Bead p[N_beads];
    double z0;
    double K_tot;
    double U_tot;
    int i;
    int j;
    void getValues();
    void computeKinEn();
    void computePotEn();
    void initialize();
    void evolve();
    void printCoord();
};

void System::getValues(){
  N_springs=N_beads-1;
  N_steps=200000;
  z0=2.0;

  for (i=0; i<N_beads; i++){
    p[i].k=1;
    p[i].J=0.0;
    p[i].s0=0.1;
    p[i].b=0.0;

    if (i!=0){
      p[i].m=0.0005;
    } else{
      p[i].m=0.01;
    }

    p[i].g=0.0;
    p[i].dt=0.0001;
    for (j=0; j<3; j++){
      p[i].vel_in[j]=0;
    }
  }
};

void System::initialize(){
  getValues();
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(1.0, 0.5);
  std::normal_distribution<double> distribution2(0.0, 0.5);
  for (i=0; i<N_beads; i++){
    double a=distribution(generator);
    double b=distribution2(generator);
    double c=distribution2(generator);
    p[i].q_old[0]=p[i].s0*i+p[i].s0*a;                       //initialize positions
    p[i].q_old[1]=1.0+p[i].s0*b;
    p[i].q_old[2]=z0+p[i].s0*c;
    for (j=0; j<3; j++){
      p[i].q[j]=p[i].q_old[j];
      //cout << p[i].q_old[j] << " ";
    }
    p[i].n=i;
    //cout << "\n";
  }                                       //pass position of the bead

  for (i=0; i<N_beads; i++){
    for (j=0; j<3; j++){
      p[i].qleft[j]=&p[(i-1+N_beads)%N_beads].q_old[j];       //initialize pointers to neighbors with PBC
      p[i].qlleft[j]=&p[(i-2+N_beads)%N_beads].q_old[j];
      p[i].qright[j]=&p[(i+1)%N_beads].q[j];
      p[i].qrright[j]=&p[(i+2)%N_beads].q[j];
    }
    p[i].getDistances();
    p[i].getScalarProd();
  }

  for (i=0; i<N_beads; i++){
    p[i].r_left=&p[(i-1+N_beads)%N_beads].r;
    p[i].r_right=&p[(i+1)%N_beads].r;
    p[i].r_fright=&p[(i+2)%N_beads].r;
    p[i].sp_left=&p[(i-1+N_beads)%N_beads].sp;
    p[i].sp_right=&p[(i+1)%N_beads].sp;
    p[i].getAngles();
  }

  for (i=0; i<N_beads; i++){
    p[i].theta_neigh=&p[(i-1+N_beads)%N_beads].theta;
    p[i].phi_neigh=&p[(i-1+N_beads)%N_beads].phi;
  }

  for (i=0; i<N_beads; i++){
    p[i].firstStepProp();
  }
};

void System::computeKinEn(){
  K_tot=0;
  for (i=0; i<N_beads; i++){
    p[i].getKinEn();
    K_tot+=p[i].Kin_En;
  }
};

void System::computePotEn(){
  U_tot=0;
  for (i=1; i<N_beads; i++){
    p[i].getPotEn();
    U_tot+=p[i].Pot_En;
  }
};

void System::evolve(){
  ofstream outdata;
  outdata.open("en.txt");
  int flag=0;
  initialize();
  computeKinEn();
  computePotEn();

  for (j=0; j<N_steps; j++){
    for (i=0; i<N_beads; i++){
      p[i].oneStepProp();
      if (isnan(p[i].q[0])){
        flag=1;
      }
      //if ((i<N_beads-1) && (p[i].phi!=0)){
      //cout << j << " " << p[i].q[0] << "\n";
      //}
    }
    computeKinEn();
    computePotEn();
    outdata << j << " " << K_tot << " " << U_tot << " " << K_tot+U_tot << "\n";
    if (flag==1){
      cout << j << "\n";
      break;
    }
    //if (j>22500){
      //cout << j << " " << p[5].q[0] << "\n";
    //}
  }
};

void System::printCoord(){
  ofstream outdata;
  outdata.open("coord.txt");
  for (i=0; i<N_beads; i++){
    for (j=0; j<3; j++){
      outdata << p[i].q_old[j] << " ";
      cout << p[i].q[j] << " ";
    }
    outdata << "\n";
    cout << "\n";
  }
  outdata.close();
};

int main(){
  System S1;
  int i;
  int j;
  //S1.initialize();
  S1.evolve();
  S1.printCoord();
  return 0;
};

    int i;
    int n;
    double *qlleft[3];
    double *qleft[3];
    double *qright[3];
    double *qrright[3];
};

void Bead::getDistances(){
  r=0;
  for (i=0; i<3; i++){
    r+=pow((q[i]-*qleft[i]),2);
  }
  r=sqrt(r);
};

void Bead::getAngles(){
  getDistances();
  theta=acos((*qright[2]-q[2])/(*r_right));
  //phi=atan2((*qright[1]-q[1]),(*qright[0]-q[0]));
  phi=0.0;
};

void Bead::getFel(){
  getAngles();
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

  //cout << *sp_right << " " << (n+1)%N_beads << "\n";
};

void Bead::getFKP(){
  getDistances();
  getScalarProd();
  for (i=0; i<3; i++){
    F_KP[i]=0;
  }
  //cout << r << "\n";
  for (i=0; i<3; i++){
    if (n>1){
      F_KP[i]=-J*((*qlleft[i]-*qleft[i])/(*r_left*r)-*sp_left*(q[i]-*qleft[i])/(*r_left*pow(r,3)));
    }
    if ((n>0) && (n<N_beads-1)){
      F_KP[i]-=J*((2*q[i]-*qleft[i]-*qright[i])/(r*(*r_right))+sp*((*qleft[i]-q[i])/(pow(r,3)*(*r_right))+(*qright[i]-q[i])/(r*pow(*r_right,3))));
    }
    if (n<N_beads-2){
      F_KP[i]-=J*((*qrright[i]-*qright[i])/(*r_right*(*r_fright))-*sp_right*(q[i]-*qright[i])/(pow(*r_right,3)*(*r_fright)));
    }
  }

  //for (i=0; i<3; i++){
    //cout << theta << " " << n << "\n";
  //}
};

void Bead::getVel(){
  for (i=0; i<3; i++){
    vel[i]=(q_new[i]-q_old[i])/(2*dt);
  }
};

void Bead::getWindVel(){
  wind_vel[0]=1.0;
  wind_vel[1]=0.0;
  wind_vel[2]=0.0;
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

void Bead::firstStepProp(){
  getWindVel();
  getWeight();
  getFKP();
  getFel();
  for (i=0; i<3; i++){
    q[i]=q_old[i]+vel_in[i]*dt+(pow(dt,2)/(2*m))*(b*(wind_vel[i]-vel_in[i])+F_el[i]+F_KP[i]+F_weight[i]);
  }
};

void Bead::oneStepProp(){
  getWindVel();
  getFKP();
  getFel();
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
    Bead p[N_beads];
    double z0;
    int i;
    int j;
    void getValues();
    void initialize();
    void evolve();
    void printCoord();
};

void System::getValues(){
  N_springs=N_beads-1;
  N_steps=500000;
  z0=2.0;

  for (i=0; i<N_beads; i++){
    p[i].k=1.0;
    p[i].J=0.1;
    p[i].s0=0.1;
    p[i].b=0.5;

    if (i!=0){
      p[i].m=0.00005;
    } else{
      p[i].m=0.01;
    }

    p[i].g=9.81;
    p[i].dt=0.0001;
    for (j=0; j<3; j++){
      p[i].vel_in[j]=0;
    }
  }
};

void System::initialize(){
  getValues();
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(1.0, 0.5);
  std::normal_distribution<double> distribution2(0.0, 0.5);
  for (i=0; i<N_beads; i++){
    double a=distribution(generator);
    double b=distribution2(generator);
    double c=distribution2(generator);
    p[i].q_old[0]=p[i].s0*i+p[i].s0*a;                       //initialize positions
    p[i].q_old[1]=1.0+p[i].s0*b;
    p[i].q_old[2]=z0+p[i].s0*c;
    for (j=0; j<3; j++){
      p[i].q[j]=p[i].q_old[j];
      cout << p[i].q_old[j] << " ";
    }
    p[i].n=i;
    cout << "\n";
  }                                       //pass position of the bead

  for (i=0; i<N_beads; i++){
    for (j=0; j<3; j++){
      p[i].qleft[j]=&p[(i-1+N_beads)%N_beads].q_old[j];       //initialize pointers to neighbors with PBC
      p[i].qlleft[j]=&p[(i-2+N_beads)%N_beads].q_old[j];
      p[i].qright[j]=&p[(i+1)%N_beads].q[j];
      p[i].qrright[j]=&p[(i+2)%N_beads].q[j];
    }
    p[i].getDistances();
    p[i].getScalarProd();
  }

  for (i=0; i<N_beads; i++){
    p[i].r_left=&p[(i-1+N_beads)%N_beads].r;
    p[i].r_right=&p[(i+1)%N_beads].r;
    p[i].r_fright=&p[(i+2)%N_beads].r;
    p[i].sp_left=&p[(i-1+N_beads)%N_beads].sp;
    p[i].sp_right=&p[(i+1)%N_beads].sp;
    p[i].getAngles();
  }

  for (i=0; i<N_beads; i++){
    p[i].theta_neigh=&p[(i-1+N_beads)%N_beads].theta;
    p[i].phi_neigh=&p[(i-1+N_beads)%N_beads].phi;
  }

  for (i=0; i<N_beads; i++){
    p[i].firstStepProp();
  }
};

void System::evolve(){
  initialize();
  for (j=0; j<N_steps; j++){
    for (i=0; i<N_beads; i++){
      p[i].oneStepProp();
      if (p[i].q[2]<0){
        break;
      }
      //if ((i<N_beads-1) && (p[i].phi!=0)){
      //cout << j << " " << p[i].q[0] << "\n";
      //}
    }
    //if (j>22500){
      //cout << j << " " << p[5].q[0] << "\n";
    //}
  }
};

void System::printCoord(){
  ofstream outdata;
  outdata.open("coord.txt");
  for (i=0; i<N_beads; i++){
    for (j=0; j<3; j++){
      outdata << p[i].q[j] << " ";
      cout << p[i].q[j] << " ";
    }
    outdata << "\n";
    cout << "\n";
  }
  outdata.close();
};

int main(){
  System S1;
  int i;
  int j;
  S1.initialize();
  //S1.evolve();
  //S1.printCoord();
  return 0;
};
