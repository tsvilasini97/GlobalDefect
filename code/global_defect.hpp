#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP


#include "LATfield2.hpp"
#include <gsl/gsl_rng.h>

using namespace LATfield2;

class Global_defect
{
private:
  string runID_;
  int latticeSize_;
  int nComponents_;
  double physicalSize_;
  double courantFactor_;
  double dt_;
  double dx_;
  double t_start_;
  double t_end_;
  double t_;

  double lambda_;
  double eta2_;

  int step;


  double a_;
  double adot_overa_;

  Lattice lat_;
  Field<double> phi_;
  Field<double> pi_;
  Field<double> pi_prev_;
  Field<double> rho_;

public:
  Global_defect(){;}
  Global_defect(string settings_filename)
  {
    initialize(settings_filename);
  }
  ~Global_defect()
  {

  }

  void initialize(string settings_filename);
  void loadSettings(string settings_filename);


  void generate_initCond();

  void evolve();

  void next_cosmology(); //this fuction compute a and adot_over_a
  void field_leapfrog();

  double potential(Site & x);
  double potentialprime(Site & x, int comp);


  void computeT();


  bool output_now();
  void output();



};

void Global_defect::initialize(string settings_filename)
{
  loadSettings(settings_filename);


  dx_ = physicalSize_/latticeSize_;
  dt_ = courantFactor_ * dx_;

  lat_.initialize(3,latticeSize_,1);
  phi_.initialize(lat_,nComponents_);
  phi_.alloc();

  pi_.initialize(lat_,nComponents_);
  pi_.alloc();

  pi_prev_.initialize(lat_,nComponents_);
  pi_prev_.alloc();

  rho_.initialize(lat_);
  rho_.alloc();

  generate_initCond();
}


void Global_defect::loadSettings(string settings_filename)
{

  COUT<< "loadSettings: reading settings from : "<<settings_filename<<endl;


  SettingsFile setfile;
  setfile.open(settings_filename, SettingsFile::autoCreate);

  setfile.read("runID",runID_);
  setfile.read("nComponents",nComponents_);
  setfile.read("latticeSize",latticeSize_);
  setfile.read("physicalSize",physicalSize_);
  setfile.read("courantFactor",courantFactor_);
  setfile.read("t_start",t_start_);
  setfile.read("t_end",t_end_);

  setfile.read("lambda",lambda_);
  setfile.read("eta2",eta2_);

  setfile.close();
}

void Global_defect::generate_initCond()
{

  // TODO: understand how to change seed with gsl....
  Site x(lat_);


  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for(x.first();x.test();x.next())
  {
    double phiNorm2 = 0;
    for(int c=0;c<nComponents_;c++)
    {
      phi_(x,c) = gsl_rng_uniform (r);
      phiNorm2 += phi_(x,c)*phi_(x,c);
    }
    double ratio =  sqrt(eta2_/phiNorm2);

    for(int c=0;c<nComponents_;c++)
    {
      phi_(x,c) *= ratio;
      pi_(x,c) = 0;
    }
  }

  gsl_rng_free (r);

  phi_.saveHDF5(runID_ + "_phi_initCond.h5");
  pi_.saveHDF5(runID_ + "_pi_initCond.h5");

  COUT<< "Initial Condition generated"<<endl;

}

void Global_defect::evolve()
{
    COUT<< "Starting main loop"<<endl;
    step = 0;
    t_ = t_start_;

  while(t_ <= t_end_)
  {
    next_cosmology();
    field_leapfrog();
    t_ += dt_;

    if(output_now()){
      output();
    }

    step++;

  }
}

void Global_defect::field_leapfrog()
{
  Site x(lat_);

  for(x.first();x.test();x.next())
  {
    for(int c = 0;c<nComponents_;c++)
    {
      phi_(x,c) += dt_ * pi_(x,c);
    }
  }

  phi_.updateHalo(); //update the value of phi in the halo

  //
  double c1 = (1.0 - dt_ * adot_overa_) / (1.0 + dt_ * adot_overa_);
  double c2 = dt_ / (1.0 + dt_ * adot_overa_);
  double a2 = a_*a_;
  //


  //cout<<dt_<<" "<<adot_overa_<<" "<<a_<<endl;
  //cout<<c1<<" "<<c2<<" "<<a2<<endl;



  // put what is in pi in pi_prev
  //.... we then switch the data between pi and pi_prev:
  double * temp = pi_prev_.data_;
  pi_prev_.data_ = pi_.data_;
  pi_.data_ = temp;

  for(x.first();x.test();x.next())
  {
    for(int c = 0;c<nComponents_;c++)
    {
      double lapPhi = -6.0 * phi_(x,c) ;
      for(int i = 0 ; i<3 ; i++)lapPhi += phi_(x+i,c) + phi_(x-i,c);
      lapPhi /= dx_*dx_;
      pi_(x,c) = c1 * pi_prev_(x,c) + c2 * ( lapPhi +  a2 * potentialprime(x,c) );
    }
  }


}


double Global_defect::potential(Site & x)
{
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_(x,i)*phi_(x,i);
  return lambda_ * ( phiNorm2 - eta2_) * ( phiNorm2 - eta2_) / 2.0;
}
double Global_defect::potentialprime(Site & x, int comp){
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_(x,i)*phi_(x,i);
  return 2.0 * lambda_ * ( phiNorm2 - eta2_) *  phi_(x,comp);
}


void Global_defect::next_cosmology()
{
  a_ = t_ * t_ / (t_end_*t_end_);
  adot_overa_  = 2 / t_;
}

bool Global_defect::output_now()
{
  return step%20==0?true:false;
}

void Global_defect::output()
{
  COUT<<"outputing field at t="<<t_<<endl;
  string filename_end= int2string(step,99999)+".h5";

  phi_.saveHDF5(runID_ + "_phi_" + filename_end);
  pi_.saveHDF5(runID_ + "_pi_" + filename_end);

  computeT();
  rho_.saveHDF5(runID_ + "_rho_" + filename_end);

}

void Global_defect::computeT()
{
  Site x(lat_);


  double a2 = a_*a_;
  for(x.first();x.test();x.next())
  {
    double mpidot = 0;
    double temp;
    double gradPhi2 = 0;

    for(int c=0;c<nComponents_;c++)
    {
      temp = (pi_prev_(x,c)+pi_(x,c))/2.0;
      mpidot = temp*temp;

      for(int i = 0;i<3;i++)
      {
        temp = ( phi_(x+i,c) - phi_(x-i,c) ) / 2.0 / dx_;
        gradPhi2 += temp*temp;
      }
    }


    rho_(x) = mpidot / 2.0 + potential(x) + gradPhi2  / 2.0 / a2;
  }
}


#endif
