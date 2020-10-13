#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP


#include "LATfield2.hpp"

using namespace LATfield2;

class Global_defect
{
private:
  int size_;
  int N_;
  double physicalSize_;
  double courantFactor_;
  double dt_;
  double dx_;


  Lattice lat_;
  Field<double> phi_;
  Field<double> pi_;

public:
  Global_defect(){;}
  Global_defect(int size, int N) : size_(size), N_(N)
  {
    initialize(size_,N_);
  }
  ~Global_defect()
  {

  }

  void initialize(int size, int N);

  void evolve();

  void evolve_step();

  bool output_now();
  void output();



};

void Global_defect::initialize(int size, int N)
{
  physicalSize_ = 1;
  courantFactor_ = 10 ;
  dx_ = physicalSize_ / N_;
  dt_ =  dx_ / courantFactor_ ;

  lat_.initialize(3,size_,1);
  phi_.initialize(lat_,N_);
  phi_.alloc();
  pi_.initialize(lat_,N_);
  pi_.alloc();

}

void Global_defect::evolve()
{
  for(int i = 0 ; i<10 ; i++)
  {
    evolve_step();
    if(output_now())
    {

    }
  }
}

void Global_defect::evolve_step()
{
  Site x(lat_);

  for(x.first();x.test();x.next())
  {
    for(int i = 0 ; i < N_; i++)
    {
      phi_(x,i) = phi_(x,i) + dt_ * pi_(x,i);
    }
  }


  for(x.first();x.test();x.next())
  {
    for(int i = 0 ; i < N_; i++)
    {
      pi_(x,i) =  phi_(x,i) + dt_*phi_(x,i);
    }
  }



}



bool Global_defect::output_now()
{
  return true;
}

void Global_defect::output()
{
  cout<<"process "<<parallel.rank()<< " : Oh we did a step!!!"<<endl;
}



#endif
