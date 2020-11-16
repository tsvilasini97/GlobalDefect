#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP
#include <sys/time.h>
#include<thread>
#include<chrono>
#include <fstream>
#include <bits/stdc++.h>

#include "LATfield2.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
  double phiavg_;
  double rhoavg_;

  double H0;
  double G;
  double omega_r;
  double omega_m;
  double omega_lambda;
  double a_i;
  double rho_c0;
  
  double fric_term_1;
  double fric_term_2;
  double Fric_term_1;
  double Fric_term_2;
  double dis_fric_term_1;
  double dis_fric_term_2;
  double t_dis;
  
  double *karr;
  double *pkarr;
  double *indicesarr;
  double *kbins_edges;
  double *kbins;
  double *PK;
//  double *kx;
//  double *ky;
//  double *kz;
  double *numberinbins;
  int k_arr_size;
  int binnos;
  
  
  Lattice lat_;
  Lattice klat_;
  Field<double> phi_;
  Field<double> pi_;
  Field<double> pi_prev_;
  Field<double> rho_;
  Field<double> P_;
  
  Field<Imag> rho_k_;
  
  PlanFFT<Imag> planrho_;
  
 

public:
  Global_defect(){;}
  Global_defect(string settings_filename)
  {
    initialize(settings_filename);
  }
  ~Global_defect()
  {
    delete[] karr;
    delete[] pkarr;
    delete[] numberinbins;
    delete[] PK;
    delete[] kbins;
    delete[] kbins_edges;
//    delete[] indicesarr;
  }
  unsigned long int random_seed();
  void initialize(string settings_filename);
  void loadSettings(string settings_filename);

  
  void generate_initCond();

  void evolve();
  double Friedmann_eq(double a);
  void next_cosmology(); //this fuction compute a and adot_over_a
  void field_leapfrog();
  double potential(Site & x);
  double potentialprime(Site & x, int comp);
  
  
  template<typename T>
  void averageField(Field<T> &f, string filename,int val);
  void averagephi();
  void averagerho();
  
  double modsqphi(Site & x);
  double modsqrho(Site & x);
  
  void compute_rho_P_();
  void compute_pk_();
  
  bool output_now();
  void output();



};

void Global_defect::initialize(string settings_filename)
{
  loadSettings(settings_filename);


  dx_ = physicalSize_/latticeSize_;
  dt_ = courantFactor_ * dx_;

  lat_.initialize(3,latticeSize_,1);
  klat_.initializeRealFFT(lat_,0);
  
  phi_.initialize(lat_,nComponents_);
  phi_.alloc();

  pi_.initialize(lat_,nComponents_);
  pi_.alloc();

  pi_prev_.initialize(lat_,nComponents_);
  pi_prev_.alloc();

  rho_.initialize(lat_);
  rho_.alloc();
  
  rho_k_.initialize(klat_);
  rho_k_.alloc();
  planrho_.initialize(&rho_,&rho_k_);
//  planrho_.alloc();
  
  k_arr_size = latticeSize_*latticeSize_*latticeSize_;
 
//  kx = new double [latticeSize_];
//  ky = new double [latticeSize_];
//  kz = new double [latticeSize_];
//  
  karr = new double [k_arr_size];
  pkarr = new double [k_arr_size];
//  indicesarr = new double [k_arr_size];
  kbins_edges = new double [binnos];
  kbins = new double [binnos-1];
  PK = new double [binnos-1];
  numberinbins = new double [binnos-1];
  
  
  P_.initialize(lat_);
  P_.alloc();
  
  rho_c0 = 3*H0*H0/(8*M_PI*G);

  generate_initCond();
  
  COUT<< "The initial conditions set are:"<<endl;
  COUT<< "start time is = "<<t_start_<<endl;
  COUT<< " End time is = "<<t_end_<<endl;
  COUT<< "Time interval is = "<<dt_<<endl;
  COUT<< "Lattice interval is = "<< dx_<<endl;
  COUT<< "Omega matter is = "<< omega_m<<endl;
  COUT<< "Omega radiation is = "<< omega_r<<endl;
  COUT<< "Omega lambda is = "<< omega_lambda<<endl;
  COUT<< "H0 is = "<< H0 <<endl;
  COUT<< "initial scale factor is = "<< a_i<<endl;
  COUT<< "The G is ="<<G<<endl;
  COUT<<"The current critical density is ="<<rho_c0<<endl;
  COUT<<"The val of pi is="<<M_PI<<endl;
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

  setfile.read("friction coefficient 1",Fric_term_1);
  setfile.read("friction coefficient 2",Fric_term_2);
  setfile.read("Dissipation time end",t_dis);
  setfile.read("Dissipation friction coefficient 1",dis_fric_term_1);
  setfile.read("Dissipation friction coefficient 2",dis_fric_term_2);
  
  setfile.read("omega_r",omega_r);
  setfile.read("omega_m",omega_m);
  setfile.read("omega_lambda",omega_lambda);
  setfile.read("H0",H0);
  setfile.read("a_i",a_i);
  setfile.read("G",G);
  setfile.read("bin_numbers",binnos);
  setfile.close();
  
}

unsigned long int Global_defect::random_seed()
  {
   struct timeval tv;
   gettimeofday(&tv,0);
   return(tv.tv_sec + tv.tv_usec);
  }
    
void Global_defect::generate_initCond()
{

  // TODO: understand how to change seed with gsl....
  Site x(lat_);

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
  
  gsl_rng_default_seed = random_seed();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for(x.first();x.test();x.next())
  {
    double phiNorm2 = 0;
    //theta = gsl_rng_uniform (r);
    //phi_(x,0) = eta2*sin(theta);
    //phi_(x,1) = eta2*cos(theta);
    for(int c=0;c<nComponents_;c++)
    {
      phi_(x,c) = gsl_ran_gaussian (r,1);
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
    averagephi();
    averagerho();
    
    ofstream phifile;
    phifile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_average_phi_1.txt",ios::trunc);
    phifile.close();
    
    ofstream rhofile;
    rhofile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_average_rho.txt",ios::trunc);
    rhofile.close();
    
    ofstream omegafile;
    omegafile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_omega_.txt",ios::trunc);
    omegafile.close();
  
  while(t_ <= t_end_)
  {
    t_ += dt_;
    next_cosmology();
    field_leapfrog();
    compute_rho_P_();
    averagephi();
    averagerho();
    
    if(output_now())
    {
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

//  if (t_<t_dis)
//  {
//    fric_term_1 = dis_fric_term_1;
//    fric_term_2 = dis_fric_term_2;
//  }
//  else
//  {
//    fric_term_1 = Fric_term_1;
//    fric_term_2 = Fric_term_2;
//  }

  fric_term_1 = Fric_term_1;
  fric_term_2 = Fric_term_2;
  double c1 = (1.0 - dt_ * (fric_term_1*adot_overa_ + fric_term_2)) / (1.0 + dt_ * (fric_term_1*adot_overa_ + fric_term_2));
  double c2 = dt_ / (1.0 + dt_ * adot_overa_);
  double a2 = a_*a_;

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
      pi_(x,c) = c1 * pi_prev_(x,c) + c2 * ( lapPhi -  a2 * potentialprime(x,c) );
    }
  }

}


void Global_defect::averagephi()
{
  Site x(lat_);
  double phisum_ = 0;
  for(x.first();x.test();x.next())
  {
    phisum_ += modsqphi(x);
  }
  parallel.sum(phisum_);
  phiavg_ = phisum_/pow(latticeSize_,3);
  
  if(parallel.rank() == 0)
	{
    ofstream phifile;
    phifile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_average_phi_1.txt",std::ios_base::app);
    phifile << t_<<" "<<phiavg_<<endl;
    phifile.close();
  }
  //averageField(phi_,"/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_average_phi.txt",0);
}

void Global_defect::averagerho()
{
  Site x(lat_);
  double rhosum_ = 0;
  for(x.first();x.test();x.next())
  {
    rhosum_ += rho_(x);
  }
  parallel.sum(rhosum_);
  rhoavg_ = rhosum_/pow(latticeSize_,3);
  if(parallel.rank() == 0)
	{
    ofstream rhofile;
    rhofile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_average_rho.txt",std::ios_base::app);
    rhofile << t_<<" "<<rhoavg_<<endl;
    rhofile.close();
  }
}

//template<typename T>
//void Global_defect::averageField(Field<T> &f, string filename,int val)
//{
//	Site x(f.lattice());
//	T sum = 0;
//	
//	if (val =0)
//	{
//	  for(x.first();x.test();x.next())
//    	{
//      	sum += modsqphi(x);
//    	}
//  }
//  else
//	{
//	  for(x.first();x.test();x.next())
//    	{
//      	sum += modsqrho(x);
//    	}
//  }
//  
//  parallel.sum(sum);
//  
//	T ave = sum/ (latticeSize_*latticeSize_*latticeSize_);
//	
//	if(parallel.rank() == 0)
//	{
//		ofstream file;
//		file.open(filename,std::ios_base::app);
//		file<< t_<<" "<<ave<<endl;
//		file.close();
//	}
//}

double Global_defect::potential(Site & x)
{
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_(x,i)*phi_(x,i);
  return lambda_ * ( phiNorm2 - eta2_) * ( phiNorm2 - eta2_) / 2.0;
}
double Global_defect::potentialprime(Site & x, int comp)
{
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_(x,i)*phi_(x,i);
  return 2.0 * lambda_ * ( phiNorm2 - eta2_) *  phi_(x,comp);
}


double Global_defect::modsqphi(Site &x)
{
  double phiNorm2 = 0;
  for(int i =0;i<nComponents_;i++)phiNorm2 += phi_(x,i)*phi_(x,i);
  return pow(phiNorm2,0.5);
}

//double Global_defect::modsqrho(Site &x)
//{
//  double phiNorm2 = 0;
//  for(int i =0;i<nComponents_;i++)phiNorm2 += rho_(x,i)*rho_(x,i);
//  return pow(phiNorm2,0.5);
//}

double Global_defect::Friedmann_eq(double a)
{
  return H0*pow((omega_r + omega_m*a + omega_lambda*(a*a*a*a)),0.5);
}

void Global_defect::next_cosmology()
{
  a_ = t_ * t_ / (t_end_*t_end_);
  adot_overa_  = 2 / t_;

//  double k1;
//  double k2;
//  double k3;
//  double k4;
//  
//  double rho_c;

//  k1 = dt_*Friedmann_eq(a_i);
//  k2 = dt_*Friedmann_eq(a_i + k1/2);
//  k3 = dt_*Friedmann_eq(a_i + k2/2);
//  k4 = dt_*Friedmann_eq(a_i + k3);
//  a_i += k1/6 + k2/3 +k3/3 +k4/6;
//  a_ = a_i;
//  adot_overa_ = Friedmann_eq(a_)/a_;

//  rho_c = 3*(Friedmann_eq(a_)*Friedmann_eq(a_)/(a_*a_*a_*a_))/(8*M_PI*G);
//  
//  double omegam=0;
//  double omegar=0;
//  double omegal=0;
//  
//  omegam = omega_m*rho_c0/ (a_*a_*a_*rho_c);
//  omegar = omega_r*rho_c0/ (a_*a_*a_*a_*rho_c);
//  omegal = 1 - omegam - omegar;
//  
//  parallel.sum(omegam);
//  parallel.sum(omegar);
//  parallel.sum(omegal);
//  
//  if(parallel.rank() == 0)
//	{
//	  ofstream omegafile;
//	  omegafile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_omega_.txt",std::ios_base::app);
//	  omegafile << a_<<" "<<omegam<<" "<<omegar<<" "<<omegal<<endl;
//	  omegafile.close();
//	}
}

bool Global_defect::output_now()
{
  return step%100==0?true:false;
}

void Global_defect::output()
{
  COUT<<"outputing field at t="<<t_<<endl;
  string filename_end= int2string(t_,99999)+".h5";
  

  //phi_.saveHDF5(runID_ + "_phi_" + filename_end);
  //pi_.saveHDF5(runID_ + "_pi_" + filename_end);

  //computeT();
  rho_.saveHDF5("/media/vilasini/DATA/UNIGE/Thesis/plots/" + runID_ + "_rho_" + filename_end);
//  rho_k_.saveHDF5(runID_ + "_rho_k_" + filename_end);
  compute_pk_();
}

void Global_defect::compute_rho_P_()
{
  Site x(lat_);

  double a2 = a_*a_;
  rhoavg_ = 0;
  for(x.first();x.test();x.next())
  {
    double mpidot = 0;
    double temp;
    double gradPhi2 = 0;
    //double phinorm2 = 0;
    for(int c=0;c<nComponents_;c++)
    {
      temp = (pi_prev_(x,c)+pi_(x,c))/2.0;
      mpidot = temp*temp;
      //phinorm2 = phi_(x,c)*phi_(x,c);
      for(int i = 0;i<3;i++)
      {
        temp = ( phi_(x+i,c) - phi_(x-i,c) ) / 2.0 / dx_;
        gradPhi2 += temp*temp;
      }
    }  
    rho_(x) = mpidot / 2.0 / a2 + potential(x) + gradPhi2  / 2.0 / a2;
    P_(x) = mpidot / 2.0 /a2 - potential(x) - gradPhi2 / 6.0 / a2;
  }
}

void Global_defect::compute_pk_()
{
  planrho_.execute(FFT_FORWARD); 
  rKSite k(klat_);
  
  cout<<"The k array size is:"<<k_arr_size<<endl;
  
  for(int i =0;i<k_arr_size;i++)
  {
    karr[i] = 0;
    pkarr[i] = 0;
  }

  double knrm;
  int index = 0;
  int index_negative = 0;
  
  string filename_end= int2string(t_,99999)+".txt";

  ofstream kfile;
  kfile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/powerspectrum/"+ runID_+"_k_"+filename_end,std::ios_base::app);
    
  for(k.first();k.test();k.next())
  {
    int kx;
    int ky;
    int kz;
    
    double k_x;
    double k_y;
    double k_z;
    
    kx = k.coord(0);
    ky = k.coord(1);
    kz = k.coord(2);
      
//    cout<<kx<<" "<<ky<<" "<<kz<<" "<<k<<" "<<endl;

    k_x = 2 *M_PI*kx*dx_/latticeSize_;
    
    if(ky<=latticeSize_/2)
    {
        k_y = 2*M_PI*ky*dx_/latticeSize_;
    }
    else
    {
        k_y = 2*M_PI*(latticeSize_ - ky)*dx_/latticeSize_;
    }
    
    if(kz<=latticeSize_/2)
    {
        k_z = 2*M_PI*kz*dx_/latticeSize_;
    }
    else
    {
        k_z = 2*M_PI*(latticeSize_ - kz)*dx_/latticeSize_;
    }
    
    knrm = pow(k_x*k_x + k_y*k_y + k_z*k_z,0.5);
    
    if(kx != 0)
    {
      index_negative = k_arr_size - index;
      karr[index_negative] = knrm;
//    cout<<karr[index]<<" "<<"knrm"<<endl;
      pkarr[index_negative] = (rho_k_(k)*rho_k_(k).conj()).real();
    }

    karr[index] = knrm;
//    cout<<karr[index]<<" "<<"knrm"<<endl;
    pkarr[index] = (rho_k_(k)*rho_k_(k).conj()).real();
    kfile<<kx<<" "<<ky<<" "<<kz<<" "<<knrm<<" "<<pkarr[index]<<endl;
    index = index + 1;
  }
  
  parallel.sum(karr,k_arr_size);
  cout<<index<<endl;
  kfile.close();
  
  if(parallel.rank() == 0)
  {
    double minimum_knrm = *min_element(karr, karr + k_arr_size);
    double maximum_knrm = *max_element(karr, karr + k_arr_size)+0.1;
    
    double delta_k = (maximum_knrm - minimum_knrm) / (binnos-1);
    double a =0;
//    cout<<"min"<<" "<<minimum_knrm<<" "<<"max"<<" "<<maximum_knrm<<endl;

    for(int i=0;i<binnos;i++)
    {
      kbins_edges[i] = 0;
    }

    for(int i=0;i<binnos-1;i++)
    {
      kbins[i] = 0;
      PK[i] = 0;
      numberinbins[i] = 0;
    }

    for(int i =0;i<binnos;i++)
    {
      kbins_edges[i] = minimum_knrm + a;
//      cout<<kbins_edges[i]<<" "<<"kedges"<<endl;
      a += delta_k;
    }

    for(int i =0;i<binnos-1;i++)
    {
      kbins[i] = (kbins_edges[i] + kbins_edges[i+1])/2.0;
//      cout<<kbins[i]<<" "<<"kbins"<<endl;
    }

    for(int i=0;i<k_arr_size;i++)
    {
      for(int j=0;j<binnos-1;j++)
      {
        if(kbins_edges[j] <= karr[i])
        {
          if(karr[i] < kbins_edges[j+1])
          {
//            cout<<kbins_edges[j] <<" "<< karr[i] <<" "<< kbins_edges[j+1]<<" "<<"ks"<<endl;
            PK[j] += pkarr[i];
            numberinbins[j] += 1;
          }
        }
      }
    }
    
    string filename_end= int2string(t_,99999)+".txt";

    ofstream pkfile;
    pkfile.open ("/media/vilasini/DATA/UNIGE/Thesis/plots/powerspectrum/"+ runID_+"_powerspectrum_"+filename_end,std::ios_base::app);
    
    for(int i=0;i<binnos;i++)
    {
//      cout<<PK[i]<<" "<<numberinbins[i]<<" "<<kbins[i]<<" "<<kbins_edges[i+1]<<" "<<pkarr[i]<<endl;
      PK[i] /= numberinbins[i];
      pkfile << kbins[i]<<" "<<PK[i]<<endl;
    }
    pkfile.close();
  }
  

}






/* The following loop computes the ks. Here the same algorithm used for fftfreq in numpy is used.*/

//  for(int i=0;i<latticeSize_/2;i++)
//  {
//    double pos_k;
//    double neg_k;
//    pos_k = i*dx_ / latticeSize_;
//    neg_k = -((latticeSize_/2) - i)*dx_ / latticeSize_;
//    kx[i] = pos_k;
//    ky[i] = pos_k;
//    kz[i] = pos_k;
//    kx[latticeSize_/2 + i] = neg_k;
//    ky[latticeSize_/2 + i] = neg_k;
//    kz[latticeSize_/2 + i] = neg_k;
//  }
//  
//  int ind = 0;
//  
//  for(int i=0;i<latticeSize_;i++)
//  {
//    for(int j=0;j<latticeSize_;j++)
//    {
//      for(int p=0;p<latticeSize_;p++)
//      {
//        knrm = pow(kx[i]*kx[i] + ky[j]*ky[j] + kz[p]*kz[p],0.5);
//        karr[ind] = knrm;
//        ind +=1;
//      }
//    }
//  }
//  

/* The below code computes the kbin edges computes the power spectrum by checking the bin in which knrm falls in. For now it is being computed in a single core*/

  
//  parallel.sum(karr,k_arr_size);
//  parallel.sum(indicesarr,k_arr_size);
//  
//  for(int i = 0;i<k_arr_size;i++)
//  {
//	  if (indicesarr[i] != 0)
//	  {
//	  	karr[i] /= indicesarr[i];
//	  	pkarr[i] /= indicesarr[i];
//	  }
//  }

//    if(kx>latticeSize_/2)
//    {
//        k_x = (-latticeSize_/2 + (kx % (latticeSize_/2)))*dx_ / latticeSize_;
//    }
//    if(kx<=latticeSize_/2)
//    {
//        k_x = kx*dx_/latticeSize_;
//    }
//    if(ky>latticeSize_/2)
//    {
//        k_y = (-latticeSize_/2 + (ky % (latticeSize_/2)))*dx_ / latticeSize_;
//    }
//    if(ky<=latticeSize_/2)
//    {
//        k_y = ky*dx_/latticeSize_;
//    }
//    if(kz>latticeSize_/2)
//    {
//        k_z = (-latticeSize_/2 + (kz % (latticeSize_/2)))*dx_ / latticeSize_;
//    }
//    if(kz<=latticeSize_)
//    {
//        k_z = kz*dx_/latticeSize_;
//    }

#endif
