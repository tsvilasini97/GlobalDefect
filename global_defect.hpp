#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP

#include <math.h>
#include <sys/time.h>
#include<thread>
#include<chrono>
#include <fstream>
#include <bits/stdc++.h>

#include "LATfield2.hpp"
using namespace LATfield2;
#include "powerSpectra.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <bits/stdc++.h> 
#include <sys/stat.h> 
#include <sys/types.h> 

class Global_defect
{
private:
	string runID_;
	string path_;

	int latticeSize_;
	int nComponents_;
	float physicalSize_;
	float courantFactor_;
	double dt_;
	float dx_;
	float t_start_;
	float t_end_;
	double t_;

	float lambda_;
	float lambda0_;
	float eta2_;

	int step;

	double a_;
	double adot_overa_;
	double phiavg_;
	double rhoavg_;
	string univ;
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

	char *snap_times;
	vector<string> snap_times_;
	vector<float> Snap_times_;
	string ic;

//  double *karr;
//  double *pkarr;
//  double *indicesarr;
//  double *kbins_edges;
//  double *kbins;
//  double *PK;
//  double *kx;
//  double *ky;
//  double *kz;
//  double *numberinbins;
//  int k_arr_size;
	int binnos;

	Lattice lat_;
	Lattice klat_;

	Field<float> phi_defect_;
	Field<float> pi_defect_;

	Field<float> pi_defect_prev_;

	Field<double> rho_;
	Field<float> P_;

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
	//    delete[] karr;
	//    delete[] pkarr;
	//    delete[] numberinbins;
	//    delete[] PK;
	//    delete[] kbins;
	//    delete[] kbins_edges;
	//    delete[] indicesarr;
	}
	unsigned long int random_seed();
	void initialize(string settings_filename);
	void loadSettings(string settings_filename);

	void generate_initCond();
	void create_directory();

	void evolve();
	double Friedmann_eq(float a);
	void next_cosmology(float t); //this fuction compute a and adot_over_a

	void field_leapfrog();
	void update_phi();
	void update_pi();

	double potential(Site & x);
	double potentialprime(Site & x, int comp);
	double potentialprime_t(Site & x, int comp);
	
	template<typename T>
	void averageField(Field<T> &f, string filename,int val);
	void averagephidefect();
	void averagerhodefect();

	double modsqphi(Site & x);
	double modsqrho(Site & x);

	void compute_rho_P_();
	void compute_pk_();

	bool output_now();
	void output();
};


void Global_defect::initialize(string settings_filename)
{
	snap_times = new char[100];

	loadSettings(settings_filename);

	dx_ = physicalSize_/latticeSize_;
	dt_ = courantFactor_ * dx_;
//	dt_ = 0.1000000;
//	dt_ = ceil(dt_ * 10.0) / 10.0;
	lat_.initialize(3,latticeSize_,1);
	klat_.initializeRealFFT(lat_,0);

	phi_defect_.initialize(lat_,nComponents_);
	phi_defect_.alloc();

	pi_defect_.initialize(lat_,nComponents_);
	pi_defect_.alloc();

	pi_defect_prev_.initialize(lat_,nComponents_);
	pi_defect_prev_.alloc();

	rho_.initialize(lat_);
	rho_.alloc();

	rho_k_.initialize(klat_);
	rho_k_.alloc();

	//  planrho_.initialize(&rho_,&rho_k_);
	//  planrho_.alloc();

	//  k_arr_size = latticeSize_*latticeSize_*latticeSize_;

	//  kx = new double [latticeSize_];
	//  ky = new double [latticeSize_];
	//  kz = new double [latticeSize_];
	//  karr = new double [k_arr_size];
	//  pkarr = new double [k_arr_size];
	//  indicesarr = new double [k_arr_size];
	//  kbins_edges = new double [binnos];
	//  kbins = new double [binnos-1];
	//  PK = new double [binnos-1];
	//  numberinbins = new double [binnos-1];

	P_.initialize(lat_);
	P_.alloc();

	rho_c0 = 3*H0*H0/(8*M_PI*G);

	generate_initCond();

	COUT<<"Saving the files to:"<<path_<<endl<<endl;
	COUT<< "The initial conditions set are:"<<endl;
	COUT<< "start time is = "<<t_start_<<endl;
	COUT<< " End time is = "<<t_end_<<endl;
	COUT<< "Time interval is = "<<dt_<<endl;
	COUT<< "Lattice interval is = "<< dx_<<endl<<endl;
	COUT<< "The value set for univ is: "<<univ<<endl;
	COUT<< "The output times are: " << snap_times << endl;

	stringstream ss(snap_times);
	int i = 0;
	while (ss.good())
	{
		string substr;
		getline(ss, substr, ',');
		snap_times_.push_back(substr);
	}
	std::transform(snap_times_.begin(), snap_times_.end(), back_inserter(Snap_times_), [](const string & astr){ return stod( astr) ; } ) ;

//  COUT<< "Omega matter is = "<< omega_m<<endl;
//  COUT<< "Omega radiation is = "<< omega_r<<endl;
//  COUT<< "Omega lambda is = "<< omega_lambda<<endl;
//  COUT<< "H0 is = "<< H0 <<endl;
//  COUT<< "initial scale factor is = "<< a_i<<endl;
//  COUT<< "The G is ="<<G<<endl;
//  COUT<<"The current critical density is ="<<rho_c0<<endl;
//  COUT<<"The val of pi is="<<M_PI<<endl<<endl;
}


void Global_defect::loadSettings(string settings_filename)
{

	COUT<< "loadSettings: reading settings from : "<<settings_filename<<endl<<endl;


	SettingsFile setfile;
	setfile.open(settings_filename, SettingsFile::autoCreate);

	setfile.read("runID",runID_);
	setfile.read("path",path_);
	setfile.read("nComponents",nComponents_);
	setfile.read("latticeSize",latticeSize_);
	setfile.read("physicalSize",physicalSize_);
	setfile.read("courantFactor",courantFactor_);
	setfile.read("t_start",t_start_);
	setfile.read("t_end",t_end_);
	setfile.read("lambda0",lambda0_);
	setfile.read("eta2",eta2_);
	setfile.read("friction coefficient 1",Fric_term_1);
	setfile.read("Dissipation time end",t_dis);
	setfile.read("Dissipation friction coefficient 1",dis_fric_term_1);
	setfile.read("univ",univ);
	setfile.read("omega_r",omega_r);
	setfile.read("omega_m",omega_m);
	setfile.read("omega_lambda",omega_lambda);
	setfile.read("H0",H0);
	setfile.read("a_i",a_i);
	setfile.read("G",G);
	setfile.read("bin_numbers",binnos);
	setfile.read("snap_times", snap_times);
	setfile.read("ic_gen",ic);
	setfile.close();

	if(parallel.rank() == 0)
	{ 
	if (mkdir(path_.c_str (), 0777) == -1) 
	cerr << "There is an error in creating directory! The directory already exists:  " << strerror(errno) << endl; 
	else
	cout << "Directory for" <<" "<<runID_<<" created"<<endl<<endl ;
	}

	ofstream settingsfile;
	settingsfile.open (path_  + runID_+"_settingsfile.txt",ios::trunc);
	settingsfile << "runID = " <<runID_ <<endl << "nComponents =" << nComponents_ <<endl << "Lattice Size ="<<latticeSize_ <<endl<<"Physical Size ="<<physicalSize_<<endl<<"Courant Factor ="<< courantFactor_ <<endl<<"t_start ="<<t_start_<<endl<<"t_end ="<<t_end_<<endl<<"lambda0 ="<<lambda0_<<endl<<"eta2 ="<<eta2_<<endl<<"omega_r ="<<omega_r<<endl<<"omega_m ="<<omega_m<<endl<<"omega_r ="<<omega_r<<endl<<"omega_lambda ="<<omega_lambda<<endl<<"bin_numbers ="<<binnos<<"ic_gen"<<ic<<endl;
	settingsfile.close();
}


unsigned long int Global_defect::random_seed()
{
	struct timeval tv;
	gettimeofday(&tv,0);
	return(tv.tv_sec + tv.tv_usec);
}


void Global_defect::generate_initCond()
{
  // Todo: understand how to change seed with gsl....
	if(ic == "default_gen")
	{
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
			//phi_defect_(x,0) = eta2*sin(theta);
			//phi_defect_(x,1) = eta2*cos(theta);
			for(int c=0;c<nComponents_;c++)
			{
				phi_defect_(x,c) = gsl_ran_gaussian (r,1);
				phiNorm2 += phi_defect_(x,c)*phi_defect_(x,c);
			}
			double ratio =  sqrt(eta2_/phiNorm2);
			for(int c=0;c<nComponents_;c++)
			{
				phi_defect_(x,c) *= ratio;
				pi_defect_(x,c) = 0;
			}
		}
		gsl_rng_free (r);
	}
	else if(ic == "read_file")
	{
		pi_defect_.loadHDF5("fieldPi_test_t10.h5");
		phi_defect_.loadHDF5("fieldPhi_test_t10.h5");
		phi_defect_.updateHalo();
	}

	phi_defect_.saveHDF5(path_ + runID_ + "_phi_defect_initCond.h5");
	pi_defect_.saveHDF5(path_  + runID_ + "_pi_defect_initCond.h5");
	COUT<< "Initial Condition generated"<<endl<<endl;
}

void Global_defect::evolve()
{
	COUT<< "Starting main loop"<<endl;
	step = 0;
	t_ = t_start_;

//	averagephidefect();
//	averagerhodefect();
//	ofstream phifile;
//	phifile.open (path_  + runID_+"_average_phi_defect.txt",ios::trunc);
//	phifile.close();
//	ofstream rhofile;
//	rhofile.open (path_  + runID_+"_average_rho_defect.txt",ios::trunc);
//	rhofile.close();

	next_cosmology(t_);
	compute_rho_P_();
	rho_.saveHDF5(path_ + runID_ + "_rho_defect_initCond.h5");

	while(t_ <= t_end_)
	{
		field_leapfrog();
//		averagephidefect();
//		averagerhodefect();
		COUT << setprecision(8) << "Current time is at: " << t_ << " and a is: " << a_ << endl;
		if(output_now())
		{
			compute_rho_P_();
			output();
		}
		step++;
	}
}


void Global_defect::field_leapfrog()
{
	update_phi();
	t_ += dt_;
	next_cosmology(t_);
	update_pi();
//	t_ += dt_/2.0;
//	t_ += dt_/2.0;
//	next_cosmology(t_);
}

void Global_defect::update_phi()
{
	Site x(lat_);
	for(x.first();x.test();x.next())
	{
		for(int c = 0;c<nComponents_;c++)
		{
			phi_defect_(x,c) += dt_ * pi_defect_(x,c);
		}
	}
	phi_defect_.updateHalo(); //update the value of phi in the halo
}


void Global_defect::update_pi()
{
	Site x(lat_);
	double c1 = (1.0 - dt_ * adot_overa_)/(1.0 + dt_ * adot_overa_);
	double c2 = dt_ / (1.0 + dt_ * adot_overa_);
	double a2 = a_*a_;

  // put what is in pi in pi_defect_prev
  //.we then switch the data between pi and pi_defect_prev:
//	float * temp = pi_defect_prev_.data_;
//	pi_defect_prev_.data_ = pi_defect_.data_;
//	pi_defect_.data_ = temp;

	for(x.first();x.test();x.next())
	{
		double lapPhi = 0;
		for(int c = 0;c<nComponents_;c++)
		{
			lapPhi = -6.0 * phi_defect_(x,c);
			for(int i = 0 ; i<3 ; i++)lapPhi += (phi_defect_(x+i,c) + phi_defect_(x-i,c));
			lapPhi = lapPhi / dx_ / dx_;
			pi_defect_(x,c) = c1 * pi_defect_(x,c) + c2 * ( lapPhi -  a2 * potentialprime(x,c) );
		}
	}
}


double Global_defect::potentialprime(Site & x, int comp)
{
	double phiNorm2 = 0;
	for(int i =0;i<nComponents_;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
	return 1.0 * lambda_ * ( phiNorm2 - eta2_) *  phi_defect_(x,comp);
}

double Global_defect::modsqphi(Site &x)
{
	double phiNorm2 = 0;
	for(int i =0;i<nComponents_;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
	return pow(phiNorm2,0.5);
}

void Global_defect::averagephidefect()
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
		phifile.open (path_ + runID_+"_average_phi_defect.txt",std::ios_base::app);
		phifile << t_<<" "<<phiavg_<<endl;
		phifile.close();
	}
}

void Global_defect::averagerhodefect()
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
		rhofile.open (path_ + runID_ + "_average_rho_defect.txt",std::ios_base::app);
		rhofile << t_<<" "<<rhoavg_<<endl;
		rhofile.close();
	}
}

double Global_defect::potential(Site & x)
{
	double phiNorm2 = 0;
	for(int i =0;i<nComponents_;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
	return lambda_ * ( phiNorm2 - eta2_) * ( phiNorm2 - eta2_) / 4.0;
}

double Global_defect::Friedmann_eq(float a)
{
	return H0*pow((omega_r + omega_m*a + omega_lambda*(a*a*a*a)),0.5);
}

void Global_defect::next_cosmology(float t)
{
	if(univ== "matter")
	{
		a_ = t * t / (t_end_*t_end_);
		adot_overa_  = 2 / t;
	}
	else if(univ == "radiation")
	{
		a_ = t/t_end_;
		adot_overa_ = 1/t;
		lambda_ = lambda0_ / a_ / a_;
	}
	else if(univ=="lcdm")
	{
		double k1;
		double k2;
		double k3;
		double k4;

		double rho_c;

		k1 = dt_*Friedmann_eq(a_i);
		k2 = dt_*Friedmann_eq(a_i + k1/2);
		k3 = dt_*Friedmann_eq(a_i + k2/2);
		k4 = dt_*Friedmann_eq(a_i + k3);
		a_i += k1/6 + k2/3 +k3/3 +k4/6;
		a_ = a_i;
		adot_overa_ = Friedmann_eq(a_)/a_;

		rho_c = 3*(Friedmann_eq(a_)*Friedmann_eq(a_)/(a_*a_*a_*a_))/(8*M_PI*G);

		double omegam=0;
		double omegar=0;
		double omegal=0;

		omegam = omega_m*rho_c0/ (a_*a_*a_*rho_c);
		omegar = omega_r*rho_c0/ (a_*a_*a_*a_*rho_c);
		omegal = 1 - omegam - omegar;

		parallel.sum(omegam);
		parallel.sum(omegar);
		parallel.sum(omegal);

		if(parallel.rank() == 0)
		{
			ofstream omegafile;
			omegafile.open (path_ + runID_+"_omega_.txt",std::ios_base::app);
			omegafile << a_<<" "<<omegam<<" "<<omegar<<" "<<omegal<<endl;
			omegafile.close();
		}
	}
}

bool Global_defect::output_now()
{
	//  return step%100==0?true:false;
	int val = 1;
	double temp = t_;
	for(int i =0; i<snap_times_.size(); i++)
	{
		double temp0 = Snap_times_[i];
//		double temp1 = t_;
		temp0 = ceil(Snap_times_[i] * 100.0) / 100.0;
//		temp1 = ceil(t_ * 100.0) / 100.0;
		if(temp - dt_/2.0 < temp0 && temp0 < temp + dt_/2.0 )
		{
			COUT << Snap_times_[i] << " " << temp << endl;
			val = 2;
		}
	}

	return val%2==0?true:false;
}


////////////////////////////////////////////////
//  output
//
////////////////////////////////////////////////
// This function writws the output 
//
///////////////////////////////////////////////

void Global_defect::output()
{
	COUT<<"outputing field at t="<<t_<<endl;
	string filename_end= int2string(t_,99999)+".h5";

	phi_defect_.saveHDF5(path_ + runID_ + "_phi_defect_" + filename_end);
	pi_defect_.saveHDF5(path_ + runID_ + "_pi_defect_" + filename_end);
	rho_.saveHDF5(path_ + runID_ + "_rho_defect_" + filename_end);
	//  compute_pk_();
}

////////////////////////////////////////////////
//  compute_rho_P
//
////////////////////////////////////////////////
// This function computes the T00 and Tii 
//
///////////////////////////////////////////////

void Global_defect::compute_rho_P_()
{
	Site x(lat_);
	double a2 = a_*a_;
	for(x.first();x.test();x.next())
	{
		double mpidot = 0;
		double temp;
		double gradPhi2 = 0;
		double phiNorm2 = 0;
		for(int c=0;c<nComponents_;c++)
		{
			mpidot += 0.5*pi_defect_(x,c)*pi_defect_(x,c);
			phiNorm2 += phi_defect_(x,c)*phi_defect_(x,c);

			for(int i = 0;i<3;i++)
			{
				temp = ( phi_defect_(x+i,c) - phi_defect_(x,c) ) / dx_;
				gradPhi2 += 0.5*temp*temp ;
			}
		}
		temp = mpidot - gradPhi2;
		temp -= a2 * lambda_ * (phiNorm2 - eta2_)*(phiNorm2 - eta2_) / 4.0; 
		rho_(x) = 2.0 * mpidot - temp;
		if(rho_(x)>50)
		{
		COUT << a2 * lambda_ * (phiNorm2 - eta2_)*(phiNorm2 - eta2_) / 4.0 << " " << mpidot << " " << gradPhi2 << endl;
		}
//		rho_(x) *= 2;
//		COUT << mpidot /2.0 << " " << gradPhi2 / 2.0 << endl;
//		rho_(x) = mpidot/2.0  + a2*potential(x) + gradPhi2 / 2.0 ;
//		rho_(x) *= a2;
//		P_(x) = mpidot / 2.0  - potential(x) - gradPhi2 / 6.0 ;
	}
}


////////////////////////////////////////////////
//  compute_pk
//
////////////////////////////////////////////////
// This function computes the powerspectrum 
//
///////////////////////////////////////////////
//  Parameters:
//
//    rho_k_ = the field used to compute the pk
//    binnos = contains the no of bins in k: taken from settings file
//    physicalSize_ = the box size: taken from settings file 
///////////////////////////////////////////////

//void Global_defect::compute_pk_()
//{
//  planrho_.execute(FFT_FORWARD); 
//  string filename_end= path_ + runID_ + "_powerspectrum" + int2string(t_,99999)+".txt";
//  output_powerSpectrum(rho_k_,
//                          filename_end,
//                          binnos,
//                          physicalSize_,
//                          false,
//                          false,
//                          true,
//                          false);
//}

#endif

//}
