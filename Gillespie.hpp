//LIBRARY, gsim.hpp

#ifndef GILLESPIE
#define GILLESPIE

#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>
#include <chrono>
#include <tuple>

using precision = double; //use to control precision of time values
using rand_eng = std::mt19937; //use to control random engine
using state_t = std::vector<int>;
using rate_list = std::vector<precision>;

/*
 * The first set of templates and functions are useful utilities
 * used throughout my own research, which I've enclosed here for
 * anyone who might also find them useful
 */


//Following template taken from Marco Merlini for vector streaming
template <typename T>
std::ostream& operator<< (std::ostream &o, std::vector<T> const& v) {  
	o << "[";  
	auto delim = "";  
	for (auto const& i : v) {    
		o << delim << i;    
		delim = ", ";  
		}  
	return o << "]";
	}

template <typename T,typename T2>
std::ostream& operator<< (std::ostream &o, std::pair<T,T2> const& v) {  
	o << "(";  
	o << v.first << "," << v.second;
	return o << ")";
	}

using namespace std;
using namespace std::chrono;

class timer
{
	clock_t start;
public:
	timer();
	~timer();
};

timer::timer(){start=clock();}
timer::~timer()
{
	int hour,min;
	float sec;
	clock_t end;

	end=clock();

	sec=(float)(end-start)/CLOCKS_PER_SEC;
	min=(int)sec/60;
        sec=sec-60*min;
	hour=(int)min/60;
        min=min-60*hour;

	cout<<"Time Elapsed : "<<hour<<" hour: "<<min<<" min: "<<sec<<" sec. "<<endl;
}


int gen_seed(){
	//return 10; //this lets you fix a seed for the whole program (deterministic)
	//Clock based seeding. This is okay.
	return static_cast<int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	//To do: add a check on std::random_device().entropy. If its 0, use clock based. Otherwise use random device.
}

/*Relative Error Calculation
 */

double rel_error(double RHS, double LHS){
	double val;
	if(RHS == 0){
		if(LHS == 0){
			val = 0;
		}
		else{
			val = abs(LHS);
		}
	}
	else if(LHS==0){val = abs(RHS);}
	else{val = abs(RHS - LHS)/abs((RHS+LHS)/2);}
	return val;
}



/*
 * Kahan Sum for stable average calculations
 */
 
void KahanSum(precision& current_total,precision& compensator,precision added_val){
	precision y = added_val - compensator;
	precision t = current_total + y;
	compensator = (t-current_total) - y;
	current_total = t;
}

/*
 * Struct for stable online weighted covariance/variance calculations
 */

struct online_weighted_avg_cov_var{
    double meanx,meany,wsum,wsum2,C,Vx,Vy;
  online_weighted_avg_cov_var(){
    meanx = 0;
    meany = 0;
    wsum = 0; //Total weights
    wsum2 = 0; //Total squared weights
    C = 0; //The covariance
    Vx = 0; //The variance in x
	Vy = 0; //The variance in y
	}
    //one loop of stable calculation. X is first data, Y is second data, w is weight
    void append(double const& x,double const& y,double const& w){
        wsum += w;
        wsum2 += w * w;
        double dx = x - meanx;
        double dy = y - meany;
        meanx += (w / wsum) * dx;
        meany += (w / wsum) * dy;
        C += w * dx * dy;
        Vx += w * pow(dx,2);
        Vy += w * pow(dy,2);
	}
	//Function to return meanx,meany,varx,vary,covxy
	std::vector<double> genMoments(){
		double covar = C / wsum;
		double varx = Vx / wsum;
		double vary = Vy / wsum;
		//# Bessel's correction for sample variance
		//# Frequency weights
		//sample_frequency_covar = C / (wsum - 1)
		//# Reliability weights
		//sample_reliability_covar = C / (wsum - wsum2 / wsum)
		std::vector<double> moms = {meanx,meany,varx,vary,covar};
		return moms;
	}

};

/*
 */
template<typename system>
struct balance_eqs{
	//Calculated for an ND system:
	//An ND vector of system averages (+ deltas for gamma)
	//An ND vector of birth rate averages
	//An ND vector of death rate averages
	//An ND vector of death-birth rate averages (+deltas for gamma)
	//An NxN matrix with element i,j = Cov(R(-)[i] - R(+)[i],x[j])
	// [I'll call gamma
	//
	//Provided: 
	//the s(i,j) matrix [stored in system object as member function]
	//Assumes that s.getRates returns {birth0,birth1...birthN,death0,death1...deathN}
    std::vector<precision> means,mean_deltas;
    std::vector<precision> birthrates,deathrates,birthcomps,deathcomps;
    std::vector<precision> netrates,net_deltas;
    std::vector<std::vector<precision>> gamma;
    precision tsum,tsum2;
    bool calculated = false;
    int N;
    system s;
    
    balance_eqs(system s_obj){
		s = s_obj;
		N = s.dimensions;
		//Initialise ND vectors
		means = std::vector<precision> (N,0.0);
		mean_deltas = std::vector<precision> (N,0.0);
		birthrates = std::vector<precision> (N,0.0);
		deathrates = std::vector<precision> (N,0.0);
		birthcomps = std::vector<precision> (N,0);
		deathcomps = std::vector<precision>(N,0);
		
		netrates = std::vector<precision> (N,0.0);
		net_deltas = std::vector<precision> (N,0.0);
		//Initialise gamma
		gamma = std::vector<std::vector<precision>> (N);
		for(int i = 0;i < N;i++){
			gamma.at(i).resize(N);
		}
		//total weighting, and sum of squared weights
		tsum = 0;
		tsum2 = 0;
	}
	balance_eqs(){}
	void recordBalance(state_t currentstate,rate_list rates,precision time_jump){
		tsum += time_jump;
		tsum2 += time_jump*time_jump;
		//update system means and get deltas
        for(int i = 0; i < N;i++){
			double d = currentstate.at(i) - means.at(i);
			mean_deltas.at(i) = d;
			means.at(i) += (time_jump/tsum)*d;
			
			//For birth and death rate averages, I don't
			//need to save the deltas for any covariances

			//double b_d = rates.at(i) - birthrates.at(i);
			//birthrates.at(i) += (time_jump/tsum)*b_d;
			//double d_d = rates.at(i+N) - deathrates.at(i);
			//deathrates.at(i) += (time_jump/tsum)*d_d;

			KahanSum(birthrates.at(i),birthcomps.at(i),(time_jump)*rates.at(i));
			KahanSum(deathrates.at(i),deathcomps.at(i),(time_jump)*rates.at(i+N));

			//For net rate, need average for delta calculation
			//i+N is where the death rate is, i is where the birth rate is, for system var i
			double n_d = (rates.at(i+N)-rates.at(i)) - (netrates.at(i));
			net_deltas.at(i) = n_d;
			netrates.at(i) += (time_jump/tsum)*n_d;
		}
		//Now calculate the gammas for covariance
		//Note gammas are NOT symmetric: i is the rate variable, and j is the system var
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				gamma.at(i).at(j) += time_jump*net_deltas.at(i)*(currentstate.at(j)-means.at(j));
			}
		}
	}
	
	void normalise(){
		if(!calculated){
			for(int i = 0; i < N; i++){
				birthrates.at(i)  /= tsum;
				deathrates.at(i) /= tsum;
				for(int j = 0;j < N;j++){
					gamma.at(i).at(j) /= tsum;
				}
				calculated = true;
			}
		}
	}
	
	void unnormalise(){
		if(calculated){
			for(int i = 0; i < N; i++){
				birthrates.at(i)  *= tsum;
				deathrates.at(i) *= tsum;
				for(int j = 0;j < N;j++){
					gamma.at(i).at(j) *= tsum;
				}
				calculated = false;
			}
		}
	}
	
	//To calculate covariance balance eqs, the covariances need to be normalised
	//First, Flux balance check (did avg birth rate in i = average deathrate in i
	double checkFluxBalance(int i){
		normalise();
		return rel_error(birthrates.at(i),deathrates.at(i));
	}
	double checkCovBalance(int i,int j){
		//The i,j covariance balance equation is:
		/*
		 * Cov(netrate_j,x_i) + Cov(netrate_i,x_j) = (birth/death)rate_i*s_ij + (birth/death)rate_j*s_ji
		 * You can divide both sides by <x_i><x_j> to reduce numerical errors: makes values relative to system size
		 */
		 
		 //If RHS and LHS were not already normalised (like above) gammas also do not need to be.
		 //Rearrange the equation to minimize the chance either side is precisely 0
		 normalise();
		 double LHS = gamma.at(i).at(j) + gamma.at(j).at(i) ;
		 double RHS = birthrates.at(i)*s.stepmatrix(i,j)  + deathrates.at(j)*s.stepmatrix(j,i);

	return rel_error(RHS,LHS);
	}
	
	std::tuple<double,double> resCovBalance(int i,int j){
		//The i,j covariance balance equation is:
		/*
		 * Cov(netrate_j,x_i) + Cov(netrate_i,x_j) = (birth/death)rate_i*s_ij + (birth/death)rate_j*s_ji
		 * You can divide both sides by <x_i><x_j> to reduce numerical errors: makes values relative to system size
		 */
		 
		 //If RHS and LHS were not already normalised (like above) gammas also do not need to be.
		 //Rearrange the equation to minimize the chance either side is precisely 0
		 normalise();
		 double RHS = gamma.at(i).at(j) + gamma.at(j).at(i) ;
		 double LHS = birthrates.at(i)*s.stepmatrix(i,j)  + deathrates.at(j)*s.stepmatrix(j,i);
	return {RHS,LHS};
	}
	
};



/*
 * The following functions and templates are the actual Gillespie 
 * simulation utilities.
 */

//given list of rates find jumptime
precision findJumptime(precision const& r_tot,rand_eng& engine){
	std::exponential_distribution<precision> distribution (r_tot);
	precision t_jump = distribution(engine);
	return t_jump;
}

//Jumptime function that doesn't need exponential dist
precision findJumptime_uniform(precision const& r_tot, rand_eng& engine,std::uniform_real_distribution<precision> distribution){
	precision rand_uni = distribution(engine);
	precision t_jump = -1.0*std::log(1 - rand_uni) / r_tot;
	return t_jump;
}

unsigned int chooseProcess(rate_list const& rates, precision const& r_tot, rand_eng& engine,std::uniform_real_distribution<precision> distribution){
	precision p = distribution(engine);
	precision pointed_val = p*r_tot;
	precision ratesum = 0.0;
	for(unsigned int proc = 0; proc < rates.size()-1;proc++){
		ratesum += rates[proc]; //Keeps running sum of rates up until now
		if(pointed_val <= ratesum){
			return proc;
		}
	}
	//Triggers if you didn't find the value in the partial summation. 
	//Needs error check in case you get p*r_tot > r_tot
	return rates.size()-1; 
	
}
 /*
 * Gillespie Simulation Driver function. 
 * type system needs to have member functions getRates and getStates,
 * appropriately formatted. Examples can be found in Systems.hpp.
 * 
 */
template<typename system,typename datastore>
void gsim(system& s,datastore& d, state_t& initialstate,long int Nsteps,long int Nmax,rand_eng& engine){

	state_t currentstate = initialstate;
	std::vector<int> counts(s.Numprocs,Nsteps);
	long int steps = 0; //Check for maximum number of steps
	long int globalcount = Nsteps*s.Numprocs;
	
	std::uniform_real_distribution<precision> uni_distribution(0.0,1.0);

	while( globalcount > 0 && steps < Nmax){
		//loop: while each process has run less than Nsteps times, or system hasnt reached Nmx steps
		rate_list rates = s.getRates(currentstate);

		//Get total rate
		precision r_tot = std::accumulate(rates.begin(),rates.end(),0.0);
		if(r_tot <= 0.0){ //change to within an epsilon
			std::cout << "No processes are going to occur. \n";
			std::cout << "The absorbing state is " << currentstate <<"\n";
			std::cout << "The remaining rates are " <<rates << "\n";
			std::cout << "This had a total rate of " <<r_tot <<"\n";
			break;
		}
        
		//Get time in state before jump
		precision t_jump = findJumptime_uniform(r_tot,engine,uni_distribution);
		//Record time spent in current state
		d.recordData(currentstate,t_jump);
		//Choose the process that occured and update state
		int proc = chooseProcess(rates,r_tot,engine,uni_distribution);
		currentstate = s.getState(proc,currentstate);
		//Add step, and remove count as appropriate
		steps++;
		if(counts.at(proc) > 0){
			globalcount -= 1;
			counts.at(proc) -= 1;
		}
	}
	std::cout << "Steps ran: " << steps << "\n";
	
}

/*
 * The following keeps track of the covariance balance equations
 * automatically. Requires system has a step matrix
 */

template<typename system,typename datastore>
balance_eqs<system> gsim_cb(system& s,datastore& d, state_t& initialstate,long int Nsteps,long int Nmax,rand_eng& engine){

	state_t currentstate = initialstate;
	std::vector<int> counts(s.Numprocs,Nsteps);
	long int steps = 0; //Check for maximum number of steps
	long int globalcount = Nsteps*s.Numprocs;
	balance_eqs<system> balances(s);
	std::uniform_real_distribution<precision> uni_distribution(0.0,1.0);

	while( globalcount > 0 && steps < Nmax){
		//loop: while each process has run less than Nsteps times, or system hasnt reached Nmx steps
		rate_list rates = s.getRates(currentstate);

		//Get total rate
		precision r_tot = std::accumulate(rates.begin(),rates.end(),0.0);
		if(r_tot <= 0.0){ //change to within an epsilon
			std::cout << "No processes are going to occur. \n";
			std::cout << "The absorbing state is " << currentstate <<"\n";
			break;
		}
        
		//Get time in state before jump
		precision t_jump = findJumptime_uniform(r_tot,engine,uni_distribution);
		//Record time spent in current state
		d.recordData(currentstate,t_jump);
		//Record the RHS and LHS of each of the covariance balance equations
		//To use this version of gsim, systems must have a birthrates and deathrates
		//member which should be a vector of which processes are the birth rates and
		//which processes are the death rates
		balances.recordBalance(currentstate,rates,t_jump);
		//Choose the process that occured and update state
		int proc = chooseProcess(rates,r_tot,engine,uni_distribution);
		currentstate = s.getState(proc,currentstate);
		//Add step, and remove count as appropriate
		steps++;
		if(counts.at(proc) > 0){
			globalcount -= 1;
			counts.at(proc) -= 1;
		}
	}
	std::cout << "Steps ran: " << steps << "\n";
	return balances;
}

/*
 * Version which runs until CB is met under threshold
 */
 

template<typename system,typename datastore>
balance_eqs<system>  gsim_cb_thres(system& s,datastore& d, state_t& initialstate,double threshold,rand_eng& engine){

	state_t currentstate = initialstate;
	balance_eqs<system> balances(s);
	std::uniform_real_distribution<precision> uni_distribution(0.0,1.0);
	bool run = true;
	long int steps = 0;
	
	while(run){
		//loop: while each process has run less than Nsteps times, or system hasnt reached Nmx steps
		rate_list rates = s.getRates(currentstate);

		//Get total rate
		precision r_tot = std::accumulate(rates.begin(),rates.end(),0);
		if(r_tot <= 0.0){ //change to within an epsilon
			std::cout << "No processes are going to occur. \n";
			std::cout << "The absorbing state is " << currentstate <<"\n";
			break;
		}
        
		//Get time in state before jump
		precision t_jump = findJumptime_uniform(r_tot,engine,uni_distribution);
		//Record time spent in current state
		d.recordData(currentstate,t_jump);
		//Record the RHS and LHS of each of the covariance balance equations
		//To use this version of gsim, systems must have a birthrates and deathrates
		//member which should be a vector of which processes are the birth rates and
		//which processes are the death rates
		balances.recordBalance(currentstate,rates,t_jump);
		//Choose the process that occured and update state
		int proc = chooseProcess(rates,r_tot,engine,uni_distribution);
		currentstate = s.getState(proc,currentstate);
		
		run = false;
		for(int i = 0;i < s.dimensions;i++){
			if(balances.checkFluxBalance(i) >= threshold){
				run = true;
				break;
			}
			for(int j = 0; j < i+1;j++){
				if(balances.checkCovBalance(i,j) >= threshold){
					run = true;
					break;
				}
			}
		}
		if(run){
			steps++;
			balances.unnormalise();
		}
	}
	std::cout << "Steps ran: " << steps << "\n";
	return balances;
}

#endif
