//Library, Data.hpp

/*
 * Library containing structure to store data. All data storage 
 * structures must have a recordData() method that takes a state and 
 * a time spent in that state. Some might need to know dimension
 * (can be retrieved by system) while others (ex running avgs)
 * might need to store information in member variables. This offers 
 * a great deal of flexibility in how you want your information 
 * presented. You can also add output methods for sending data to Python
 * for example
 * 
 * Structures must have a
 * -A data storage mechanism and void recordData function which takes the 
 * 		current state and jumptime. Could be time traces, map of points
 * 		and times, array for prob dist, or combo.
 * 
 * If pre-existing structs have the functionality needed, you can
 * merge them by building a new struct that makes the old ones
 * act as members. Then, the new structs recordData() can call the
 * recordData() of its members.
 * 
 */
#ifndef DATA
#define DATA
#include "Gillespie.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <map>
#include <numeric>
#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <deque>

/*
 * Utility Functions and structs for working with the following data types
 */
//Given a vector with the size of an array in each dimension,
//and a tuple-set of indices, return the lienar index of the array
int linear_index(std::vector<int> dimension_sizes,std::vector<int> index){
	int idx = 0;
	for(unsigned int x = 0; x < dimension_sizes.size();x++){
		int temp_idx = index.at(x);
		for(unsigned int y = 0; y < x; y++){
			temp_idx *= dimension_sizes.at(y);
		}
		idx += temp_idx;
	}
	return idx;
}
//Struct for holding an array with its ND sizes. Usually floats, but you never know.
template <typename T>
struct sized_array{
	int totalelements;
	std::vector<int> dim_sizes;
	std::vector<T> array;
	sized_array(std::vector<int>& dimension_sizes){
		totalelements = std::accumulate(dimension_sizes.begin(), dimension_sizes.end(), 1, std::multiplies<int>());
		dim_sizes = dimension_sizes;
		array.resize(totalelements);
	}
	
};

////Data storage structs. Unordered maps, time traces, arrays
//struct pdist_unordered{
	//int dimensions;
	
	//void recordData(){};
//}

//Stores time traces of data. Inappropriate for long runtimes
//as this is quite memory intensive
struct time_traces{
	//Nsteps in the Gillespie Sim*num procs is the best case scenario 
	//size. time traces will be slow, but you could start with this 
	//size anyway. quicker but more memory intensive: you know that you'll only
	//ever have Nmax points in a Gillespie sim, so feed that in as initsize.
	//No reallocation so things are fast. Not terrible assuming you're
	//only ever taking traces of say a few thousand points.
	//Running counter of time. By default start at 0
	precision currenttime = 0.0;
	//vectors with the actual time data
	std::vector<state_t> statehistory;
	std::vector<precision> timehistory;
	time_traces(int const& initsize){
		statehistory.resize(initsize);
		timehistory.resize(initsize);
		
	}
	
	void recordData(state_t const& currentstate,precision time_jump){
			statehistory.push_back(currentstate);
			timehistory.push_back(time_jump);
	}
};

struct online_averages{
	std::vector<precision> averages;
	precision totaltime = 0.0;
    std::vector<precision> avgcompensators;
    std::vector<std::vector<precision>> covariances;
    int dimension;
    
	online_averages(int const& dimensionality){
		dimension = dimensionality;
		averages = std::vector<precision> (dimension,0.0);
		avgcompensators =  std::vector<precision> (dimension,0.0);
		}
	void recordData(state_t const& currentstate,precision time_jump){
			totaltime += time_jump;
			for(unsigned int i = 0;i < currentstate.size();i++){
				KahanSum(averages[i],avgcompensators[i],time_jump*currentstate[i]);
			}
	}
	std::vector<float> normaliseAverages(){
		int dimension = averages.size();
		std::vector<float> true_avgs(dimension);
		for(int i = 0; i < dimension;i++){
			true_avgs[i] = averages[i]/totaltime;
		}
		return true_avgs;
	}
};

//While slowest, a vector_map allows the system to be unbounded
//while collecting the probability distribution information itself
//Because arrays are a bit unqieldly in C++, you can stream this
//to a file and read it into Python for better analysis.
struct vector_map{
	precision totaltime = 0.0;
	std::map<state_t,precision> time_hist;
    
	void recordData(state_t const& currentstate,precision time_jump){
		time_hist[currentstate] += time_jump;
		totaltime += time_jump;
	}
    
	std::map<state_t,precision> normaliseHist(){
		std::map<state_t,precision> pdist;
			for(const auto& x : time_hist){
				pdist[x.first] = x.second/totaltime;
			}
		return pdist;
	}
	
};

//From https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
//and then edited to go into a template specialization for std::hash
namespace std {
    template<> struct hash<state_t> {
        std::size_t operator()(state_t const& vec) const {
          std::size_t seed = vec.size();
          for(auto& i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
          }
          return seed;
        }
    };
}

struct vector_map_unordered{
	precision totaltime = 0.0;
	std::unordered_map<state_t,precision> time_hist;
	
	void recordData(state_t const& currentstate, precision const& time_jump){
		time_hist[currentstate] += time_jump;
		totaltime += time_jump;
	}
	
	std::vector<std::pair<state_t,precision> > sortedVector(std::unordered_map<state_t,precision> hist){
		std::vector<std::pair<state_t,precision> > vec(hist.begin(),hist.end());
		sort(vec.begin(),vec.end());
		return vec;
	}
	
	std::unordered_map<state_t,precision> normaliseHist(){
		std::unordered_map<state_t,precision> pdist;
			for(const auto& x : time_hist){
				pdist[x.first] = x.second/totaltime;
			}
		return pdist;
	}
	
	void stream_to_file(std::ofstream& targetfile){
		std::vector<std::pair<state_t,precision> > outputvector = sortedVector(normaliseHist());
		for(const auto& x : outputvector){
			targetfile << x << std::endl;
		}
	}
	
	sized_array<precision> gen_array(){
		std::unordered_map<state_t,precision> pdist = normaliseHist();
		std::vector<int> dims((pdist.begin()->first).size(),0);
		for(const auto& x : pdist){
			for(unsigned int i = 0; i < dims.size();i++){
				if(x.first.at(i)+1 > dims.at(i)){
					dims.at(i) = x.first.at(i)+1;
				}
			}
		}
		sized_array<precision> data_array (dims);
		for(const auto& x : pdist){
			data_array.array.at(linear_index(dims,x.first)) = x.second;
		}
		return data_array;	
	}
	
};

struct donothing_storage {
    void recordData(state_t ignore, precision alsoignore) {}
};


//See: https://youtu.be/YA-nB2wjVcI?t=88
//I'm assuming that we don't actually need to do the second pass that is 
//normally needed to make sure our heaviest hitters really are the correct
//ones. This is because the simulation is stochastic, and because things
//"settle down" after a while.
//A little note about terminology: a counter is said to be occupied if it
//has any value in it, and it is said to be available (i.e. a candidate for 
//eviction) if its value is <= 0.0.
template <int k, int threshold = 500>
struct k_heaviest {
    //vals will contain the k heaviest hitters
    state_t vals[k];
    //In the THERMALIZE state, counters[i] contains the tally as would
    //normally be done in the k-heavy hitters algorithm. In the RECORD 
    //state, it contains the total time seen for vals[i] and is NOT 
    //decreased when an unknown state is seen.
    precision counters[k];
    //Says how many counters are occupied
    int occupied;
    
    //We check if this exceeds the threshold template parameter to see if we
    //can exit the THERMALIZE state.
    int therm_time;
    
    //Keeps track of total time jumps seen in the RECORD state
    precision total_time;
    
    int state;
    
    static int const THERMALIZE = 0;
    static int const RECORD = 1;
    
    //Initialize state
    k_heaviest() {
        state = THERMALIZE;
        
        for (int i = 0; i < k; i++) {
            counters[i] = 0.0;
        }
        occupied = 0;
        
        therm_time = 0;
        
        total_time = 0.0;
    }
    
    //Little helper function that increments our internal counter, and
    //transitions to the RECORD state if we exceed the threshold. Resets all
    //counters to 0, and for speed, sorts the array of heavy hitters.
    void incr_and_try_transition() {
        therm_time++;
        if (therm_time > threshold) {
            state = RECORD;
            for (int i = 0; i < occupied; i++) {
                counters[i] = 0.0;
            }
            
            std::sort(vals, vals+occupied);
        }
    }
    
    //Main idea: we run the k-heavy hitters algorithm as usual when we are
    //in state THERMALIZE. Then, once the therm_time counter reaches the
    //threshold template parameter, we reset all the times transition to 0.0
    //In the RECORD state, we add the time jump to counters[i] if sim_state
    //is equal to vals[i], otherwise we ignore it.
    void recordData(state_t const& sim_state, precision jump_time) {
        if (state == THERMALIZE) {
            //Check if this state is present in our heavy hitters list. 
            //While we're at it, we can search to see if there's a free slot
            int avail_slot_idx = -1;
            for (int i = 0; i < occupied && i < k; i++) {
                //If this state is present in the values, increase its count
                //and then quit
                if (vals[i] == sim_state) {
                    counters[i] += jump_time;
                    incr_and_try_transition();
                    return;
                }
                //Check if this slot is available. If there is more than one
                //available slot, this idx will find the last one.
                if (counters[i] <= 0.0) avail_slot_idx = i;
            }
            
            if (occupied < k) {
                //If there is an unoccupied slot, occupy it. 
                counters[occupied] = jump_time;
                vals[occupied] = sim_state;
                occupied++;
            } else if (avail_slot_idx != -1) {
                //There was an occupied (but available) slot, evict it
                vals[avail_slot_idx] = sim_state;
                counters[avail_slot_idx] = jump_time;
            } else {
                //No free slot. Decrement all counters
                for (int i = 0; i < occupied; i++) {
                    counters[i] -= jump_time/occupied;
                }
            }
            
            //See if we reached the thermalization time
            incr_and_try_transition();
        } else {
            //RECORD state. 
            total_time += jump_time;
            
            //Search for sim_state in vals. Unintuitive: std::binary_search
            //doesn't actually return an index, so we need to use 
            //std::upper_bound instead
            auto it = std::upper_bound(vals, vals + occupied, sim_state);
            int idx = it - vals - 1; //This -1 is because upper_bound 
            //returns an iterator to the element after the potential match
            
            //Check for equality. This is an ugly workaround because we were
            //forced to use std::upper_bound
            if (idx >= 0 && vals[idx] == sim_state) {
                counters[idx] += jump_time;
            }
        }
    }
  
    void normalize() {
        for (int i = 0; i < occupied; i++) {
            counters[i] /= total_time;
        }
    }
    
    //////////////////////////////////////////////////////////////
    //Type conversions to standard container types//
    //////////////////////////////////////////////////////////////
    
    operator std::unordered_map<state_t,precision>() const {
        std::unordered_map<state_t,precision> ret;
        for (int i = 0; i < occupied; i++) {
            ret[vals[0]] = counters[i]/total_time;
        }
        
        return ret;
    }
    
    operator std::vector<std::pair<state_t, precision> >() const {
        std::vector<std::pair<state_t, precision> > ret;
        ret.reserve(occupied);
        
        for (int i = 0; i < occupied; i++) {
            ret.push_back({vals[i], counters[i]/total_time});
        }
        
        return ret;
    }
};

struct avg_var_cov_2D{
	//Name indicates that with the saved moments you can calculate
	//the average, variance, and covariance of the system steady state
	vector<precision> firstmoments;
    vector<precision> firstcompensators;
    vector<precision> secondmoments;
    vector<precision> secondcompensators;
    precision mixedmoment = 0.0;
    precision mixedcomp = 0.0;
    bool divided = false;
    
    precision totaltime = 0.0;

	avg_var_cov_2D(){
		firstmoments = vector<precision> (2,0.0);
		firstcompensators = vector<precision> (2,0.0);
		secondmoments = vector<precision> (2,0.0);
		secondcompensators = vector<precision> (2,0.0);
		
		}


    //recordData will do the summation, but not the final division by total time to get the 'average'
	void recordData(state_t const& currentstate,precision time_jump){
			if(divided){divided = false;}
			totaltime += time_jump;
			for(unsigned int i = 0;i < 2;i++){
				KahanSum(firstmoments[i],firstcompensators[i],time_jump*currentstate[i]);
				KahanSum(secondmoments[i],secondcompensators[i],time_jump*pow(currentstate[i],2));
			}
			KahanSum(mixedmoment,mixedcomp,time_jump*currentstate[0]*currentstate[1]);
	}
	
	//effects division to create averages
	void generateAverages(){
		if(!divided){
			for(int i = 0; i < 2;i++){
				firstmoments[i] /= totaltime;
				secondmoments[i] /= totaltime;	
			}
			mixedmoment /= totaltime;
			divided = true;
		}
	}
	
	//Averages are the first moments, but to get variances and the cov
	//you need to specifically calculate them
	
	precision calccovariance(){
		generateAverages();
		return mixedmoment - firstmoments[0]*firstmoments[1];
	}
	
	vector<float> calcvariance(){
		generateAverages();
		vector<float> variances(2,0.0);
		for(int i = 0; i <2;i++){
			variances[i] = secondmoments[i] - pow(firstmoments[i],2.0);
		}
		return variances;
	}
	
};

/*Structure to record the averages, variances, and covariances
 * of ND systems. To accompany, linearisation for upper triangular
 * matrices is included.
 */

unsigned int n_from_ij(unsigned int ii,unsigned int jj,unsigned int N){
	//return linear index from i,j pair with N values
	int smaller = std::min(ii,jj);
	int bigger = std::max(ii,jj);
	unsigned int n = N*(N-1)/2 - (N-smaller)*(N-smaller-1)/2 + bigger;
	return n;
} 


struct ND_avg_var_cov{
    std::vector<double> means,vars,covs,deltas;
    double tsum,tsum2;
    bool calculated = false;
    int dim;
    ND_avg_var_cov(int const & dimension){
		//there will be N choose 2 covariances = (N-1)(N)/2 + N variances
		int numcovs = dimension + (dimension-1)*dimension/2;
		covs = std::vector<double> (numcovs,0.0);
		//There will be ND vectors for averages, deltas, and variances
		means = std::vector<double> (dimension,0.0);
		//vars = std::vector<double> (dimension,0.0);
		deltas = std::vector<double> (dimension,0.0);
		//total weighting, and sum of squared weights
		tsum = 0;
		tsum2 = 0;
		//Dimension
		dim = dimension;
	}
	ND_avg_var_cov(){}
    //one loop of stable calculation.
    //Here, we have to go through whole state for all the deltas. Once
    //we have all the deltas, we can do covariances
    void recordData(state_t const& currentstate,precision time_jump){
        if(calculated){calculated = false;}

        tsum += time_jump;
        tsum2 += pow(time_jump,2);
        for(int i = 0; i < dim;i++){
			double d = currentstate.at(i) - means.at(i);
			deltas.at(i) = d;
			means.at(i) += (time_jump/tsum) * d;
			//vars.at(i) += time_jump * pow(d,2);
		}
		        
		//Calculate covariances
		for(int i = 0;i < dim;i++){
			for(int j=i;j<dim;j++){
				covs.at(n_from_ij(i,j,dim)) += time_jump*deltas.at(i)*(currentstate.at(j)-means.at(j));
			}
		}
	}
	//Function to effect divisions on variances and covariances
	void normCovs(){
		if(!calculated){
			for(unsigned int i=0;i<covs.size();i++){
				covs.at(i) /= tsum;
			}
			//for(unsigned int i=0;i<vars.size();i++){
				//vars.at(i) /= tsum;
			//}
		calculated = true;
		}
	}
	
	double covariance(int i,int j){
		normCovs();
		return covs.at(n_from_ij(i,j,dim));
	}
	//Get the correlation
	double correlation(int i,int j){
		double num = covariance(i,j);
		double dem = sqrt(covariance(i,i)*covariance(j,j));
		return num/dem;
	}
	
};

struct store_N_ordered{
	//Stores N last states visited in the simulation.
		std::deque<state_t> states;
		std::deque<precision> times;
	unsigned int N_states;
	int oldest_loc = 0;
	store_N_ordered(int const& N){
		states = std::deque<state_t>(N);
		times = std::deque<precision>(N);
		N_states = N;
	}
	
	void recordData(state_t const& currentstate,precision time_jump){
			if(states.size()>=N_states){states.pop_front();}
			if(times.size()>=N_states){times.pop_front();}
			states.push_back(currentstate);
			times.push_back(time_jump);
	}
	
	
	void stream_to_file(std::ofstream& targetfile){
		targetfile << "[" ;
		for (unsigned int i = 0; i < N_states; i++){
			if(i == 0){targetfile << times[i];}
			else{targetfile << "," << times[i];}
		}
		targetfile << "]" << std::endl;
		targetfile << "[" ;
		for (unsigned int i = 0; i < N_states; i++){
			if(i == 0){targetfile << states[i];}
			else{targetfile << "," << states[i];}
		}
		targetfile << "]";
	}
	
	void stream_to_np_bin(std::fstream& fstates,std::fstream& ftimes){

		for (unsigned int i = 0; i < N_states; i++){
			fstates.write((char*) states.at(i).data(),
                     states.at(i).size() * sizeof(int));
			ftimes.write(reinterpret_cast<const char*>( &times.at(i) ), sizeof(precision));
		}
	}
	
};


struct N_ordered_pdist{
	//Stores N last states visited in the simulation and probability distribution of whole simulation
		vector_map_unordered pdist;
		store_N_ordered trajectory;
	N_ordered_pdist(int const& N):trajectory(N){}
	
	void recordData(state_t const& currentstate,precision time_jump){
			pdist.recordData(currentstate,time_jump);
			trajectory.recordData(currentstate,time_jump);
	}
	
	
	void stream_to_file(std::ofstream& probabilityfile,std::ofstream& trajectoryfile){
		pdist.stream_to_file(probabilityfile);
		trajectory.stream_to_file(trajectoryfile);
	}
	
};


//Nice pretty-printing operator
template <int a, int b>
std::ostream& operator<< (std::ostream& o, k_heaviest<a,b> const& kh) {
    return o << std::vector<std::pair<state_t, precision> >(kh);
}


//struct CrossCorrelation{
	////Calculate the cross correlation and means of systems.
	////Cross corr starts at t and goes to times t + tau_i where tau_i 
	////is given in a vector
    //std::vector<double> means,vars,covs,deltas;
    //double tsum,tsum2;
    //bool calculated = false;
    //int dim;
    //bool divided = false;
    
    //precision totaltime = 0.0;

	//YanHsu_FDT_Calc(){
		//firstmoments = vector<precision> (2,0.0);
		//firstcompensators = vector<precision> (2,0.0);
		//secondmoments = vector<precision> (2,0.0);
		//secondcompensators = vector<precision> (2,0.0);
		
		//}


    ////recordData will do the summation, but not the final division by total time to get the 'average'
	//void recordData(state_t const& currentstate,precision time_jump){
			//if(divided){divided = false;}
			//totaltime += time_jump;
			//for(unsigned int i = 0;i < 2;i++){
				//KahanSum(firstmoments[i],firstcompensators[i],time_jump*currentstate[i]);
				//KahanSum(secondmoments[i],secondcompensators[i],time_jump*pow(currentstate[i],2));
			//}
			//KahanSum(mixedmoment,mixedcomp,time_jump*currentstate[0]*currentstate[1]);
	//}
	
	////effects division to create averages
	//void generateAverages(){
		//if(!divided){
			//for(int i = 0; i < 2;i++){
				//firstmoments[i] /= totaltime;
				//secondmoments[i] /= totaltime;	
			//}
			//mixedmoment /= totaltime;
			//divided = true;
		//}
	//}
	
	////Averages are the first moments, but to get variances and the cov
	////you need to specifically calculate them
	
	//precision calccovariance(){
		//generateAverages();
		//return mixedmoment - firstmoments[0]*firstmoments[1];
	//}
	
	//vector<float> calcvariance(){
		//generateAverages();
		//vector<float> variances(2,0.0);
		//for(int i = 0; i <2;i++){
			//variances[i] = secondmoments[i] - pow(firstmoments[i],2.0);
		//}
		//return variances;
	//}
	
//};

#endif
