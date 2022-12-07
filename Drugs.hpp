#ifndef DRUGS
#define DRUGS

#include "Gillespie.hpp"
#include "Data.hpp"
#include <random>
#include <vector>
#include <iostream>
#include <cmath>
#include <tuple>

//Data storage type for drug perturbation. Stores last Nstore states.
//And has method for averaging them.

struct store_N{
	//Stores N last states visited in the simulation.
	//Contains member function to average them later
	std::vector<state_t> states;
	std::vector<precision> times;
	int N_states;
	int oldest_loc = 0;
	store_N(int const& N){
		states = std::vector<state_t> (N);
		times = std::vector<precision> (N);
		N_states = N;
	}
	
	void recordData(state_t const& currentstate,precision time_jump){
			states.at(oldest_loc) = currentstate;
			times.at(oldest_loc) = time_jump;
			oldest_loc = (oldest_loc+1)%N_states;
	}
	
	std::vector<precision> averageState(){
		unsigned int dim = states.at(0).size();
		std::vector<precision> float_avg(dim,0.0);
		std::vector<precision> compensators(dim,0.0);
		precision totaltime = 0.0;
		precision timecomp = 0.0;
		std::vector<precision> avgstate(dim,0);
		for(int i = 0; i < N_states; i++){
			if(states.at(i).size() == 0){
				std::cout << "This run only averaged " << i+1 << " states. \n";
				states.resize(i);
				times.resize(i);
				break;
			};
			for(unsigned int j = 0; j < dim; j++){
				KahanSum(float_avg.at(j),compensators.at(j),times.at(i)*states.at(i).at(j));
			}
			KahanSum(totaltime,timecomp,times.at(i));
		}

		for(unsigned int i = 0; i < dim; i++){
			avgstate.at(i) = float_avg.at(i)/totaltime;
		}
		return avgstate;
	}
	
	void stream_to_file(std::ofstream& targetfile){
		targetfile << times << std::endl;
		targetfile << states;
	}
	
};

struct store_corrs_Ns{
	store_N avgstore;
	ND_avg_var_cov covstore;
	store_corrs_Ns(int Navg,int dimension):avgstore(Navg),covstore(dimension){}
	
	void recordData(state_t const& currentstate,precision time_jump){
		avgstore.recordData(currentstate,time_jump);
		covstore.recordData(currentstate,time_jump);
	}
};

struct averagestate{
		std::vector<precision> states;
		std::vector<precision> compensators;
		double totaltime = 0;
		averagestate(){}
		averagestate(int dimension){
			states = std::vector<precision> (dimension);
			compensators = std::vector<precision> (dimension);
		}
		
		void recordData(state_t const& currentstate,precision time_jump){
			for(unsigned int i = 0; i< currentstate.size();i++){
				KahanSum(states.at(i),compensators.at(i),time_jump*currentstate.at(i));
			}
			totaltime += time_jump;
		}
	std::vector<precision> averageState(){
		std::vector<precision> avgstate (states.size());
		for(unsigned int i = 0;i<states.size();i++){
			double s = states.at(i) / totaltime;
			avgstate.at(i) = s;
		}
		return avgstate;
	}
};

struct store_corrs_avgs{
	averagestate avgstore;
	ND_avg_var_cov covstore;
	store_corrs_avgs(int dimension):avgstore(dimension),covstore(dimension){}
	
	void recordData(state_t const& currentstate,precision time_jump){
		avgstore.recordData(currentstate,time_jump);
		covstore.recordData(currentstate,time_jump);
	}
};

template<typename system>
std::tuple<std::vector<precision>,ND_avg_var_cov> average_and_corrs_from_system(system s, state_t i_state, long int Nsteps,long int Nmax,int Navg,rand_eng eng){
	store_corrs_avgs res(i_state.size());
	gsim(s,res,i_state,Nsteps,Nmax,eng);
	std::vector<precision> uavg = res.avgstore.averageState();
	return {uavg,res.covstore};
}


template<typename system>
std::vector<precision> average_from_system(system s, state_t u_state, long int Nsteps,long int Nmax,int Navg,rand_eng eng){
	averagestate pert_res(u_state.size());
	gsim(s,pert_res,u_state,Nsteps,Nmax,eng);
	std::vector<precision> newavg = pert_res.averageState();
	return newavg;
}

//Template to run gillespie algorithm and retrieve the delta m_i values
//Uses above defined data storgae type. Deltas compared to an input
//unperturbed/initial state
template<typename system>
std::vector<precision> deltas_from_drug(system s, std::vector<precision> u_state_p, long int Nsteps,long int Nmax,int Navg,rand_eng engine){
	state_t u_init_state (u_state_p.size());
	for(unsigned int i = 0; i< u_state_p.size();i++){
		u_init_state.at(i) = round(u_state_p.at(i));
	}
	std::vector<precision> deltas = average_from_system(s,u_init_state,Nsteps,Nmax,Navg,engine);//get new averages, then subtract
	for(unsigned int i = 0; i < deltas.size(); i++){
		deltas.at(i) -= u_state_p.at(i);
	}
	return deltas;
}

//Given an unperturbed system, and a vector of other systems,collect 
//deltas for each system in a vector
template<typename system>
std::tuple<std::vector<std::vector<precision>>,ND_avg_var_cov> runExperiment(system unpert_system,std::vector<system> pert_systems,state_t initstate,long int Nsteps,long int Nmax,int Navg){
	unsigned int numperts = pert_systems.size();
	std::vector<std::vector<precision>> deltas(numperts);
	rand_eng eng(gen_seed());
	std::vector<precision> unpert_state;
	ND_avg_var_cov unpert_cov;
	std::tie(unpert_state,unpert_cov) = average_and_corrs_from_system(unpert_system, initstate,Nsteps,Nmax,Navg,eng);
	for(unsigned int n = 0;n < numperts;n++){
		rand_eng neweng(gen_seed());
		deltas.at(n) = deltas_from_drug(pert_systems.at(n),unpert_state,Nsteps,Nmax,Navg,neweng);
		std::cout << "Perturbations remaining: " << numperts-(n+1) << "\n";
	}
	return {deltas,unpert_cov};
}



/*
 * Version to run until CB meets threshold
 */
template<typename system>
std::tuple<state_t,ND_avg_var_cov> average_and_corrs_from_system_threshold(system s, state_t i_state, double threshold,int Navg,rand_eng eng){
	store_corrs_avgs res(i_state.size());
	gsim_cb_thres(s,res,i_state,threshold,eng);
	std::vector<precision> uavg_p = res.avgstore.averageState();
	state_t uavg(uavg_p.begin(), uavg_p.end());
	return {uavg,res.covstore};
}


template<typename system>
state_t average_from_system_threshold(system s, state_t u_state,double threshold,int Navg,rand_eng eng){
	averagestate pert_res(u_state.size());
	gsim_cb_thres(s,pert_res,u_state,threshold,eng);
	std::vector<precision> newavg_p = pert_res.averageState();
	state_t newavg(newavg_p.begin(), newavg_p.end());
	return newavg;
}

//Template to run gillespie algorithm and retrieve the delta m_i values
//Uses above defined data storgae type. Deltas compared to an input
//unperturbed/initial state
template<typename system>
state_t deltas_from_drug_threshold(system s, state_t u_state, double threshold,int Navg,rand_eng engine){
	state_t deltas = average_from_system_threshold(s,u_state,threshold,Navg,engine);//get new averages, then subtract
	for(unsigned int i = 0; i < deltas.size(); i++){
		deltas.at(i) -= u_state.at(i);
	}
	return deltas;
}

//Given an unperturbed system, and a vector of other systems,collect 
//deltas for each system in a vector
template<typename system>
std::tuple<std::vector<state_t>,ND_avg_var_cov> runExperiment_threshold(system unpert_system,std::vector<system> pert_systems,state_t initstate,double threshold,int Navg){
	unsigned int numperts = pert_systems.size();
	std::vector<state_t> deltas(numperts);
	rand_eng eng(gen_seed());
	state_t unpert_state;
	ND_avg_var_cov unpert_cov;
	std::tie(unpert_state,unpert_cov) = average_and_corrs_from_system_threshold(unpert_system, initstate,threshold,Navg,eng);
	for(unsigned int n = 0;n < numperts;n++){
		rand_eng neweng(gen_seed());
		deltas.at(n) = (deltas_from_drug_threshold(pert_systems.at(n),unpert_state,threshold,Navg,neweng));
		std::cout << "Perturbations remaining: " << numperts-(n+1) << "\n";
	}
	return {deltas,unpert_cov};
}

#endif
