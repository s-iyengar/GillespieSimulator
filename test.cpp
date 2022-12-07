#include "Systems.hpp"
#include "Gillespie.hpp"
#include "Drugs.hpp"
#include <random>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <ctime>

using namespace std;
namespace fs = std::filesystem;


std::string return_current_time_and_date()
{
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%m-%d %H%M");
    return ss.str();
}


int main(){
	string resultsfolder = "results/mRNAProteintrace" + return_current_time_and_date();
	long int Nmax = 1000000000;
	long int Nsteps = 1000000000;
	timer ob;

	state_t initialstate = {10,30};
	
	mRNA_Protein_Basic system{20,15.0/20,1,1};

	rand_eng engine(gen_seed());
	
	store_N_ordered data(1000);
	balance_eqs<mRNA_Protein_Basic> bals = gsim_cb_thres(system,data,initialstate,0.01,engine);
		
	
	fs::create_directories(resultsfolder);
	
	ofstream pdist;
	pdist.open(resultsfolder+"/trace.dat");
	data.stream_to_file(pdist);
	pdist.close();


	
	pdist.close();
	return 0;
}
