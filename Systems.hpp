/*
 * Define Structures that are systems. Systems are constituted
 * of a set of possible processes. Each process should indicate how to 
 * calculate its rate given a state, and what the new state after it 
 * will be. For example, we can define an 'mRNA-Protein model' as:
 * 
 * [old state] --(rate)--> [new state]
 * 
 * m,p --(lambda) --> m+1,p  (1)
 * m,p --(gamma*m) --> m,p+1 (2)
 * m,p --(beta_m*m) --> m-1,p (3)
 * m,p --(beta_p*m) --> m,p-1 (4)
 * 
 * Where [m,p] is your state and lambda, gamma, and the betas are 
 * constants. If all constants were set to 1 and the system was in
 * state [10,5], the rates of each process and the new states they lead
 * to are as follows:
 * (1): 1 [11,5]
 * (2): 10 [10,6]
 * (3): 10 [9,5]
 * (4): 5 [10,4]
 * 
 * More complicated processes are possible: for example there may be a 
 * feedback loop where the rate of (1) is dependent on variable p. Or
 * there may be an additional process where two mRNA and a protein, if 
 * they interact spontaneously degrade.
 * 
 * (ex 1) m,p --(lambda**n/K**n + p**n) --> m+1,p 
 * (ex 2) m,p -- (delta*p*m**2) --> m-2,p-1 
 * 
 * The template for a system needs:
 *  -All of the relevant constants (almost always floats)
 *  - The dimensionality and number of processes in a system
 * 	-A getRates() function which given a state outputs the rates of
 * 		N procedures as {rate_1...rate_N}
 *  -A getState() function which given an ID number (index 
 * 		correspending to getRates output) tells you what state to move
 * 		to.
 * 
 * New to July 29th: (not implemented in all systems)
 * 	-A stepmatrix() function which given two integers less than the dimension
 * 		s_ij = sum_rates[fractionofxiinrate*|changeinjduringrate|*sign(changeiniduringrate*changeinjduringrate)]
 *  _Automatic covariance balance checking assumes rates are ordered {births...deaths...}
 *    This might be changed in the future but no guarantees it is.
 */
 #ifndef SYSTEMS
 #define SYSTEMS
 #include "Gillespie.hpp"
 #include <cmath>

//many systems have the identity matrix as their step matrix. Defined here for simplicity
double identity_mat(int i,int j){
	double val = 0;
	if(i==j){
		val = 1;
	}
	return val;
}

 //Example system from top comment
struct mRNA_Protein_Basic{
	 //Members: four constants chosen at initialisation, two
	 //constants of the simulation
	 float lambda;
	 float gamma;
	 float beta_m;
	 float beta_p;
	 static int const dimensions = 2;
	 //Numprocs must match the length of rates output
	 static int const Numprocs = 4;
	 
	 //Given a state, calculate the four rates.
	 rate_list getRates(state_t const& state){
		 //retrieve state vars
		 int m = state.at(0);
		 int p = state.at(1);
		 //calculate rates
		 float mrnabirth = lambda;
		 float mrnadeath = beta_m*m;
		 float proteinbirth = gamma*m;
		 float proteindeath = beta_p*p;
		 //join in list
		 rate_list rates = {mrnabirth,proteinbirth,mrnadeath,proteindeath};
		 return rates;
	 }
	 
	 //Given an integer correspending to an index of the outputs of
	 //getRates() and the current state, yield the new state.
	 state_t getState(int const& proc,state_t const& currentstate){
		 state_t newstate = currentstate;
		 switch(proc){
			case 0: //mRNA birth
				newstate.at(0) += 1;
				break;
			case 2: //mRNA death
				newstate.at(0) -= 1;
				break;
			case 1: //protein birth
				newstate.at(1) += 1;
				break;
			case 3: //protein death
				newstate.at(1) -= 1;
				break;
		}
		return newstate;
	 } 
	 double stepmatrix(int i,int j){
		return identity_mat(i,j);
	 }
};

struct mRNA_Protein_HillFeedback{
	 //Members: four constants chosen at initialisation, two
	 //constants of the simulation
	 float lambda;
	 int n;
	 float K;
	 float gamma;
	 float beta_m;
	 float beta_p;
	 static int const dimensions = 2;
	 //Numprocs must match the length of rates output
	 static int const Numprocs = 4;
	 
	 //Given a state, calculate the four rates.
	 rate_list getRates(state_t const& state){
		 //retrieve state vars
		 int m = state.at(0);
		 int p = state.at(1);
		 //calculate rates
		 float mrnabirth = lambda/(std::pow(p/K,n) + 1);
		 float mrnadeath = beta_m*m;
		 float proteinbirth = gamma*m;
		 float proteindeath = beta_p*p;
		 //join in list
		 rate_list rates = {mrnabirth,mrnadeath,proteinbirth,proteindeath};
		 return rates;
	 }
	 
	 //Given an integer correspending to an index of the outputs of
	 //getRates() and the current state, yield the new state.
	 state_t getState(int& proc,state_t const& currentstate){
		 state_t newstate = currentstate;
		 switch(proc){
			case 0: //mRNA birth
				newstate.at(0) += 1;
				break;
			case 1: //mRNA death
				newstate.at(0) -= 1;
				break;
			case 2: //protein birth
				newstate.at(1) += 1;
				break;
			case 3: //protein death
				newstate.at(1) -= 1;
				break;
		}
		return newstate;
	 }
	 double stepmatrix(int i,int j){
		return identity_mat(i,j);
	 }
};

struct association_dissociation{
	float associationconst;
	float dissociationconst;
	
	int m1order;
	int m2order;
	
	static int const dimensions = 3;
	static int const Numprocs = 2;
	rate_list getRates(state_t const& state){
		int m1 = state.at(0);
		int m2 = state.at(1);
		int complex = state.at(2);
		
		float associate = associationconst*pow(m1,m1order)*pow(m2,m2order);
		float dissociate = dissociationconst*complex;
		
		rate_list rates = {associate,dissociate};
		return rates;
	}
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		switch(proc){
		case 0:
			newstate.at(0) -= m1order;
			newstate.at(1) -= m2order;
			newstate.at(2) += 1;
			break;
		case 1:
			newstate.at(0) += m1order;
			newstate.at(1) += m2order;
			newstate.at(2) -= 1;
			break;
		}
		return newstate;
	}
};



/*
 * The following systems have constructors which take their parameters as vectors of numbers
 * This is important for my work simulating experiments on perturbed systems:
 * it lets you initiate a state with an object you can mutate at random
 */

struct dual_reporting_cascade{
	/*
	 * System is the following
	 * 		    --->m1--->m2--->out
	 * 		  m0
	 * 			---|3--->m4--->out
	 * m0 will have constant birth, and all have linear deaths.
	 * Others have births linear to upstream molecule.
	 * Births denoted lambda,deaths denoted beta
	 * 
	 * The catch that's different from a true dual reporter:
	 * branch 3 is suppressed by the presence of 0. For now in the form
	 * of a reduction in birth rate (imagine 0 binding a needed
	 * component for production of 3, like a transcription factor)
	 */
	 std::vector<float> births;
	 std::vector<float> deaths;
	 static int const dimensions = 2;
	 static int const Numprocs = 10;
	 dual_reporting_cascade(std::vector<float> lambdas,std::vector<float> betas){
		births = lambdas;
		deaths = betas;
	 }
	 dual_reporting_cascade(){};
	 rate_list getRates(state_t const& state){
		
		float birth0 = births.at(0);
		float birth1 = births.at(1)*state.at(0);
		float birth2 = births.at(2)*state.at(1);
		float birth3 = births.at(3); 
		float birth4 = births.at(4)*state.at(3);
		float death0 = deaths.at(0)*state.at(0);
		float death1 = deaths.at(1)*state.at(1);
		float death2 = deaths.at(2)*state.at(2);
		float death3 = deaths.at(3)*state.at(3)*state.at(0);
		float death4 = deaths.at(4)*state.at(4);
		
		rate_list rates = {birth0,birth1,birth2,birth3,birth4,
							death0,death1,death2,death3,death4};
		return rates;
	}
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		if(proc < 5){
			newstate.at(proc) += 1;
		}
		else{
			newstate.at(proc%5) -= 1;
		};
		return newstate;
	}
	double stepmatrix(int i,int j){return identity_mat(i,j);}
};


struct dual_reporting_cascade_hill{
	/*
	 * System is the following
	 * 		    --->m1--->m2--->out
	 * 		  m0
	 * 			---|3--->m4--->out
	 * m0 will have constant birth, and all have linear deaths.
	 * Others have births linear to upstream molecule.
	 * Births denoted lambda,deaths denoted beta
	 * 
	 * The catch that's different from a true dual reporter:
	 * branch 3 is suppressed by the presence of 0. For now in the form
	 * of a reduction in birth rate (imagine 0 binding a needed
	 * component for production of 3, like a transcription factor)
	 */
	 std::vector<float> births;
	 std::vector<float> deaths;
	 double K;
	 double n;
	 static int const dimensions = 5;
	 static int const Numprocs = 10;
	 dual_reporting_cascade_hill(std::vector<float> lambdas,std::vector<float> betas,double k,double N){
		births = lambdas;
		deaths = betas;
		K = k;
		n = N;
	 }
	 dual_reporting_cascade_hill(){};
	 rate_list getRates(state_t const& state){
		
		float birth0 = births.at(0);
		float birth1 = births.at(1)*state.at(0);
		float birth2 = births.at(2)*state.at(1);
		float birth3 = births.at(3)/(1+pow(state.at(0)/K,n)); 
		float birth4 = births.at(4)*state.at(3);
		float death0 = deaths.at(0)*state.at(0);
		float death1 = deaths.at(1)*state.at(1);
		float death2 = deaths.at(2)*state.at(2);
		float death3 = deaths.at(3)*state.at(3);
		float death4 = deaths.at(4)*state.at(4);
		
		rate_list rates = {birth0,birth1,birth2,birth3,birth4,
							death0,death1,death2,death3,death4};
		return rates;
	}
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		if(proc < 5){
			newstate.at(proc) += 1;
		}
		else{
			newstate.at(proc%5) -= 1;
		};
		return newstate;
	}
	double stepmatrix(int i,int j){return identity_mat(i,j);}
};

struct dual_feedback{
		 std::vector<float> births;
	 std::vector<float> deaths;
	 static int const dimensions = 5;
	 static int const Numprocs = 10;
	 dual_feedback(std::vector<float> lambdas,std::vector<float> betas){
		births = lambdas;
		deaths = betas;
	 }
	 
	 rate_list getRates(state_t const& state){
		
		float birth0 = births.at(0)/pow(state.at(3)+1,0.5);
		float birth1 = births.at(1)*state.at(0);
		float birth2 = births.at(2)*state.at(1);
		float birth3 = births.at(3)/pow((state.at(0)+1),0.5); //The division is the birth rate suppression. Offset by 1 to avoid divide by zero and square root to make less severe
		float birth4 = births.at(4)*state.at(3);
		float death0 = deaths.at(0)*state.at(0);
		float death1 = deaths.at(1)*state.at(1);
		float death2 = deaths.at(2)*state.at(2);
		float death3 = deaths.at(3)*state.at(3);
		float death4 = deaths.at(4)*state.at(4);
		
		rate_list rates = {birth0,birth1,birth2,birth3,birth4,
							death0,death1,death2,death3,death4};
		return rates;
	}
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		if(proc < 5){
			newstate.at(proc) += 1;
		}
		else{
			newstate.at(proc%5) -= 1;
		};
		return newstate;
	}
};


struct dual_feedback_hill{
		 std::vector<float> births;
	 std::vector<float> deaths;
	 static int const dimensions = 5;
	 static int const Numprocs = 10;
	 float K;
	 float n;
	 char type;
	 dual_feedback_hill(std::vector<float> lambdas,std::vector<float> betas,float k,float N){
		births = lambdas;
		deaths = betas;
		K = k;
		n = N;
	 }
	 dual_feedback_hill(){}
	 rate_list getRates(state_t const& state){
		float birth3 = births.at(3)*state.at(0);
		float birth1 = births.at(1)*state.at(0);
		float birth2 = births.at(2)*state.at(1);
		//C++ pow does not take kindly to being asked by zero, which is really a limit. Code it by hand!
		float birth0;
		if(state.at(3)==0){birth0 = births.at(0);}
		else{birth0 = births.at(0)*pow(state.at(3),n)/(pow(K,n)+pow(state.at(3),n));}
		float birth4 = births.at(4)*state.at(3);
		float death0 = deaths.at(0)*state.at(0);
		float death1 = deaths.at(1)*state.at(1);
		float death2 = deaths.at(2)*state.at(2);
		float death3 = deaths.at(3)*state.at(3);
		float death4 = deaths.at(4)*state.at(4);
		
		rate_list rates = {birth0,birth1,birth2,birth3,birth4,
							death0,death1,death2,death3,death4};
		return rates;
	}
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		if(proc < 5){
			newstate.at(proc) += 1;
		}
		else{
			newstate.at(proc%5) -= 1;
		};
		return newstate;
	}
};

struct repressilator{
	float birth0,birth1,birth2;
	float death0,death1,death2;
	float K0,K1,K2;
	int n0,n1,n2;
	static int const dimensions = 3;
	static int const Numprocs = 6;
	rate_list getRates(state_t const& state){
		
		float birthrate0;
		if(state.at(2)==0){birthrate0 = birth0;}
		else{birthrate0 = birth0/(std::pow(state.at(2)/K0,n0) + 1);}
		float birthrate1;
		if(state.at(0)==0){birthrate1 = birth1;}
		else{birthrate1 = birth1/(std::pow(state.at(0)/K1,n1) + 1);}
		float birthrate2;
		if(state.at(1)==0){birthrate2 = birth2;}
		else{birthrate2 = birth2/(std::pow(state.at(1)/K2,n2) + 1);}
		
		float deathrate0 = death0*state.at(0);
		float deathrate1 = death1*state.at(1);
		float deathrate2 = death2*state.at(2);
		
		rate_list rates = {birthrate0,birthrate1,birthrate2,
							deathrate0,deathrate1,deathrate2};
		return rates;		
	}
	
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		if(proc < 3){
			newstate.at(proc) += 1;
		}
		else{
			newstate.at(proc%3) -= 1;
		};
		return newstate;
	}
};


struct poisson_repressor{
	float birth0,birth1;
	float death0,death1;
	float K;
	int n;
	static int const dimensions = 2;
	static int const Numprocs = 4;
	rate_list getRates(state_t const& state){
		
		float birthrate0 = birth0;
		float birthrate1;
		if(state.at(0)==0){birthrate1 = birth1;}
		else{birthrate1 = birth1/(std::pow(state.at(0)/K,n) + 1);}
		
		float deathrate0 = death0*state.at(0);
		float deathrate1 = death1*state.at(1);
		
		rate_list rates = {birthrate0,birthrate1,
							deathrate0,deathrate1};
		return rates;		
	}
	
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		if(proc < 2){
			newstate.at(proc) += 1;
		}
		else{
			newstate.at(proc%2) -= 1;
		};
		return newstate;
	}
};

struct AIF{
	float lambda;
	float theta;
	float mu;
	float beta;
	float C;

	AIF(){}
	AIF(float l,float t,float m,float b, float c){
		lambda = l;
		theta = t;
		mu = m;
		beta = b;
		C = c;
	}
	static int const dimensino = 3;
	static int const Numprocs = 5;
	rate_list getRates(state_t const& state){
		rate_list rates = {lambda*state.at(2),theta*state.at(0),mu,beta*state.at(0),C*state.at(1)*state.at(2)};
		return rates;
	}
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		
		switch(proc){
			case 0:
				newstate.at(0) += 1;
				break;
			case 1:
				newstate.at(1) += 1;
				break;
			case 2:
				newstate.at(2) += 1;
				break;
			case 3:
				newstate.at(0) -= 1;
				break;
			case 4:
				newstate.at(2) -= 1;
				newstate.at(1) -= 1;
				break;
		}
		
		return newstate;
	}
};

struct AIF_ext{
	std::vector<float> lambda; //
	float theta;
	float mu;
	std::vector<float> beta;
	float C;
	AIF_ext(){}
	AIF_ext(std::vector<float> l,float t,float m,std::vector<float> b, float c){
		lambda = l;
		theta = t;
		mu = m;
		beta = b;
		C = c;
	}
	static int const dimensions = 5;
	static int const Numprocs = 9;
	rate_list getRates(state_t const& state){
		rate_list rates = {lambda.at(0)*state.at(1),mu,lambda.at(1)*state.at(1),theta*state.at(0),lambda.at(2)*state.at(3),
							beta.at(0)*state.at(0),beta.at(1)*state.at(2),beta.at(2)*state.at(4),C*state.at(1)*state.at(3)};
		return rates;
	}
	state_t getState(int& proc,state_t const& currentstate){
		state_t newstate = currentstate;
		
		switch(proc){
			case 0:
				newstate.at(0) += 1;
				break;
			case 1:
				newstate.at(1) += 1;
				break;
			case 2:
				newstate.at(2) += 1;
				break;
			case 3:
				newstate.at(3) += 1;
				break;
			case 4:
				newstate.at(4) += 1;
				break;
			case 5:
				newstate.at(0) -= 1;
				break;
			case 6:
				newstate.at(2) -= 1;
				break;
			case 7:
				newstate.at(4) -= 1;
				break;
			case 8:
				newstate.at(3) -= 1;
				newstate.at(1) -= 1;
				break;
		}
		
		return newstate;
	}
};

#endif
