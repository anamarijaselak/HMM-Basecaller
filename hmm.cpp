#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <math.h>

#include <chrono>
#include <utility>
#include "logsumset.hpp"


typedef std::chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
//#define durationMili(a) std::chrono::duration_cast<std::chrono::miliseconds>(a).count()

#define timeNow() std::chrono::high_resolution_clock::now()

using namespace std;

//g++ -std=c++11 hmm.cpp -o hmm
///home/mare/minion-basecaller/notebook

//lambda = NC_001416.fa

int getdir (string dir, vector<string> &files)
{   
	unsigned char isFile =0x8;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if ( dirp->d_type == isFile && '.' != dirp->d_name[0]){
            files.push_back(string(dirp->d_name));
        }
        
    }
    closedir(dp);
    return 0;
}

void parseFile(vector<int> &v, string file){
	ifstream infile;
	string tmpstring;

	infile.open((file).c_str()); 

	int j = 0;
	while(getline(infile,tmpstring)) {
		std::vector<string> tokens;
		boost::algorithm::split(tokens, tmpstring, boost::is_any_of(" "));
		v[j] = (int) stod(tokens[0]);
		j++;
	}
	infile.close();
}

void getNeighbours(vector<int> &neighbours, string neighboursFile){
	parseFile(neighbours, neighboursFile);
}

void getNeighboursMask(vector<int> &neighboursMask, string neighboursMaskFile){
	parseFile(neighboursMask, neighboursMaskFile);
}

void getBNeighbours(vector<int> &bNeighbours, string bNeighboursFile){
	parseFile(bNeighbours, bNeighboursFile);
}

void getBNeighboursMask(vector<int> &bNeighboursMask, string bNeighboursMaskFile){
	parseFile(bNeighboursMask, bNeighboursMaskFile);
}

void getSixMers(vector<string> &sixmers, string poreModelFile){
	ifstream infile;
	string tmpstring;

	infile.open((poreModelFile).c_str()); 

	int j = 0;
	while(getline(infile,tmpstring)) {
		std::vector<string> tokens;
		boost::algorithm::split(tokens, tmpstring, boost::is_any_of("\t "));
		sixmers[j] = tokens[0];
		j++;
	}
	infile.close();
}

void getIndexesForBW(vector<int> &K, string indexesForBWFile){
	ifstream infile;
	string tmpstring;

	infile.open((indexesForBWFile).c_str()); 

	int j = 0;
	while(getline(infile,tmpstring)) {
		K[j] = stoi(tmpstring);
		j++;
	}
	infile.close();	
}

double getTransProbBetween(string sixmerX, string sixmerY, double pStay, double pSkip){
	double prob = 0;
	double p_skip_1 = pSkip / (1+pSkip);

	if(sixmerX.compare(sixmerY) == 0){
		prob += pStay;
	}

	if (sixmerX.compare(1,5, sixmerY, 0, 5) == 0){
		 prob += (0.25)*(1-pStay-pSkip);
	}

	//4 equal
	if(sixmerX.compare(2, 4, sixmerY, 0, 4) == 0){
		prob += 0.0625*p_skip_1;
	}

	//3 equal
	if(sixmerX.compare(3, 3, sixmerY, 0, 3) == 0){
		prob += pow(p_skip_1, 2) / (1u << (2 * 3));
	}

	//2 equal
	if(sixmerX.compare(4, 2, sixmerY, 0, 2) == 0){
		prob += pow(p_skip_1, 3) / (1u << (2 * 4));
	}

	//1 equal
	if(sixmerX.compare(4, 1, sixmerY, 0, 1) == 0){
		prob += pow(p_skip_1, 4) / (1u << (2 * 5));
	}

	prob += (pow(p_skip_1, 5) / (1.0 - p_skip_1)) / 4096;



	return prob;
}

void defineTransitionProbabilities(vector<double> &transProb, double pStay, double pSkip, vector<string> &sixMers){
	TimeVar t1=timeNow();
	for(int i=0; i<4096; i++){
		for(int j=0; j<4096; j++){
			//vjer prijelaza iz sixmer_x u sixmer_y
			transProb[i*4096 + j] = getTransProbBetween(sixMers[i], sixMers[j], pStay, pSkip);
		}
	}
	cout << duration(timeNow()-t1) << " seconds for transition probs.\n";
}

void scalePoreModelsMOM(vector<double> &emissionModel, vector<double> &events){
	//Method of moments applied:
    double eventsMean = 0.0;
    for(int i = 1; i< events.size(); i+=4){
    	eventsMean += events[i];
    }
    eventsMean = eventsMean / (events.size()/4.0);
    for(int i=0; i<4096*5; i+=5){
    	emissionModel[i] = eventsMean;
    }
}

void scalePoreModelsEM(vector<double> &emissionModel, double shift, double scale, double drift, 
		double var, double scale_sd, double var_sd){
	//[(mean dev meanInv devInv lambda) x 4096] = emission model

	for(int i=0; i < 4096; i++){
		emissionModel[i*5] = emissionModel[i*5] * scale + shift;
		emissionModel[i*5 + 1] = emissionModel[i*5 + 1]*var;
		emissionModel[i*5 + 2] = emissionModel[i*5 + 2] * scale_sd;
		emissionModel[i*5 + 4] = emissionModel[i*5 + 4] * var_sd;
		emissionModel[i*5 + 3] = pow(pow(emissionModel[i*5 + 2], 3.0) / emissionModel[i*5 + 4], .5);
	}

}

vector<double> parsePoreModels(string poreModelFile, vector<double> &events ){
	//[(mean dev meanInv devInv lambda) x 4096]
	vector<double> emissionModel(4096*5);
	ifstream infile;
	string tmpstring;
	infile.open((poreModelFile).c_str()); 
	int j = 0;
	while(getline(infile,tmpstring)) {
		std::vector<string> tokens;
		boost::algorithm::split(tokens, tmpstring, boost::is_any_of("\t "));
		for(int k=1; k< 6; k++){
			emissionModel[j*5 + k - 1] = stod(tokens[k]);
		}
			
		j++;
	}
	infile.close();
	scalePoreModelsMOM(emissionModel, events);
	return emissionModel;
}

void getEmissionProbability(double mean, double variance, vector<double> &emissionModel, vector<double> &emissionProbsForT){
	static const double inv_sqrt_2pi = 0.3989422804014327;
	int j = 0;
	for(int i = 0; i < 4096*5; i+=5){
		double meanGaussian = emissionModel[i];
		double stDevGaussian = emissionModel[i+1];
		double a = (mean - meanGaussian) / stDevGaussian;
		double gaussianPdf = (inv_sqrt_2pi / stDevGaussian) * exp(-0.5f * a * a);

		double meanInvGaussian = emissionModel[i+2];
		double lamda = emissionModel[i+4];
		double b = (variance - meanInvGaussian) / meanInvGaussian;
		double inverseGaussianPdf = (inv_sqrt_2pi * sqrt(lamda /(variance*variance*variance))) * exp(-0.5f * b * b * (lamda / variance));
    	
    	emissionProbsForT[j] = gaussianPdf * inverseGaussianPdf;
    	j++;
	}  
}

void getEmissionProbabilityLog(double mean, double variance, vector<double> &emissionModel, vector<double> &emissionProbsForT){
	static const double log_2pi = std::log(2.0 * M_PI);
	int j = 0;
	for(int i = 0; i < 4096*5; i+=5){
		double meanGaussian = emissionModel[i];
		double stDevGaussian = emissionModel[i+1];

		if (stDevGaussian == 0.0)
		{
		    stDevGaussian = 0.01;
		}

		double a = (mean - meanGaussian) / stDevGaussian;
		double gaussianPdf = - log(stDevGaussian) - (log_2pi + a * a) / 2.0;


		double meanInvGaussian = emissionModel[i+2];
		double lambda = emissionModel[i+4];
		double b = (variance - meanInvGaussian) / meanInvGaussian;
		double inverseGaussianPdf = (log(lambda) - log_2pi - 3.0 * log(variance) - lambda * b * b / variance) / 2.0;	
    	
    	emissionProbsForT[j] = gaussianPdf + inverseGaussianPdf;
    	j++;
	}   
}

void runViterbi(vector<double> &events, vector<int> &states, vector<double> &transProb, 
		vector<double> &initialProb, vector<double> &emissionModel, vector<int> &neighbours,
		int start, int end){

	int numberOfEvents = end - start;
	//events = M x 4 

	//za svako stanje u svakom trenutku racunamo maks vjer
	vector<double> trellis;
	vector<int> backpt;

	trellis.assign(4096*numberOfEvents, 1.0f);
	//for every state only one kmer and each state has 4 properties
	backpt.assign(4096*numberOfEvents, -1);

	vector<double> emissionProbsForT(4096);

	//build trellis and backtracking diagrams

	//trellis[0] = initialProb * b_j0
	getEmissionProbability(events[start*4+1], events[start*4+2], emissionModel, emissionProbsForT);
	for(int i=0; i<4096; i++){
		trellis[i] = (1.0f/4096.0f) * emissionProbsForT[i];
		
	}

	int eventIndex = 1;
	//for each event
	for(int i=start+1; i< end; i++){
		double sum = 0.0f;
		getEmissionProbability(events[i*4 + 1], events[i*4 + 2], emissionModel, emissionProbsForT);
		//for each state
		for(int j=0; j<4096; j++){
			double max = 0.0f;
			int backptIndex = -1;
			//for each possible neighbour
			for(int k=0; k<21; k++){
				int neighbour = neighbours[21*j + k];
				//max(k) {transProb(neighbour, j)*trellis(neighbour, eventIndex-1)} - neighbour of j
				double x = transProb[neighbour*4096 + j] * trellis[(eventIndex-1)*4096 + neighbour];
				if(x > max){
					max = x;
					backptIndex = neighbour;
				}
			}
			trellis[eventIndex*4096 + j] = max * emissionProbsForT[j];
			backpt[eventIndex*4096 + j] = backptIndex;

			sum += trellis[eventIndex*4096 + j];
		}
		for(int j=0; j<4096; j++){
			trellis[eventIndex*4096 + j] = trellis[eventIndex*4096 + j] / sum;
		}
		eventIndex++;
		
	}

	//backtrack to get states
	//last state with highest probability
	int lastState = 0;
	double maxProb = trellis[eventIndex*4095];
	for(int i=1; i<4096; i++){
		if(transProb[eventIndex*4095 +i] < maxProb){
			maxProb = trellis[eventIndex*4095 + i];
			lastState = i;
		}
	}

	states[start + eventIndex-1] = lastState;

	for(int i = eventIndex -1; i > 0 ; i--){
		states[start + i-1] = backpt[4096 * i + lastState];
		lastState = states[start + i-1];
	}


}

void decodeEvents(vector<double> &events, vector<int> &states, vector<double> &transProb, 
		vector<double> &initialProb, vector<double> &emissionModel, vector<int> &neighbours){

	TimeVar t1=timeNow();
	int numberOfEvents = events.size() / 4;
	int chunks = numberOfEvents/1000;

	if((chunks * 1000) < numberOfEvents){
		chunks+=1;
	}

	//1000 by 1000 events - each event has 4 parameters -> multiply by 4
	for(int i = 0; i < chunks; i++){
		int start = i*1000;
		int end = start + 1000;
		if(end > numberOfEvents){
			end = numberOfEvents;
		}
		
		runViterbi(events, states, transProb, initialProb, emissionModel, neighbours, start, end);
		cout <<"Viterbi, chunk " << i <<". \n";
	}
	
	cout << duration(timeNow()-t1) << " seconds for Viterbi algorithm.\n";
}

string decodeStates(vector<int> &states, vector<string> &sixmers){
	ostringstream bases;
	if(states[0] != -1){
		bases << sixmers[states[0]];
	}
	
	for(int i=1; i<states.size(); i++){
		if(states[i] == -1){
			continue;
		}
		if(states[i-1] == -1){
			string currentSixmer = sixmers[states[i]];
			bases << currentSixmer;
			continue;
		}
		string previousSixmer = sixmers[states[i-1]];
		string currentSixmer = sixmers[states[i]];
		if(previousSixmer.compare(currentSixmer) != 0){
			if (previousSixmer.compare(1,5, currentSixmer, 0, 5) == 0){
		 		bases << currentSixmer.at(5);
			}else if(previousSixmer.compare(2,4, currentSixmer, 0, 4) == 0){
				bases << currentSixmer.at(4);
				bases << currentSixmer.at(5);
			}
		}
	}
	return bases.str();
}

int gE(vector<double> &A, vector<double> &X, vector<double> &B){
	if(A[0] == 0){
		return -1;
	}
	double r = B[1] - (A[3]*B[0])/A[0];
	double l = A[4] - (A[3]*A[1])/A[0];
	if( r == 0){
		return -1;
	}
	X[1] = l / r;
	X[2] = 0;
	X[0] = (B[0] - A[1]*X[1])/A[0];

	for(int i=0; i< 3; i++){
    	cout << "Short ge : X[ "<<i<<"]" << " = "<<X[i] <<"\n";
    }

	return 1;
}

//solving linear equations for A = 3x3
//return 1 if there is a solution, return -1 if there is no solution
int gaussEliminations(vector<double> &A, vector<double> &X, vector<double> &B){
	//compute scaling factors for scaled partial pivoting
	TimeVar t1 = timeNow();
	vector<double> scalingFactors(3);
	for(int i=0; i<9; i+=3){
		double max = A[i];
		for(int j=i; j<i+3;j++){
			if(A[j] > max){
				max = A[j];
			}
		}
		scalingFactors[i/3] = max;
	}

	//gaussEliminations with scaled partial pivoting for Ax = B
	for(int i=0; i<3; i++){
		int p = i;
		//pivot - diagonal element divided by scaling factor
		double p_val = abs(A[i*3 + i] / scalingFactors[i]); 

		//check if we need to interchange rows
		for(int j=i+1; j<3; j++){
			double p_val2 = abs(A[j*3 + i] / scalingFactors[j]);
			if(p_val2 > p_val){
				p_val = p_val2;
				p = j;
			}
		}

		//check if the system has a solution
		if (p_val == 0)
        {
            return -1;
        }

        //exchange rows i and p
        if( p > i){
        	cout << "Changeam " << p << " sa " << i << "\n";
        	double temp = A[i*3];
        	A[i*3] = A[p*3];
        	A[p*3] = temp;

        	temp = A[i*3 + 1];
        	A[i*3 + 1] = A[p*3 + 1];
        	A[p*3 + 1] = temp;

        	temp = A[i*3 + 2];
        	A[i*3 + 2] = A[p*3 + 2];
        	A[p*3 + 2] = temp;

        	temp = B[i];
        	B[i] = B[p];
        	B[p] = temp;

        	temp = scalingFactors[i];
        	scalingFactors[i] = scalingFactors[p];
        	scalingFactors[p] = temp;

        }

        //eliminate variables
        for(int j = i + 1; j < 3; j++){
        	double m = A[j*3 + i] / A[i*3 + i];
        	A[j*3 + i] = 0;
	        for (int k = i + 1; k < 3; k++)
	        {	
	            A[j*3 + k] =  A[j*3 + k] - m * A[i*3 + k];
	        }
	        B[j] = B[j] - m * B[i];
        }
	}

	//upper triangular system
    X[2] = B[2] / A[8];
    X[1] = (B[1] - A[5]*X[2]) / A[4];
    X[0] = (B[0] - A[1]*X[1] - A[2]*X[2]) / A[0];

    for(int i=0; i< 3; i++){
    	cout << "X[ "<<i<<"]" << " = "<<X[i] <<"\n";
    }

    cout << duration(timeNow()-t1) << " seconds for Gaussian eliminations algorithm.\n";
	return 1;
}

void getBackwardValuesLog(vector<double> &backwardMatrix, vector<double> &transProb,
			vector<double> &events, vector<double> &emissionModel, vector<int> &backwardN, vector<int> &bMask){

	TimeVar t1=timeNow();

	typedef logsum::logsumset<double> LogSumSet_Type;
	LogSumSet_Type s(false);

	vector<double> emissionProbsForT(4096);

	for(int i=0; i<4096; i++){
		backwardMatrix[(events.size()/4-1)*4096 +i] = 0;
	}

	//for each event backwards
	for(int i = (events.size() / 4 -1); i > 0; i--){
		if(i % 1000 == 0){
			cout << "i = " << i;
			cout << " Bwm : " << backwardMatrix[i*4096+1] << "\n";
		}
		int index = i -1;
		getEmissionProbabilityLog(events[i*4 +1], events[i*4 +2], emissionModel, emissionProbsForT);
		for(int j=0; j<4096; j++){
			s.clear();
			for(int k=0; k < 21; k++){
				if(bMask[k] == 0){
					continue;
				}
				int bNeighbour = backwardN[j*21 +k];
				s.add(log(transProb[j*4096 + bNeighbour])+ emissionProbsForT[bNeighbour] + backwardMatrix[i*4096 + bNeighbour]);
			}

			backwardMatrix[index*4096 + j] += s.val();
		}
	}

	cout << duration(timeNow()-t1) << " seconds for Backward algorithm.\n";
}

void getBackwardValues(vector<double> &backwardMatrix, vector<double> &transProb,
			vector<double> &events, vector<double> &emissionModel, vector<int> &backwardN, vector<int> &bMask){

	TimeVar t1=timeNow();

	vector<double> emissionProbsForT(4096);

	for(int i=0; i<4096; i++){
		backwardMatrix[(events.size()/4-1)*4096 +i] = 1;
	}

	//for each event backwards
	for(int i = (events.size() / 4 -1); i > 0; i--){
		double sumScaling = 0;
		if(i % 1000 == 0){
			//cout << "i = " << i;
			//cout << " Bwm : " << backwardMatrix[i*4096+1] << "\n";
		}
		int index = i -1;
		getEmissionProbability(events[i*4 +1], events[i*4 +2], emissionModel, emissionProbsForT);
		for(int j=0; j<4096; j++){
			double sum = 0;
			for(int k=0; k < 21; k++){
				if(bMask[k] == 0){
					continue;
				}
				int bNeighbour = backwardN[j*21 +k];
				sum += transProb[j*4096 + bNeighbour]*emissionProbsForT[bNeighbour]*backwardMatrix[i*4096 + bNeighbour];
			}

			backwardMatrix[index*4096 + j] += sum;
			sumScaling += sum;
		}
		for(int j=0; j<4096; j++){
			backwardMatrix[index*4096 + j] = backwardMatrix[index*4096 + j] / sumScaling;
		}
	}

	cout << duration(timeNow()-t1) << " seconds for Backward algorithm.\n";
}

//returnlog(P(X | HMM))
double getForwardValuesLog(vector<double> &forwardMatrix, vector<double> &transProb, vector<double> &events, 
			vector<double> &emissionModel, vector<int> &neighbours, vector<int> &mask){
	
	TimeVar t1=timeNow();

	typedef logsum::logsumset<double> LogSumSet_Type;
	LogSumSet_Type s(false);

	vector<double> emissionProbsForT(4096);
	getEmissionProbabilityLog(events[1], events[2], emissionModel, emissionProbsForT);

	double logNmbStates = log(4096);
	for(int i = 0; i < 4096; i++){
		forwardMatrix[i] = emissionProbsForT[i] - logNmbStates;
	}

	//for each event
	for(int i = 1; i < events.size() / 4; ++i){
		if(i % 1000 == 0){
			//cout << "i = " << i;
			//cout << " Fwm : " << forwardMatrix[i*4095+1] << "\n";
		}
		getEmissionProbabilityLog(events[i*4 +1], events[i*4 +2], emissionModel, emissionProbsForT);	
		for(int j=0; j< 4096; ++j){
			//for each possible neighbour
			s.clear();
			for(int k=0; k<21; k++){
				if(mask[21*j + k] == 0){
					continue;
				}
				int neighbour = neighbours[21*j + k];				
				s.add(forwardMatrix[(i-1)*4096 + neighbour]+log(transProb[neighbour*4096 + j]));
			}
			forwardMatrix[i*4096 + j] = s.val() + emissionProbsForT[j];			
		}
	}

	s.clear();

	for(int i = ((events.size() / 4)-1)*4096; i < (events.size()/4)*4096; i++){
		s.add(forwardMatrix[i]);
	}
	cout << duration(timeNow()-t1) << " seconds for Forward algorithm.\n";
	cout << "Val = " << s.val() << "\n";
	return s.val();
}

//returnlog(P(X | HMM))
double getForwardValues(vector<double> &forwardMatrix, vector<double> &transProb, vector<double> &events, 
			vector<double> &emissionModel, vector<int> &neighbours, vector<int> &mask){
	
	TimeVar t1=timeNow();

	vector<double> emissionProbsForT(4096);
	getEmissionProbability(events[1], events[2], emissionModel, emissionProbsForT);
	for(int i = 0; i < 4096; i++){
		forwardMatrix[i] = emissionProbsForT[i]*(1.0/4096.0);
	}

	//for each event
	for(int i = 1; i < events.size() / 4; ++i){
		double sumScaling = 0;
		if(i % 1000 == 0){
			//cout << "i = " << i;
			//cout << " Fwm : " << forwardMatrix[i*4095+1] << "\n";
		}
		getEmissionProbability(events[i*4 +1], events[i*4 +2], emissionModel, emissionProbsForT);	
		for(int j=0; j< 4096; ++j){
			//for each possible neighbour
			double sum = 0;
			for(int k=0; k<21; k++){
				if(mask[21*j + k] == 0){
					continue;
				}
				int neighbour = neighbours[21*j + k];				
				sum+=forwardMatrix[(i-1)*4096 + neighbour]*transProb[neighbour*4096 + j];
			}
			forwardMatrix[i*4096 + j] =sum * emissionProbsForT[j];	
			sumScaling +=forwardMatrix[i*4096 + j];	
		}
		for(int j=0; j<4096;j++){
			forwardMatrix[i*4096 + j] = forwardMatrix[i*4096 + j] /sumScaling;
		}
	}

	double posteriori = 0;

	for(int i = ((events.size() / 4)-1)*4096; i < (events.size()/4)*4096; i++){
		posteriori += forwardMatrix[i];
	}



	cout << duration(timeNow()-t1) << " seconds for Forward algorithm.\n";
	cout << "Val = " << posteriori << "\n";
	return posteriori;
}


//return -1 if we can scale and 1 if we can't
int trainScalingParametersLog(vector<double> &forwardMatrix, vector<double> &backwardMatrix, 
		vector<double> &events, vector<double> &emissionModel, double logPrData, vector<double> &X, vector<double> &scalingParams){

	//X[0] = a, X[1] = b, X[2] = c, scalingParams[0] = d, scalingParams[1] = v scalingParams[2] = u
	//a = shift
	//b=scale
	//c=drift
	//d=var^2
	//v=scale_sd
	//u=var_sd

	TimeVar t1=timeNow();

	std::vector<double> A(9);
	std::vector<double> B(3);

	A.assign(9, 0.0);
	B.assign(3, 0.0);

	double D       = 0.0; 
    double V_numer = 0.0; 
    double V_denom = 0.0;
	double U_pos = 0.0; 

	for(int i=0; i<events.size()/4; i++){
		std::vector<double> s;
		s.assign(3, 0.0f);

		std::vector<double> l;
		l.assign(3, 0.0f);

		double sum_jp_ij = 0;
		for(int j=0; j<4096;j++){
			double p_ij = std::exp(forwardMatrix[i*4096+j] + backwardMatrix[i*4096 +j] - logPrData);
            double term_s0 = p_ij / (emissionModel[j*5 +1] * emissionModel[j*5 +1]);
            double term_s1 = term_s0 * emissionModel[j*5];
            double term_s2 = term_s1 * emissionModel[j*5];
            double term_l0 = p_ij * emissionModel[j*5+4];
            double term_l1 = term_l0 / emissionModel[j*5];
            double term_l2 = term_l1 / emissionModel[j*5];
            s[0] += term_s0;
            s[1] += term_s1;
            s[2] += term_s2;
            l[0] += term_l0;
            l[1] += term_l1;
			l[2] += term_l2;
			sum_jp_ij += p_ij;
		} // states 

		A[0] += s[0];
		A[1] += s[1];
        A[4] += s[2];
        //A[2] += s[0] * events[i*4 + 3];
        A[2] = 0.0;
        //A[5] += s[1] * events[i*4 + 3];
        A[5] = 0.0;
        //A[8] += s[0] * events[i*4 + 3] * events[i*4 + 3];
        A[8] = 1.0;

        //events[i*4 + 3] = start_i events[i*4 + 1] = mean_i events[i*4 + 2] = st_dev_i
        B[0] += s[0] * events[i*4 + 1];
        B[1] += s[1] * events[i*4 + 1];
        //B[2] += s[0] * events[i*4 + 1] * events[i*4 + 3];
        B[2] = 0;

        D += s[0] * events[i*4 + 1] * events[i*4 + 1];
        V_numer += l[2] * events[i*4 + 2];
        V_denom += l[1];
		U_pos += l[0] / events[i*4 + 2];

	} // events 

	A[3] = A[1];
    A[6] = A[2];
	A[7] = A[5];

	cout << duration(timeNow()-t1) << " seconds for training scaling parameters - before Gauss eliminations.\n";

	vector<double> ACopy(9);
	vector<double> BCopy(3);

	for(int i=0; i<9; i++){
		ACopy[i] = A[i];
	}

	for(int i=0; i< 3; i++){
		BCopy[0] = B[0];
	}

	if(gaussEliminations(ACopy, X, BCopy) == -1){
		cout << "Matrix was singular!\n";
		return -1;
	}

	//if(gE(ACopy, X, BCopy) == -1){
	//	cout << "Matrix was singular!\n";
	//	return -1;
	//}

	//sqrt(var) = scaling factor for standard deviation, in files we have that written
	double d = sqrt((D + X[0] * X[0] * A[0]
                          + X[1] * X[1] * A[4]
                          + X[2] * X[2] * A[8]
                          + 2.0 * X[0] * X[1] * A[1]
                          + 2.0 * X[0] * X[2] * A[2]
                          + 2.0 * X[1] * X[2] * A[5]
                          - 2.0 * (X[0] * B[0] + X[1] * B[1]+ X[2] * B[2]))/(double)4096.0);

	double v = V_numer / V_denom;
	double u = (double)4096 / (U_pos - V_denom / v);

	scalingParams[0] = d;
	scalingParams[1] = v;
	scalingParams[2] = u;

	cout << duration(timeNow()-t1) << " seconds for training scaling parameters including Gauss eliminations.\n";
	return 1;
}


void getNeighboursForK(vector<int> &K, vector<int> &neighbours, vector<string> &sixMers){

	for(int i=0; i< K.size();i++){
		string k = sixMers[i].substr(1,5);
		string n1 = k + "A";
		string n2 = k + "T";
		string n3 = k + "C";
		string n4 = k + "G";
		for(int j=0; j<4096; j++){
			if(sixMers[j] == n1){
				neighbours[i*4] = j;
			}
			if(sixMers[j] == n2){
				neighbours[i*4+1] = j;
			}
			if(sixMers[j] == n3){
				neighbours[i*4+2] = j;
			}
			if(sixMers[j] == n4){
				neighbours[i*4+3] = j;
			}
		}

	}
}


vector<double> trainTransitionParametersLog(vector<double> &forwardMatrix, vector<double> &backwardMatrix, 
		vector<double> &events, vector<double> &emissionModel, 
		vector<int> K, vector<string> &sixMers, double logPrData, double pStep, double pStay){

	cout << "Tu1\n";

	typedef logsum::logsumset<double> LogSumSet_Type;
	LogSumSet_Type denom(false);
	LogSumSet_Type pStayNum(false);
	LogSumSet_Type pSkipNum(false);

	vector<double> emmisionProbsForT(4096);
	vector<int> neighboursForK(K.size()*4);

	double logStep = log(0.25*pStep);
	double logStay = log(pStay);

	cout << "Tu2\n";
	getNeighboursForK(K, neighboursForK, sixMers);

	cout << "Tu3\n";

	vector<double> result(2);

	for(int i=0; i< events.size()/4 -1; i++){
		getEmissionProbabilityLog(events[(i+1)*4 +1], events[(i+1)*4 +2], emissionModel, emmisionProbsForT);
		for(int j=0; j<K.size(); j++){
			LogSumSet_Type v(false);
			int k = K[j];
			double logPosteriori = forwardMatrix[i*4096 + k]+backwardMatrix[i*4096 + k]-logPrData;
			denom.add(logPosteriori);

			double jointProbLogStay = forwardMatrix[i*4096 + k]+logStay + emmisionProbsForT[k]+ backwardMatrix[(i+1)*4096 + k]-logPrData;
			pStayNum.add(jointProbLogStay);

			v.add(jointProbLogStay);
			for(int l = 0; l<4; l++){
				if(sixMers[k] == sixMers[neighboursForK[k*4+l]]){
					continue;
				}
				int neighK = neighboursForK[k*4+l];
				double jointProbLogSkip = forwardMatrix[i*4096 + k]+logStep + emmisionProbsForT[neighK]+ backwardMatrix[(i+1)*4096 + neighK]-logPrData;
				v.add(jointProbLogSkip);
			}
			double probSkip = exp(logPosteriori) - exp(v.val());

			pSkipNum.add(log(probSkip));

		}
	}

	cout << "Tu4\n";
	result[0] = exp(pStayNum.val() - denom.val());
	result[1] = exp(pSkipNum.val() - denom.val());
	return result;
}

vector<double> trainTransitionParameters(vector<double> &forwardMatrix, vector<double> &backwardMatrix, 
		vector<double> &events, vector<double> &emissionModel, 
		vector<int> K, vector<string> &sixMers, double pStep, double pStay, double pSkip){
	double denom = 0.0;
	double pStayNum = 0.0;
	double pSkipNum = 0.0;

	double prData= 0;

	for(int i=0; i < K.size(); i++){
		//cout << (events.size()/4-1)*4096 + K[i] << "-----------------" <<  ;
		prData += forwardMatrix[(events.size()/4-1)*4096 + K[i]];
	}

	cout << "Pr data " << prData << "\n";
	vector<double> emmisionProbsForT(4096);
	vector<int> neighboursForK(K.size()*4);

	getNeighboursForK(K, neighboursForK, sixMers);

	vector<double> result(2);

	for(int i=0; i< events.size()/4 -1; i++){
		getEmissionProbability(events[(i+1)*4 +1], events[(i+1)*4 +2], emissionModel, emmisionProbsForT);
		for(int j=0; j<K.size(); j++){
			int k = K[j];
			double posteriori = forwardMatrix[i*4096 + k]*backwardMatrix[i*4096 + k]/prData;
			denom+=posteriori;

			double jointProbStay = forwardMatrix[i*4096 + k]*pStay*emmisionProbsForT[k]*backwardMatrix[(i+1)*4096 + k]/prData;
			pStayNum+=jointProbStay;

			pSkipNum -=jointProbStay;

			for(int l = 0; l<4; l++){
				if(sixMers[k] == sixMers[neighboursForK[k*4+l]]){
					continue;
				}
				int neighK = neighboursForK[k*4+l];
				double jointProbSkip = forwardMatrix[i*4096 + k]*0.25*pStep*emmisionProbsForT[neighK]*backwardMatrix[(i+1)*4096 + neighK]/prData;
				pSkipNum-=jointProbSkip;
			}
		}
	}

	cout << "Denom = " << denom << " pStayNum " << pStayNum << " pSkipNum " << pSkipNum << "\n";
	pSkipNum += K.size();
	result[0] = pStayNum / denom;
	result[1] = pSkipNum / denom;
	if(result[0] <= 0 || result[0] >= 1){
		result[0] = pStay;
	}
	if(result[0] <= 0 || result[0] >= 1){
		result[0] = pSkip;
	}
	return result;
}

int main(int argc, char const *argv[])
{	
	string readsDirectory = argv[1];
	string poreModelFile = argv[2];
	string neighboursFile = argv[3];
	string bNeighboursFile = argv[4];
	string neighboursMaskFile = argv[5];
	string bNeighboursMaskFile = argv[6];
	string indexesForBWFile = argv[7];
	string resultsFile = argv[8];

	//string readsDirectory = "/home/mare/Radna površina/lambda_R9_6_7_16/ProcessedFiles";
	//string poreModelFile = "/home/mare/Radna površina/Projekt/r9_250bps.nucleotide.6mer.template.model";
	//string neighboursFile = "/home/mare/Radna površina/Projekt/Neighbours.txt";
	//string bNeighboursFile = "/home/mare/Radna površina/Projekt/BackwardNeighbours.txt";
	//string neighboursMaskFile = "/home/mare/Radna površina/Projekt/Mask.txt";
	//string bNeighboursMaskFile = "/home/mare/Radna površina/Projekt/BackwardMask.txt";
	//string indexesForBWFile = "/home/mare/Radna površina/Projekt/indexesForBW";
	std::vector<double> initialProb;
	std::vector<double> transProb(4096*4096);
	std::vector<string> sixmers(4096);

	std::vector<int> neighbours(86016);
	std::vector<int> bNeighbours(86016);
	std::vector<int> nMask(86016);
	std::vector<int> bNMask(86016);
	std::vector<int> indexesForBW(2160);

	initialProb.assign(4096, 1.0/4096.0);

	double pStay = 0.1;
	double pSkip = 0.3; 

	getSixMers(sixmers, poreModelFile);
	getNeighbours(neighbours, neighboursFile);
	getBNeighbours(bNeighbours, bNeighboursFile);
	getNeighboursMask(nMask, neighboursMaskFile);
	getBNeighboursMask(bNMask, bNeighboursMaskFile);
	getIndexesForBW(indexesForBW, indexesForBWFile);
	defineTransitionProbabilities(transProb, pStay, pSkip, sixmers);

	ifstream infile;
	string tmpstring;

	vector<string> readsNames = vector<string>();
    getdir(readsDirectory,readsNames);

    //cout << readsNames.size() << "\n";

    for (unsigned int i = 0;i < readsNames.size();i++) {
    	//1 read
       	infile.open((readsDirectory+"/"+readsNames[i]).c_str()); 
       	cout << "File : " << (readsDirectory+"/"+readsNames[i]).c_str() << "\n";
       	cout << "Template : " << poreModelFile << "\n";
       //read first line to see how many events there are
       	//infile.open((readsDirectory+"/"+"Andrej_HC_20160706_FNFAD12873_MN17271_sequencing_run_lambda_R9_6_7_16_21932_ch251_read23752_strand.fast5").c_str());
       	getline(infile,tmpstring);
       	int numberOfEvents = stoi(tmpstring);

    	//[(lengt mean stdv start) x numberOfEvents]
    	vector<double> events(numberOfEvents*4);
    	int j = 0;
		while(getline(infile,tmpstring)) {
			std::vector<string> tokens;
			boost::algorithm::split(tokens, tmpstring, boost::is_any_of(" "));
			for(int k=0; k< 4; k++)
				events[j*4 + k] = stod(tokens[k]);
			j++;
		}
		infile.close();	

		vector<int> states(events.size() / 4);
		vector<double> emissionModel = parsePoreModels(poreModelFile, events);

		double oldLogPrData = 0.0;

		//TRAINING
		for(int trainingRound = 0; trainingRound < 1; trainingRound++){
			cout << "---------TRAINING_ROUND "<< trainingRound << "-----------------\n";
			cout << "Forward running ...\n";
			vector<double> forwardMatrixLog(events.size() * 1024);
			double logPrData = getForwardValuesLog(forwardMatrixLog, transProb, events, emissionModel, neighbours, nMask);		
			
			vector<double> forwardMatrix(events.size() * 1024);
			double prData = getForwardValues(forwardMatrix, transProb, events, emissionModel, neighbours, nMask);
			
			cout << "Backward running ...\n";
			
			vector<double> backwardMatrix(events.size() * 1024);			
			getBackwardValues(backwardMatrix, transProb, events, emissionModel, bNeighbours, bNMask);	

			cout << "Log pr data = " << logPrData << "\n";	

			if(oldLogPrData < logPrData && logPrData > 0.0){
				vector<double> backwardMatrixLog(events.size() * 1024);
				getBackwardValuesLog(backwardMatrixLog, transProb, events, emissionModel, bNeighbours, bNMask);
				cout << "Scaling parameters training...\n";
				vector<double> X(3);
				vector<double> scalingParams(3);

				int res = trainScalingParametersLog(forwardMatrixLog, backwardMatrixLog, events, emissionModel, logPrData, X, scalingParams);
				if( res == -1){
					continue;
				}
				cout << "Scaling pore models...\n";
				scalePoreModelsEM(emissionModel, X[0], X[1], X[2], scalingParams[0], scalingParams[1], scalingParams[2]);
				oldLogPrData = logPrData;

			}
			cout << "Running Baum Welch...\n";
			vector<double> results = trainTransitionParameters(forwardMatrix, backwardMatrix, events, emissionModel,  indexesForBW, 
																sixmers, 1-pStay-pSkip, pStay, pSkip);

			defineTransitionProbabilities(transProb, results[0], results[1], sixmers);
			cout << "Result of BW: " << "pStay = " << results[0] << " pSkip = " << results[1] << "\n";
		}		
		//FINISHED TRAINING

		//DECODING
		cout << "Viterbi running ...\n";
		decodeEvents(events, states, transProb, initialProb, emissionModel, neighbours);

		cout << "Decoding states...\n";
		
		string resultingRead = decodeStates(states, sixmers);

		cout << "Writing results...\n";
		ofstream out(resultsFile + to_string(i)+string(".fa"));
		out << readsNames[i].c_str() << "\n";
		out << resultingRead;

		out.close();
		//break;
    }
	return 0;
}