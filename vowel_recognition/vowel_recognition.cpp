// vowel_recognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <cstring>
#include <limits>
#include "stdlib.h"
#include "math.h"
#include "ctype.h"

#define N 320 //total sample in one frame
#define P 12 
#define F 5 // count of stable frames
#define PI 3.142857142857

//constants
//threshold multiplier value is taken after tuning
const double limit = 5000.0;

//global variables
double mx; //to find max value sample
double dcShift, nFactor; //dc shift and normalization factor

//array for samples, energy, framewise samples, storing tokhura distance 
double x[100000], energy[100000], steadyFrames[F][N], tokhuraDist[5];

//array for storing the values of ri, ai, ci, avg ci, allci, reference ci framewise
double R[F][P+1], A[F][P+1], C[F][P+1], avgCi[25][P+1], Ci[F][P+1], allCi[50][F][P+1], restoreCi[F][P+1];

//tokhura weights
double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//index helper to store the values of ci to allCi
int u = 0, v = 0, files = 0;

//vowels
char vowels[5] = {'a', 'e', 'i', 'o', 'u'};

//xSize - size of array x
//enSize - size of array energy
//start_i & end_i - start and end marker of steady frames
long int xSize, enSize, start_i, end_i;

//variable to find accuracy
int totalCorrect = 0, individualCorrect= 0;

//function to apply the Raised Sine Window in Ci of each frame
void raisedSinWindow(){
	long double sum=0;
	for(int f = 0; f<F; f++){
		for(int m=1;m<=P;m++){
			sum = (P/2)*sin((PI*m)/P);
			C[f][m]*=sum;	
		}	
	}
}

//function for applying hamming window to all the stable frames
void applyHammingWindow(){
	for(int i=0; i<F; ++i){
		for(int j=0; j<N; ++j){
			steadyFrames[i][j] *= 0.54-0.46*cos(2*PI*steadyFrames[i][j]/N-1);
		}
	}
}

//storing the values of C to 3D matrix Ci
void storeCito3D(){
	//saving ci values to 3D matrix
	for(int f=0; f<F; f++){
		for(v=0;v<P;v++){
			allCi[files][f][v+1]=C[f][v+1];
		}
	}
	files++;
}

// storing the avg ci values to file
void dumpAvgCi(){
	FILE *indi;
	char filename[30];
	int index = 0;
	
	for(int vl=0; vl<5; vl++){//loop over vowel to save file/vowel 
		sprintf(filename, "Reference/reference_ci_%c.txt", vowels[vl]);
		indi = fopen(filename, "w");
		//looping over frames
		for(int f=0; f<F; f++){
			//looping over p
			for(int p=0; p<P; p++){
				double sum = 0;
				//summing up all the file wise values and taking average
				for(int file=vl*10; file<(vl+1)*10; file++){
					sum += allCi[file][f][p+1];
				}
				sum /= 10.0;
				avgCi[index][p+1] = sum;
				fprintf(indi, "%lf ", sum);
			}
			index++;
			fprintf(indi, "\n");
		}
		printf("%s is created\n", filename);
		fclose(indi);
	}
}

//This function calulate the cepstral coeff Ci's
void calculate_Cis(){
	double sum=0;
	
	for(int f = 0; f<F; f++){
		C[f][0]=log(R[f][0]*R[f][0]);

		for(int m=1;m<=P;m++){
			sum=0;
			for(int k=1;k<m;k++){
				sum += (k*C[f][k]*A[f][m-k])/(m*1.0);
			}
			C[f][m]=A[f][m]+sum;
		}
	}
	
	//applying raised sin window on cis
	raisedSinWindow();
	//storing the ci values to 3D matrix
	storeCito3D();
}

// Function to apply Durbin Algorithm And Find The value of ai's 
void durbinAlgo(){
	double Alpha[13][13],E[13],K[13];
	double sum=0;


	for(int f = 0; f<F; f++){
		E[0] = R[f][0];
		for(int i=1;i<=P;i++){
			sum=0;
			for(int j=1;j<=i-1;j++){
				sum += Alpha[i-1][j]*R[f][i-j];	
			}
			
			K[i]=(R[f][i]-sum)/E[i-1];
				
			Alpha[i][i]=K[i];
		
			for(int j=1;j<=i-1;j++){
				Alpha[i][j]=Alpha[i-1][j] - K[i]*Alpha[i-1][i-j];
			}
		
			E[i]=(1-(K[i]*K[i]))*E[i-1];
		}

		//storing the ai values
		for(int i=1;i<=P;i++){
			A[f][i]= Alpha[P][i];
		}
	}
	
	//finding cepstral constants
	calculate_Cis();
}

//calculating the Ris values
void calculate_Ris(){
	
    //calculating Ris
	for(int f = 0; f<5; f++){
		for(int m =0; m<=P; m++){
			R[f][m] = 0;
			for(int k=0; k<N-m; k++){
				R[f][m] += steadyFrames[f][k]*steadyFrames[f][k+m];
			}
		}
	}
	//calling durbinAlgo to find ai values
	durbinAlgo();
}

//function to get dcshift value and set in global variable
double getDCShift(char *filename){

    long int sample_count = 0;
    FILE *fp;
    char line[80];

    //reading dc_shift.txt file
    fp = fopen(filename, "r");
    
    if(fp == NULL){
        printf("File not found\n");
        exit(1);
    }
    
	dcShift = 0;
    while(!feof(fp)){
        fgets(line, 80, fp);
        dcShift += atof(line);
        sample_count++;
    }
    dcShift /= sample_count;
    
    fclose(fp);
}

//function to setup the global variable like, max and nFactor
//max and nFactor depends on the vowel recording file and are used to do the normalization
void setupGlobal(char *filename){
    FILE *fp;
    long int totalSample = 0;
    char line[100];

    fp = fopen(filename, "r");
    if(fp == NULL){
        printf("Error opening file\n");
    }

    //get max value
    mx = 0;
    while(!feof(fp)){
        fgets(line, 100, fp);
        if(!isalpha(line[0])){
            totalSample++;
            if(mx < abs(atoi(line)))
                mx = abs(atoi(line));
        }
    }
    
    nFactor = (double)limit/mx;
   
    //setup dcShift
    getDCShift(filename);
    fclose(fp);
}

//marking the stable frames using STE
void framesMarker(){
	long int totalSample = 0, max_i = 0;
	int n = 0;
	double en = 0, m = 0;
	enSize = 0;
	//calculating the ste and marking the highest energy index
	while(totalSample < xSize){
		//if the n is equal to N than we have to store energy to the array
		if(n == N){
			//taking average
			en /= N;
			
			//checking whether it is max energy frame?
			if(m < en) m = en, max_i = enSize;

			energy[enSize++] = en; //store the energy value
			
			//resetting en and n to count for the next frame
			en = 0;
			n = 0;
		}

		n++;
		en += x[totalSample] * x[totalSample];
		totalSample++;
	}

	//marking start
	start_i = max_i > 2 ? (max_i-2)*N : 0;
	//marking end
	end_i = max_i <enSize - 3 ? (max_i+3)*N : enSize*N;
	
	int f = 0;
	//storing the frames into steadyFrames array
	for(int i = start_i, j=0; i<end_i; i++){
		steadyFrames[f][j++] = x[i];
		if(j == N) f++, j=0;
	}
}

//utility function which reads the files and normalizes it and store the values in array x
void read_file(char *filename){
	char line[70];
    FILE *ip;

	setupGlobal(filename);

	ip = fopen(filename, "r");

    if(ip == NULL) printf("Error in opening file %s\n", filename);

	//start storing from 0
	xSize = 0;
	//reading the values from input file, normalizing it and writing it to output file
    while(!feof(ip)){
        fgets(line, 100, ip);
        
        //input file may contain header, so we skip it
        if(!isalpha(line[0])){
            int y = atof(line);
            double normalizedX = floor((y-dcShift)*nFactor);
            //if(abs(normalizedX) > 1)
             x[xSize++] = normalizedX;
        }
    }
	
	fclose(ip);
	//find stable frames
	framesMarker();
	//apply hamming window on stable frames
	applyHammingWindow();
	//count Ris using stable frames
    calculate_Ris();
}

//driver function to execute training
void training(){
	printf("training started\n");
	char filenames[5][30] = {"Training/214101023_a_0.txt", "Training/214101023_e_0.txt", "Training/214101023_i_0.txt", "Training/214101023_o_0.txt", "Training/214101023_u_0.txt"};//

	for(int i=0; i<5; i++){//loop over vowels
		for(int j = 0; j<10; j++){//loop over files
			//read file
			read_file(filenames[i]);
			//next file
			filenames[i][21]++;
		}
	}
}

//fucntion which calculates the distance using dump Ci values of training set
double tokhuraDistance(FILE *ip){
	char line[3000];
	int f = 0;
	//read the reference file and store in restoredCi
	while(!feof(ip) && f<F){
		fgets(line, sizeof(line), ip);
		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &restoreCi[f][1], &restoreCi[f][2], &restoreCi[f][3], &restoreCi[f][4], &restoreCi[f][5], &restoreCi[f][6], &restoreCi[f][7], &restoreCi[f][8], &restoreCi[f][9], &restoreCi[f][10], &restoreCi[f][11], &restoreCi[f][12]);
		f++;
	}

	double finalDist = 0;
	for(int i=0; i<F; i++){
		double dist = 0;
		for(int p=1; p<=P; p++){
			double d = (C[i][p]- restoreCi[i][p]);
			dist += tokhuraWeights[p-1]*d*d;
		}
		finalDist += dist/(P*1.0); //taking average
	}
	return finalDist/(F*1.0); // taking average
}

//function to calculate the distance and making prediction
char calculateTokhura(){
	char filename[30];
	FILE *ip;
	double minDist = DBL_MAX;
	char predictedVowel;

	for(int i=0; i<5; i++){
		sprintf(filename, "Reference/reference_ci_%c.txt", vowels[i]);
		ip = fopen(filename, "r");
		if(ip == NULL) printf("Error in opening file %s\n", filename);
		//calculating distance from vowels[i]'s reference file
		double distance = tokhuraDistance(ip);
		//storing distance to print
		tokhuraDist[i] = distance;
		//checking whether it is minimum than the minimum distance we found?
		if(minDist > distance){
			//update incase new distance is minimum
			minDist = distance;
			//update the vowel to be returned
			predictedVowel = vowels[i];
		}
	}
	return predictedVowel;
}

void testing(){
	//testing files
	char filename[5][30] = {"Testing/214101023_a_10.txt", "Testing/214101023_e_10.txt", "Testing/214101023_i_10.txt", "Testing/214101023_o_10.txt", "Testing/214101023_u_10.txt"};
	//making files variable 0 to store the ci values in 3D array
	files = 0;

	//whether user wants to see the tokhura distance or not
	int choice;
	printf("testing of text files has started\n");
	printf("Do you want to see the report of all the tokhura weights? (1- yes | 0-no) - ");
	scanf("%d", &choice);

	for(int vl = 0; vl<5; vl++){
		individualCorrect = 0;
		for(int file=0; file<10; file++){
			//reading test file
			read_file(filename[vl]);
			//output
			printf("\nfile: %s recognized as--> ", filename[vl]);
			//tokhura calculation will return the predicted character
			char prediction = calculateTokhura();
			printf("%c\n", prediction);
			if(choice){
				printf("Tokhura weights distance from (a, e, i, o, u) = (%lf, %lf, %lf, %lf, %lf)\n", tokhuraDist[0], tokhuraDist[1], tokhuraDist[2], tokhuraDist[3], tokhuraDist[4]);
			}
			//if prediction is correct than increasing counts
			if(prediction == vowels[vl]) totalCorrect++, individualCorrect++;
			//going for next file
			filename[vl][21]++;
		}
		//printing the vowel recognition accuracy / vowel
		printf("********** Accuracy of vowel %c: %.2lf %% **********\n\n", vowels[vl], (individualCorrect/10.0)*100);
	}
	//printing overall accuracy
	printf("***************************************************\n");
	printf("***** Overall Accuracy of the system: %.2lf %% *****\n", (totalCorrect/50.0)*100);
	printf("***************************************************\n\n");
}

int _tmain(int argc, _TCHAR* argv[]){
    //training using 50 recordings of 10 recording/vowel
	training();

	//generating reference file
	dumpAvgCi();

	//testing of the text files
	testing();
	
	system("pause");
	return 0;
}