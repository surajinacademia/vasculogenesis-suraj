#!/bin/bash

#******************************************************************************************************
#This program creates folders with different simulation parameters and runs the program. It also checks 
#if the code has executed the simulation code and moves to different parameter.
#******************************************************************************************************



#Reads the simulation parameter range file and checks if it exists in the folder or not
INPUT=parameters.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

#Reads the inputs from the parameters file 
while read IF IR spring AA density angle
do

d='IF_'$IF'_IR_'$IR'_spring_'$spring'_AA_'$AA'_density_'$density'_angle_'$angle

if [ -d "${d}" ]; then
echo "Code Running or Data Exits. Runing the Next parameters in line"
elif [ ! -d "${d}" ]; then
echo "Data DOES NOT Exits. Runing the Next parameters in queue"
mkdir $d                  #Make the simulation folder 
cd $d

for i in {1..6}           #Samples for ensembles
do
mkdir ${i}
cd ${i}

mkdir data                #make a data directory

cat>"main.cpp"<<END
/******************************************************************************************************
Which was originally written by Katie Copehgen during her PhD at UC Merced Physics dept. She is now at Princeton University as Postdoc. 

Then it was passed down to Farnaz Golnaraghi during her PhD at UC Merced, Physics. Now she works at Univeristy of Chicago. 

The major contribution on editing the code is by Mikahl Banwarth-Kuhn who is a Post-doc at Suzzane Sindi's Lab.

Now it's Suraj Sahu's turn to work on this code for his PhD. He is just running this code for different simulation parameters.

    This simulates a group of agents that are connected to neighboring agents by springs 
    which can be fairly long range. Agents can intercalate other agents by interrupting 
    the space inbetween the agents connected by the spring. The agents have angular noise. 
\******************************************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include <vector>
#include <random>
using namespace std;

#define N 225        //Number of agents. [0-10000]
#define AA 0        //Alignment interaction magnitude. [0-10]
#define V 0.01        //Agent velocity for self propulsion. [0-0.1]
#define K ${spring}       //Spring interaction strength. [0-0.1]
#define NOISE 0.01    //Noise (variance of gaussian distribution for random angular noise). [0-2*pi]
#define RS 0.0        //Variance in agent sizes, for gaussian distribution. [0-RA (defined below)]
#define GAP ${IF}        //The amount of overlap to break springs. [0-1]

#define L 25.0     //Size of the circle which the agents are initiallized within. [sqrt(N) ish]
#define RA 1        //Average size of the agents, sets lengthscale. [1]
#define T 1        //Number of timesteps the simulation is run for. [100-100000]
#define vid_steps 1
#define SAMPLES 1    //Number of samples run. Must be 1 to output files for making videos. [1-100+]
#define MAX_X 45.0
#define MIN_X -45.0
#define MAX_Y 45.0
#define MIN_Y -45.0

#define Neighbor_region ${IR}
#define anisotropic_adhesion_on 1          //control whether we see iso or aniso adhesion
#define anisotropic_adhesion_angle M_PI/${angle}      //size of the adhesion window rangee [0-pi/2]
#define lower_window cosf(anisotropic_adhesion_angle/2+M_PI)
#define upper_window cosf(anisotropic_adhesion_angle/2)
#define testing_aniso   0          //if this is 1 initial conditions of cell will be a line
#define _USE_MATH_DEFINES

#define lattice_nonrandom 1
#define lattice_start 1 // if 1 then cells will start in lattice formation
#define lattice_step_size 6.0*RA
#define lattice_width 84
#define lattice_height 84
#define lattice_neg_x -lattice_width/2 //left lower bound of lattice
#define lattice_pos_x lattice_width/2  //right upper bound of lattice
#define lattice_neg_y -lattice_height/2 //bottom lower bound of lattice
#define lattice_pos_y lattice_height/2 //top upper bound of lattice

#include "particle.cpp"		//Particle class.
#include "clustering.cpp"	//Groups agents into clusters.
#include "opramcl.cpp"		//Calculate the properties of the agnet clusters.

int main(){
	//Seed the random number generator.
	long int letseed=time(NULL);
	srand (letseed);

	float t1=clock();	//Save starting time to calculate the program runtime with.

	FILE *f00;	//For making videos.
	FILE *f01;	//Store parameters.
	f01 = fopen("params.txt", "w");
	FILE *f02;	//Cluster data at simulation finish.
	f02=fopen("data/clusters.dat","w");
	FILE *f03;	//Average order vs. t.
	f03=fopen("data/ordervst.dat","w");
	FILE *f04;	//Agent data at simulation finish.
	f04=fopen("data/agents.dat","w");

	vector<particle> agents(N, particle());		//Vector for storing all the agents.
	for (int s=0; s<SAMPLES; s++){			//Loop through the samples.
		//Initialize agents, random position and velocity direction.
        //float theta_test = 0;
        //float new_x = cosf(theta_test);
        //float new_y = sinf(theta_test);
        //int newParticleID = 0;
        //for (int i=0;i<N;i++){
            //cout << newParticleID << endl;
            //if(testing_aniso){//if this is turned on then particles will start in a straight line
               // agents[i].init(new_x,new_y,newParticleID);
            //}else{//otherwise particles will be randomly seeded
                //agents[i].init(newParticleID);
            //}
            //newParticleID = newParticleID + 1;
            //theta_test = theta_test + 2*M_PI/10;
		//}
        //Initialize agents, random position and velocity direction.
        if(lattice_start == 1){
            float new_x = lattice_neg_x;
            float new_y = lattice_neg_y;
            int newParticleID = 0;
            for (int i=0;i<N;i++){
                //cout << newParticleID << endl;
                agents[i].init(new_x,new_y,newParticleID);
                newParticleID = newParticleID + 1;
                if(new_x < lattice_pos_x){
                    new_x = new_x + lattice_step_size;
                }else{
                    new_y = new_y + lattice_step_size;
                    new_x = lattice_neg_x;
                }
		    }
        }else if (lattice_nonrandom == 1){
            float new_x = lattice_neg_x;
            float new_y = lattice_neg_y;
            int newParticleID = 0;
            std::default_random_engine e1(32767);
            std::uniform_real_distribution<float> u(0,1);
            float new_ang;
            for (int i=0;i<N;i++){
                new_ang =  2*M_PI*u(e1);
                //cout << newParticleID << endl;
                cout << new_ang << endl;
                agents[i].init(new_x,new_y,new_ang,newParticleID);
                newParticleID = newParticleID + 1;
                if(new_x < lattice_pos_x){
                    new_x = new_x + lattice_step_size;
                }else{
                    new_y = new_y + lattice_step_size;
                    new_x = lattice_neg_x;
                }
		    }

        }

		for (int t=0; t<T; t++){		//Loop through time.

			//Once enough time has passed to stablize the cluster, spread the sizes, introducing the polydispersity gradually over 100 timesteps.
			if (t>(int)3*(sqrtf((float)N)/((float)(2.0*K))*RA)&&t<(int)3*(sqrtf(N)/((float)(2.0*K))*RA)+100){
				for (int i=0;i<N;i++)
					agents[i].spreadSizes();
					//printf("I am here\n");
			}

			//Create neighbor list for all agents.
			for (int i=0; i<N; i++){
				vector<int> temp;	//Temporarily stores the neighbor list.
				temp.reserve(N);	//Pre reserve space to speed up program.

				for (int j=0; j<N; j++){	//Loop through all other particles to check if they should be added to neighborlist.
					float dx=agents[j].getX()-agents[i].getX();
					float dy=agents[j].getY()-agents[i].getY();
					float d=sqrtf(dx*dx+dy*dy);
					if (d<=agents[i].getSize()+agents[j].getSize()+Neighbor_region*RA){
						temp.push_back(j);  //If i and j are close enough together add j to the neighborlist,right now they have to be 2*RA apart
					}
				}
				agents[i].setNeigh(temp); //Set the neighborlist to temp.
				temp.clear();
				agents[i].setOldAng();  //Also save their angles for use in interactions.
			}


			//Loop throgh every particle to calculate interactions.
			for (int i=0; i<N; i++){
                vector<int> aniso_temp;
                aniso_temp.reserve(N);

                int my_id = agents[i].getParticleID();
				//cout << my_id << endl;
                float Ax=0, Ay=0, A, Sx=0, Sy=0, newx, newy, dx, dy, d, newang;	//For storing all the interactions, to calculate new directions for agent i.
				for (int j=0; j<agents[i].getNeigh().size(); j++){		//Loop through j: all of agent i's neighbors.
                    int neighbor_id = agents[agents[i].getNeigh()[j]].getParticleID();
					//Calculate separation vector between agent i and i's j'th neighbor to use for spring.
					float dx=agents[agents[i].getNeigh()[j]].getX()-agents[i].getX();
					float dy=agents[agents[i].getNeigh()[j]].getY()-agents[i].getY();
					float d=sqrtf(dx*dx+dy*dy);
					//Make sure it doesnt blow up when they are too close.
					if (d<0.001){
						d=0.001;
                    }
					//Spring interaction.
					int kspringtemp=1;	//kspringtemp is for storing whether the spring is interrupted by other cells or not.
					//this loop is for anisotropic adhesion
                    if(anisotropic_adhesion_on==1){
                            //We determine if adhesion interaction will happen based on orientation of cell i and cell j.
                            //We do this by taking the dot product of the cell-to-cell vector (center_j - center_i) and
                            //the mean orientation vector of cell i and cell j
                            //Case 1 :When the cell-to-cell vector and mean orientation vector are parallel the 
                            //dot product will equal 1. We want to allow a little wiggle room, i.e. there is a window of ahesion
                            //points on the cells, so we accept any angle greater than min_window_size.
                            //Case 2: if the two cell orientation vectors are anti-parellel you will be dotting with the zero vector. 
                            //To fix this we do the folllowing. We have an if statement to catch when the mean orientation is zero and
                            //we pass this connection along since this case is allowed.
                            //Otherwise, the  cell-to-cell vector and mean orientation will be close to anti-parallel but not exactly
                            //and we will get that the dot product will be close to 0. To keep a reasonable window range
                            //we accept all values less than max_window_size.
                            float orientation_x_curr_particle = 0;
                            float orientation_y_curr_particle =  0;
                            float orientation_x_neighbor = 0;
                            float orientation_y_neighbor = 0;
                            float ang_diff_x = 0;
                            float ang_diff_y = 0;
                            orientation_x_curr_particle = cosf(agents[i].getOldAng());
                            orientation_y_curr_particle = sinf(agents[i].getOldAng());
                            orientation_x_neighbor = cosf(agents[agents[i].getNeigh()[j]].getOldAng());
                            orientation_y_neighbor = sinf(agents[agents[i].getNeigh()[j]].getOldAng());
                            ang_diff_x = orientation_x_neighbor - orientation_x_curr_particle;
                            ang_diff_y = orientation_y_neighbor - orientation_y_curr_particle;
                            //antiparallel catch 
                            if((ang_diff_x == 0)&&(ang_diff_y == 0)){
                                kspringtemp = 1;
                                aniso_temp.push_back(neighbor_id);
                            }else{
                                //compute center to center vec
                                //dx computed above is x value of center to center vec
                                //dy computed above is y value of center to center vec
                                //d computed above is the magnitude of the center to center vec
                                //mean  orientation vec
                                ang_diff_x = ang_diff_x*.5;
                                ang_diff_y = ang_diff_y*.5;
                                //get magnitude of mean orientation vec
                                float mag_ang_diff_vec = sqrt(ang_diff_x*ang_diff_x + ang_diff_y*ang_diff_y);
                                //dot product is (a_x*b_x + a_y*b_y)/(mag_a*mag_b)
                                float dot_product = (dx*ang_diff_x + dy*ang_diff_y)/(d*mag_ang_diff_vec);
                                if((lower_window < dot_product < upper_window)){
                                    kspringtemp = 0;
                                }else{
                                    aniso_temp.push_back(neighbor_id);
                                }
                            }
                    }//end anisotropic adhesion loop
                    float dx2, dy2, d2;	//For storning separation vector to other neighbors than j.
					if (my_id!=neighbor_id){		
						float avgr=(agents[i].getSize()+agents[agents[i].getNeigh()[j]].getSize());
						//this loop is for GAP interaction which we probably won't consider
                        for (int k=0; k<agents[i].getNeigh().size(); k++){	//Loop through all other neighbors of i that aren't j, they shallt be called k. 
							int second_neighbor_id = agents[agents[i].getNeigh()[k]].getParticleID();
                            if (neighbor_id!=second_neighbor_id && my_id!=second_neighbor_id){				//Only do the thing if they are 3 different agents.
								//Separation vector between agent i and it's k'th neighbor.
								dx2=agents[agents[i].getNeigh()[k]].getX()-agents[i].getX();	
								dy2=agents[agents[i].getNeigh()[k]].getY()-agents[i].getY();
								d2=sqrtf(dx2*dx2+dy2*dy2);
								//b basically records the size of the space between agent i and the j'th neighbor, at the relative distance to the k'th neighbor.
								float b=agents[i].getSize()+(agents[agents[i].getNeigh()[j]].getSize()-agents[i].getSize())/(d*d)*(dx2*dx+dy2*dy);
								//Then check if the k'th neighbor is within the GAP parameter fraction of the space between agent i and the j'th neighbor. If it is then interrupt the spring by setting kspringtemp to 0.
								if (sqrtf(d2*d2-powf(dx*dx2/d+dy*dy2/d,2))<agents[agents[i].getNeigh()[k]].getSize()+GAP*b&&(dx*dx2+dy*dy2)>0){
									kspringtemp=0;
								}
							}
						}//end gap interaction loop
						
                        //Alignment interaction.
						Ax+=cosf(agents[agents[i].getNeigh()[j]].getOldAng());
						Ay+=sinf(agents[agents[i].getNeigh()[j]].getOldAng());
						// Ax += 0;
						// Ay += 0;
						if (d<avgr){	//Always have the repulsive part of the spring interaction.
							kspringtemp=1;
                        }
						//Spring interaction.
						Sx+=-kspringtemp*(avgr-d)*(dx/d);
						Sy+=-kspringtemp*(avgr-d)*(dy/d);
					}
				}
                float total_adh_neighs = agents[i].getNeigh().size();
                if(total_adh_neighs > 0){
                    Ax = Ax/total_adh_neighs;
                    Ay = Ay/total_adh_neighs;
                    //cout << Ax << Ay << endl;
                }
                agents[i].setAnisoNeighs(aniso_temp);
                aniso_temp.clear();
				//Parameters u and v just for calculating the gaussion random variable for noise.
				float u=rand()/(float)RAND_MAX;
				float v=rand()/(float)RAND_MAX;
				//New direction of propulsion for agent i is a gaussian random variable with variance NOISE, plus the angle that comes from the alignment interaction, and it's own travel direction.
				agents[i].setAng(fmodf((NOISE*sqrtf(-2*log(u))*cosf(2*M_PI*v)+atan2(AA*Ay+sinf(agents[i].getAng()),AA*Ax+cosf(agents[i].getAng()))),2*M_PI));
				//Set the velocity of the agent to the propulsion with magnitude V in the propulsion direction stored in Ang, and the spring interaction with magnitude K.
				agents[i].setVx(V*cosf(agents[i].getAng())+K*Sx);
				agents[i].setVy(V*sinf(agents[i].getAng())+K*Sy);

			}

			//Move agents according to their new directions.
			for (int i=0; i<N; i++){
				agents[i].updatePos();
                agents[i].boundaryCheck();
            }
			int vidtime=vid_steps;		//How many timesteps to skip between each frame for the video.
			if (SAMPLES==1 && t%vidtime==0){
				char buff[32];		//buff saves the name of the file to store this frame data in.
				snprintf(buff, sizeof(char)*32, "video/frames/fr%06d", t);
				f00=fopen(buff, "w");	//Files for saving video frame data to make videos (vid.py).
				//Loop through each agent to put data into video frame file.
				for (int i=0;i<N;i++){
					//Find the relative velocity of each agent with all of it's neighbors (will set the color of agents in the video).
					float relv=0;
					for (int j=0;j<agents[i].getNeigh().size();j++){
						//Find how far agent i has moved since the last video frame.
						float mxi=agents[i].getX()-agents[i].getOldX();
						float myi=agents[i].getY()-agents[i].getOldY();
						//Find how far the j'th neighbor has moved since last video frame.
						float mxj=agents[agents[i].getNeigh()[j]].getX()-agents[agents[i].getNeigh()[j]].getOldX();
						float myj=agents[agents[i].getNeigh()[j]].getY()-agents[agents[i].getNeigh()[j]].getOldY();
						relv+=sqrtf(powf((mxj-mxi),2)+powf((myj-myi),2));	//Relative movement of agent i and all neighbors.
					}
					relv=relv/(agents[i].getNeigh().size()*V*vidtime);	//Normalize relative velocity..
					//Write agent data into file for video frame structure 4 columns: relative velocity	X position	Y position	agent size.
					fprintf(f00, "%f	%f	%f	%f	%f\n",relv,agents[i].getX(),agents[i].getY(),agents[i].getAng(), agents[i].getSize());
				}
				//Save the current agent positions to use for calculating relative velocities of the next frame.
				for (int i=0; i<N; i++)
					agents[i].storePos();
				fclose(f00);

			}
            //if(t%100 == 0){
            //    cout << "Tim step is: " << t << endl;
           // }
            for(int i=0; i<N;i++){
                for(int j=0;j<agents[i].getAnisoNeighs().size();j++){
                    //cout.flush() << "I am particle " << agents[i].getParticleID() << " and my adhesion neighbor is " << agents[agents[i].getNeigh()[j]].getParticleID() << endl;
                }
            }
		}	//End of time loop.

		//Calculate properties of each cluster.
		int labels[N], clusternum;
		clusternum = clustering(agents, labels);		//Function to assign a cluster label to each agent, and returns the total number of clusters.
		for (int i=0; i<N; i++)
			agents[i].setLabel(labels[i]);

		float Order[N],MedianR[N],AvgR[N],StdR[N];		//Arrays for storing the properties of each cluster, only clusternum elements will be saved into these arrays.
		int Size[N];
		opramcl(agents, Size, Order, MedianR, AvgR, StdR, clusternum);	//Calculate and store the size, order, median, mean, and standard deviation of the size of agents within each cluster.

		for (int i=0; i<clusternum; i++){	//Loop through the clusters and store the data about them.
			fprintf(f02,"%d	%f	%f	%f	%f\n", Size[i], Order[i], MedianR[i], AvgR[i], StdR[i]);
		}
		for (int i=0; i<N; i++){	//Save the individual position, direction, size, cluster, and neighbor lists for each agent. 
			fprintf(f04,"%f	%f	%f	%f	%d	%d	",agents[i].getX(),agents[i].getY(),agents[i].getAng(),agents[i].getSize(),agents[i].getLabel(),(int)agents[i].getNeigh().size());
			for (int j=0; j<agents[i].getNeigh().size(); j++)
				fprintf(f04,"%d	",agents[i].getNeigh()[j]);
			fprintf(f04,"\n");
		}
	}

	float t2=clock();
	//Write all of the system parameters and run time to params.txt file.
	fprintf(f01, "N = 	%d \nv = 	%f \nra = 	%f\nrs =	%f\na =	%f\nk =	%f\nnoise = 	%f\nsamples = 	%d\nmaxtime =	%d\nseed = 	%ld \n runtime = %f", N, V, (float)RA, (float)RS, (float)AA, (float)K, (float)NOISE, SAMPLES, T, letseed, (t2-t1)/(float)CLOCKS_PER_SEC);
	fclose(f01);
	return 0;
}


END

wait

cat>"opramcl.cpp"<<END
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

void opramcl(vector<particle> agents, int Size[], float Order[], float MedianR[], float AvgR[], float StdR[], int clusternum){
    float txvel[N], tyvel[N], avgr[N], stdr[N];        //For storing the total velocity, average and standard deviation of size to calculate the order parameters.
    int i, size[N];                        //Store the size (number of particles in) cluster.

    vector< vector<float> > sizesformed(clusternum); //Sizes of the agents in each cluster for calculating median of the agent size for each cluster.

    //Initialize things to zero for however many clusters there are.
    for (i=0; i<clusternum; i++){
        txvel[i]=0;
        tyvel[i]=0;
        size[i]=0;
        avgr[i]=0;
        stdr[i]=0;
    }

    for (i=0; i<N; i++){    //agents[i].getLable() returns the index for the cluster that agent i is in.
        //Calculate total velocities for finding order parameter.
        txvel[agents[i].getLabel()]+=cosf(agents[i].getAng());
        tyvel[agents[i].getLabel()]+=sinf(agents[i].getAng());
        //Total number of particles in the cluster.
        size[agents[i].getLabel()]++;
        //Average size of the particles in the cluster, and a list of all the sizes for finding the median.
        avgr[agents[i].getLabel()]+=agents[i].getSize();
        sizesformed[agents[i].getLabel()].push_back(agents[i].getSize());
    }

    for (i=0; i<N; i++){
        //Standard deviation of particle size.
        stdr[agents[i].getLabel()]+=powf(agents[i].getSize()-avgr[agents[i].getLabel()]/(float)size[agents[i].getLabel()],2.0);
    }

    //Loop through all the clusters to finalize calculations of parameters and store them in the arrays to be accessed in main.cpp.
    for (i=0; i<clusternum; i++){
        if (txvel[i]*txvel[i]+tyvel[i]*tyvel[i]>0)
            Order[i]=(sqrtf(txvel[i]*txvel[i]+tyvel[i]*tyvel[i])/(float)size[i]);        //Order paramter is magnitude of the sum of the velocities divided by the size of the cluster.
        Size[i]=(size[i]);        //The size is the number of particles in the cluster.

        size_t n = sizesformed[i].size();        //Last index of the cluster, so that when the vector of particle sizes for this cluster is sorted, the n/2 element will be the median.
        sort(sizesformed[i].begin(), sizesformed[i].end());
        //If there are an even number of particles in the cluster then median is the after of the middle 2.
        if (n%2==0)
            MedianR[i]=float(sizesformed[i][n/2 -1] + sizesformed[i][n/2])/2.0;
        else        //If there are an odd number then it is the middle particle.
            MedianR[i]=sizesformed[i][n/2];

        AvgR[i]=(avgr[i]/(float)size[i]);        //Average size of particles in the cluster.
        StdR[i]=(sqrtf(stdr[i]/(float)size[i]));    //Standard deviation of particle sizes in the cluster.
    }
    return;
}


END

wait

cat>"particle.cpp"<<END
#include <iostream>

using namespace std;
class particle {
    private:
        //particle ID
        int particleID;
        //Position of particle.
        float X;
        float OldX;
        float Y;
        float OldY;
        //Velocity of particle.
        float Vx;
        float Vy;
        //Direction of self propulsion of particle.
        float Ang;
        float OldAng;
        //Current size of particle, and equilibrium size after spread has happened.
        float Size;
        float Sizeeq;
        //Lable of the cluster that particle is in.
        int Label;
        //Neighborlist.
        vector<int> Neigh;
        //Anisotropic Adhesion Neighbors
        vector<int> AnisoNeighs;
        
    public:
        //Set the variables.
        void setAng(float ang){Ang=ang;}
        void setOldAng(void){OldAng=Ang;}
        void setLabel(int label){Label=label;}
        void setNeigh(vector<int> neigh){Neigh=neigh;}
        void setAnisoNeighs(vector<int> aniso_neighs){AnisoNeighs = aniso_neighs;}
        void setVx(float vx){Vx=vx;}
        void setVy(float vy){Vy=vy;}
        
        //Calculate and run some internal things.
        void spreadSizes();
        void init(int newParticleID);
        void init(float new_x,float new_y,int newParticleID);
        void init(float new_x, float new_y, float newang, int newParticleID);
        void updatePos(void);
        void boundaryCheck(void);
        void storePos(void);

        //Return the variables for use.
        int   getParticleID(void){return particleID;}
        float getX(void){return X;}
        float getOldX(void){return OldX;}
        float getY(void){return Y;}
        float getOldY(void){return OldY;}
        float getVx(void){return Vx;}
        float getVy(void){return Vy;}
        float getAng(void){return Ang;}
        float getOldAng(void){return OldAng;}
        float getSize(void){return Size;}
        vector<int> getNeigh(void){return Neigh;}
        vector<int> getAnisoNeighs(void){return AnisoNeighs;}
        int getLabel(void){return Label;}


        particle();
};
particle::particle(void){
    particleID = 0;
    X=0;
    Y=0;
    Ang=0;
    OldAng=0;
    Size=0;
    Label=0;
}
void particle::init(float new_x, float new_y, int newParticleID){
    X = new_x;
    Y = new_y;
    particleID = newParticleID;
    Ang=2*M_PI*(float)rand()/(float)RAND_MAX;        //Random starting direction between 0 and 2*pi.
    OldAng=Ang;
    float rmu=log(RA/sqrtf((RS*RA*RS*RA/((float)(RA*RA)))+1));    //Mu and Sigma for generating size of the particle from a log normal distribution.
    float rsi=sqrtf(log(((RS*RA*RS*RA)/((float)(RA*RA)))+1));
    float u=(float)rand()/(float)RAND_MAX;                //u and v are random variables for generating a normal distribution.
    float v=(float)rand()/(float)RAND_MAX;
    Sizeeq=exp(rmu+rsi*sqrtf(-2*log(u))*cosf(2*M_PI*v));        //Generate an equilibrium size for particle from a log normal distribution.
    Size=RA;                            //Set the initial size to the average of the system.
    Label=-1;                            //Label everything -1 initially for later assignment to clusters.
}

void particle::init(float new_x, float new_y, float newang, int newParticleID){
    X = new_x;
    Y = new_y;
    particleID = newParticleID;
    Ang=newang;        //Start with a fixed direction between 0 and 2*pi.
    OldAng=Ang;
    float rmu=log(RA/sqrtf((RS*RA*RS*RA/((float)(RA*RA)))+1));    //Mu and Sigma for generating size of the particle from a log normal distribution.
    float rsi=sqrtf(log(((RS*RA*RS*RA)/((float)(RA*RA)))+1));
    float u=(float)rand()/(float)RAND_MAX;                //u and v are random variables for generating a normal distribution.
    float v=(float)rand()/(float)RAND_MAX;
    Sizeeq=exp(rmu+rsi*sqrtf(-2*log(u))*cosf(2*M_PI*v));        //Generate an equilibrium size for particle from a log normal distribution.
    Size=RA;                            //Set the initial size to the average of the system.
    Label=-1;                            //Label everything -1 initially for later assignment to clusters.
}

void particle::init(int newParticleID){ //For initialization.
    float dO=0;    //dO stores the distance from the particle to the origin.
    do {        //Assign a random X and Y to the particle within a box of size 2*L centered on the origin.
        X=2*L*(float)rand()/(float)RAND_MAX-L;
        Y=2*L*(float)rand()/(float)RAND_MAX-L;
        dO=sqrtf(X*X+Y*Y);    //Calculate distance to the origin.
    } while (dO>L);            //If it is outside of the circle regenerate a random point untill you get one inside of the circle with radius L.
    particleID = newParticleID;
    Ang=2*M_PI*(float)rand()/(float)RAND_MAX;        //Random starting direction between 0 and 2*pi.
    OldAng=Ang;
    float rmu=log(RA/sqrtf((RS*RA*RS*RA/((float)(RA*RA)))+1));    //Mu and Sigma for generating size of the particle from a log normal distribution.
    float rsi=sqrtf(log(((RS*RA*RS*RA)/((float)(RA*RA)))+1));
    float u=(float)rand()/(float)RAND_MAX;                //u and v are random variables for generating a normal distribution.
    float v=(float)rand()/(float)RAND_MAX;
    Sizeeq=exp(rmu+rsi*sqrtf(-2*log(u))*cosf(2*M_PI*v));        //Generate an equilibrium size for particle from a log normal distribution.
    Size=RA;                            //Set the initial size to the average of the system.
    Label=-1;                            //Label everything -1 initially for later assignment to clusters.
}
void particle::updatePos(void){         //Move the position by the velocity amount.
    X=X+Vx;
    Y=Y+Vy;
}
void particle::boundaryCheck(void){
    if(X > MAX_X){
        X = MIN_X;
    }
    if(X < MIN_X){
        X = MAX_X;
    }
    if(Y > MAX_Y){
        Y = MIN_Y;
    }
    if(Y < MIN_Y){
        Y = MAX_Y;
    }

}
void particle::storePos(void){            //Save the current position, for use in calculating relative velocities for videos.
    OldX=X;
    OldY=Y;
}
void particle::spreadSizes(){            //Change the size by a small amount towards the equilibrium size, so after 100 changes it will be the right size.
    Size+=(Sizeeq-RA)/(100.0);
}

END
wait 

cat>"clustering.cpp"<<END
#include <iostream>

void addNewNeighbor(int i, int currentLabel, vector<particle> agents, int labels[]){
    int j;
    if (labels[i]==-1){        //If particle i isnt in a cluster yet add it to the current cluster.
        labels[i]=(currentLabel);
        for (j=0; j<agents[i].getNeigh().size(); j++){        //Go through all of particle i's neighbors and add them to the same cluster as particle i.
            if (agents[i].getNeigh()[j]!=i)            //Also add all of particle j's neighbors to the cluster, do this recursively until all of the particles in this cluster have a label that is not -1 anymore.
                addNewNeighbor(agents[i].getNeigh()[j], currentLabel, agents, labels);
        }
    }
    return;        //Once all of the neighbors and neighbors neighbors etc... are added to the current cluster return to start the next cluster.
}

int clustering(vector<particle> agents, int labels[]){        //Assign a label for each particle by cluster.
    int j, currentLabel=0;                    //Start with label 0.
    for (j=0; j<N; j++)
        labels[j]=-1;                    //Set all the labels for all particles to -1 to begin with.

    for (j=0; j<N; j++){                    //Loops through all the particles to assign them to a cluster
        if (labels[j]==-1){                //If particle j is not yet in a cluster, add it to the current cluster.
            addNewNeighbor(j, currentLabel, agents, labels);        //Add particle j and it's neighbors to the current cluster.
            currentLabel++;                //Increase the current cluster label by 1 to start labeling the particles in the next cluster in the next loop through.
        }
    }
    return currentLabel;        //currentLabel at the end is the total number of cluster, so return that.
}

END
wait


mkdir video                   #make a video directory
cd video                      #go to the video directory
mkdir frames                  #make a frames directory

cat>"vid.py"<<END

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.cm as cm
import math
import os


files=[]
for (dirpath, dirnames, filenames) in os.walk("frames/"):
	files.extend(["frames/"+filenamest for filenamest in filenames if "DS_Store" not in filenamest])
	break
files = sorted(files)
print(files)
fig = plt.figure(figsize=(20,20))
# ax = plt.figure(figsize=(20,20))
ax = plt.gca()
plt.clf()

def readData(i):
	#print(i)
	plt.clf() ## CLEARS THE PROJECTIONS
	data=[]
	data=np.loadtxt(files[i])

	data = list(zip(*data))
	#print (data)
	x=data[1]
	y=data[2]
	colors=data[0]
	size=[100*i*i for i in data[4]]
	# ## BACKGROUND GRID -> BLACK
	# COMx=np.mean(np.array(x))
	# COMy=np.mean(np.array(y))
	# crossx=int(math.ceil(COMx/50.0))*50
	# crossy=int(math.ceil(COMy/50.0))*50
	# plt.plot([crossx,crossx],[crossy-200,crossy+200],"k--",alpha=0.5)
	# plt.plot([crossx-200,crossx+200],[crossy,crossy],"k--",alpha=0.5)
	# plt.plot([crossx+50,crossx+50],[crossy-200,crossy+200],"k--",alpha=0.5)
	# plt.plot([crossx-200,crossx+200],[crossy+50,crossy+50],"k--",alpha=0.5)
	# plt.plot([crossx-50,crossx-50],[crossy-200,crossy+200],"k--",alpha=0.5)
	# plt.plot([crossx-200,crossx+200],[crossy-50,crossy-50],"k--",alpha=0.5)
	# plt.plot([crossx-100,crossx-100],[crossy-200,crossy+200],"k--",alpha=0.5)
	# plt.plot([crossx-200,crossx+200],[crossy-100,crossy-100],"k--",alpha=0.5)
	## PLOT PARTICLES
	plt.scatter(x, y, c='black', s=size, vmin=-0.25, vmax=2.5, edgecolors="none", cmap = cm.get_cmap("Spectral"))
	## DRAW PBC BOX -> RED
	plt.plot([-45,45],[-45,-45],"r-", linewidth = 2.5) ## xmin to xmax at ymin
	plt.plot([-45,45],[45,45],"r-", linewidth = 2.5) ## xmin to xmax at ymax
	plt.plot([-45,-45],[-45,45],"r-", linewidth = 2.5) ## ymin to ymax at xmin
	plt.plot([45,45],[-45,45],"r-", linewidth = 2.5) ## ymin to ymax at xmax
	## INTERIOR BOX GRIC -> BLACK
	#plt.plot([-40,40],[0,0],"k--", linewidth = 0.5) ## xmin to xmax at y = 0 Horizontal line
	#plt.plot([0,0],[-40,40],"k--", linewidth = 0.5) ## ymin to ymax at x = 0 Vertical line
	## AXES LIMITS -> STATIC
	plt.xlim([-90,90])
	plt.ylim([-90,90])
	# Move left y-axis and bottim x-axis to centre, passing through (0,0)
	# xax = np.array(range(-40,40))
	# yax = np.array(range(-40,40))
	# plt.plot(xax,yax)
	ax.spines["left"].set_position(("axes", 1))
	ax.spines["bottom"].set_position(("axes",1))
	ax.spines["right"].set_color("none")
	ax.spines["top"].set_color("none")
	#plt.xlabel("X")
	## THESE will cause axis to move, makes axis dynamic not static
	# plt.xlim([COMx-100,COMx+100])
	# plt.ylim([COMy-100,COMy+100])

	plt.axis("on")
	plt.subplots_adjust(left=-0.205,bottom=-0.205,right=1.05,top=1.05)


anim = ani.FuncAnimation(fig,readData,frames=len(files), blit=False)
anim.save("animation.gif", fps=10)#, extra_args=["-vcodec", "libx264"])

END

cd ..                       #Come out of video directory

g++ main.cpp                #Execute the code
./a.out

cd ..                       #Come out of sample directory
done

cd ..                       #Come out of simulation folder
fi

done < $INPUT
IFS=$OLDIFS


