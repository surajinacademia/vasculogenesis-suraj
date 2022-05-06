/******************************************************************************************************The code was originally written by Katie Copehgen during her PhD at UC Merced Physics dept. She is now at Princeton University as Postdoc. Then it was passed down to Farnaz Golnaraghi during her PhD at UC Merced, Physics. Now she works at Univeristy of Chicago. The major contribution on editing the code is by Mikahl Banwarth-Kuhn who is a Post-doc at Suzzane Sindi's Lab.Now it's Suraj Sahu's turn to work on this code for his PhD. He is just running this code for different simulation parameters.    This simulates a group of agents that are connected to neighboring agents by springs     which can be fairly long range. Agents can intercalate other agents by interrupting     the space inbetween the agents connected by the spring. The agents have angular noise. \******************************************************************************************************/#include <iostream>#include <stdlib.h>#include <math.h>#include <fstream>#include <time.h>#include <vector>#include <random>using namespace std;#define N 225          //Number of agents. #define V 0.01        //Agent velocity for self propulsion. Low velocity choosen#define NOISE 0.01    //Noise (variance of gaussian distribution for random angular noise). [0-2*pi]#define RS 0.0        //Variance in agent sizes, for gaussian distribution. [0-RA (defined below)]#define RA 1        //Average size of the agents, sets lengthscale. [1] Diameter or radius???#define K 1       //Spring interaction strength. [0-0.1]???#define anisotropic_adhesion_on 1          //control whether we see iso or aniso adhesion#define testing_aniso   0          //if this is 1 initial conditions of cell will be a line#define lattice_nonrandom 1#define lattice_start 1 // if 1 then cells will start in lattice formation#define Neighbor_region 2#define anisotropic_adhesion_angle M_PI/4      //size of the adhesion window rangee [0-pi]#define AA 1.0          //Alignment interaction magnitude. [0-1 or 10]?#define GAP 1        //The amount of overlap to break springs. [0-1]#define T 2000        //Number of timesteps the simulation is run for. [100-100000]#define vid_steps 10#define SAMPLES 1    //Number of samples run. Must be 1 to output files for making videos. [1-100+]//INITIAL CONDITIONS//Uniform circle#define L 25.0     //Size of the circle which the agents are initiallized within. [sqrt(N) ish]//Lattice Configuration#define lattice_step_size 4.0#define lattice_width (sqrtf(N)-1)*lattice_step_size#define lattice_height (sqrtf(N)-1)*lattice_step_size#define lattice_neg_x -lattice_width/2  //left lower bound of lattice#define lattice_pos_x lattice_width/2  //right upper bound of lattice#define lattice_neg_y -lattice_width/2 //bottom lower bound of lattice#define lattice_pos_y lattice_width/2 //top upper bound of lattice#define MAX_X (sqrtf(N)*lattice_step_size)/2#define MIN_X -MAX_X#define MAX_Y MAX_X#define MIN_Y -MAX_Y#define _USE_MATH_DEFINES #include "particle.cpp"		//Particle class.#include "clustering.cpp"	//Groups agents into clusters.#include "opramcl.cpp"		//Calculate the properties of the agnet clusters.int main(){	//Seed the random number generator.	long int letseed=time(NULL);	srand (letseed);	float t1=clock();	//Save starting time to calculate the program runtime with.	FILE *f00;	//For making videos.	FILE *f01;	//Store parameters.	f01 = fopen("params.txt", "w");	FILE *f02;	//Cluster data at simulation finish.	f02=fopen("data/clusters.dat","w");	FILE *f03;	//Average order vs. t.	f03=fopen("data/ordervst.dat","w");	FILE *f04;	//Agent data at simulation finish.	f04=fopen("data/agents.dat","w");	vector<particle> agents(N, particle());		//Vector for storing all the agents.	for (int s=0; s<SAMPLES; s++){			//Loop through the samples.		//Initialize agents, random position and velocity direction.        //float theta_test = 0;        //float new_x = cosf(theta_test);        //float new_y = sinf(theta_test);        //int newParticleID = 0;        //for (int i=0;i<N;i++){            //cout << newParticleID << endl;            //if(testing_aniso){//if this is turned on then particles will start in a straight line               // agents[i].init(new_x,new_y,newParticleID);            //}else{//otherwise particles will be randomly seeded                //agents[i].init(newParticleID);            //}            //newParticleID = newParticleID + 1;            //theta_test = theta_test + 2*M_PI/10;		//}        //Initialize agents, random position and velocity direction.        if(lattice_start == 1){            float new_x = lattice_neg_x;            float new_y = lattice_neg_y;            int newParticleID = 0;            for (int i=0;i<N;i++){                //cout << newParticleID << endl;                agents[i].init(new_x,new_y,newParticleID);                newParticleID = newParticleID + 1;                if(new_x < lattice_pos_x){                    new_x = new_x + lattice_step_size;                }else{                    new_y = new_y + lattice_step_size;                    new_x = lattice_neg_x;                }		    }        }else if (lattice_nonrandom == 1){            float new_x = lattice_neg_x;            float new_y = lattice_neg_y;            int newParticleID = 0;            std::default_random_engine e1(32767);            std::uniform_real_distribution<float> u(0,1);            float new_ang;            for (int i=0;i<N;i++){                new_ang =  2*M_PI*u(e1);                //cout << newParticleID << endl;                cout << new_ang << endl;                agents[i].init(new_x,new_y,new_ang,newParticleID);                newParticleID = newParticleID + 1;                if(new_x < lattice_pos_x){                    new_x = new_x + lattice_step_size;                }else{                    new_y = new_y + lattice_step_size;                    new_x = lattice_neg_x;                }		    }        }		for (int t=0; t<T; t++){		//Loop through time.			//Once enough time has passed to stablize the cluster, spread the sizes, introducing the polydispersity gradually over 100 timesteps.			if (t>(int)3*(sqrtf((float)N)/((float)(2.0*K))*RA)&&t<(int)3*(sqrtf(N)/((float)(2.0*K))*RA)+100){				for (int i=0;i<N;i++)					agents[i].spreadSizes();					//printf("I am here\n");			}			//Create neighbor list for all agents.			for (int i=0; i<N; i++){				vector<int> temp;	//Temporarily stores the neighbor list.				temp.reserve(N);	//Pre reserve space to speed up program.				for (int j=0; j<N; j++){	//Loop through all other particles to check if they should be added to neighborlist.					float dx=agents[j].getX()-agents[i].getX();					float dy=agents[j].getY()-agents[i].getY();					float d=sqrtf(dx*dx+dy*dy);					if (d<=agents[i].getSize()+agents[j].getSize()+Neighbor_region*RA){						temp.push_back(j);  //If i and j are close enough together add j to the neighborlist,right now they have to be 2*RA apart					}				}				agents[i].setNeigh(temp); //Set the neighborlist to temp.				temp.clear();				agents[i].setOldAng();  //Also save their angles for use in interactions.			}			//Loop throgh every particle to calculate interactions.            int kaligntemp = 0; //kaligntemp is for storing whether agents will be aligned or not.			for (int i=0; i<N; i++){                vector<int> aniso_temp;                aniso_temp.reserve(N);                int my_id = agents[i].getParticleID();				//cout << my_id << endl;                float Ax=0, Ay=0, A, Sx=0, Sy=0, newx, newy, dx, dy, d, newang;	//For storing all the interactions, to calculate new directions for agent i.				for (int j=0; j<agents[i].getNeigh().size(); j++){		//Loop through j: all of agent i's neighbors.                    int neighbor_id = agents[agents[i].getNeigh()[j]].getParticleID();					//Calculate separation vector between agent i and i's j'th neighbor to use for spring.					float dx=agents[agents[i].getNeigh()[j]].getX()-agents[i].getX();					float dy=agents[agents[i].getNeigh()[j]].getY()-agents[i].getY();					float d=sqrtf(dx*dx+dy*dy);					//Make sure it doesnt blow up when they are too close.					if (d<0.001){						d=0.001;                    }					//Spring interaction.					int kspringtemp=0;	//kspringtemp is for storing whether the spring is interrupted by other cells or not.                    float dx2, dy2, d2;	//For storning separation vector to other neighbors than j.					if (my_id!=neighbor_id){								float avgr=(agents[i].getSize()+agents[agents[i].getNeigh()[j]].getSize());                        //this loop is for anisotropic adhesion                        if(anisotropic_adhesion_on==1){                            float orientation_curr_particle = 0; //Oritentation of current particle                            float orientation_neighbor = 0;      //Orientation of neighbour                            float cos_ang_diff = 0;               //Difference of orientation angle between current particle and neighbour. Always take the positive value                            float position_x_curr_particle = 0;  //x position of current particle                            float position_y_curr_particle = 0;  //y position of current particle                            float position_x_neigh_particle = 0; //x position of neighbour particle                            float position_y_neigh_particle = 0;  //y position of neighbour particle                            float center_to_center_orientation = 0;       //center to center vector                             float cos_window = 0;                            float cos_window_half = 0;                            float cos_center_curr = 0;                            float cos_center_neigh = 0;                                                        position_x_curr_particle = agents[i].getX();                            position_y_curr_particle = agents[i].getY();                            position_x_neigh_particle = agents[agents[i].getNeigh()[j]].getX();                            position_y_neigh_particle = agents[agents[i].getNeigh()[j]].getY();                            orientation_curr_particle = agents[i].getAng();                            orientation_neighbor = agents[agents[i].getNeigh()[j]].getAng();                            cos_ang_diff = cosf(orientation_neighbor - orientation_curr_particle); //cos of i and j                                                        if (position_x_neigh_particle >= position_x_curr_particle){                                center_to_center_orientation = atan2((position_y_neigh_particle-position_y_curr_particle), (position_x_neigh_particle-position_x_curr_particle));                            } else if (position_x_neigh_particle<position_x_curr_particle){                                center_to_center_orientation = atan2((position_y_neigh_particle-position_y_curr_particle),(position_x_neigh_particle-position_x_curr_particle)) + M_PI;                                                            }                                                                                    cos_center_curr = cosf(center_to_center_orientation - orientation_curr_particle); //Angle between Rij and i                            cos_center_neigh = cosf(center_to_center_orientation - orientation_neighbor); //Angle between Rij and j                            cos_window = cosf(anisotropic_adhesion_angle); //cos(alpha)                            cos_window_half = cos_window/2.0;   //cos(alpha/2)                                                        if(cos_ang_diff >= cos_window){ //Parallel condition for i and j                                //if(case a)                                if((cos_center_curr >= cos_window_half)&&(cos_center_neigh>=cos_window_half)){                                    kspringtemp = 1;                                    kaligntemp = 1;                                    Ax+=cosf(agents[agents[i].getNeigh()[j]].getAng());                                    Ay+=sinf(agents[agents[i].getNeigh()[j]].getAng());                                    aniso_temp.push_back(neighbor_id);                                //else if(case d)                                }else if((cos_center_curr<=-1.0*cos_window_half)&&(cos_center_neigh<=-1.0*cos_window_half)){                                     kspringtemp = 1;                                    kaligntemp = 1;                                    Ax+=cosf(agents[agents[i].getNeigh()[j]].getAng());                                    Ay+=sinf(agents[agents[i].getNeigh()[j]].getAng());                                      aniso_temp.push_back(neighbor_id);                                //no attach                                }else{                                    kspringtemp = 0;                                    kaligntemp = 0;                                }                            }else if (cos_ang_diff <= -1.0*cos_window){ //Antiparallel condition                                //if(case b)                                if((cos_center_curr>=cos_window_half)&&(cos_window>=-0.1*cos_window_half)){                                    kspringtemp = 1;                                    kaligntemp = 1;                                    Ax+=cosf(agents[agents[i].getNeigh()[j]].getAng());                                    Ay+=sinf(agents[agents[i].getNeigh()[j]].getAng());                                    aniso_temp.push_back(neighbor_id);                                                                //else if(case c)                                }else if((cos_center_curr<=-0.1*cos_window_half)&&(cos_center_neigh>=cos_window_half)){                                    kspringtemp = 1;                                    kaligntemp = 1;                                    Ax+=cosf(agents[agents[i].getNeigh()[j]].getAng());                                    Ay+=sinf(agents[agents[i].getNeigh()[j]].getAng());                                    aniso_temp.push_back(neighbor_id);                                                                      //no attach                                }else{                                    kspringtemp = 0;                                    kaligntemp = 0;                                }                                                    }else{                                kspringtemp = 0;                                kaligntemp = 0;                            }                        }//end anisotropic adhesion loop						//this loop is for GAP interaction which we probably won't consider                        for (int k=0; k<agents[i].getNeigh().size(); k++){	//Loop through all other neighbors of i that aren't j, they shallt be called k. 							int second_neighbor_id = agents[agents[i].getNeigh()[k]].getParticleID();                            if (neighbor_id!=second_neighbor_id && my_id!=second_neighbor_id){				//Only do the thing if they are 3 different agents.								//Separation vector between agent i and it's k'th neighbor.								dx2=agents[agents[i].getNeigh()[k]].getX()-agents[i].getX();									dy2=agents[agents[i].getNeigh()[k]].getY()-agents[i].getY();								d2=sqrtf(dx2*dx2+dy2*dy2);								//b basically records the size of the space between agent i and the j'th neighbor, at the relative distance to the k'th neighbor.								float b=agents[i].getSize()+(agents[agents[i].getNeigh()[j]].getSize()-agents[i].getSize())/(d*d)*(dx2*dx+dy2*dy);								//Then check if the k'th neighbor is within the GAP parameter fraction of the space between agent i and the j'th neighbor. If it is then interrupt the spring by setting kspringtemp to 0.								if (sqrtf(d2*d2-powf(dx*dx2/d+dy*dy2/d,2))<agents[agents[i].getNeigh()[k]].getSize()+GAP*b&&(dx*dx2+dy*dy2)>0){									kspringtemp=0;								}							}						}//end gap interaction loop						                        //Alignment interaction.						// Ax += 0;						// Ay += 0;						if (d<avgr){	//Always have the repulsive part of the spring interaction.							kspringtemp=1;                        }                        //Spring interaction.                        Sx+=-kspringtemp*(avgr-d)*(dx/d);                        Sy+=-kspringtemp*(avgr-d)*(dy/d);					}				}                float total_adh_neighs = agents[i].getAnisoNeighs().size();                if(total_adh_neighs > 0){                    Ax = kaligntemp*Ax/total_adh_neighs;                    Ay = kaligntemp*Ay/total_adh_neighs;                    //cout << Ax << Ay << endl;                }                agents[i].setAnisoNeighs(aniso_temp);                aniso_temp.clear();				//Parameters u and v just for calculating the gaussion random variable for noise.				float u=rand()/(float)RAND_MAX;				float v=rand()/(float)RAND_MAX;				//New direction of propulsion for agent i is a gaussian random variable with variance NOISE, plus the angle that comes from the alignment interaction, and it's own travel direction.				agents[i].setAng(fmodf((NOISE*sqrtf(-2*log(u))*cosf(2*M_PI*v)+atan2(Ay+sinf(agents[i].getAng()),Ax+cosf(agents[i].getAng()))),2*M_PI));				//Set the velocity of the agent to the propulsion with magnitude V in the propulsion direction stored in Ang, and the spring interaction with magnitude K.				agents[i].setVx(V*cosf(agents[i].getAng())+K*Sx);				agents[i].setVy(V*sinf(agents[i].getAng())+K*Sy);			}			//Move agents according to their new directions.			for (int i=0; i<N; i++){				agents[i].updatePos();                agents[i].boundaryCheck();            }			int vidtime=vid_steps;		//How many timesteps to skip between each frame for the video.			if (SAMPLES==1 && t%vidtime==0){				char buff[32];		//buff saves the name of the file to store this frame data in.				snprintf(buff, sizeof(char)*32, "video/frames/fr%06d", t);				f00=fopen(buff, "w");	//Files for saving video frame data to make videos (vid.py).				//Loop through each agent to put data into video frame file.				for (int i=0;i<N;i++){                    int ID = agents[i].getParticleID();					//Find the relative velocity of each agent with all of it's neighbors (will set the color of agents in the video).					float relv=0;					for (int j=0;j<agents[i].getAnisoNeighs().size();j++){						//Find how far agent i has moved since the last video frame.						float mxi=agents[i].getX()-agents[i].getOldX();						float myi=agents[i].getY()-agents[i].getOldY();						//Find how far the j'th neighbor has moved since last video frame.						float mxj=agents[agents[i].getAnisoNeighs()[j]].getX()-agents[agents[i].getAnisoNeighs()[j]].getOldX();						float myj=agents[agents[i].getAnisoNeighs()[j]].getY()-agents[agents[i].getAnisoNeighs()[j]].getOldY();						relv+=sqrtf(powf((mxj-mxi),2)+powf((myj-myi),2));	//Relative movement of agent i and all neighbors.					}					relv=relv/(agents[i].getAnisoNeighs().size()*V*vidtime);	//Normalize relative velocity..					//Write agent data into file for video frame structure 4 columns: relative velocity	X position	Y position	agent size.					fprintf(f00, "%d %f	%f	%f	%f	%f",ID,relv,agents[i].getX(),agents[i].getY(),agents[i].getAng(), agents[i].getSize());                    for (int j=0; j<agents[i].getAnisoNeighs().size(); j++){                        fprintf(f00,"%d	",agents[i].getAnisoNeighs()[j]);                    }                    fprintf(f00,"\n");				}				//Save the current agent positions to use for calculating relative velocities of the next frame.				for (int i=0; i<N; i++)					agents[i].storePos();				fclose(f00);			}            //if(t%100 == 0){            //    cout << "Tim step is: " << t << endl;           // }            for(int i=0; i<N;i++){                for(int j=0;j<agents[i].getAnisoNeighs().size();j++){                    //cout.flush() << "I am particle " << agents[i].getParticleID() << " and my adhesion neighbor is " << agents[agents[i].getNeigh()[j]].getParticleID() << endl;                }            }		}	//End of time loop.		//Calculate properties of each cluster.		int labels[N], clusternum;		clusternum = clustering(agents, labels);		//Function to assign a cluster label to each agent, and returns the total number of clusters.		for (int i=0; i<N; i++)			agents[i].setLabel(labels[i]);		float Order[N],MedianR[N],AvgR[N],StdR[N];		//Arrays for storing the properties of each cluster, only clusternum elements will be saved into these arrays.		int Size[N];		opramcl(agents, Size, Order, MedianR, AvgR, StdR, clusternum);	//Calculate and store the size, order, median, mean, and standard deviation of the size of agents within each cluster.		for (int i=0; i<clusternum; i++){	//Loop through the clusters and store the data about them.			fprintf(f02,"%d	%f	%f	%f	%f\n", Size[i], Order[i], MedianR[i], AvgR[i], StdR[i]);		}		for (int i=0; i<N; i++){	//Save the individual position, direction, size, cluster, and neighbor lists for each agent. 			fprintf(f04,"%f	%f	%f	%f	%d	%d	",agents[i].getX(),agents[i].getY(),agents[i].getAng(),agents[i].getSize(),agents[i].getLabel(),(int)agents[i].getNeigh().size());			for (int j=0; j<agents[i].getNeigh().size(); j++)				fprintf(f04,"%d	",agents[i].getNeigh()[j]);			fprintf(f04,"\n");		}	}	float t2=clock();//	//Write all of the system parameters and run time to params.txt file.//	fprintf(f01, "N = 	%d \nv = 	%f \nra = 	%f\nrs =	%f\na =	%f\nk =	%f\nnoise = 	%f\nsamples = 	%d\nmaxtime =	%d\nseed = 	%ld \n runtime = %f", N, V, (float)RA, (float)RS, (float)AA, (float)K, (float)NOISE, SAMPLES, T, letseed, (t2-t1)/(float)CLOCKS_PER_SEC);	fclose(f01);	return 0;}