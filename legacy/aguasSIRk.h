
#include <math.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include "distearth.h"
#include "aguasData.h"

//#include <boost/random.hpp>
//#include <boost/random/lognormal_distribution.hpp>

#define MUTMAX 200
#define NPAR_PER_DEME_FIX 150
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MAX_FACT_FOR_DOUBLE 170 //maximo que double suporta para calcular factoriais [0,170]


#define ARTESUNATE       10
#define PIPERAQUINE      11
#define PIPARTE          12
#define NODRUG           13
#define LUMEFANTRINE     14
#define AL               15
#define PRIMAQUINE       16

#define MDA      20
#define MONO     21
#define MST      22
#define NOSTRAT  23
#define ACTSWITCH  24

#define NOITN  30
#define ITN    31
#define ITNeffective 32

#define MOBILE  80
#define STATIC  81
#define TEMP    82

#define NO     70
#define YES    71
#define ACTIVE 72
#define WAIT   73
#define CHECK  74

#define PREGNANT  331
#define NOTPREGNANT 332

#define MALE  130
#define FEMALE  131

#define RESISTRO  41
#define RESISTRA  42
#define RESISTRB  43
#define RESISTRAB 44

#define SUSCEPTIBLE 1
#define INFECTIOUS 2
#define BLOODSTAGE 3
#define LIVERSTAGE 4

#define CLINICAL 5
#define ASYMPTOMATIC 6
#define NOINFECTION 7

#define IMMUNE 55
#define NONIMMUNE 56

#define TIMELIMIT 36500000
#define TIMESTART 1
#define TENYEARS 365*10

using namespace std;

std::ofstream fout_resistance;
std::ofstream fout_incidence;

std::ofstream fout_village;
std::ostringstream strfile_village;

std::ofstream fout_taus;
std::ostringstream strfile_taus;

std::ofstream fout_nei;
std::ostringstream strfile_nei;

std::ofstream fout_age;
std::ostringstream strfile_age;

std::ofstream fout_regions;
std::ostringstream strfile_regions;

std::ofstream fout_incmps;
std::ostringstream strfile_incmps;

std::ofstream fout_clinical;
std::ostringstream strfile_clinical;

std::ofstream fout_regions2;
std::ostringstream strfile_regions2;

std::ofstream fout_loglike;
std::ostringstream strfile_loglike;

int G_N;
int G_tmax;
int G_tmed;
int G_runs;
int G_set;
double G_I0;
int G_yearsToDrugs;
int G_yearsToDrugsI1;
int G_maxa;
float G_move;
float G_migrate;
float G_migrate_mobile;

float G_beta;
float G_beta0;
float G_gamma;
float G_sigma;
float G_delta;
float G_mu;

float G_xa0;
float G_xai;
float G_xb;
float G_xitn;
float G_xab;
float G_xai2;
float G_xprim;
float G_xala;
float G_xall;
float G_xlum;
float G_xc;

float G_itneffect;
float G_radius;
float G_itnprop;
float G_amp;
float G_phase;
float G_alpha;
float G_c;
float G_artprophylaxis;
float G_eta;
float G_covmsat;
float G_phic;
float G_phia;
float G_mellow;

float G_cBroda;
float G_cBrada;
float G_cBrbda;
float G_clroda;
float G_clrada;
float G_clrbda;

float G_cBrodo;
float G_cBrado;
float G_cBrbdo;
float G_clrodo;
float G_clrado;
float G_clrbdo;

float G_cBradlf;
float G_clradlf;
float G_cBrbdlf;
float G_clrbdlf;
float G_cBrodlf;
float G_clrodlf;

float G_cBradal;
float G_clradal;
float G_cBrbdal;
float G_clrbdal;
float G_cBrodal;
float G_clrodal;

float G_clrdpmc;
float G_clrdpma;

float G_cBrodab;
float G_cBradab;
float G_cBrbdab;
float G_clrodab;
float G_clradab;
float G_clrbdab;

float G_cBrodb;
float G_cBradb;
float G_cBrbdb;
float G_clrodb;
float G_clradb;
float G_clrbdb;

//float G_cLp;
//float G_clp;
int mdastot;
float G_initclin;
float G_trigger;
int G_trecACT;
int G_treca;
int G_trecb;
int G_trecbma;

float G_erada;
float G_erbdb;
float G_elvsBdb;

float G_pa;
float G_pb;
float G_po;
float G_pinf;

float G_propRxa;
float G_propRxab;
float G_propRxi1;
float G_propRxai1;
float G_pab;
float G_covab;
float G_kab;
float G_lam;

int G_dura;
int G_mteams;
int G_dayspervill;
float G_timecomp;
float G_fullcourse;
float G_nomp;
float G_propasymtreat;
float G_immunity;
float G_import;
int G_importevents;

float G_seas1;
float G_seas2;
float G_seas3;
float G_seas4;

float G_tau;
float G_tau1;
float G_tauab;
float G_tau2;
float G_sensitivity;
float G_prop_immune;

float G_mobility_static;
float G_static_return;
float G_mobility_mobile;
float G_mobile_return;
float G_temp_mobile; // this is not initialised in karenjan17.cpp

float G_mosqdeath1;
float G_mosqdeath2;
float G_mosqinc;

float G_mda_start1;
float G_mda_end1;
float G_mda_start2;
float G_mda_end2;
float G_mda_start3;
float G_mda_end3;
bool G_switch_mda1;
bool G_switch_mda2;
bool G_switch_mda3;

float G_taumagi1;
float G_duri1;
float G_taumagab;
float G_durab;
float G_startab;

float G_taumagabi1;
float G_durabi1;

float G_taumagai1;
float G_durai1;

float G_bn;
float G_bn_start;
float G_bn_end;
float G_bn_value;
float G_bnmag;
float G_bndur;

float G_age1;
float G_age2;
float G_age3;
float G_age4;
float G_age5;
float G_age6;
float G_age7;
float G_age8;

float G_k;
float G_r;
float G_ke;
float G_mu1;
float G_mu2;
float G_sig;
float G_pclin;
float G_pclinimm;
float G_natality;

int G_dimsurveys;
int G_first;
int G_delay;
float betaf;

int G_xregion;
int G_yregion;
int G_nstartpoints;

float seasons;
float G_poormp;
float G_daysoff;
float G_antroph;
float G_indoorpref;


vector<vill*> G_villages;
vector<Deme*> G_demes;
vector<listv*> G_list;
vector<mosq*> G_mosquitoes;

StateDeme *G_groups;
vector<villdata*> G_villdata;

vector<int> pop;
vector<int> inf;
vector<int> inf2;
vector<int> tot;
vector<int> aux;
vector<int> auxage;
vector<int> aux2;
vector<int> number_pop_area;
vector<int> i_init;
vector<float> lat;
vector<float> lon;
vector<float> betas;
vector<int> mp;
vector<int> mptime;
vector<int> mdatime;
vector<int> elevation;
vector<int> malpost;
vector<float> betaseason;
vector<float> treatab;
vector<float> treatasympt;
vector<float> treatclinical;
vector<int>triggs;
vector<int>surveys;
vector<float>mob;
vector<float>pa;
vector<float>pb;
vector<float>seasonvec;
vector<int>agev;
vector<int>moiv;
vector<int>clinstat;
vector<int>levelimm;
vector<int>cumminf;
vector<float>guesspars;

void infection();
void livertoblood();
void tosymptoms();
void recovery();
void birthdeath();
void printResistance(int timestep);
void printTreat();
void treatment(int timestep);
void inline squarepulse(int timestep);
void printTaus();
void printbednet();
void drugeffect();
void drugloss();
void itnloss();
void printage();
void printinfo(int timestep);
void printinfo2(int timestep);
void ageing();
void setclinicalprob();
void setsusceptibility();
void setinfectiousness();
void setprobdeath();
void seasonality(int timestep);
void openposts(int timestep);
void openposts2(int timestep);
void settreatrate();
void immunityloss();
void immunitycount();
void clinicalresolution();
void startvillages();
void startvillages2(int timestep);
void startvillages3(int timestep);
void updatevillages(int timestep);
//int updatevillages(int timestep, int mdastot);
void updatevillages2(int timestep);
void updatevillages3(int timestep);
void updatevilldata(int timestep);
void triggerMDA(int timestep);
void treattrigger(int timestep);
void startvillagestrigger();
void neighbouring();
void surveytimes();
int  countmda();
void countposts();
void infectmosquitoes();
void survivemosqs();
void mosqincubation();
void printmosqs();
void mobility();
void village_here_now();
void migration(int timestep);
void temp_migration(int timestep);
void seasonal_network();
void migration_network();
void static_network();
void printagestruct(int timestep);
void printregion(int timestep);
void superinfectiontoblood();
void superinfectiontoinfectious();
void superrecovery();
void superbirth();
void superdrugeffect();
void naturalhistoryinfection(int timestep);
void villages(int timestep, float seasons);
void population_movement(int timestep);
void mosquitodynamics();
void updateLL(std::vector<float> vec_guesses, double loglike, long thisrun);
double calculateincidence(int G_xregion,int G_tmax);
double PoissonLikelihood(std::vector<double> thedata,std::vector<double> themodel);
vector<float> guessparams( std::vector<float> bestparams);
float gasdev(long g_seed);
double square(float);
double NegativeBinomialLikelihood(std::vector<double> thedata, std::vector<double> themodel,double ascale);
double LogGamma(double x);
void importation(int timestep);
void printimmunity(int timestep);

void inline fitnessinit(){
     int deme, ia, ib, io, iall;
     float randomEvent;
     float randomEvent2;
     float randomEvent3;
     float randomEvent4;
     int sumi = 0;
     mosq* mosquito;
     auxage=vector<int>(G_N,0);

        int k=0;
        float agep[G_maxa];
        int np=0;
        int begp=0;

//        for (k=0; k<G_maxa; k++){
//
//            agep[k]= G_age1 + G_age2*k + G_age3*pow(k,2) + G_age4*pow(k,3) + G_age5*pow(k,4) + G_age6*pow(k,5) + G_age7*pow(k,6) + G_age8*pow(k,7);
//            agep[k] = agep[k]/0.94;
//            begp = np + 0;
//            np = np + int(agep[k]*G_N);
//
//            for (deme=begp; deme<np; deme ++){
//                auxage[deme]=k;
//                //cout<<k<<endl;
//                if (k>75){
//                    int dage = int(ran2()*35)+40;
//                    auxage[deme]=dage;
//                }
//            }
//        }
//        random_shuffle(auxage.begin(),auxage.end());
//        cout<<"agedistdone"<<endl;



//        for (deme=0; deme<G_N; deme ++) {
//            //cout<<"age"<<auxage[deme]<<endl;
//            int rage = int(ran2()*G_N);
////            G_demes[deme]->age = auxage[rage];
//            if (G_demes[deme]->age == 0){
//                G_demes[deme]->cummulative_exposures=0;
//                G_demes[deme]->immunity_level=0;
//                G_demes[deme]->immunity=NONIMMUNE;
//                G_demes[deme]->immunity_days=0;
//            }
//            else{
//                G_demes[deme]->cummulative_exposures=int(exp(0.05*float(G_demes[deme]->age)));
//                G_demes[deme]->immunity_level=int(G_demes[deme]->cummulative_exposures/10);
//                if (G_demes[deme]->immunity_level>0){
//                    G_demes[deme]->immunity=IMMUNE;
//                    G_demes[deme]->immunity_days=40;
//                }
//                else{G_demes[deme]->immunity=NONIMMUNE;}
//            }
//        }
//        cout<<"assign deme age"<<endl;


        for (int i=0;i<G_xregion;i++){

           int people_region = pop[i];
           float mobile_pop = mob[i];
           //cout<<"i: "<<i<<"popi: "<<pop[i]<<endl;
           //cout<<"HOLA1"<<endl;

           ia = int (pa[i] * i_init[i]);
           ib = int (pb[i] * i_init[i]);
           io = i_init[i]-(ia+ib);
//           cout<<ia<<" "<<ib<<" "<<io<<" INIT: "<<i_init[i]<<endl;

           iall = io+ia+ib;

           char Region;
           char buffer [50];
           int n;

           mosquito=G_mosquitoes[i];
           mosquito->infected=0;
           mosquito->infectious=0;

           if (i==0){

              aux=vector<int>(people_region,0);
              for (int u=0; u<people_region;u++){
                aux[u]=u;
                //cout<<"aux is: "<<aux[u]<<endl;
              }
              random_shuffle(aux.begin(),aux.end());


              for(deme = 0; deme<io; deme++){
              int ind=aux[deme];
              //cout<<"ind"<<ind<<endl;
              n=1000+i;
              G_demes[ind]->latitude = lat[i];
              G_demes[ind]->longitude = lon[i];
              G_demes[ind]->mda_days = 0;
              G_demes[ind]->mda_times = 0;
              G_demes[ind]->home = n;
              G_demes[ind]->oldhome = n;
              G_demes[ind]->now = n;
              G_demes[ind]->region = n;
              G_demes[ind]->resistance = RESISTRO;
              G_demes[ind]->dominant = RESISTRO;
              G_demes[ind]->transmission = INFECTIOUS;
              G_demes[ind]->drug = NODRUG;
              G_demes[ind]->prim = NO;
              G_demes[ind]->moi = 1;
              G_demes[ind]->infections.push_back(RESISTRO);
              G_demes[ind]->strategy = NOSTRAT;
//              G_groups->SUM_IRO = G_groups->SUM_IRO + 1;
              G_villdata[i]->SUM_IRO = G_villdata[i]->SUM_IRO + 1;

//
                    float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[aux[deme]]->mobile = TEMP;}
                            else{G_demes[aux[deme]]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}

                    float randomEvent3 = ran2();
                    if(randomEvent3 <= 0.02){
                       mosquito->infectious += 1;
                       mosquito->resistancei.push_back(RESISTRO);
                    }
                    if(randomEvent3 <= 0.04){
                       mosquito->infected += 1;
                       mosquito->resistance.push_back(RESISTRO);
                    }
              }

              for(deme = io; deme<io+ia; deme++){
              int ind=aux[deme];
              G_demes[ind]->latitude = lat[i];
              G_demes[ind]->longitude = lon[i];
              n=1000+i;
              G_demes[ind]->region = n;
              G_demes[ind]->home = n;
              G_demes[ind]->oldhome = n;
              G_demes[ind]->now = n;
              G_demes[ind]->resistance = RESISTRA;
              G_demes[ind]->dominant = RESISTRA;
              G_demes[ind]->transmission = INFECTIOUS;
              G_demes[ind]->moi = 1;
              G_demes[ind]->infections.push_back(RESISTRA);
              G_demes[ind]->mda_days = 0;
              G_demes[ind]->mda_times = 0;
              G_demes[ind]->drug = NODRUG;
              G_demes[ind]->prim = NO;
              G_demes[ind]->strategy = NOSTRAT;
//              G_groups->SUM_IRA = G_groups->SUM_IRA + 1;
              G_villdata[i]->SUM_IRA = G_villdata[i]->SUM_IRA + 1;
//
              float randomEvent = ran2();
                    if(randomEvent <= G_initclin){
                             G_demes[ind]->clinical = CLINICAL;
//                             G_groups->SUM_CLIN = G_groups->SUM_CLIN+ 1;
                    }
                    else { G_demes[ind]->clinical = ASYMPTOMATIC;
//                    G_groups->SUM_ASYM = G_groups->SUM_ASYM+ 1;
                    }

              float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[ind]->mobile = TEMP;}
                            else{G_demes[ind]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}

              float randomEvent3 = ran2();
                    if(randomEvent3 <= 0.02){
                       mosquito->infectious += 1;
                       mosquito->resistancei.push_back(RESISTRA);
                    }
                    if(randomEvent3 <= 0.04){
                       mosquito->infected += 1;
                       mosquito->resistance.push_back(RESISTRA);
                    }
              }


              for(deme = io+ia; deme<io+ia+ib; deme++){
              int ind = aux[deme];
              G_demes[ind]->latitude = lat[i];
              G_demes[ind]->longitude = lon[i];
              n=1000+i;
              G_demes[ind]->region = n;
              G_demes[ind]->home = n;
              G_demes[ind]->oldhome = n;
              G_demes[ind]->now = n;
              G_demes[ind]->resistance = RESISTRB;
              G_demes[ind]->dominant = RESISTRB;
              G_demes[ind]->transmission = INFECTIOUS;
              G_demes[ind]->moi = 1;
              G_demes[ind]->infections.push_back(RESISTRB);
              G_demes[ind]->drug = NODRUG;
              G_demes[ind]->prim = NO;
              G_demes[ind]->strategy = NOSTRAT;
//              G_groups->SUM_IRB = G_groups->SUM_IRB + 1;
              G_villdata[i]->SUM_IRB = G_villdata[i]->SUM_IRB + 1;
              G_demes[ind]->mda_days = 0;
              G_demes[ind]->mda_times = 0;

              float randomEvent = ran2();
                    if(randomEvent <= G_initclin){
                             G_demes[ind]->clinical = CLINICAL;
//                             G_groups->SUM_CLIN = G_groups->SUM_CLIN + 1;
                             G_demes[ind]->clintoday=1;
                             }
                    else { G_demes[ind]->clinical = ASYMPTOMATIC;
//                    G_groups->SUM_ASYM = G_groups->SUM_ASYM + 1;
                    G_demes[ind]->clintoday=0;}

              float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[ind]->mobile = TEMP;}
                            else{G_demes[ind]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}

              float randomEvent3 = ran2();
                    if(randomEvent3 <= 0.02){
                       mosquito->infectious += 1;
                       mosquito->resistancei.push_back(RESISTRB);
                    }
                    if(randomEvent3 <= 0.04){
                       mosquito->infected += 1;
                       mosquito->resistance.push_back(RESISTRB);
                    }
              }

              for(deme = io+ia+ib; deme<people_region; deme++){
              int ind=aux[deme];
              G_demes[ind]->latitude = lat[i];
              G_demes[ind]->longitude = lon[i];
              n=1000+i;
              G_demes[ind]->region = n;
              G_demes[ind]->home = n;
              G_demes[ind]->oldhome = n;
              G_demes[ind]->now = n;
              G_demes[ind]->resistance = RESISTRO;
              G_demes[ind]->dominant = RESISTRO;
              G_demes[ind]->transmission = SUSCEPTIBLE;
              G_demes[ind]->moi = 0;
              G_demes[ind]->infections.clear();
              G_demes[ind]->infectious_qeue.clear();
              G_demes[ind]->drug = NODRUG;
              G_demes[ind]->prim = NO;
              G_demes[ind]->mda_days = 0;
              G_demes[ind]->mda_times = 0;
              G_demes[ind]->strategy = NOSTRAT;
//              G_groups->SUM_S = G_groups->SUM_S+ 1;
              G_villdata[i]->SUM_S = G_villdata[i]->SUM_S+ 1;
              G_demes[ind]->clinical = NOINFECTION;
//              G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;

              float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[ind]->mobile = TEMP;}
                            else{G_demes[ind]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}
              }
           }
     //I
           else{
                //cout<<"entrei aki"<<i<<endl;
//                int k=0;
//                float agep[G_maxa];
                int np=0;
                //int begp=0;
                int people_prev = pop[i-1];
                //cout<<"which village: "<<i<<" count till here: "<<people_prev+people_region<<endl;
                //cout<<"pop"<<pop[i]<<"/tpeople region:"<<people_prev<<endl;
//
//                    for (k=0; k<G_maxa; k++){
//
//                        agep[k]= G_age1 + G_age2*k + G_age3*pow(k,2) + G_age4*pow(k,3) + G_age5*pow(k,4) + G_age6*pow(k,5) + G_age7*pow(k,6) + G_age8*pow(k,7);
//                        agep[k] = agep[k]/0.88;
//
//                        begp = np + people_prev;
//                        np = np + int(agep[k]*(people_region-people_prev));
//                        //cout<<"entrei aki tb"<< "beggining"<<begp<<"popsize"<<np<<endl;
//
//                        for (deme=begp; deme<people_region; deme ++) {
//                        G_demes[deme]->age = k;
//                        if (G_demes[deme]->age ==0){
//                            G_demes[deme]->cummulative_exposures=0;
//                            G_demes[deme]->immunity_level=0;
//                            G_demes[deme]->immunity=NONIMMUNE;
//                            G_demes[deme]->immunity_days=0;
//                        }
//                        else{
//                            G_demes[deme]->cummulative_exposures=int(4-4*(exp(-0.24*k)));
//                            G_demes[deme]->immunity_level=int(ran2()*4)+2;
//                            if (G_demes[deme]->immunity_level>0){
//                                G_demes[deme]->immunity=IMMUNE;
//                                G_demes[deme]->immunity_days=40;
//                            }
//                            else{G_demes[deme]->immunity=NONIMMUNE;}
//                        }
//                        //cout<<"e aki!!!"<<endl;
//                        //cout<<"deme age:"<< k <<endl;
//                        }
//                    }

                aux2=vector<int>(people_region-people_prev,0);
                for (int u=people_prev; u<people_region;u++){
                    aux2[u-people_prev]=u;
                    //cout<<"aux2 is: "<<aux2[u-people_prev]<<endl;
                }
                random_shuffle(aux2.begin(),aux2.end());


                for(deme = people_prev; deme<people_prev+io; deme++){
                    int ind=aux2[deme-people_prev];
                    //cout<<"ind is :"<<ind<<endl;

                    G_demes[ind]->latitude = lat[i];
                    G_demes[ind]->longitude = lon[i];
                    n=1000+i;
                    G_demes[ind]->mda_days = 0;
                    G_demes[ind]->mda_times = 0;
                    G_demes[ind]->region = n;
                    G_demes[ind]->home = n;
                    G_demes[ind]->oldhome = n;
                    G_demes[ind]->now = n;
                    G_demes[ind]->resistance = RESISTRO;
                    G_demes[ind]->dominant = RESISTRO;
                    G_demes[ind]->transmission = INFECTIOUS;
                    G_demes[ind]->moi = 1;
                    G_demes[ind]->infections.push_back(RESISTRO);
                    G_demes[ind]->drug = NODRUG;
                    G_demes[ind]->prim = NO;
                    G_demes[ind]->strategy = NOSTRAT;
//                    G_groups->SUM_IRO = G_groups->SUM_IRO + 1;
                    G_villdata[i]->SUM_IRO = G_villdata[i]->SUM_IRO + 1;

                    float randomEvent = ran2();
                    if(randomEvent <= G_initclin){
                             G_demes[ind]->clinical = CLINICAL;
//                             G_groups->SUM_CLIN = G_groups->SUM_CLIN + 1;
                             }
                    else { G_demes[ind]->clinical = ASYMPTOMATIC;
//                    G_groups->SUM_ASYM = G_groups->SUM_ASYM + 1;
                    }

                    float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[ind]->mobile = TEMP;}
                            else{G_demes[ind]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}

                    float randomEvent3 = ran2();
                    if(randomEvent3 <= 0.02){
                       mosquito->infectious += 1;
                       mosquito->resistancei.push_back(RESISTRO);
                    }
                    if(randomEvent3 <= 0.04){
                       mosquito->infected += 1;
                       mosquito->resistance.push_back(RESISTRO);
                    }
                }

                for(deme = people_prev+io; deme<people_prev+io+ia; deme++){
                int ind = aux2[deme-people_prev];
                G_demes[ind]->latitude = lat[i];
                G_demes[ind]->longitude = lon[i];
                n=1000+i;
                G_demes[ind]->region = n;
                G_demes[ind]->home = n;
                G_demes[ind]->oldhome = n;
                G_demes[ind]->now = n;
                G_demes[ind]->resistance = RESISTRA;
                G_demes[ind]->dominant = RESISTRA;
                G_demes[ind]->transmission = INFECTIOUS;
                G_demes[ind]->moi = 1;
                G_demes[ind]->infections.push_back(RESISTRA);
                G_demes[ind]->drug = NODRUG;
                G_demes[ind]->prim = NO;
                G_demes[ind]->mda_days = 0;
                G_demes[ind]->mda_times = 0;
                G_demes[ind]->strategy = NOSTRAT;
//                G_groups->SUM_IRA = G_groups->SUM_IRA + 1;
                G_villdata[i]->SUM_IRA = G_villdata[i]->SUM_IRA + 1;

                float randomEvent = ran2();
                    if(randomEvent <= G_initclin){
                             G_demes[ind]->clinical = CLINICAL;
//                             G_groups->SUM_CLIN = G_groups->SUM_CLIN + 1;
                             }
                    else { G_demes[ind]->clinical = ASYMPTOMATIC;
//                    G_groups->SUM_ASYM = G_groups->SUM_ASYM + 1;
                    }

                float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[ind]->mobile = TEMP;}
                            else{G_demes[ind]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}

                float randomEvent3 = ran2();
                    if(randomEvent3 <= 0.02){
                       mosquito->infectious += 1;
                       mosquito->resistancei.push_back(RESISTRA);
                    }
                    if(randomEvent3 <= 0.04){
                       mosquito->infected += 1;
                       mosquito->resistance.push_back(RESISTRA);
                    }
                }

                for(deme = people_prev+io+ia; deme<people_prev+io+ia+ib; deme++){
                //cout<<"HOLA7777"<<endl;
                int ind = aux2[deme-people_prev];
                G_demes[ind]->latitude = lat[i];
                G_demes[ind]->longitude = lon[i];
                n=1000+i;
                G_demes[ind]->region = n;
                G_demes[ind]->home = n;
                G_demes[ind]->oldhome = n;
                G_demes[ind]->now = n;
                G_demes[ind]->resistance = RESISTRB;
                G_demes[ind]->dominant = RESISTRB;
                G_demes[ind]->transmission = INFECTIOUS;
                G_demes[ind]->moi = 1;
                G_demes[ind]->infections.push_back(RESISTRB);
                G_demes[ind]->drug = NODRUG;
                G_demes[ind]->prim = NO;
                G_demes[ind]->strategy = NOSTRAT;
                G_demes[ind]->mda_days = 0;
                G_demes[ind]->mda_times = 0;
//                G_groups->SUM_IRB = G_groups->SUM_IRB + 1;
                G_villdata[i]->SUM_IRB = G_villdata[i]->SUM_IRB + 1;

                float randomEvent = ran2();
                    if(randomEvent <= G_initclin){
                             G_demes[ind]->clinical = CLINICAL;
//                             G_groups->SUM_CLIN = G_groups->SUM_CLIN + 1;
                             G_demes[ind]->clintoday=1;
                             }
                    else { G_demes[ind]->clinical = ASYMPTOMATIC;
//                    G_groups->SUM_ASYM = G_groups->SUM_ASYM + 1;
                    G_demes[ind]->clintoday=0;
                    }

                float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[ind]->mobile = TEMP;}
                            else{G_demes[ind]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}

                float randomEvent3 = ran2();
                    if(randomEvent3 <= 0.02){
                       mosquito->infectious += 1;
                       mosquito->resistancei.push_back(RESISTRB);
                    }
                    if(randomEvent3 <= 0.04){
                       mosquito->infected += 1;
                       mosquito->resistance.push_back(RESISTRB);
                    }
                }

                for(deme = people_prev+io+ia+ib; deme<people_region; deme++){
                int ind = aux2[deme-people_prev];
                G_demes[ind]->latitude = lat[i];
                G_demes[ind]->longitude = lon[i];
                n=1000+i;
                G_demes[ind]->region = n;
                G_demes[ind]->home = n;
                G_demes[ind]->oldhome = n;
                G_demes[ind]->now = n;
                G_demes[ind]->resistance = RESISTRO;
                G_demes[ind]->dominant = RESISTRO;
                G_demes[ind]->transmission = SUSCEPTIBLE;
                G_demes[ind]->moi = 0;
                G_demes[ind]->infections.clear();
                G_demes[ind]->infectious_qeue.clear();
                G_demes[ind]->drug = NODRUG;
                G_demes[ind]->prim = NO;
                G_demes[ind]->mda_days = 0;
                G_demes[ind]->mda_times = 0;
                G_demes[ind]->strategy = NOSTRAT;
//                G_groups->SUM_S = G_groups->SUM_S + 1;
                G_villdata[i]->SUM_S =  G_villdata[i]->SUM_S + 1;
                G_demes[ind]->clinical = NOINFECTION;
//                G_groups->SUM_NOINF = G_groups->SUM_NOINF+ 1;

                float randomEvent2 = ran2();
                    if(randomEvent2 <= mobile_pop){
                        float randomEvent3 = ran2();
                            if(randomEvent3 <= G_temp_mobile){
                                G_demes[ind]->mobile = TEMP;}
                            else{G_demes[ind]->mobile = MOBILE;}
                    }
                    else { G_demes[ind]->mobile = STATIC;}
                }
           }// close else
           //cout<<"HOLA111"<<endl;

          sumi = sumi+iall;
     } // for each region
     //cout<<"prev: "<<float (float(sumi)/float(G_N))<<endl;

     //initial bednet distribution
     for(deme = 0; deme < G_N; deme++){
        randomEvent = ran2();
        if (randomEvent < G_itnprop){
            G_demes[deme]->itn = ITN;
//            G_groups->SUM_ITN = G_groups->SUM_ITN + 1;
        }
        else{
            G_demes[deme]->itn = NOITN;
//            G_groups->SUM_NOITN = G_groups->SUM_NOITN + 1;
        }
    }
    //cout<<"HOLA2222"<<endl;


    //initial gender distribution
    for(deme = 0; deme < G_N; deme++){
        randomEvent = ran2();
        if (randomEvent <= 0.482){ //in accordance with 1998 census
            G_demes[deme]->gender = MALE;
        }
        else{
            G_demes[deme]->gender = FEMALE;
        }
    }
    //cout<<"finishgender"<<endl;

//    for(deme = 0; deme < G_N; deme++){
//        float randomEvent2 = ran2();
//          if(randomEvent2 <= G_prop_immune){
//                    G_demes[deme]->immunity = NONIMMUNE;
//                    G_groups->SUM_NIMM = G_groups->SUM_NIMM + 1;
//                    }
//          else {
//               G_demes[deme]->immunity = IMMUNE;
//               G_groups->SUM_IMM = G_groups->SUM_IMM + 1;
//               }
//    }

} //close fitnessinit


void ageing(int timestep){
     Deme* deme;

     for (int d=0; d<G_N; d++){
         deme=G_demes[d];
         if (timestep%365==0){
            deme->age = deme->age + 1;
         }
     }
}

void startvillagestrigger(){
     vill* village;

     for (int i=0; i<G_xregion; i++){
         village=G_villages[i];
         village->mda = NO;
         village->mda_days = 0;
         village->mda_times = 0;
         village->counter = 0;
         village->post = 0;
         village->delay = 0;
         village->mda_start = 0;
     }
}


void startvillages(){
     vill* village;
     listv* listt;

    for (int i=0; i<G_xregion; i++){
        int region = i;
        village=G_villages[region];
        village->mda = NO;
        village->mda_days = 0;
        village->mda_times = 0;
        village->counter = 0;
        village->post = 0;
        village->clinical=0;
    }

    for (int j=0; j<G_nstartpoints; j++){
        listt = G_list[j];
        //cout<<"reg: "<<listt->listvill[0]<<endl;
        //int g=getchar();
        int region = listt->listvill[0];
        //cout<<"reg2: "<<region<<endl;
        village=G_villages[region];
        village->mda = YES;
        village->mda_days = 1;
        village->mda_times = 1;
        village->counter = 1;
        //cout<<"tt: "<<G_nstartpoints<<endl;
    }
}

void startvillages2(int timestep){
     vill* village;
     listv* listt;
     Deme* deme;

     if (timestep==G_mda_start2){

         for (int i=0; i<G_xregion; i++){
             int region = i;
             village=G_villages[region];
             village->mda = NO;
             village->mda_days = 0;
             village->mda_times = 0;
             village->counter = 0;
             village->post = 0;
         }

         for (int j=0; j<G_nstartpoints; j++){
    	     listt = G_list[j];
             listt->listorder = 0;
             //cout<<"reg: "<<listt->listvill[0]<<endl;
             //int g=getchar();
             int region = listt->listvill[0];
             //cout<<"reg2: "<<region<<endl;
             village=G_villages[region];
             village->mda = YES;
             village->mda_days = 1;
             village->mda_times = 1;
             village->counter = 1;
          }

         for (int w=0; w<G_N; w++){
             deme=G_demes[w];
             deme->mda_times = 0;
             deme->mda_days = 0;
         }
    }
}

void startvillages3(int timestep){
     vill* village;
     listv* listt;
     Deme* deme;

     if (timestep==G_mda_start3){

         for (int i=0; i<G_xregion; i++){
             int region = i;
             village=G_villages[region];
             village->mda = NO;
             village->mda_days = 0;
             village->mda_times = 0;
             village->counter = 0;
             village->post = 0;
         }

         for (int j=0; j<G_nstartpoints; j++){
    	     listt = G_list[j];
             listt->listorder = 0;
             //cout<<"reg: "<<listt->listvill[0]<<endl;
             //int g=getchar();
             int region = listt->listvill[0];
             //cout<<"reg2: "<<region<<endl;
             village=G_villages[region];
             village->mda = YES;
             village->mda_days = 1;
             village->mda_times = 1;
             village->counter = 1;
          }

         for (int w=0; w<G_N; w++){
             deme=G_demes[w];
             deme->mda_times = 0;
             deme->mda_days = 0;
         }
    }
}

//int updatevillages(int timestep, int mdastot){
//     vill* village;
//
//    for(int i =0; i<G_xregion; i++){
//        village=G_villages[i];
//        int times = mdatime[i];
//
//        if(village->mda==ACTIVE){
//            if (village->mda_days>G_dayspervill){
//                village->mda==NO;
//            }
//        }
//
//        if (times==timestep){
//            //cout<< "MDA NOW"<<endl;
//            village->mda = YES;
//            village->mda_days=1;
//            village->delay=0;
//            mdastot+=1;
//        }
//     }
//
////     cout<<"MDAS: "<<mdastot<<endl;
//     return mdastot;
//}

void updatevillages(int timestep){
     vill* village;
     listv* listt;
     vill* village2;

     if (timestep>G_mda_start1 && timestep<G_mda_start2){

     for(int i =0; i<G_xregion; i++){
             village=G_villages[i];
             if (village->mda==YES){
             //cout<<"which: "<< i <<endl;
             village->mda_days = village->mda_days+1;}
     }


     for (int j=0; j<G_nstartpoints; j++){
         listt = G_list[j];
         int k=listt->listorder;
         int region=listt->listvill[k];
         village2=G_villages[region];

         if (village2->mda_days>G_dayspervill){
            village2->mda=NO;

            if (k<(G_xregion-1)){
                k=k+1;
                //cout<<"regsum:"<<k<<endl;
                int k2=k;
                int region2=listt->listvill[k2];
                //cout<<"list: "<<region2<<"\tk2: "<<k2<< "\tlist num: "<< j <<endl;
                int qq = G_villages[region2]->mda_times;

                //cout<<"before: "<<region2<<"\tmda times: "<<qq<<endl;
                if (qq > 0 ){
                   int tt=0;
                   while(tt==0){
                   if (k<(G_xregion-1)){
                      //cout<<"here: "<< endl;
                      k=k+1;
                      int k2=k;
                      int region2=listt->listvill[k2];

                      qq = G_villages[region2]->mda_times;
                      //cout<<"change: "<< region2 <<"\tmda times: "<<qq<< "\tlist num: "<< j <<endl;
                      if (qq == 0){
                             tt=tt+1;
                             listt->listorder=k2;
                             int region3=listt->listvill[k2];
                             village2=G_villages[region3];
                             village2->mda=YES;
                             village2->mda_times=village->mda_times+1;
                             listt->listorder=k2;
                             village2->mda_days=1;
                             }
                   }
                   else{tt=tt+1;}
                   }
                }
                else{
                //cout<<"after: "<<region2<<"\tmda times: "<<qq<<endl;
                //cout<<"after: "<<region2<<"\tmda times: "<<qq<<endl;
                if(qq==0){
                village2=G_villages[region2];
                village2->mda=YES;
                village2->mda_times=village->mda_times+1;
                listt->listorder=k2;
                village2->mda_days=1;}
                }
            }
        }
        }
    }
}

void updatevillages2(int timestep){
     vill* village;
     listv* listt;
     vill* village2;

    if (timestep>G_mda_start2){
         for(int i =0; i<G_xregion; i++){
                 village=G_villages[i];
                 if (village->mda==YES){
                 //cout<<"which: "<< i <<endl;
                 village->mda_days = village->mda_days+1;}
         }


         for (int j=0; j<G_nstartpoints; j++){
             listt = G_list[j];
             int k=listt->listorder;
             int region=listt->listvill[k];
             village2=G_villages[region];

             if (village2->mda_days>G_dayspervill){
                village2->mda=NO;

                if (k<(G_xregion-1)){
                    k=k+1;
                    //cout<<"regsum:"<<k<<endl;
                    int k2=k;
                    int region2=listt->listvill[k2];
                    //cout<<"list: "<<region2<<"\tk2: "<<k2<< "\tlist num: "<< j <<endl;
                    int qq = G_villages[region2]->mda_times;

                    //cout<<"before: "<<region2<<"\tmda times: "<<qq<<endl;
                    if (qq > 0 ){
                       int tt=0;
                       while(tt==0){
                       if (k<(G_xregion-1)){
                          //cout<<"here: "<< endl;
                          k=k+1;
                          int k2=k;
                          int region2=listt->listvill[k2];

                          qq = G_villages[region2]->mda_times;
                          //cout<<"change: "<< region2 <<"\tmda times: "<<qq<< "\tlist num: "<< j <<endl;
                          if (qq == 0){
                                 tt=tt+1;
                                 listt->listorder=k2;
                                 int region3=listt->listvill[k2];
                                 village2=G_villages[region3];
                                 village2->mda=YES;
                                 village2->mda_times=village->mda_times+1;
                                 listt->listorder=k2;
                                 village2->mda_days=1;
                                 }
                       }
                       else{tt=tt+1;}
                       }
                    }
                    else{
                    //cout<<"after: "<<region2<<"\tmda times: "<<qq<<endl;
                        if(qq==0){
                            village2=G_villages[region2];
                            village2->mda=YES;
                            village2->mda_times=village->mda_times+1;
                            listt->listorder=k2;
                            village2->mda_days=1;
                        }
                    }
                }
             }
         }
    }
}

void updatevillages3(int timestep){
     vill* village;
     listv* listt;
     vill* village2;

    if (timestep>G_mda_start3){
        for(int i =0; i<G_xregion; i++){
             village=G_villages[i];
             if (village->mda==YES){
             //cout<<"which: "<< i <<endl;
             village->mda_days = village->mda_days+1;}
        }


        for (int j=0; j<G_nstartpoints; j++){
             listt = G_list[j];
             int k=listt->listorder;
             int region=listt->listvill[k];
             village2=G_villages[region];

             if (village2->mda_days>G_dayspervill){
                village2->mda=NO;

                if (k<(G_xregion-1)){
                    k=k+1;
                    //cout<<"regsum:"<<k<<endl;
                    int k2=k;
                    int region2=listt->listvill[k2];
                    //cout<<"list: "<<region2<<"\tk2: "<<k2<< "\tlist num: "<< j <<endl;
                    int qq = G_villages[region2]->mda_times;

                    //cout<<"before: "<<region2<<"\tmda times: "<<qq<<endl;
                    if (qq > 0 ){
                       int tt=0;
                       while(tt==0){
                       if (k<(G_xregion-1)){
                          //cout<<"here: "<< endl;
                          k=k+1;
                          int k2=k;
                          int region2=listt->listvill[k2];

                          qq = G_villages[region2]->mda_times;
                          //cout<<"change: "<< region2 <<"\tmda times: "<<qq<< "\tlist num: "<< j <<endl;
                          if (qq == 0){
                                 tt=tt+1;
                                 listt->listorder=k2;
                                 int region3=listt->listvill[k2];
                                 village2=G_villages[region3];
                                 village2->mda=YES;
                                 village2->mda_times=village->mda_times+1;
                                 listt->listorder=k2;
                                 village2->mda_days=1;
                                 }
                       }
                       else{tt=tt+1;}
                       }
                    }
                    else{
                    //cout<<"after: "<<region2<<"\tmda times: "<<qq<<endl;
                        if(qq==0){
                            village2=G_villages[region2];
                            village2->mda=YES;
                            village2->mda_times=village->mda_times+1;
                            listt->listorder=k2;
                            village2->mda_days=1;
                        }
                    }
                }
            }
        }
    }
}




void neighbouring(){
     vill* village;
     vill* village2;
     float lat;
     float longit;
     float lat2;
     float longit2;
     int d;
     int i;
     int rows = G_xregion, cols = G_xregion;
     float** prob = new float*[rows];
     int** nprob = new int*[rows];
     int tot = 5000;

     for (d=0; d<G_xregion; d++){
         village=G_villages[d];
         lat=G_villages[d]->latitude;
         longit=G_villages[d]->longitude;
         //cout<<"tvill: "<<G_villages[d]->latitude<<endl;
         prob[d]=new float[cols];
         float sump[G_xregion];
         float eps =0.05;
         for (i=0; i<G_xregion; i++){
            village2=G_villages[i];
            lat2=G_villages[i]->latitude;
            longit2=G_villages[i]->longitude;

                  if (d==i){
                  prob[d][i]=0;
                  }
                  else{
                  double dist = distanceEarth(lat,longit,lat2,longit2);
                  prob[d][i]=float(village2->npeople)*float(G_villages[i]->npeople)/(float(pow(dist,1.5))+eps);
                  //cout<<"hey2 :\t"<<prob[d][i]<<endl;
                  }
                  if (i==0){
                     sump[i]=prob[d][i];
                  }
                  else{
                    sump[i]=sump[i-1]+prob[d][i];
                  }
                  //cout<<"sum :\t"<<sump[i]<<endl;
         }
         //cout<<"sum :\t"<<sump[G_xregion-1]<<endl;
         prob[d][d]=(sump[G_xregion-1]*G_move);
         //cout<<"new :\t"<<prob[d][d]<<"\tbefore: "<<sump[G_xregion-1]<<endl;
         float nsump=sump[G_xregion-1]+prob[d][d];
         //cout<<"sum :\t"<<nsump<<"\tbefore: "<<sump[G_xregion-1]<<endl;

         int k=0;
         nprob[d]=new int[cols];
         for (int w=0; w<G_xregion; w++){
              prob[d][w]=prob[d][w]/nsump;
              nprob[d][w]=int(prob[d][w]*tot);
              //cout<<"hey2 :\t"<<nprob[d][w]<<endl;
              k=k+nprob[d][w];
              //cout<<"k :\t"<<k<<endl;
              if (k>tot){
                 nprob[d][w]=tot-nprob[d][w-1];
              }
              for (int q=0;q<nprob[d][w];q++){
                  village->nearest.push_back(w);
              }
         }
         //cout<<"village: "<< d<<"\tfirst: "<<village->nearest[1]<<endl;
         //cout<<"\t size"<<village->nearest.size()<<endl;

     //    ofstream fout_nei;
//         ostringstream strfile_nei;
//

     }
     strfile_nei << "outputs/nei-vill.txt" ;
     fout_nei.open(strfile_nei.str().c_str());

     for (int r=0;r<G_xregion;r++){
         village=G_villages[r];
         fout_nei << village->nearest.size() << endl;
     }
     fout_nei.close();
}

void migration_network(){
    vill* village;
    int d;
    int i;
    int rows = G_xregion, cols = G_xregion;
    float** flow = new float*[rows];
    int** nflow = new int*[rows];
    int tot = 5000;

    for (d=0; d<G_xregion; d++){
        village=G_villages[d];
        flow[d]=new float[rows];
        float sumf[G_xregion];

        for (i=0; i<G_xregion; i++){
            if (d==i){
                flow[d][i]=0;
            }
            else{
                flow[d][i] = village->migrantflows[i];
            }
            if (i==0){
                sumf[i]=flow[d][i];
            }
            else{
                sumf[i]=sumf[i-1]+flow[d][i];
            }
        }

        flow[d][d]=(sumf[G_xregion-1]*G_move);
        float nsump=sumf[G_xregion-1]+flow[d][d];

        int k=0;
        nflow[d]=new int[rows];
        for (int w=0; w<G_xregion; w++){
            flow[d][w]=flow[d][w]/nsump;
            nflow[d][w]=int(flow[d][w]*tot);
            k=k+nflow[d][w];
            if (k>tot){
                nflow[d][w]=tot-nflow[d][w-1];
            }
            for (int q=0;q<nflow[d][w];q++){
                village->migrantnet.push_back(w);
            }
        }
        //cout<< "region: "<<d<<" size: "<<  village->migrantnet.size()<<endl;
    }
}

void seasonal_network(){
     vill* village;
     vill* village2;
     float lat;
     float longit;
     float lat2;
     float longit2;
     int d;
     int i;
     int rows = G_xregion, cols = G_xregion;
     float** prob = new float*[rows];
     int** nprob = new int*[rows];
     int tot = 5000;

     for (d=0; d<G_xregion; d++){
         village=G_villages[d];
         lat=G_villages[d]->latitude;
         longit=G_villages[d]->longitude;
         //cout<<"tvill: "<<G_villages[d]->latitude<<endl;
         prob[d]=new float[rows];
         float sump[G_xregion];
         float eps =0.05;
         for (i=0; i<G_xregion; i++){
            village2=G_villages[i];
            lat2=G_villages[i]->latitude;
            longit2=G_villages[i]->longitude;

                  if (d==i){
                    prob[d][i]=0;
                  }
                  else{
                      double dist = distanceEarth(lat,longit,lat2,longit2);
                      prob[d][i]=float(village2->npeople)*float(G_villages[i]->npeople)/(float(pow(dist,.5))+eps);
                  }
                  if (i==0){
                     sump[i]=prob[d][i];
                  }
                  else{
                    sump[i]=sump[i-1]+prob[d][i];
                  }
         }
         prob[d][d]=(sump[G_xregion-1]*G_move);
         float nsump=sump[G_xregion-1]+prob[d][d];

         int k=0;
         nprob[d]=new int[rows];
         for (int w=0; w<G_xregion; w++){
              prob[d][w]=prob[d][w]/nsump;
              nprob[d][w]=int(prob[d][w]*tot);
              k=k+nprob[d][w];
              if (k>tot){
                 nprob[d][w]=tot-nprob[d][w-1];
              }
              for (int q=0;q<nprob[d][w];q++){
                  village->seasonalnet.push_back(w);
              }
         }
        //cout<< "region: "<<d<<" size: "<<  village->seasonalnet.size()<<endl;
     }
}
//
//void static_network(){
//     vill* village;
//    int d;
//    int i;
//    int rows = G_xregion, cols = G_xregion;
//    float** flow = new float*[rows];
//    int** nflow = new int*[rows];
//    int tot = 5000;
//
//    for (d=0; d<G_xregion; d++){
//        village=G_villages[d];
//        flow[d]=new float[rows];
//        float sumf[G_xregion];
//
//        for (i=0; i<G_xregion; i++){
//            if (d==i){
//                flow[d][i]=0;
//            }
//            else{
//                flow[d][i] = village->staticflows[i];
//            }
//            if (i==0){
//                sumf[i]=flow[d][i];
//            }
//            else{
//                sumf[i]=sumf[i-1]+flow[d][i];
//            }
//        }
//
//        flow[d][d]=(sumf[G_xregion-1]*G_move);
//        float nsump=sumf[G_xregion-1]+flow[d][d];
//
//        int k=0;
//        nflow[d]=new int[rows];
//        for (int w=0; w<G_xregion; w++){
//            flow[d][w]=flow[d][w]/nsump;
//            nflow[d][w]=int(flow[d][w]*tot);
//            k=k+nflow[d][w];
//            if (k>tot){
//                nflow[d][w]=tot-nflow[d][w-1];
//            }
//            for (int q=0;q<nflow[d][w];q++){
//                village->staticnet.push_back(w);
//            }
//        }
//    }
//}

void static_network(){
     vill* village;
     vill* village2;
     float lat;
     float longit;
     float lat2;
     float longit2;
     int d, i;
     float crit = 1.5;
     int rows = G_xregion;
     int cols = G_xregion;
     float** prob = new float*[rows];
     int** nprob = new int*[rows];
     int tot = 5000;

     for (d=0; d<G_xregion; d++){
         village=G_villages[d];
         lat=G_villages[d]->latitude;
         longit=G_villages[d]->longitude;
         //cout<<"tvill: "<<G_villages[d]->latitude<<endl;
         prob[d]=new float[rows];
         float sump[G_xregion];
         float eps =0.05;
         for (i=0; i<G_xregion; i++){
            village2=G_villages[i];
            lat2=G_villages[i]->latitude;
            longit2=G_villages[i]->longitude;

                  if (d==i){
                    prob[d][i]=0;
                  }
                  else{
                      double dist = distanceEarth(lat,longit,lat2,longit2);
                      prob[d][i]=float(village2->npeople)*float(G_villages[i]->npeople)/((1.0+(1.0/(1.0+float(exp(10*(float(dist-crit)))))))*(float(pow(dist,.5))+eps));
                  }
                  if (i==0){
                     sump[i]=prob[d][i];
                  }
                  else{
                    sump[i]=sump[i-1]+prob[d][i];
                  }
         }
         prob[d][d]=(sump[G_xregion-1]*G_move);
         float nsump=sump[G_xregion-1]+prob[d][d];

         int k=0;
         nprob[d]=new int[rows];
         for (int w=0; w<G_xregion; w++){
              prob[d][w]=prob[d][w]/nsump;
              nprob[d][w]=int(prob[d][w]*tot);
              k=k+nprob[d][w];
              if (k>tot){
                 nprob[d][w]=tot-nprob[d][w-1];
              }
              for (int q=0;q<nprob[d][w];q++){
                  village->staticnet.push_back(w);
              }
         }
        //cout<< "region: "<<d<<" size: "<<  village->staticnet.size()<<endl;
     }
}


void setsusceptibility(){
     Deme* deme;
     int age;
     float susceptibility;

     for (int d=0; d<G_N; d++){
         deme=G_demes[d];
         age=deme->age;
         //susceptibility = (1-G_r*exp(-G_k*age));
         susceptibility=1.0;
         deme->susceptibility = susceptibility;
         deme->infectiousness = susceptibility;

         if (deme->itn==ITN){
            if (deme->clinical == CLINICAL){
            deme->susceptibility =  deme->susceptibility*(1.0-G_itneffect);
            deme->infectiousness =  deme->infectiousness*(1.0-G_itneffect)*G_phic;
            }
            else if (deme->clinical == ASYMPTOMATIC){
            deme->susceptibility =  deme->susceptibility*(1.0-G_itneffect);
            deme->infectiousness =  deme->infectiousness*(1.0-G_itneffect)*G_phia;
            }
         }
         else if (deme->itn==NOITN){
            if (deme->clinical == CLINICAL){
            deme->susceptibility =  deme->susceptibility;
            deme->infectiousness =  deme->infectiousness*(1.0-G_itneffect)*G_phic;
            }
            else if (deme->clinical == ASYMPTOMATIC){
            deme->susceptibility =  deme->susceptibility;
            deme->infectiousness =  deme->infectiousness*(1.0-G_itneffect)*G_phia;
            }
         }

         // prophylaxis from artesunate uptake
         if (deme->drug==ARTESUNATE || deme->drug==PIPARTE){
            deme->susceptibility =  deme->susceptibility* G_artprophylaxis;
         }
     }
}

void setinfectiousness(){
     Deme* deme;
     int age;
     float infectiousness;

     for (int d=0; d<G_N; d++){
         deme=G_demes[d];
         age=deme->age;
         //susceptibility = (1-G_r*exp(-G_k*age));
         infectiousness=1.0;
         deme->infectiousness = infectiousness;

         if (deme->itn==ITN){
            if (deme->clinical == CLINICAL){
            deme->infectiousness =  deme->infectiousness*G_phic*(1.0-G_itneffect);
            }
            else if (deme->clinical == ASYMPTOMATIC){
            deme->infectiousness =  deme->infectiousness*G_phia*(1.0-G_itneffect);
            }
         }
         else if (deme->itn==NOITN){
            if (deme->clinical == CLINICAL){
            deme->infectiousness =  deme->infectiousness*G_phic;
            }
            else if (deme->clinical == ASYMPTOMATIC){
            deme->infectiousness =  deme->infectiousness*G_phia;
            }
         }
     }
}

void setclinicalprob(){
    Deme* deme;
    for (int d=0; d<G_N; d++){
        deme=G_demes[d];
        int level=deme->immunity_level;
        int moi=deme->moi;
        int cumm=deme->cummulative_exposures;

        deme->probclinicallb=(0.25*exp(-(cumm-1))+exp(-0.6*cumm))/level;
        deme->probclinicalbb=exp(-0.8*(moi-1))*deme->probclinicallb;
    }
}


void setDemes(int maxDemes){
     Deme* newDeme;
	 for(int i=0; i<maxDemes ;i++) {
        G_demes[i] = new Deme();
	 }
}

void setMosqs(int G_xregion){
     mosq* newMosq;
	 for(int i=0; i<G_xregion ;i++) {
             G_mosquitoes[i] = new mosq();
	 }
}

void setvillages(int G_xregion){
     vill* newvillage;
	 for(int i=0; i<G_xregion ;i++) {
        G_villages[i] = new vill();
	 }
}

void setvillagesums(int G_xregion){
     villdata* newvillagedata;
	 for(int i=0; i<G_xregion ;i++) {
        G_villdata[i] = new villdata();
	 }
}

void setlists(int G_nstartpoints){
     listv* newlist;
     for(int i=0; i<G_nstartpoints ;i++) {
             G_list[i] = new listv();
     }
}

void setprobdeath(){
     Deme* deme;
     int age;
     float deathprob;
     float k1;

     for (int d=0; d<G_N; d++){
         deme=G_demes[d];
         age=deme->age;
         k1 = age/G_mu2;
         deathprob = (1-(exp(-k1)))+(1/(G_sig*sqrt(2*3.14)))*exp(-(pow(age-G_mu1,2))/(2*(pow(G_sig,2))));
         deme->death = deathprob*G_mu/0.8;
     }
}

void seasonality (int timestep){
     betaseason=vector<float>(G_xregion,0.0);
     for (int i=0;i<G_xregion;i++){
        betaseason[i] = betas[i]+(G_amp*betas[i]*cos(2.0*3.1416*((double(timestep)-90.0)/365.0)));
        //betaseason[i] = betas[i]+(G_amp*betas[i]*cos(2.0*3.1416*(timestep*G_phase-G_phase)));
        //cout<<"\t:"<<betas[i]<<endl;
    }
}

void openposts (int timestep){
     vill* village;
     //malpost=vector<int>(G_xregion,0);
     for (int i=0;i<G_xregion;i++){
        village= G_villages[i];
        //cout<<"timetoopen:\t"<<mptime[i]<<endl;
        if (timestep>mptime[i]){
           village->post = mp[i];

        //cout<<"open:\t"<<i<<"\ttime:\t"<<timestep<<"\t"<<mp[i]<<endl;
 //       if (timestep>365*5){
//           malpost[i] = 1;
        }
    }
}

void openposts2 (int timestep){
     vill* village;
     malpost=vector<int>(G_xregion,0);
     for (int i=0;i<G_xregion;i++){
        //cout<<"timetoopen:\t"<<mptime[i]<<endl;
        if (timestep>mptime[i]){
           village->post = 1;
        }

        if (timestep>G_mda_start1){
           village= G_villages[i];
           village->post= 1;
        }
    }
}

void countposts(){
   int t=0;
   vill* village;
   for (int i=0;i<G_xregion;i++){
       village = G_villages[i];
       if (village->post== 1){
          t=t+1;
       }
   }
//   cout << "NUMBER POSTS: "<<t<<endl;
}


void settreatrate(){
     //allows for different treatment rates in different regions - aimed at targetting
     //the regions with the highest transmission - could be defined by bitting rates or
     //initial prevalence in pilot studies e.g.
     Deme* deme;
     vill* village;
     int d;
     int t=0;
     for (d=0;d<G_N;d++){
         deme=G_demes[d];
         int region=deme->now;
         int ind=region-1000;
         village = G_villages[ind];
         float b=betas[ind];
         int malp=village->post;
         int mpst=malp;

      //   if (b <=0.02215){
//            treatclinical[d] = (1.0/G_timecomp)*G_fullcourse*G_covab;
//            treatasympt[d] = (1.0/G_timecomp)*G_fullcourse*G_covab;
//            }
//         else {
//              treatclinical[d] = (1.0/G_timecomp)*G_fullcourse*G_covab;
//              treatasympt[d] = (1.0/G_timecomp)*G_fullcourse*G_covab;
//         }

           //cout<<"malpost:\t"<<mpst<<"\tregion\t:"<<ind<<endl;

         if (mpst == 1){
            treatclinical[d] = (1.0/G_timecomp)*G_fullcourse*G_covab;
            //cout<<"postratec\t"<<treatclinical[d]<<endl;
            treatasympt[d] = (1.0/G_timecomp)*G_fullcourse*G_covab;
            }
         else{
              treatclinical[d] = (1.0/G_timecomp)*G_fullcourse*G_covab*G_nomp;
              //cout<<"treatclin: "<<treatclinical[d]<<endl;
              treatasympt[d] = (1.0/G_timecomp)*G_fullcourse*G_covab*G_nomp;
              //cout<<"treatssym: "<<treatasympt[d]<<endl;
         }
     }
}

int countmda(){

     Deme* deme;
     int d;
     int c = 0;

     for (d=0;d<G_N;d++){
         deme=G_demes[d];

         if (deme->mda_times>0){
            c=c+1;
         }
     }
     return c;
}

double PoissonLikelihood(std::vector<double> thedata, std::vector<double> themodel){

    double like = 1.0e9;  // default value
    if (thedata.size()==themodel.size()){
        like = 0.0;
        for (int i=0;i<thedata.size();i++){
            double mu = themodel[i];
            if (mu<=0.0) mu = 1.0e-9;
            double x = thedata[i];
            like = like-mu+x*log(mu);
         }
    }else{
        cout << "(UsefulUtils::PoissonLikelihood): vectors need to be the same size! "
             << thedata.size() << " "
             << themodel.size() << " "
             << endl;
    } // end check that the vectors are appropriate sizes

//    like = -like; // return the *negative* log likelihood
    return like;
}

double NegativeBinomialLikelihood(std::vector<double> thedata,std::vector<double> themodel,double ascale){

      double like = 1.0e9;  // default value
      //http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0000180
      if (thedata.size()==themodel.size()){
         double k = 1.0/ascale;
         like = 0.0;
         for (int i=0;i<thedata.size();i++){
            double mu = themodel[i];
            if (mu<=0.0) mu = 1.0e-9;
            double x = thedata[i];
            double alike = LogGamma(k+x)
                         - LogGamma(x+1.0)
                         - LogGamma(k)
                         + x*log(mu/(mu+k))
                         - k*log(1.0+mu/k);
            //cout << i  << " "
                 //<< mu << " "
                 //<< x  << " "
                 //<< k  << " "
                 //<< LogGamma(k+x) << " "
                 //<< LogGamma(x+1.0) << " "
                 //<< LogGamma(k) << " "
                 //<< x*log(mu/(mu+k)) << " "
                 //<< k*log(1.0+mu/(k)) << " "
                 //<< alike << " "
                 //<< like << endl;
            like = like + alike;
         }

         /*
         double datasum = Sum(thedata);
         like = 0.0;
         for (int i=0;i<thedata.size();i++){
            double atemp=0.0;
            if ((thedata[i]-1)>=1){
               //vector<double> rvec(thedata[i]-1) ;
               //std::generate (rvec.begin(), rvec.end(), Generator()); // Fill with 0, 1, ..., 99.
               //rvec = Scale(rvec,ascale);
               //vector<double> x = rvec;
               //transform(rvec.begin()
               //         ,rvec.end()
               //         ,x.begin()
               //         ,bind2nd(plus<double>()
               //         ,double(1.0)));
               //atemp = Sum(x);

               for (int r=1;r<=(thedata[i]-1);r++){
                 atemp = atemp + log(1.0 + ascale*double(r));
               }
            }
            like = like
                 + atemp
                 - thedata[i]*log(ascale)
                 + thedata[i]*log(ascale*themodel[i]/datasum)
                 - (thedata[i]+1.0/ascale)*log(1.0+ascale*themodel[i]/datasum);
         }
         */
      }else{
        cout << "(UsefulUtils::NegativeBinomialLikelihood): vectors need to be the same size! "
             << thedata.size() << " "
             << themodel.size() << " "
             << endl;
      } // end check that the vectors are appropriate sizes

      like = -like; // return the *negative* log likelihood
      return like;
}


float gasdev(long g_seed)
{
	static int iset=0;
	static float gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran2()-1.0;
			v2=2.0*ran2()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

double square(float x)
{
   double square_of_x;
   square_of_x = x * x;
   return square_of_x;
}

double LogGamma(double x){
    double a;
    double b;
    double c;
    double p;
    double q;
    double u;
    double w;
    double z;
    double logpi;
    double ls2pi;
    double tmp;
    double result;

    logpi = log(acos(-1.0));
    ls2pi = 0.91893853320467274178;

    if (x<(-34.0)){
        q = -x;
        w = LogGamma(q);
        p = double(int(q));
        //i = ae_round(p, _state);

        z = q-p;
        if (z>0.5){
            p = p+1;
            z = p-q;
        }
        z = q*sin(acos(-1.0)*z);
        result = logpi-log(z)-w;
        return result;
    }

    if (x<13.0){
        z = 1;
        p = 0;
        u = x;
        while (u>=3.0){
            p = p-1;
            u = x+p;
            z = z*u;
        }

        while(u<2){
            z = z/u;
            p = p+1;
            u = x+p;
        }

        if (z<0){
            z = -z;
        }

        if(u==2){
            result = log(z);
            return result;
        }
        p = p-2;
        x = x+p;
        b = -1378.25152569120859100;
        b = -38801.6315134637840924+x*b;
        b = -331612.992738871184744+x*b;
        b = -1162370.97492762307383+x*b;
        b = -1721737.00820839662146+x*b;
        b = -853555.664245765465627+x*b;

        c = 1;
        c = -351.815701436523470549+x*c;
        c = -17064.2106651881159223+x*c;
        c = -220528.590553854454839+x*c;
        c = -1139334.44367982507207+x*c;
        c = -2532523.07177582951285+x*c;
        c = -2018891.41433532773231+x*c;
        p = x*b/c;

        result = log(z)+p;
        return result;
    }

    q = (x-0.5)*log(x)-x+ls2pi;

    if (x>100000000){
        result = q;
        return result;
    }

    p = 1/(x*x);

    if (x>=1000.0){
        q = q+((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
    }else{
        a = 8.11614167470508450300*0.0001;
        a = -5.95061904284301438324*0.0001+p*a;
        a = 7.93650340457716943945*0.0001+p*a;
        a = -2.77777777730099687205*0.001+p*a;
        a = 8.33333333333331927722*0.01+p*a;
        q = q+a/x;
    }
    result = q;
    return result;
 }


