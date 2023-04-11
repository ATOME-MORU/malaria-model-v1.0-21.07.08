#include <time.h>
#include "poisson.h"
#include "aguasSIRk.h"
#include <vector>
#include <omp.h>
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <cstdlib>

// ciclo de vida da metapopulacao : infection- migracao-selecao(sd=0)-amostragem.
//caso o numero de infectados seja 0 o program para

int main(int ac, char **av){
    int cl=clock();
    G_xregion=962;
    G_yregion=6;
    G_nstartpoints=50;
    G_runs=10;
    G_runs=1;
    float G_prev = 0.06;
    float betaf=0.83;
    int nummdas = 1;
    float initresa = 1.0;
    G_trigger = 4.0*G_prev;
    G_tauab = 1.0/1.5;
    G_phia = 0.5;
    G_delta = 1.0/80.0;   // natural recovery rate
    G_alpha = 1.0/60.0;    // immunity loss

    int thisrun = 0;
    while( (++thisrun)<=G_runs ){
    int timestep = 0;
    int i,d;

    VectorSpace* vregion[G_xregion];
    listv* ls[G_nstartpoints];
    listv* ordl[G_nstartpoints];
    G_list = vector<listv*>(G_nstartpoints, new listv);
    setlists(G_nstartpoints);

    for (int i =0;i<G_nstartpoints;i++){
        listv* lis = new listv();
        lis->listvill = vector<int>(G_xregion,0);
        ls[i] = lis;
    }

    ifstream ifs2( "data/order100last.txt"  );
    for(int i= 0;i<G_nstartpoints;i++){
        vector<int> space_aux2; //criar novo auxiliar
        space_vector2( ifs2, space_aux2 ); //ler para o auxiliar
        G_list[i]->listvill = space_aux2; //guardar o vector por apontador em vregion[i]
        G_list[i]->listorder = 0;
    }


    VectorSpace2* vmigrate[G_xregion];
    for(int i=0; i<G_xregion; i++){
        VectorSpace2* vm = new VectorSpace2();
        vm->v2 = vector<float>(G_xregion,0.0);
        vmigrate[i] = vm;
    }
    ifstream ifs3( "data/movmigratejuly.txt"  );
    for(int d=0; d<G_xregion; d++){
        vector<float> space_aux3; //criar novo auxiliar
        migrate_vector( ifs3, space_aux3 ); //ler para o auxiliar
        vmigrate[d]->v2 = space_aux3; //guardar o vector por apontador em vregion[i]
    }

    VectorSpace3* vstatic[G_xregion];
    for(int i=0; i<G_xregion; i++){
        VectorSpace3* vst = new VectorSpace3();
        vst->v22 = vector<float>(G_xregion,0.0);
        vstatic[i] = vst;
    }
    ifstream ifs4( "data/movmigratejuly.txt"  );
    for(int d=0; d<G_xregion; d++){
        vector<float> space_aux4; //criar novo auxiliar
        static_vector( ifs4, space_aux4 ); //ler para o auxiliar
        vstatic[d]->v22 = space_aux4; //guardar o vector por apontador em vregion[i]
    }

//    for(int i= 0;i<G_nstartpoints;i++){
//            listv* listt[G_nstartpoints];
//            listv* lis = new listv();
//          lis->listvill = G_list[i]->listvill;
//          listt[i] = lis;
//            cout << "vector: " << G_list[i]->listvill[0] << endl;
////               for(int j= 0;j<G_xregion;j++){
////                       cout<<" l:"<< listt[i]->listvill[j]<<" ";
////               }
//               //cout<< "order\t"<<G_list[i]->listorder<< endl;
//                //cout<<"first: "<< vregion[i]->v[0]<<" "<< endl;
//    }
//
//    int gg=getchar();

    G_villages = vector<vill*>(G_xregion, new vill);
    setvillages(G_xregion);
    G_villdata = vector<villdata*>(G_xregion, new villdata);
    setvillagesums(G_xregion);

    for(int i =0;i<G_xregion;i++){
        VectorSpace* vs = new VectorSpace();
        vs->v = vector<float>(G_yregion,0.0);
        vregion[i] = vs;
    }

    ifstream ifs( "data/regions_datajuly2.txt"  );
    for(int i=0; i<G_xregion; i++){
        vector<float> space_aux; //criar novo auxiliar
        space_vector( ifs, space_aux ); //ler para o auxiliar
        vregion[i]->v = space_aux; //guardar o vector por apontador em vregion[i]
    }

    //for(int i= 0;i<G_xregion;i++){
               //cout << "vector: " << i << endl;
      //         for(int j= 0;j<G_yregion;j++){
        //               cout<<" conteudo: "<< vregion[i]->v[j]<<" ";
          //     }
            //    cout<< "\t"<<endl;
              //  cout<<"first: "<< vregion[i]->v[0]<<" "<< endl;
    //}



    G_N=0;
    pop=vector<int>(G_xregion,0);
    mob=vector<float>(G_xregion,0.0);
    number_pop_area=vector<int>(G_xregion,0);
    i_init=vector<int>(G_xregion,0);
    lon=vector<float>(G_xregion,0.0);
    lat=vector<float>(G_xregion,0.0);
    betas=vector<float>(G_xregion,0.0);
    mp=vector<int>(G_xregion,0);
    mptime=vector<int>(G_xregion,0);
    malpost=vector<int>(G_xregion,0);
    pa=vector<float>(G_xregion,0.0);
    pb=vector<float>(G_xregion,0.0);

    for (int i=0;i<G_xregion;i++){
        G_N=G_N+int(vregion[i]->v[0]);
        int number=int(vregion[i]->v[0]);
        int I=int((vregion[i]->v[3]*5.1/0.75)*vregion[i]->v[0]);
        I = ( I > number ? number : I );
        float lo=vregion[i]->v[1];
        float la=vregion[i]->v[2];
        float beta_values=betaf*vregion[i]->v[3];
        int malpst=int(vregion[i]->v[4]);
        int malpstime=int(vregion[i]->v[5]);
        float mobpop = float(vregion[i]->v[6]);
        vector<float>flows = vmigrate[i]->v2;
        vector<float>sflows = vstatic[i]->v22;
        float presa=float(initresa*(vregion[i]->v[7]));

        pop[i]=G_N;
        mob[i]=mobpop;
        number_pop_area[i]=number;
        i_init[i]=I;
        lon[i]=lo;
        lat[i]=la;
        betas[i]=beta_values;
        mp[i]=malpst;
        mptime[i]=malpstime;
        pa[i]=presa;

        //cout<< "betas:"<< mptime[i]<<"\t"<<endl;
        villdata* villd;
        vill* village;
        village=G_villages[i];
        villd=G_villdata[i];
        village->npeople = number;
        villd->popsize = number;
        village->mda = NO;
        village->mda_days=0;
        village->mda_times=0;
        G_villages[i]->post=malpst;
        G_villages[i]->latitude=lat[i];
        G_villages[i]->longitude=lon[i];
        G_villages[i]->migrantflows=flows;
        G_villages[i]->staticflows=sflows;
        villd->lat=la;
        villd->lon=lo;
        G_villdata[i]->init();
        //cout<<"lineN: "<<vregion[i]->v[0]<<endl;
        //cout<<"pop:"<<pop[i]<< endl;
        //cout<<"i: "<<i_init[i]<<endl;
        //cout<<"longmax: "<<longmax[i]<<"\tlongmin: "<<longmin[i]<<endl;
    }
    cout<<"N:\t"<<G_N<<endl;

    //int f=getchar();

    /***********************************/
    //cout << "setting output files..." << endl;
    G_move = 0.95;

    G_sensitivity = 0.8;
    G_mteams = G_nstartpoints;
    G_dayspervill = 4;
    G_delay = 60;
    G_first = 365*5+G_delay;

    /******** mobility  ***********/
    G_mobility_static = 125.0/365.0;    // static pop probability of moving to another village
    G_static_return = 1.0/3.0;         // static pop average length stay in another village
    G_mobility_mobile = 125.0/365.0;     // seasonal migrant prob of moving to another village
    G_mobile_return = 1.0/3.0;         // seasonal migrant average length stay in another village
    G_migrate = 0.0;                   // seasonal migrant prob to migrate to another village other than home place - redefined in function
    G_migrate_mobile = 1.0/90.0;


//    strfile_village << "C:\\051115\\IBM\\mds_mapping\\village" ;
//    //strfile_age << "D:/Documents and Settings/ricaguas/Desktop/Individual_Based/age" ;
//    strfile_village << "-D-" <<  G_N << "-move-" << G_move << "-trigger-"<< G_trigger ;
//    strfile_village << "-2.txt";
//    fout_village.open(strfile_village.str().c_str());

//    strfile_taus << "C:\\051115\\IBM\\mds_mapping\\taus" ;
//    strfile_taus << "-D-" <<  G_N << "-move-" << G_move ;
//    strfile_taus << ".txt";
//    fout_taus.open(strfile_taus.str().c_str());
//
//    strfile_age << "C:\\051115\\IBM\\mds_mapping\\agestruct100.txt" ;
//    fout_age.open(strfile_age.str().c_str());
//
    strfile_regions << "outputs/regionresults"<< "-delta-"<< G_delta<<"-trigger-"<< G_trigger<< "-mobility-"<<G_mobility_static<< "-G_phia-"<< G_phia <<".txt" ;
    fout_regions.open(strfile_regions.str().c_str());

    /***********************************/
    allocFactorial();
    G_groups = new StateDeme();
    G_groups->init();
    printTreat();
    /***********************************/

    // G_tmax = 365*50;
    G_tmax = 365*5;
    G_tmax = 40;

    G_maxa = 85;

    G_amp = 0.9;
    G_phase = 90.0/365.0;

    G_pa = initresa;   // initial artemisinin resistance level
    G_pb = 0.01;
    G_po = 1.0-G_pa-G_pb;

    G_age1 = 6.802e-03;
    G_age2 = 6.755e-03;
    G_age3 = -7.683e-04;
    G_age4 = 3.663e-05;
    G_age5 = -8.972e-07;
    G_age6 = 1.181e-08;
    G_age7 = -7.965e-11;
    G_age8 = 2.164e-13;


    G_mu = 0.0222/365.0;   // death rate 0.0222 - 45 yrs life expectancy
    G_gamma = 1.0/5.0;     // liver to blood
    G_sigma = 1.0/15.0;    // blood to gametocytes

    G_natality = 0.062/365.0; // for constant pop size

    G_artprophylaxis = 1.0; //change
    G_c = 1.0;              // mosquito infectiousness to humans - change this
    G_phic = 1.0;

    G_eta = 0.1;       //change
    G_covmsat = 0.1;   //change

    G_radius = 2.0; //radius squared

    G_k = 0.14;
    G_r = 0.99;
    G_ke = 1.0;
    G_sig = 1.0/15.0;
    G_mu1 = 80.0;
    G_mu2 = 50.0;

    G_propRxa = 0.052;
    G_pab = 0.62;       //MST 0.8- value dependes on intervention
    G_covab = 0.9;
    G_kab = 900.0;
    G_lam = 300.0;
    G_propRxi1 = 0.8;
    G_propRxai1 = 0.8;

    G_timecomp = 1.5;
    G_fullcourse = 0.9;
    G_nomp = 0.1;
    G_pclin = 0.9;
    G_pclinimm = 0.1;
    G_mellow = 1.0/10.0;
    G_initclin = 0.02;
    //G_prop_immune = 0.2;

    G_tau = (1.0/7.0)*G_propRxa;
    //G_tau1 = G_propRxi1*G_taumagi1;
    G_tau2 = 14.0/365.0; //change

    G_mosqdeath1 = 1.0/7.0;   //life expectancy of an infectious mosquito
    G_mosqdeath2 = 1.0/20.0;  //life expectancy of an exposed mosquito
    G_mosqinc = 1.0/14.0;     // extrinsic incubation period

    G_mda_start1=365*5;
    G_mda_end1=TIMELIMIT;
    G_mda_start2=365*6;
    G_mda_end2=TIMELIMIT;
    G_mda_start3=TIMELIMIT;
    G_mda_end3=TIMELIMIT;
    G_switch_mda1=false;
    G_switch_mda2=false;
    G_switch_mda3=false;

        // average length of stay of mobile migrants in each village

    /******** treatment in transmission  ***********/
     G_cBroda = 1.0/5.0;      // richard - 1/7
     G_cBrada = G_cBroda*0.27;
     G_cBrbda = G_cBroda;

     G_clroda = 1.0/3.0;       // richard - 1/4
     G_clrada = 0.27*G_clroda;
     G_clrbda = G_clroda;

     G_cBrodab = (1.0/5.0);    // richard - 1/7
     G_cBradab = 0.27*G_cBrodab + (1.0-0.27)*G_cBrodb; // Uninitialised variable?
     G_cBrbdab = 0.8*G_cBrodab + (1.0-0.8)*G_cBroda;

     G_clrodab = 1.0/3.0;      // richard - 1/3
     G_clradab = 0.27*G_clrodab + (1.0-0.27)*G_clrodb; // Uninitialised variable?
     G_clrbdab = 0.8*G_clrodab + (1.0-0.8)*G_clroda ;

     G_cBrodb = 1.0/3.0;
     G_cBradb = 1.0/3.0;
     G_cBrbdb = (1.0/3.0)*0.8;

     G_clrodb = (1.0/15.0);
     G_clradb = G_clrodb;
     G_clrbdb = G_clrodb*0.8;

     G_cBradlf = 1.0/3.0;
     G_cBrbdlf = 1.0/3.0;
     G_cBrodlf = 1.0/3.0;
     G_clradlf = 1.0/7.0;
     G_clrodlf = 1.0/7.0;
     G_clrbdlf = 1.0/7.0;

     G_cBradal = 1.0/5.0;
     G_cBrbdal = 1.0/5.0;
     G_cBrodal = 1.0/5.0;
     G_clradal = 1.0/3.0;
     G_clrbdal = 1.0/3.0;
     G_clrodal = 1.0/3.0;

     G_clrdpmc = 1.0/1.4;
     G_clrdpma = 1.0/1.1;

     G_xa0 = (1.0/7.0);     // lose ART
     G_xai = 1.0/3.0;       // lose DHA in ACT
     G_xb = 1.0/30.0;
     G_xab = (1.0/30.0);    // lose PIP

     G_xprim = (1.0/2.0);   // lose primaquine
     G_xala = 1.0/3.0;      // lose ART in AL
     G_xall = 1.0/15.0;     // lose lumefantrine in AL
     G_xlum = 1.0/15.0;     // lose LUM

     G_antroph = 0.95;
     G_indoorpref = 0.95;
     G_xitn =(0.0)/365.0;   //change
     G_itneffect = 0.5*G_antroph*G_indoorpref;
     G_itnprop = 0.2;       // change or use data to inform

    /***********************************/
    //cout << "setting demes and parasites..." << endl;
    G_demes = vector<Deme*>(G_N, new Deme);
    setDemes(G_N);
    G_mosquitoes = vector<mosq*>(G_xregion, new mosq);
    setMosqs(G_xregion);
    //cout << "done." << endl;
    treatab=vector<float>(G_N,0.0);
    treatclinical=vector<float>(G_N,0.0);
    treatasympt=vector<float>(G_N,0.0);
    seasonvec=vector<float>(G_tmax,0.0);
    /***********************************/

        // g_seed=-251978884+thisrun;
        std::ostringstream strfile_resistance;

        strfile_resistance << "outputs/logmda-run-"<<thisrun << "-prev-"<<G_prev<<"-G_phia-"<< G_phia ;
        strfile_resistance << "-mdacov-" <<  G_tauab << "-delta-"<< G_delta << "-teams-" << G_mteams << "-trigger-"<< G_trigger<< "-resistinit-"<< pa[0];
        strfile_resistance << "-mobility-"<< G_mobility_static <<".txt";
        fout_resistance.open(strfile_resistance.str().c_str(),ios::out);

//        //fout_resistance.open(strfile_resistance.str().c_str());
//         if (!fout_resistance.is_open()) {
//              cerr << "error: open file for output failed!" << endl;
//              abort();
//           }

        VectorSpaceage* vage[G_N];
        for(int i=0; i<G_N; i++){
            VectorSpaceage* vag = new VectorSpaceage();
            vag->v2a = vector<int>(4,0);
            vage[i] = vag;
        }
        ifstream ifage( "data/agestruct100.txt"  );
        for(int i=0; i<G_N; i++){
            vector<int> space_aux; //criar novo auxiliar
            age_vector( ifage, space_aux ); //ler para o auxiliar
            vage[i]->v2a = space_aux; //guardar o vector por apontador em vregion[i]
        }
//        for(int i= 0;i<G_N;i++){
//               cout << "vector: " << i << endl;
//               for(int j=0;j<4;j++){
//                    cout<<" conteudo: "<< vage[i]->v2a[j]<<" ";
//               }
//                cout<< "\t"<<endl;
//                cout<<"first: "<< vage[i]->v2a[0]<<" "<< endl;
//        }


        agev=std::vector<int>(G_N,0);
        moiv=std::vector<int>(G_N,0);
        clinstat=std::vector<int>(G_N,0);
        levelimm=std::vector<int>(G_N,0);
        cumminf=std::vector<int>(G_N,0);

        for (int i=0;i<G_N;i++){
            int agev=int(vage[i]->v2a[0]);
            int moiv=int(vage[i]->v2a[1]);
            int clinstat=int(vage[i]->v2a[2]);
            int levelimm=int(vage[i]->v2a[3]);
            int cumminf=int(vage[i]->v2a[4]);

            G_demes[i]->age=agev;
            G_demes[i]->moi=moiv;
            G_demes[i]->clinical=clinstat;
            G_demes[i]->immunity_level=levelimm;
            G_demes[i]->cummulative_exposures=cumminf;

            if (G_demes[i]->immunity_level>0){
                G_demes[i]->immunity=IMMUNE;
                G_demes[i]->immunity_days=40;
            }
            else{G_demes[i]->immunity=NONIMMUNE;}
        }


        VectorSpace5* seasonvec[G_tmax];
        for(int i=0; i<G_tmax; i++){
            VectorSpace5* vss = new VectorSpace5();
            vss->vseas= vector<float>(G_xregion,0.0);
            seasonvec[i] = vss;
        }
        ifstream ifseason( "data/seasonality_100yrs.txt"  );
        for(int d=0; d<G_tmax; d++){
            vector<float> space_aux4; //criar novo auxiliar
            season_vector( ifseason, space_aux4 ); //ler para o auxiliar
            seasonvec[d]->vseas = space_aux4; //guardar o vector por apontador em vregion[i]
        }
//        for(int i= 0;i<G_tmax;i++){
//            cout << "time: " << i << endl;
//            cout<<"beta: "<< seasonvec[i]->vseas[0]<< endl;
//        }

            /***********************************/
            cout << "\trunning network " << thisrun << endl;
            //cout << "\tinitializing... " << endl;
            /***********************************/

            fitnessinit();
            cout<<"before"<<endl;
            printResistance(timestep);
            cout<<"res"<<endl;
            printmosqs();
            village_here_now();
            surveytimes();
            //agedistribution();
            setsusceptibility();

            /** interventions **/
            //startvillagestrigger();
            startvillages();
            //cout<<"startvillages complete"<<endl;

            seasonality(timestep);
            openposts(timestep);
            countposts();

            /** connectivity networks **/
            //cout<<"before"<<endl;
            seasonal_network();
            migration_network();
            static_network();
            neighbouring();
            //cout<<"after"<<endl;
            //cout << "\tdone neighbouring " << endl;

            infectmosquitoes();
            settreatrate(); // this uses ->now but is commented in the simulation loop
            G_groups->printTransmission();
            printResistance(timestep);
            //printTreat();
            setprobdeath();
            //printbednet();
            //printage();
            //int c = getchar();
            //cout<<"12"<<endl;
            /***********************************/
            //cout << "\tgoing for " << G_tmax << " generations..." << endl;
            /***********************************/


            while( timestep<G_tmax ){
            timestep++;
                    //if(timestep%(G_tmax/20) == 0)
                   cout << "going for timestep: " << timestep << endl;
                   //printagestruct(timestep);

                   seasons=seasonvec[timestep-1]->vseas[0];
//                   cout << "seasonality done " << endl;

                   squarepulse(timestep);
                   village_here_now(); /** creates a vector for each village with the index of the individuals currently in that village **/
                   //birthdeath(); /** given their death probability evaluate if each person is to die. If so, replace that person with a newborn **/
                   //superbirth();
                   //settreatrate(); /** set the treatment probabilities depending on presentation of symptoms and whether there is a malaria post **/
                   //cout<<"1"<<endl;
                   //ageing(timestep); /** increase individual's age by one year every 365 days **/
                   //seasonality(timestep); /** explicit formulation for the biting rate fluctuation over time **/
                   //openposts(timestep);
                   //openposts2(timestep); /** make all villages malaria posts upon MDA start **/
                   //setsusceptibility(); /** sets the susceptibility to infection - depends on ITN usage and maybe age **/

                   /** interventions **/
                   //updatevillages(timestep);
                   //triggerMDA(timestep); /** checks if at the time of the surveys, if any village is over the TME define threshold. creates a list of
                                            //all villages over the threshold and assigns them an MDA **/
                   /** 2nd MDA **/
                   startvillages2(timestep);
                   updatevillages2(timestep);

                   /** 3rd MDA **/
                   //startvillages3(timestep);
                   //updatevillages3(timestep);

                   //livertoblood(); /** evaluates the transition from liverstage to bloodstage **/
                   //drugeffect(); /** checks the parasite clearance given the circulating drugs and their effectiveness given the parasite resistance **/
                   //tosymptoms(); /** evaluates the transition from bloodstage to gametocyte carriage #infectious **/
                   //recovery(); /** evaluates the transition from infectious to susceptible, granting clinical immunity **/

//                   immunitycount();  /** counts days since recovery **/
//                   setprobdeath(); /** sets the probability of dying depending on age **/
//                   clinicalresolution(); /** evaluates the transition from clinical to asymptomatic **/
//                   immunityloss(); /** evaluates the loss of clinical immunity **/
//                   superinfectiontoblood();
//                   superinfectiontoinfectious();
//                   superrecovery();
//                   superdrugeffect();

                   villages(timestep, seasons);
//                   cout << "villages done " << endl;
                   naturalhistoryinfection(timestep);
//                   cout << "humans done " << endl;
                   mosquitodynamics();
//                   cout << "mosquitoes done " << endl;
                   population_movement(timestep);
//                   cout << "movement done " << endl;

                   //infectmosquitoes(); /** generates a list of infected mosquitoes given the biting rates and the infected people in each village **/
                   //survivemosqs(); /** evaluates the survival of infected and infectious mosquitoes **/
                   //mosqincubation(); /** checks if infected mosquitoes turn infectious **/
//                   printmosqs();

                   treatment(timestep);
                   //treattrigger(timestep); /** deploys drugs according to the village MDA status and hands out baseline treatment based on symptomatology and whether the village is a malaria post **/
//                   cout << "treatment done " << endl;
                   infection(); /** given the infectious mosquitoes list for each village generate new human infections **/
//                   cout << "infection done " << endl;
                   //drugloss();  /** clears drugs from circulation according to PK-PD parameters **/
                   //itnloss();
                   //countposts(); /** check the number of malaria posts **/

                   //mobility();      /** short term mobility **/
                   //cout<<"before"<<endl;
                   //migration(timestep); /** long term migration to a different village **/
                   //cout<<"middle"<<endl;
                   //temp_migration(timestep); /** transient seasonal migration **/
                   //cout<<"after"<<endl;

                   /** print outputs **/
                   //printregion(timestep);
//                   cout << "print regions done " << endl;
                   printResistance(timestep);

                   //printinfo2(timestep);
                   //printTreat();
                   //cout<<"44"<<endl;
                   //printbednet();
//                   printagestruct(timestep);
//                   cout << "print age done " << endl;

            cout << endl;
            }
            updatevilldata(timestep);
            fout_resistance.flush();
            fout_resistance.close();
            fout_age.flush();
            fout_age.close();
     }

     //cout << "cleaning memory...";
     printf("model ran in %lg seconds\n", ((double) (clock()-cl))/CLOCKS_PER_SEC);
     cout << endl;
     int c = getchar();
     exit(1);
}

/** this makes the program wait for the user to press any key **/
void badStateProblem(char * message){
     cout << message << endl;
     int c = getchar();
}

void drugloss(){
     Deme* deme;
     float randomEvent, randomEvent2;
     int d;

     for (d=0; d<G_N; d++){
         deme=G_demes[d];
         int region=G_demes[d]->home;
         int ind=region-1000;

        if (G_demes[d]->drug == NODRUG ){
            continue;
        }

        else if(G_demes[d]->drug == ARTESUNATE ){
            randomEvent=ran2();
            if (randomEvent<=G_xa0){
            G_demes[d]->drug = NODRUG;
                G_groups->SUM_ART = G_groups->SUM_ART - 1;
                G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART - 1;
                if (G_demes[d]->transmission ==  BLOODSTAGE ||  G_demes[d]->transmission ==  INFECTIOUS){
                   G_groups->SUM_FAIL = G_groups->SUM_FAIL + 1;
                }
                continue;
            }
         }

         else if(G_demes[d]->drug == PIPARTE ){
                 randomEvent=ran2();
                 randomEvent2=ran2();
                 if (randomEvent<=G_xai){
                    G_demes[d]->drug = PIPERAQUINE;
                    G_groups->SUM_PIP = G_groups->SUM_PIP + 1;
                    G_villdata[ind]->SUM_PIP = G_villdata[ind]->SUM_PIP + 1;
                    G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE - 1;
                    G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE - 1;
                    continue;
                 }
                 else if (randomEvent2<=G_xab){
                    G_demes[d]->drug = ARTESUNATE;
                    G_groups->SUM_ART = G_groups->SUM_ART + 1;
                    G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE - 1;
                    G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART + 1;
                    G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE - 1;
                    continue;
                 }
         }

         else if(deme->drug == PIPERAQUINE ){
                 randomEvent=ran2();
                 if (randomEvent<=G_xab){
                    G_demes[d]->drug = NODRUG;
                    G_groups->SUM_PIP = G_groups->SUM_PIP - 1;
                    G_villdata[ind]->SUM_PIP = G_villdata[ind]->SUM_PIP - 1;
                    if (G_demes[d]->transmission !=  SUSCEPTIBLE){
                       G_groups->SUM_FAIL = G_groups->SUM_FAIL + 1;
                    }
                    continue;
                 }
         }


         else if(G_demes[d]->drug == AL ){
                 randomEvent=ran2();
                 randomEvent2=ran2();
                 if (randomEvent<=G_xala){
                    G_demes[d]->drug = LUMEFANTRINE;
                    G_groups->SUM_LUM = G_groups->SUM_LUM + 1;
                    G_villdata[ind]->SUM_LUM = G_villdata[ind]->SUM_LUM + 1;
                    G_groups->SUM_AL = G_groups->SUM_AL - 1;
                    G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL - 1;
                    continue;
                 }
                 else if (randomEvent2<=G_xall){
                    G_demes[d]->drug = ARTESUNATE;
                    G_groups->SUM_ART = G_groups->SUM_ART + 1;
                    G_groups->SUM_AL = G_groups->SUM_AL - 1;
                    G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART + 1;
                    G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL - 1;
                    continue;
                 }
         }

         else if(G_demes[d]->drug == LUMEFANTRINE ){
                 randomEvent=ran2();
                 if (randomEvent<=G_xlum){
                    G_demes[d]->drug = NODRUG;
                    G_groups->SUM_LUM = G_groups->SUM_LUM - 1;
                    G_villdata[ind]->SUM_LUM = G_villdata[ind]->SUM_LUM - 1;
                    if (G_demes[d]->transmission !=  SUSCEPTIBLE){
                       G_groups->SUM_FAIL = G_groups->SUM_FAIL + 1;
                    }
                    continue;
                 }
         }


         if (G_demes[d]->prim == YES ){
             randomEvent=ran2();
                 if (randomEvent<=G_xprim){
                    G_demes[d]->prim = NO;
                    G_groups->SUM_PRIM = G_groups->SUM_PRIM - 1;
                    G_villdata[ind]->SUM_PRIM = G_villdata[ind]->SUM_PRIM - 1;
                 }
         }

    } // close for loop
}

void itnloss(){
     Deme* deme;
     float randomEvent;

     for (int d=0; d<G_N; d++){
          deme=G_demes[d];
          if (G_demes[d]->itn == NOITN){
               continue;}

          else if(G_demes[d]->itn == ITN ){
                 randomEvent=ran2();
                 if (randomEvent<=G_xitn){
                    G_demes[d]->itn = NOITN;
                    G_groups->SUM_ITN = G_groups->SUM_ITN - 1;
                    //G_groups->SUM_ITN = G_groups->SUM_ITN + 1;
                    continue;}
          }
     }
}


void drugeffect(){
     Deme* deme;
     float randomEvent;
     float randomEvent2;
     float paramToTest = 0.0;
     int cc = 0;

     for (int d=0; d<G_N; d++){
         deme=G_demes[d];
         int region=G_demes[d]->home;
         int ind=region-1000;
              //look for the classes of the deme
              //choose the correct probability and save it in paramToTest
            if (G_demes[d]->drug == NODRUG ){
               continue;}

            else if(G_demes[d]->drug == ARTESUNATE){
                     if(G_demes[d]->resistance==RESISTRA){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBrada){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          }    //close bloodstage

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrada){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRA =  G_villdata[ind]->SUM_IRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close Infectious

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;
                               }
                     } //close resistA

                     else if(G_demes[d]->resistance==RESISTRB){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                          /*****/
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBrbda){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRB = G_groups->SUM_BRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close bloodstage

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrbda){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close Infectious

                          else if (G_demes[d]->transmission==LIVERSTAGE || G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;
                               }
                     } //close resistB

                     else if(G_demes[d]->resistance==RESISTRO){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBroda){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRO = G_groups->SUM_BRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }

                          } //close bloodstage

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clroda){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } // close Infectious

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     } //close resistO

                     else {cout << " dstate " << G_demes[d]->drug << " rstate "<< G_demes[d]->resistance << " tstate " << G_demes[d]->transmission <<endl;
                     cout<<"drugeffect: bad resitance state "<<endl;}
             } //close artesunate


             else if(G_demes[d]->drug == PIPERAQUINE){
                     if(G_demes[d]->resistance==RESISTRA){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBradb){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } // close bloodstage

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clradb){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infectious

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     } //close resistA

                     else if(G_demes[d]->resistance==RESISTRB){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                          /*****/
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBrbdb){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRB = G_groups->SUM_BRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrbdb){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     } //close resistB

                     else if(G_demes[d]->resistance==RESISTRO){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBrodb){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRO = G_groups->SUM_BRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrodb){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     } //close resistO

                     else {cout << " dstate " << G_demes[d]->drug << " rstate "<< G_demes[d]->resistance << " tstate " << G_demes[d]->transmission <<endl;
                     cout<<"drugeffect: bad resitance state "<<endl;}

             } //close piperaquine

             else if(G_demes[d]->drug == PIPARTE){
                     if(G_demes[d]->resistance==RESISTRA){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                      if(randomEvent <= G_cBradab){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                      }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clradab){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;
                          }
                     } //close resistA

                     else if(G_demes[d]->resistance==RESISTRB){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBrbdab){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRB = G_groups->SUM_BRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrbdab){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } // close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     } //close resistB

                     else if(G_demes[d]->resistance==RESISTRO){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                      if(randomEvent <= G_cBrodab){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRO = G_groups->SUM_BRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                      }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrodab){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     }//close resistO

                     else {
                     cout << " dstate " << G_demes[d]->drug << " rstate "<< G_demes[d]->resistance << " tstate " << G_demes[d]->transmission <<endl;
                     cout<<"drugeffect: bad resitance state "<<endl;}
             } //close piparte


              else if(G_demes[d]->drug == LUMEFANTRINE){
                   if(G_demes[d]->resistance==RESISTRA){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                      if(randomEvent <= G_cBradlf){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                      }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clradlf){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     } //close resistA

                     else if(G_demes[d]->resistance==RESISTRB){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBrbdlf){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRB = G_groups->SUM_BRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrbdlf){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } // close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;}
                     } //close resistB

                     else if(G_demes[d]->resistance==RESISTRO){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                      if(randomEvent <= G_cBrodlf){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRO = G_groups->SUM_BRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                      }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrodlf){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;
                               }
                     }//close resistO

                     else {
                     cout << " dstate " << G_demes[d]->drug << " rstate "<< G_demes[d]->resistance << " tstate " << G_demes[d]->transmission <<endl;
                     cout<<"drugeffect: bad resitance state "<<endl;}
             } //close lumefantrine


             else if(G_demes[d]->drug == AL){
                     if(G_demes[d]->resistance==RESISTRA){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                      if(randomEvent <= G_cBradal){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                      }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clradal){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;
                          }
                     } //close resistA

                     else if(G_demes[d]->resistance==RESISTRB){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_cBrbdal){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRB = G_groups->SUM_BRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrbdal){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } // close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;
                               }
                     } //close resistB

                     else if(G_demes[d]->resistance==RESISTRO){
                          if(G_demes[d]->transmission==BLOODSTAGE){
                                      randomEvent = ran2();
                                      if(randomEvent <= G_cBrodal){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_BRO = G_groups->SUM_BRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                      }
                          } //close blood

                          else if (G_demes[d]->transmission==INFECTIOUS){
                                      randomEvent = ran2();
                                        if(randomEvent <= G_clrodal){
                                                       cc=cc+1;
                                                       G_demes[d]->resistance=RESISTRO;
                                                       G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                                                       G_groups->SUM_S = G_groups->SUM_S + 1;
                                                       G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                                                       G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                                                       G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                                                       G_demes[d]->transmission=SUSCEPTIBLE;
                                                       G_demes[d]->immunity=IMMUNE;
                                                       G_demes[d]->immunity_days=1;
                                                       if (G_demes[d]->clinical==CLINICAL){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       if (G_demes[d]->clinical==ASYMPTOMATIC){
                                                           G_demes[d]->clinical=NOINFECTION;
                                                           G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                                                           G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                                                       }
                                                       continue;
                                        }
                          } //close infect

                          else if (G_demes[d]->transmission==LIVERSTAGE|| G_demes[d]->transmission==SUSCEPTIBLE){
                               continue;
                               }
                     }//close resistO

                     else {
                     cout << " dstate " << G_demes[d]->drug << " rstate "<< G_demes[d]->resistance << " tstate " << G_demes[d]->transmission <<endl;
                     cout<<"drugeffect: bad resitance state "<<endl;}

             } //close coartem

             else {
             cout << " dstate " << G_demes[d]->drug << " rstate "<< G_demes[d]->resistance << " tstate " << G_demes[d]->transmission <<endl;
             cout<<"drug effect: unknown drug "<<endl;}

             if(G_demes[d]->prim == YES){
                   if(G_demes[d]->transmission==INFECTIOUS){
                         if (G_demes[d]->clinical==CLINICAL){
                            randomEvent = ran2();
                            if(randomEvent <= G_clrdpmc){
                               cc=cc+1;
                               G_groups->SUM_S = G_groups->SUM_S + 1;
                               G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                               G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                               G_demes[d]->transmission=SUSCEPTIBLE;
                               G_demes[d]->immunity=IMMUNE;
                               G_demes[d]->immunity_days=1;

                               if(G_demes[d]->resistance==RESISTRO){
                                   G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                                   G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                               }
                               else if(G_demes[d]->resistance==RESISTRA){
                                    G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                                    G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                               }
                               else if(G_demes[d]->resistance==RESISTRB){
                                    G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                                    G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                               }

                               G_demes[d]->clinical=NOINFECTION;
                               G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                               G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                            }
                         }

                         if (G_demes[d]->clinical==ASYMPTOMATIC){
                            randomEvent2 = ran2();
                            if(randomEvent2 <= G_clrdpma){
                               cc=cc+1;
                               G_groups->SUM_S = G_groups->SUM_S + 1;
                               G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;
                               G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
                               G_demes[d]->transmission=SUSCEPTIBLE;
                               G_demes[d]->immunity=IMMUNE;
                               G_demes[d]->immunity_days=1;

                               G_demes[d]->clinical=NOINFECTION;
                               G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                               G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;

                               if(G_demes[d]->resistance==RESISTRO){
                               G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                               G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                               }
                               else if(G_demes[d]->resistance==RESISTRA){
                               G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                               G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                               }
                               else if(G_demes[d]->resistance==RESISTRB){
                               G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                               G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                               }
                            }
                         }
                   }
             }     // close primaquine
         } //close for cycle
         cout << " sumcleared " << cc <<endl;
} //close function



/***********************************/
void treatment(int timestep){ // timestep not used
     vill* village;
     Deme* deme;
     float randomEvent;
     float randomEvent2;
     float randomEvent3;
     float randomEvent4;
     int g,t,k;
     g=t=k=0;

     for (int d=0; d<G_N; d++){
         //t=t+1;
         deme=G_demes[d];
         float G_trtasym=treatasympt[d]*0.0001; // *0.0001 in init?
         float G_trtclin=treatclinical[d]; // this is ->now village
         //cout<<"MDAtreatrate"<<G_tauab<<endl;
         //cout<<"tclin: "<<G_trtclin<<endl;
         //cout<<"tasym: "<<G_trtasym<<endl;

         int region=deme->home;
         int ind=region-1000;
         village=G_villages[ind]; // home village for mda? check ->now village?

         if (deme->mda_times>0) {deme->mda_days=deme->mda_days+1;}

         if (G_switch_mda1 || G_switch_mda2|| G_switch_mda3){
            //cout<<"MDA ON" <<endl;
            if(village->mda==YES){
                //cout<<"REGION:\t"<<ind<<"MDA??:\t"<<village->mda<<"\tTREAT THIS GUY" <<endl;
                //if (deme->drug!=ARTESUNATE && deme->drug!=PIPARTE && deme->drug!=PIPERAQUINE ){
                randomEvent = ran2();
                            if(randomEvent <= G_tauab){
                            // Any need to check for current drug? e.g.check ->mda_days > a value
                            // G_tauab likely to vary between villages?
                                           k=k+1; // k not used
                                           deme->drug = PIPARTE;
                                           deme->prim = YES;
                                           deme->strategy = MDA;
                                           deme->mda_days=0;
                                           deme->mda_times=deme->mda_times+1;
                                           // Is it possible a host can be treated multiple times during one mda visit?
                                           // Should this be avoided?

//                                           G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE + 1;
//                                           G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE + 1;
                                           //continue;
                            }
            }
            else if (deme->mda_days>30 && deme->mda_days<=30+7 &&  deme->mda_times<=3){
                 //if (deme->drug!=ARTESUNATE && deme->drug!=PIPARTE && deme->drug!=PIPERAQUINE ){
                 randomEvent4 = ran2();
                            if(randomEvent4 <= G_tauab){
                                           k=k+1;
                                           deme->drug = PIPARTE;
                                           deme->prim = YES;
                                           deme->strategy = MDA;
                                           deme->mda_days = 0;
                                           deme->mda_times=deme->mda_times+1;

//                                           G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE + 1;
//                                           G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE + 1;
                                           //continue;
                            }
                        //}
            }
         } //close mda

        if (deme->drug!=PIPERAQUINE && deme->drug!=PIPARTE && deme->drug!=ARTESUNATE ){
        // what about other drugs? should this be == NODRUG? -yes
        // Is this checking if there is any drug in the system, or new drugs only i.e. mda-ed drugs, drug_days = 0?
            if (deme->clinical==CLINICAL){
               randomEvent2 = ran2();
               if(randomEvent2 <= G_trtclin){
                               deme->drug = PIPARTE;
                               deme->prim = YES;
                               deme->strategy = ACTSWITCH;
//                               G_groups->SUM_AL = G_groups->SUM_AL + 1;
//                               G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL + 1;
                               //continue;
               }
            }
            else if (deme->clinical==ASYMPTOMATIC){
               randomEvent3 = ran2();
               if(randomEvent3 <= G_trtasym){
                               deme->drug = PIPARTE;
                               deme->prim = YES;
                               deme->strategy = ACTSWITCH;
//                               G_groups->SUM_AL = G_groups->SUM_AL + 1;
//                               G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL + 1;
                               //continue;
               }
            }
         }
     }
     //cout<<"t:\t"<<t<<"\tg\t:"<<g<<"\tk:\t"<<k<<endl;   // close for cycle
} //close treat function

void surveytimes(){
//     for (int i=G_delay; i<G_tmax; i+=365){
//         if (i>=G_first){
//         surveys.push_back(i);
//         cout<<"time: "<<i;
//         }
//     }

     surveys.push_back(1860);
     surveys.push_back(2065);
     G_dimsurveys=surveys.size();
     cout<<" dimvector: "<<G_dimsurveys<<endl;
}

//turns treatment on and off
void triggerMDA(int timestep){
     if (timestep>=G_mda_start1){

        int d;
        vill* village;
        Deme* deme;
        int k = 0;

        for (d=0;d<G_xregion;d++){
           village= G_villages[d];
           if (village->mda_times>0){
                village->mda_days=village->mda_days+1;
           }
        }

        for (int w=0;w<G_dimsurveys;w++){
            int timeref=surveys[w];
            if (timestep == timeref){
            //cout<<"1...2...3...BEGIN TESTING"<<endl;
            for (d=0;d<G_xregion;d++){
                village= G_villages[d];
                if (village->mda==NO){
                     int POP;
                     int RO;
                     int RA;
                     int RB;
                     float propinf;
                     int SAMPLES = 50;
                     float detecprev;

                     POP = G_villdata[d]->popsize;
                     RO=G_villdata[d]->SUM_IRO;
                     RA=G_villdata[d]->SUM_IRA;
                     RB=G_villdata[d]->SUM_IRB;
                     int inf=RO+RA+RB;
                     propinf=float(inf)/float(POP);
                     int detectcases = int(poissonHighPrecision(propinf*SAMPLES));
                     int seennumber = int((detectcases*G_sensitivity));
                     detecprev = float(seennumber)/float(SAMPLES);

                     if (detecprev >= G_trigger){
                        village->delay=1;
                        village->mda=YES;
                        if (village->post==0){
                            village->post = 1;
                            village->mptime = timestep;
                        }
                        village->mda_days = 0;
                        k=k+1;
                        triggs.push_back(d);
                     }
               }
            }
            cout<< "size: "<<triggs.size()<< "\tkcount: " << k<<endl;
            random_shuffle(triggs.begin(),triggs.end());
            }
        }

        int sz = triggs.size();
        for (int r=0; r<sz; r++){
             if (triggs.size()==0){continue;}
             int vv = triggs[r];
             //cout<< "which: "<< vv <<endl;
             village=G_villages[vv];
             if (village->mda_days > G_dayspervill){
                //int g =getchar();
                triggs.erase(triggs.begin()+r);
                r-=1;
                village->mda=NO;
                village->mda_days=0;
                village->mda_times=0;
                village->delay=0;
             }
             //cout<< "size: "<<triggs.size()<<endl;
        }

        if (triggs.size()>=G_mteams){
            for (int h=0; h<G_mteams; h++){
                int indv = triggs[h];
                village= G_villages[indv];
                village->mda=YES;
            }
        }
        else{
             for (int hh=0; hh<triggs.size(); hh++){
                int indv2 = triggs[hh];
                village= G_villages[indv2];
                village->mda=YES;
            }
        }
     }
}


void treattrigger(int timestep){
     vill* village;
     Deme* deme;
     float randomEvent;
     float randomEvent2;
     float randomEvent3;
     float randomEvent4;
     int g,t,k,kk;
     g=t=k=kk=0;

     for (int d=0; d<G_N; d++){
         t=t+1;
         deme=G_demes[d];
         float G_trtasym=treatasympt[d]*0.001; // 0.0001?
         float G_trtclin=treatclinical[d];

         int region=deme->now;
         int ind=region-1000;
         village=G_villages[ind];

         if (village->mda==YES){
         kk=kk+1;
         village->delay=village->delay+1;
               if (village->delay > 30 ){
                   village->mda = ACTIVE;
                   village->mda_times = 1;
                   village->mda_start = timestep;
                   deme->mda_times=0;
                   deme->mda_days=0;
               }
         }

         if (deme->mda_times>0) {deme->mda_days=deme->mda_days+1;}

         if(village->mda==ACTIVE){
                //if (deme->drug!=ARTESUNATE && deme->drug!=PIPARTE && deme->drug!=PIPERAQUINE ){
                g=g+1;
                randomEvent = ran2();
                            if(randomEvent <= G_tauab){
                                           k=k+1;
                                           if (deme->drug==ARTESUNATE) {
                                              G_groups->SUM_ART = G_groups->SUM_ART - 1;
                                              G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART - 1;
                                              }
                                           else if (deme->drug==PIPARTE) {
                                                G_groups->SUM_PIPARTE= G_groups->SUM_PIPARTE - 1;
                                                G_villdata[ind]->SUM_PIPARTE= G_villdata[ind]->SUM_PIPARTE - 1;
                                                }
                                           else if (deme->drug==PIPERAQUINE) {
                                                G_groups->SUM_PIP= G_groups->SUM_PIP - 1;
                                                G_villdata[ind]->SUM_PIP= G_villdata[ind]->SUM_PIP - 1;
                                                }

                                           deme->drug = PIPARTE;
                                           deme->strategy = MDA;
                                           deme->prim = YES;
                                           deme->mda_days=1;
                                           deme->mda_times=deme->mda_times+1;

                                           G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE + 1;
                                           G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE + 1;
                                           continue;
                            }
                        //}
            }
            else if (deme->mda_days>30 && deme->mda_days<=30+7 &&  deme->mda_times<3){
                 //if (deme->drug!=ARTESUNATE && deme->drug!=PIPARTE && deme->drug!=PIPERAQUINE ){
                 //cout<<"REGION:\t"<<ind<<"MDA??:\t"<<village->mda<<"\tTREAT THIS GUY" <<endl;
                 randomEvent4 = ran2();
                            if(randomEvent <= G_tauab){
                                           k=k+1;
                                           if (deme->drug==ARTESUNATE) {
                                              G_groups->SUM_ART = G_groups->SUM_ART - 1;
                                              G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART - 1;
                                              }
                                           else if (deme->drug==PIPARTE) {
                                                G_groups->SUM_PIPARTE= G_groups->SUM_PIPARTE - 1;
                                                G_villdata[ind]->SUM_PIPARTE= G_villdata[ind]->SUM_PIPARTE - 1;
                                                }
                                           else if (deme->drug==PIPERAQUINE) {
                                                G_groups->SUM_PIP= G_groups->SUM_PIP - 1;
                                                G_villdata[ind]->SUM_PIP= G_villdata[ind]->SUM_PIP - 1;
                                                }

                                           deme->drug = PIPARTE;
                                           deme->prim = YES;
                                           deme->strategy = MDA;
                                           deme->mda_days = 1;
                                           deme->mda_times=deme->mda_times+1;

                                           G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE + 1;
                                           G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE + 1;
                                           continue;
                            }
                        //}
            }


            if (deme->strategy!=MDA && deme->drug!=ARTESUNATE && deme->drug!=AL && deme->drug!=LUMEFANTRINE ){
            // only place where deme->strategy is used
                if (deme->clinical==CLINICAL){
                   randomEvent2 = ran2();
                   if(randomEvent2 <= G_trtclin){
                                   //cout<<"ttttttttt"<<G_trtclin<<endl;
                                   deme->drug = PIPARTE;
                                   deme->prim = YES;
                                   deme->strategy = ACTSWITCH;
                                   G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE + 1;
                                   G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE + 1;
                                   continue;
                   }
                }
                else if (deme->clinical==ASYMPTOMATIC){
                   randomEvent3 = ran2();
                   if(randomEvent3 <= G_trtasym){
                                   //cout<<"gggggggg"<<G_trtasym<<endl;
                                   deme->drug = PIPARTE;
                                   deme->prim = YES;
                                   deme->strategy = ACTSWITCH;
                                   G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE + 1;
                                   G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE + 1;
                                   continue;
                   }
                }
            }
        }
     //cout<<"t:\t"<<t<<"\tg:\t"<<g<<"\tk:\t"<<k<<"\tkk:\t"<<kk<<endl;   // close for cycle
} //close treat function

void inline squarepulse(int timestep){
     //squarepulse MDA
     if(timestep >= G_mda_start1 && timestep <= G_mda_end1)            { G_switch_mda1 = true; }
     if(timestep > G_mda_end1)                                         { G_switch_mda1 = false;}

     if(timestep >= G_mda_start2 && timestep <= G_mda_end2)            { G_switch_mda2 = true; }
     if(timestep > G_mda_end2)                                         { G_switch_mda2 = false;}

     if(timestep >= G_mda_start3 && timestep <= G_mda_end3)            { G_switch_mda3 = true; }
     if(timestep > G_mda_end3)                                         { G_switch_mda3 = false;}
}

void population_movement(int timestep){
    vill* village;
    Deme* deme;
    float migrate2 = 1.0/90.0;
    int region, ind;
    float randomEvent;
    float randomEvent2;
    float randomEvent3;
    float randomEvent4;
    float t,k;
    int g,sz,n,d,gg,g2,g3,g4;


    if (timestep % 120 == 0){
        G_migrate = 5.0/10.0;
    }
    if (timestep % 270 == 0){
        G_migrate = 0.0;
        for (d=0; d<G_N; d++){
            deme=G_demes[d];
            if (deme->home != deme-> oldhome){
                n = (deme->home)-1000;
                g = (deme->oldhome)-1000;
                deme->home = deme->oldhome;
                deme->now = deme->home;
            }
        }
    }

    for (d=0; d<G_N; d++){
        deme=G_demes[d];
        region=deme->now;
        ind=region-1000;
        village=G_villages[ind];

        /** temporary migration **/
        if (deme->mobile==TEMP){
            randomEvent = ran2();
            if (randomEvent<G_migrate){
                gg = village->seasonalnet[int(ran2()*village->seasonalnet.size())];
                //gg = int(ran2()*G_xregion);
                deme->oldhome=deme->home;
                deme->home=gg+1000;
                deme->now=deme->home;
            }
        }

        /** mobile population migration **/
        if (deme->mobile==MOBILE){
            randomEvent = ran2();
            if (randomEvent<G_migrate_mobile){
                // move them to a random place
                sz=village->migrantnet.size();
//                g2 = int(ran2()*G_xregion);
//                deme->home=g2+1000;
//                deme->now=deme->home;
                if (sz>0){
                    g2 = village->migrantnet[int(ran2()*village->migrantnet.size())];
                    deme->home=g2+1000;
                    deme->now=deme->home;
                }
            }
        }

        /** general population short term movement **/
        if (deme->mobile==STATIC){
            if (deme->now==deme->home ){
                randomEvent = ran2();
                if (randomEvent<G_mobility_static){
                    //send people to a random village

//                    g3 = int(ran2()*G_xregion);
//                    deme->now=g3+1000;
                    //int szz = village->staticnet.size();
                    //cout<<"SIZE NET"<<szz<<endl;
                    g3 = village->staticnet[int(ran2()*village->staticnet.size())];
                    deme->now=g3+1000;
                }
            }
            else{
                randomEvent2 = ran2();
                if (randomEvent2<G_static_return){
                    //return people to their home/starting village
                    deme->now=deme->home;
                }
            }
        }
        else {
            if (deme->now==deme->home ){
                randomEvent3 = ran2();
                if (randomEvent3<G_mobility_mobile){
                    //send people to a random village
//                    g4 = int(ran2()*G_xregion);
//                    deme->now=g4+1000;

                    g4 = village->staticnet[int(ran2()*village->staticnet.size())];
                    deme->now=g4+1000;
                }
            }
            else {
                randomEvent4 = ran2();
                if (randomEvent4<G_mobile_return){
                    //return people to their home/starting village
                    deme->now=deme->home;
                }
            }
        }
    }
}

void temp_migration(int timestep){
    vill* village;
    Deme* deme;
    float randomEvent;

    if (timestep % 120 == 0){
        G_migrate = 1.0/10.0;
    }
    if (timestep % 270 == 0){
        G_migrate = 0.0;
        for (int d=0; d<G_N; d++){
            deme=G_demes[d];
            if (deme->home != deme-> oldhome){
                int n = (deme->home)-1000;
                int g = (deme->oldhome)-1000;

                deme->home = deme->oldhome;
                deme->now = deme->home;
            }
        }
    }
    for (int d=0; d<G_N; d++){
        deme=G_demes[d];
        if (deme->mobile==TEMP){
            randomEvent = ran2();
            if (randomEvent<G_migrate){
                int region=deme->now;
                int ind=region-1000;
                village=G_villages[ind];
                //int g = village->seasonalnet[int(ran2()*village->seasonalnet.size())];
                int g = int(ran2()*G_xregion);
                deme->oldhome=deme->home;
                int n = (deme->now)-1000;
//
                deme->home=g+1000;
                deme->now=deme->home;
            }
        }
    }
}



void migration(int timestep){
    vill* village;
    Deme* deme;
    float randomEvent;
    float migrate2 = 1.0/90.0;
    cout<<"here"<<endl;
    for (int d=0; d<G_N; d++){
        deme=G_demes[d];
        if (deme->mobile==MOBILE){
            randomEvent = ran2();
            if (randomEvent<G_migrate_mobile){
                // move them to a random place
                int region=deme->now;
                int ind=region-1000;
                village=G_villages[ind];
                int sz=village->migrantnet.size();
                int g = int(ran2()*G_xregion);
                deme->home=g+1000;
                deme->now=deme->home;
//                if (sz>0){
//                    int g = village->migrantnet[int(ran2()*village->migrantnet.size())];
//                    int n = (deme->now)-1000;
//                    deme->home=g+1000;
//                    deme->now=deme->home;
//                }
            }
        }
    }
}


void mobility(){
    vill* village;
    Deme* deme;
    float randomEvent;
    float randomEvent2;
    float randomEvent3;
    float randomEvent4;
    float g,t,k;

    for (int d=0; d<G_N; d++){
        deme=G_demes[d];
        if (deme->mobile==STATIC){
            if (deme->now==deme->home ){
                randomEvent = ran2();
                if (randomEvent<G_mobility_static){
                    //send people to a random village
                    int g = int(ran2()*G_xregion);
                    deme->now=g+1000;
                }
            }
            else{
                randomEvent2 = ran2();
                if (randomEvent2<G_static_return){
                    //return people to their home/starting village
                    deme->now=deme->home;
                }
            }
        }
        else {
            if (deme->now==deme->home ){
                randomEvent3 = ran2();
                if (randomEvent<G_mobility_mobile){
                    //send people to a random village
                    int g = int(ran2()*G_xregion);
                    deme->now=g+1000;
                }
            }
            else {
                randomEvent4 = ran2();
                if (randomEvent2<G_mobile_return){
                    //return people to their home/starting village
                    deme->now=deme->home;
                }
            }
        }
    }
}


void village_here_now(){
    vill* village;
    Deme* deme;

    float mpp, doff;
    for (int i=0; i<G_xregion;i++){
        village=G_villages[i];
        village->herenow.clear();

//        mpp=ran2();
//        if(mpp<G_poormp){
//            village->postfail=1;
//        }
//        if (village->postfail==1){
//            doff=ran2();
//            if (doff<G_daysoff){
//                village->postfail=0;
//            }
//        }
    }

    for (int d=0; d<G_N; d++){
        deme=G_demes[d];
        int region=deme->now;
        int ind=region-1000;
        village=G_villages[ind];
        village->herenow.push_back(d);
    }

}

void mosquitodynamics(){
    float randomEvent,t1,t2;
    int db, newbites, dba, db2, i, j;
    mosq* mosquito;
    vill* village;
    Deme* deme_bite_past;
    float lambdaPoisson;
    int aux22;
    int MEXP=0;
    int MINF =0;

    for (i=0; i<G_xregion; i++){
        mosquito = G_mosquitoes[i];
        village = G_villages[i];
        lambdaPoisson = (betaseason[i]*village->herenow.size());
        newbites = int(poissonHighPrecision(lambdaPoisson));
        //cout<<"bites: "<<newbites<<endl;

        /** infect mosquitoes **/
        for ( j=0; j<newbites; j++){
            int dba = int(ran2()* village->herenow.size());
            db=village->herenow.at(dba);
//            int dba = int(ran2()* (G_N-1));
//            db=dba;
            deme_bite_past = G_demes[db];
            if (deme_bite_past->infections.empty()==false){
                randomEvent = ran2();
                if(randomEvent <= G_c* deme_bite_past->infectiousness){
//                    db2 = int(ran2()* deme_bite_past->infections.size());
//                    mosquito->resistance.push_back(deme_bite_past->infections.at(db2));
                    mosquito->infected += 1;
                    mosquito->resistance.push_back(deme_bite_past->dominant);
                    //cout<<"infect mosq dominant: "<<deme_bite_past->dominant<<endl;
                }
            }
        }


        /**  mosquito survival  **/
        if (mosquito->resistance.empty()==false){
            vector<int>::iterator it1 = mosquito->resistance.begin();
            for ( ; it1 != mosquito->resistance.end(); ){
                t2=ran2();
                if (t2<G_mosqdeath2){
                    it1 = mosquito->resistance.erase(it1);
                    mosquito->infected-=1;
                } else {
                    ++it1;
                }
            }
        }
        if (mosquito->resistancei.empty()==false){
            vector<int>::iterator it2 = mosquito->resistancei.begin();
            for ( ; it2 != mosquito->resistancei.end(); ) {
                t1=ran2();
                if (t1<G_mosqdeath1){
                    it2 = mosquito->resistancei.erase(it2);
                    mosquito->infectious-=1;
                } else {
                    ++it2;
                }
            }
        }

        /** mosquito incubation **/
        int jj=0;
        //cout<<"infected " << aux2 << "vector infected " << mosquito->resistance.size() <<endl;
        if (mosquito->resistance.empty()==false){
            vector<int>::iterator it = mosquito->resistance.begin();
            for ( ; it != mosquito->resistance.end(); ) {
                jj+=1;
                t1=ran2();
                aux22 = mosquito->resistance.at(jj-1);
                if (t1<G_mosqinc){
                    it = mosquito->resistance.erase(it);
                    mosquito->resistancei.push_back(aux22);
                    mosquito->infected-=1;
                    mosquito->infectious+=1;
                    jj-=1;
                } else {
                    ++it;
                }
            }
        }

        //cout<< "nbites"<< mosquito->infected << "resveclength" << mosquito->resistance.size() <<endl;
        /** print mosqs **/
        MEXP=MEXP+mosquito->infected;
        MINF=MINF+mosquito->infectious;
    }
    //cout << " MEXP: " << MEXP << " MINF: "<< MINF << endl;
}

void infectmosquitoes(){
    float randomEvent;
    int db;
    mosq* mosquito;
    vill* village;
    Deme* deme_bite_past;

    for ( int i=0;i<G_xregion;i++){
        mosquito=G_mosquitoes[i];
        village=G_villages[i];
        float lambdaPoisson = (betaseason[i]*village->herenow.size());
        int newbites = int(poissonHighPrecision(lambdaPoisson));

        for ( int j=0; j<newbites; j++){
            //int dba = int(ran2()* village->herenow.size());
            //db=village->herenow.at(dba);
            int dba = int(ran2()* G_N);
            db=dba;
            deme_bite_past = G_demes[db];
            if (deme_bite_past->infections.empty()==false){
                randomEvent = ran2();
                if(randomEvent <= G_c* deme_bite_past->infectiousness){
                    int db2 = int(ran2()* deme_bite_past->infections.size());
                    mosquito->infected += 1;
                    mosquito->resistance.push_back(deme_bite_past->infections.at(db2));
                }
            }
        }
        //cout<< "nbites"<< mosquito->infected << "resveclength" << mosquito->resistance.size() <<endl;
    }
}

void survivemosqs(){
    mosq* mosquito;
    int nbites;
    int aux1,aux2;
    float t1,t2;

    for (int i=0;i<G_xregion;i++){
        mosquito=G_mosquitoes[i];
        aux1 = mosquito->infectious;
        aux2 = mosquito->infected;

        if (mosquito->resistance.empty()==false){
            vector<int>::iterator it1 = mosquito->resistance.begin();
            for ( ; it1 != mosquito->resistance.end(); ) {
                t2=ran2();
                if (t2<G_mosqdeath2){
                    it1 = mosquito->resistance.erase(it1);
                    mosquito->infected-=1;
                } else {
                    ++it1;
                }
            }
        }
        if (mosquito->resistancei.empty()==false){
            vector<int>::iterator it2 = mosquito->resistancei.begin();
            for ( ; it2 != mosquito->resistancei.end(); ) {
                t1=ran2();
                if (t1<G_mosqdeath1){
                    it2 = mosquito->resistancei.erase(it2);
                    mosquito->infectious-=1;
                } else {
                    ++it2;
                }
            }
        }
     }
}

void mosqincubation(){
    mosq* mosquito;
    int aux,aux2;
    float t1;

    for (int i=0;i<G_xregion;i++){
        mosquito=G_mosquitoes[i];
        aux = mosquito->infected;
        int j=0;
        //cout<<"infected " << aux << "vector infected " << mosquito->resistance.size() <<endl;
        if (mosquito->resistance.empty()==false){
            vector<int>::iterator it = mosquito->resistance.begin();
            for ( ; it != mosquito->resistance.end(); ) {
                j+=1;
                t1=ran2();
                aux2 = mosquito->resistance.at(j-1);
                if (t1<G_mosqinc){
                    it = mosquito->resistance.erase(it);
                    mosquito->resistancei.push_back(aux2);
                    mosquito->infected-=1;
                    mosquito->infectious+=1;
                    j-=1;
                } else {
                    ++it;
                }
            }
        }
    }
}

void infection(){
//infection mechanism; number of bites determined by force of infection;
//state of deme bitten in the past determines what will happen to individuals bitten now
        //cout<<"CHEGUEI!!!!"<<endl;
        //int c = getchar();
    int i=0;
    int newbites;
    float randomEvent;
    int db, dba, dbm;
    int dl;
    vill* village;
    mosq* mosquito;
    //cout << "lamb:" << lambdaPoisson << " nb:" << newbites << " I%:" << G_groups->getI()/float(G_N) << endl;
    Deme* deme_bite;
    Deme* deme_bite_in_past;
    Deme* deme2;

    for ( i=0;i<G_xregion;i++){
        mosquito=G_mosquitoes[i];
        village = G_villages[i];
        //cout<<""<<i<<endl;
        float lambdaPoisson =(betaseason[i]*village->herenow.size());
        newbites = int(poissonHighPrecision(lambdaPoisson));
        //newbites = mosquito->infectious;
        //cout<< "infectious mosq" << mosquito->infectious <<  "mosq length" << mosquito->resistancei.size() <<endl;
        //cout<<"nb: "<<newbites<<endl;
        //int sz = G_villages[i]->herenow.size();
        //cout<< "size here now: " << sz <<endl;

        for(int ninfec=0; ninfec<newbites; ninfec++ ){
            dbm = int(ran2()* mosquito->resistancei.size());
            if (dbm>0){
            int res = mosquito->resistancei[dbm];

            dba = int(ran2()* G_villages[i]->herenow.size());
            db = G_villages[i]->herenow.at(dba);
            deme_bite = G_demes[db];
            //cout << "ninf" << ninfec << "res" << res <<endl;
             //cout<<"ff: "<<dl<<"\ti: "<<db<<"\tpop"<<i<<endl;
             //cout<<"ff: "<<deme_bite_in_past->transmission<<"\ti: "<<deme_bite->transmission<<endl;

            if (deme_bite->transmission == SUSCEPTIBLE ){

                randomEvent = ran2();
                if(randomEvent <= G_c *deme_bite->susceptibility){

                    deme_bite->transmission = LIVERSTAGE;
                    G_groups->SUM_S = G_groups->SUM_S - 1;
                    G_villdata[i]->SUM_S = G_villdata[i]->SUM_S - 1;
                    deme_bite->cummulative_exposures = deme_bite->cummulative_exposures+1;

                    if(res == RESISTRO ){
                        deme_bite->blood_qeue.push_back(RESISTRO);
//                        G_groups->SUM_LRO = G_groups->SUM_LRO + 1;
//                        G_villdata[i]->SUM_LRO = G_villdata[i]->SUM_LRO + 1;
                    }

                    else if(res == RESISTRA){
                        deme_bite->blood_qeue.push_back(RESISTRA);
//                        G_groups->SUM_LRA = G_groups->SUM_LRA + 1;
//                        G_villdata[i]->SUM_LRA = G_villdata[i]->SUM_LRA + 1;
                    }

                    else if(res == RESISTRB){
                        deme_bite->blood_qeue.push_back(RESISTRB);
//                        G_groups->SUM_LRB = G_groups->SUM_LRB + 1;
//                        G_villdata[i]->SUM_LRB = G_villdata[i]->SUM_LRB + 1;
                    }

                    else {cout << "infection: couldnt find resistance" << endl;
                    int c = getchar();
                    }

                } //close randomevent if

            } //close check infectious and susceptible
            else if (deme_bite->transmission == LIVERSTAGE || deme_bite->transmission == BLOODSTAGE){

                randomEvent = ran2();
                if(randomEvent <= deme_bite->susceptibility){
                    deme_bite->super = 1;
                    deme_bite->cummulative_exposures = deme_bite->cummulative_exposures+1;

                    if(res == RESISTRO ){
                        deme_bite->blood_qeue.push_back(RESISTRO);
                    }

                    else if(res == RESISTRA){
                        deme_bite->blood_qeue.push_back(RESISTRA);
                    }

                    else if(res == RESISTRB){
                        deme_bite->blood_qeue.push_back(RESISTRB);
                    }

                }   //close randomevent if
            }   //close check infectious and susceptible
        }
        }   //close for bites
    } //close for region
} //close function


void superdrugeffect(){
    Deme* deme;
    int aux2,aux3;
    float t1,t2;
    int cc = 0;

    for (int d=0; d<G_N; d++){
        deme=G_demes[d];
        int j=0;
        int jj=0;

        if (G_demes[d]->drug == NODRUG){
            continue;
        }
        else{
            if (G_demes[d]->infectious_qeue.empty()==false){
                //cout<<"qeue"<<deme->infectious_qeue.size()<<endl;
                vector<int>::iterator it = G_demes[d]->infectious_qeue.begin();
                for ( ; it != G_demes[d]->infectious_qeue.end(); ) {
                    j+=1;
                    aux2 = G_demes[d]->infectious_qeue.at(j-1);
                    if(G_demes[d]->drug == ARTESUNATE){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBrada){
                                    //cout<<"CLEAR 1"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbda){
                                    //cout<<"CLEAR 2"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBroda){
                                    //cout<<"CLEAR 3"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == PIPERAQUINE){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradb){
                                    //cout<<"CLEAR 4"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdb){
                                    //cout<<"CLEAR 5"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodb){
                                    //cout<<"CLEAR 6"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == PIPARTE){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradab){
                                    //cout<<"CLEAR 7"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdab){
                                    //cout<<"CLEAR 8"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodab){
                                    //cout<<"CLEAR 9"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == LUMEFANTRINE){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradlf){
                                    //cout<<"CLEAR 10"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdlf){
                                    //cout<<"CLEAR 11"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodlf){
                                    //cout<<"CLEAR 12"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == AL){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradal){
                                    //cout<<"CLEAR 13"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdal){
                                    //cout<<"CLEAR 14"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodal){
                                    //cout<<"CLEAR 15"<<endl;
                                    it=G_demes[d]->infectious_qeue.erase(it);
                                    j-=1;
                                    continue;
                                }
                                else{++it;continue;}
                            }
                    }
                    else{++it;}
                }
                if (G_demes[d]->infectious_qeue.empty()==true && G_demes[d]->infections.empty()==true){
                     //cout<<"RECOVERED 1"<<endl;
                     cc=cc+1;
                     G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                     G_groups->SUM_S = G_groups->SUM_S + 1;
                     G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                     G_demes[d]->transmission=SUSCEPTIBLE;
                     G_demes[d]->immunity=IMMUNE;
                     G_demes[d]->immunity_days=1;
                     G_demes[d]->immunity_level=G_demes[d]->immunity_level+1;
                     G_demes[d]->moi=0;

                     if (G_demes[d]->clinical==CLINICAL){
                        G_demes[d]->clinical=NOINFECTION;
                        G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                        G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                     }
                     if (G_demes[d]->clinical==ASYMPTOMATIC){
                        G_demes[d]->clinical=NOINFECTION;
                        G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                        G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                     }
                }
            } //close bloostage

            if (G_demes[d]->infections.empty()==false){
                //cout<<"infections"<<deme->infections.size()<<endl;
                vector<int>::iterator it2 = G_demes[d]->infections.begin();
                for ( ; it2 != G_demes[d]->infections.end(); ) {
                    jj+=1;
                    aux3 = G_demes[d]->infections.at(jj-1);
                    //cout<<"what value"<<aux3<<endl;
                    //cout<<"what drug"<<deme->drug<<endl;
                    if(G_demes[d]->drug == ARTESUNATE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clrada){
                                    //cout<<"BLEAR 1"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbda){
                                    //cout<<"BLEAR 2"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clroda){
                                    //cout<<"BLEAR 3"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == PIPERAQUINE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradb){
                                    //cout<<"BLEAR 4"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdb){
                                    //cout<<"BLEAR 5"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodb){
                                    //cout<<"BLEAR 6"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == PIPARTE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradab){
                                    //cout<<"BLEAR 7"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdab){
                                    //cout<<"BLEAR 8"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodab){
                                    //cout<<"BLEAR 9"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == LUMEFANTRINE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradlf){
                                    //cout<<"BLEAR 10"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdlf){
                                    //cout<<"BLEAR 11"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodlf){
                                    //cout<<"BLEAR 12"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                    }
                    else if (G_demes[d]->drug == AL){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradal){
                                    //cout<<"BLEAR 13"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdal){
                                    //cout<<"BLEAR 14"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodal){
                                    //cout<<"BLEAR 15"<<endl;
                                    it2=G_demes[d]->infections.erase(it2);
                                    jj-=1;
                                    continue;
                                }
                                else{++it2;continue;}
                            }
                            ++it2;
                    }
                    else{++it2;}
                }

                if (G_demes[d]->infectious_qeue.empty()==true && G_demes[d]->infections.empty()==true){
                     //cout<<"RECOVERED 2"<<endl;
                     cc=cc+1;
                     G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                     G_groups->SUM_S = G_groups->SUM_S + 1;
                     G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                     G_demes[d]->transmission=SUSCEPTIBLE;
                     G_demes[d]->immunity=IMMUNE;
                     G_demes[d]->immunity_days=1;
                     G_demes[d]->immunity_level=G_demes[d]->immunity_level+1;
                     G_demes[d]->moi=0;

                     if (G_demes[d]->clinical==CLINICAL){
                        G_demes[d]->clinical=NOINFECTION;
                        G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
                        G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                     }
                     if (G_demes[d]->clinical==ASYMPTOMATIC){
                        G_demes[d]->clinical=NOINFECTION;
                        G_groups->SUM_ASYM = G_groups->SUM_ASYM - 1;
                        G_groups->SUM_NOINF = G_groups->SUM_NOINF + 1;
                     }
                }
            } // close infectious
        }
    } // close for individuals
    //cout<<"HOW MANY RECOVERED FROM DRUG ACTION: "<<cc<<endl;
}


void superinfectiontoblood(){
    Deme* deme;
    int aux2;
    float t1;

    for (int i=0;i<G_N;i++){
        deme=G_demes[i];
        int j=0;
        if (deme->blood_qeue.empty()==false){
            vector<int>::iterator it = deme->blood_qeue.begin();
            for ( ; it != deme->blood_qeue.end(); ) {
                j+=1;
                t1=ran2();
                aux2 = deme->blood_qeue.at(j-1);
                if (t1<G_gamma){
                    deme->transmission=BLOODSTAGE;
                    it = deme->blood_qeue.erase(it);
                    deme->infectious_qeue.push_back(aux2);
                    j-=1;
                } else {
                    ++it;
                }
            }
        }
     }
}

void superinfectiontoinfectious(){
    Deme* deme;
    int aux2;
    float t1,pci;

    for (int i=0;i<G_N;i++){
        deme=G_demes[i];
        int j=0;
        if (deme->infectious_qeue.empty()==false){
            vector<int>::iterator it = deme->infectious_qeue.begin();
            for ( ; it != deme->infectious_qeue.end(); ) {
                j+=1;
                t1=ran2();
                aux2 = deme->infectious_qeue.at(j-1);
                if (t1<G_sigma){
                    deme->transmission=INFECTIOUS;
                    deme->infections.push_back(aux2);
                    if (deme->moi>1){
                        pci=ran2();
                        if (pci<= deme->probclinicalbb){
                            deme->clinical = CLINICAL;
                        }
                    }
                    else if (deme->moi == 1){
                        pci=ran2();
                        if (pci<= deme->probclinicallb){
                            deme->clinical = CLINICAL;
                        }
                        else{deme->clinical = ASYMPTOMATIC;}
                    }

                    it = deme->infectious_qeue.erase(it);
                    j-=1;
                } else {
                    ++it;
                }
            }
        }
    }
}

void superrecovery(){
    Deme* deme;
    int aux2,j;
    float t1,pci,delta;
    int rec=0;

    for (int d=0;d<G_N;d++){
        //deme=G_demes[i];
        int j=0;
        if (G_demes[d]->infections.empty()==false){
            vector<int>::iterator it = G_demes[d]->infections.begin();
            for ( ; it != G_demes[d]->infections.end(); ) {
                j+=1;
                delta=ran2();
                if (delta<=G_delta){
                    it = G_demes[d]->infections.erase(it);
                    G_demes[d]->moi-=1;
                    j-=1;
                }
                else{++it;}
            }
            if (G_demes[d]->infections.empty()==true){
                G_demes[d]->moi=0;
                rec+=1;
                G_demes[d]->transmission=SUSCEPTIBLE;
                G_demes[d]->immunity_level= G_demes[d]->immunity_level+1;
                G_demes[d]->immunity_days=1;
                G_groups->SUM_IMM = G_groups->SUM_IMM + 1;
                G_groups->SUM_NIMM = G_groups->SUM_NIMM - 1;
                G_demes[d]->resistance=RESISTRO;
                if (G_demes[d]->clinical==CLINICAL){
                    G_groups->SUM_CLIN = G_groups->SUM_CLIN -1 ;
                    G_groups->SUM_NOINF = G_groups->SUM_NOINF +1 ;
                    G_demes[d]->clinical = NOINFECTION;
                 }
                 else if (G_demes[d]->clinical==ASYMPTOMATIC){
                    G_groups->SUM_ASYM = G_groups->SUM_ASYM -1 ;
                    G_groups->SUM_NOINF = G_groups->SUM_NOINF +1 ;
                    G_demes[d]->clinical = NOINFECTION;
                 }
            }
        }
    }
    //cout<<"natural recovery:"<<rec<<endl;
}

void superbirth(){
    Deme* deme;
    int d;
    float mu;
    float randomEvent;
    float randomEvent2;

    for (d=0;d<G_N;d++){
         deme=G_demes[d];
         mu=ran2();
         if (mu<= deme->death){
            int region=deme->home;
            int ind=region-1000;

            if(G_demes[d]->drug == ARTESUNATE){
                G_groups->SUM_ART = G_groups->SUM_ART - 1;
                G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART - 1;
            }
            else if(G_demes[d]->drug == PIPERAQUINE){
                G_groups->SUM_PIP = G_groups->SUM_PIP - 1;
                G_villdata[ind]->SUM_PIP = G_villdata[ind]->SUM_PIP - 1;
            }
            else if(G_demes[d]->drug == PIPARTE){
                G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE - 1;
                G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE - 1;
            }
            else if(G_demes[d]->drug == AL){
                G_groups->SUM_AL = G_groups->SUM_AL - 1;
                G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL - 1;
            }
            else if(G_demes[d]->drug == LUMEFANTRINE){
                G_groups->SUM_LUM = G_groups->SUM_LUM - 1;
                G_villdata[ind]->SUM_LUM = G_villdata[ind]->SUM_LUM - 1;
            }

            G_demes[d]->transmission = SUSCEPTIBLE;
            G_demes[d]->resistance = RESISTRO;
            G_demes[d]->drug = NODRUG;
            G_demes[d]->prim = NO;
            G_demes[d]->strategy = NOSTRAT;
            G_demes[d]->age = 0;
            G_demes[d]->immunity_level = 0;
            G_demes[d]->moi = 0;
            G_demes[d]->cummulative_exposures = 0;

            randomEvent = ran2();
            if (randomEvent <= 0.482){ //in accordance with 1998 census
                G_demes[d]->gender = MALE;
            }
            else{
                G_demes[d]->gender = FEMALE;
            }

            randomEvent2 = ran2();
            if (randomEvent2 <= G_itnprop){
               if ( G_demes[d]->itn==NOITN){
                  G_groups->SUM_ITN = G_groups->SUM_ITN + 1;
                  G_groups->SUM_NOITN = G_groups->SUM_NOITN - 1;
               }
               G_demes[d]->itn = ITN;
            }
            else{
                 if ( G_demes[d]->itn==ITN){
                    G_groups->SUM_ITN = G_groups->SUM_ITN - 1;
                    G_groups->SUM_NOITN = G_groups->SUM_NOITN + 1;
                 }
                 G_demes[d]->itn = NOITN;
            }
        }
    }
}

void livertoblood(){
     Deme* deme;
     int d;
     float gamma, pci,pcni;

     for (d=0;d<G_N;d++){
         deme=G_demes[d];
         int region=deme->home;
         int ind=region-1000;

         if (deme->transmission==LIVERSTAGE){
             gamma=ran2();
             if (gamma<=G_gamma){
                 deme->transmission=BLOODSTAGE;

                  if(deme->resistance == RESISTRO ){
                        G_groups->SUM_LRO = G_groups->SUM_LRO - 1;
                        G_groups->SUM_BRO = G_groups->SUM_BRO + 1;
                        G_villdata[ind]->SUM_LRO = G_villdata[ind]->SUM_LRO - 1;
                        G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO + 1;
                  }

                  else if(deme->resistance == RESISTRA ){
                        G_groups->SUM_LRA = G_groups->SUM_LRA - 1;
                        G_groups->SUM_BRA = G_groups->SUM_BRA + 1;
                        G_villdata[ind]->SUM_LRA = G_villdata[ind]->SUM_LRA - 1;
                        G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA + 1;
                  }

                  else if(deme->resistance == RESISTRB ){
                        G_groups->SUM_LRB = G_groups->SUM_LRB - 1;
                        G_groups->SUM_BRB = G_groups->SUM_BRB + 1;
                        G_villdata[ind]->SUM_LRB = G_villdata[ind]->SUM_LRB - 1;
                        G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB + 1;
                  }
                 else { cout << "livertoblood: couldnt find res" << endl; int c = getchar();}


                 pci=ran2();
                 if (pci<= deme->probclinicallb){
                    G_groups->SUM_CLIN = G_groups->SUM_CLIN +1 ;
                    G_groups->SUM_NOINF = G_groups->SUM_NOINF -1 ;
                    deme->clinical = CLINICAL;
                 }
                 else{
                    G_groups->SUM_ASYM = G_groups->SUM_ASYM +1 ;
                    G_groups->SUM_NOINF = G_groups->SUM_NOINF -1 ;
                    deme->clinical = ASYMPTOMATIC;
                 }
             } //close gamma
         }  //close liverstage
     } //close for individuals
} //close function


void tosymptoms(){
     Deme* deme;
     int d;
     float sigma;

     for (d=0;d<G_N;d++){
         deme=G_demes[d];
         int region=deme->home;
         int ind=region-1000;
         if (deme->transmission==BLOODSTAGE){
             sigma=ran2();
             if (sigma<=G_sigma){

                 deme->transmission=INFECTIOUS;

                  if(deme->resistance == RESISTRO ){
                        G_groups->SUM_BRO = G_groups->SUM_BRO - 1;
                        G_groups->SUM_IRO = G_groups->SUM_IRO + 1;
                        G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO - 1;
                        G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO + 1;
                  }

                  else if(deme->resistance == RESISTRA ){
                        G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                        G_groups->SUM_IRA = G_groups->SUM_IRA + 1;
                        G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA - 1;
                        G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA + 1;
                  }

                  else if(deme->resistance == RESISTRB ){
                        G_groups->SUM_BRB = G_groups->SUM_BRB - 1;
                        G_groups->SUM_IRB = G_groups->SUM_IRB + 1;
                        G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB - 1;
                        G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB + 1;
                  }
                  else { cout << "bloodtosympt: couldnt find res" << endl; int c = getchar();}
            }
         }
     }
}


void recovery(){
     Deme* deme;
     int d;
     float delta;

     for (d=0;d<G_N;d++){
         deme=G_demes[d];
         int region=G_demes[d]->home;
         int ind=region-1000;
         if (G_demes[d]->transmission==INFECTIOUS){
             delta=ran2();
             if (delta<=G_delta){

                 G_demes[d]->transmission=SUSCEPTIBLE;
                 G_groups->SUM_S = G_groups->SUM_S + 1;
                 G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;

                 G_demes[d]->immunity=IMMUNE;
                 G_demes[d]->immunity_days=1;
                 G_groups->SUM_IMM = G_groups->SUM_IMM + 1;
                 G_groups->SUM_NIMM = G_groups->SUM_NIMM - 1;

                 if(G_demes[d]->resistance == RESISTRO ){
                    G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                    G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                    G_demes[d]->resistance=RESISTRO;
                 }

                 else if(G_demes[d]->resistance == RESISTRA ){
                    G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                    G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                    G_demes[d]->resistance=RESISTRO;
                 }

                 else if(G_demes[d]->resistance == RESISTRB ){
                    G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                    G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;
                    G_demes[d]->resistance=RESISTRO;
                 }
                 else { cout << "bloodtosympt: couldnt find res" << endl; int c = getchar();}

                 if (G_demes[d]->clinical==CLINICAL){
                    G_groups->SUM_CLIN = G_groups->SUM_CLIN -1 ;
                    G_groups->SUM_NOINF = G_groups->SUM_NOINF +1 ;
                    G_demes[d]->clinical = NOINFECTION;
                 }
                 else if (G_demes[d]->clinical==ASYMPTOMATIC){
                    G_groups->SUM_ASYM = G_groups->SUM_ASYM -1 ;
                    G_groups->SUM_NOINF = G_groups->SUM_NOINF +1 ;
                    G_demes[d]->clinical = NOINFECTION;
                 }
             }
         }
     }
}

void immunitycount(){
     Deme* deme;
     int d;

     for (d=0;d<G_N;d++){
         deme=G_demes[d];
         if (G_demes[d]->immunity==IMMUNE){
            G_demes[d]->immunity_days=G_demes[d]->immunity_days+1;
         }
     }
}

void immunityloss(){
     Deme* deme;
     int d;
     float alpha;

     for (d=0;d<G_N;d++){
         deme=G_demes[d];
         if (deme->immunity==IMMUNE && deme->immunity_days>40){
            alpha=ran2();
            if (alpha<=G_alpha){
                deme->immunity_level=deme->immunity_level-1;
                if (deme->immunity_level==0){
                    deme->immunity==NONIMMUNE;
                    G_groups->SUM_IMM = G_groups->SUM_IMM - 1;
                    G_groups->SUM_NIMM = G_groups->SUM_NIMM + 1;
                    deme->immunity_days=0;
                }
            }
         }
     }
}

void villages(int timestep, float seasons){
    Deme* deme;
    vill* village;
    int t=0;
    malpost=vector<int>(G_xregion,0);
    betaseason=vector<float>(G_xregion,0.0);

    for (int i=0;i<G_xregion;i++){
        village = G_villages[i];

        /** seasonality **/
//        betaseason[i] = betas[i]+(G_amp*betas[i]*cos(2.0*3.1416*((double(timestep)-90.0)/365.0)));

        betaseason[i] = (betas[i]*seasons)+(betas[i]*G_amp*0.5);
//        cout<<betaseason[i]<<endl;

        /** open posts according to data **/
        if (timestep==mptime[i] && mp[i]>0){
            if (village->post==0){
                village->post = 1;
                village->mptime=mptime[i];
            }
        }

//        /** shut down **/
//        if (timestep > 2191){
//           village->post = 0;
//        }

        /** scale up **/
        if (timestep > 2191){
            if (village->post ==0){
                float r = ran2();
                if (r<0.1){
                    village->post=1;
                    village->mptime=timestep;
                }
            }
        }

//        /** open posts with mda **/
//        if (timestep>G_mda_start1){
//           village->post= 1;
//        }

        /** count posts **/
        if (village->post== 1){
            t=t+1;
        }

    }
    cout << "NUMBER POSTS: "<<t<<endl;
}

void naturalhistoryinfection(int timestep){
    // Deme* deme;
    // vill* village;
    // int d;
    // float alpha, b, k1, mellow, mu;
    // float deathprob;
    // int aux1,aux2,aux3,auxt;
    // float t1,t2,pci,delta, randomEvent, randomEvent2;
    int rec=0;
    // int j, j1, j2;
    // int region, ind, malp, mpst;
    // int level, moi, cumm, age;
    // float susceptibility=1.0;
    int cc = 0;

    int tid = 0;

    #pragma omp parallel for private(tid)
    for (int d=0;d<G_N;d++){
    // for (d=0;d<G_N;d++){

        // tid = omp_get_thread_num();
        // cout << "tid:" << tid << ",dem:" << d << endl << std::flush;

        Deme* deme;
        vill* village;
        float alpha, b, k1, mellow, mu;
        float deathprob;
        int aux1,aux2,aux3,auxt;
        float t1,t2,pci,delta, randomEvent, randomEvent2;
        // int rec=0;
        int j, j1, j2;
        int region, ind, malp, mpst;
        int level, moi, cumm, age;
        float susceptibility=1.0;
        // int cc = 0;


        deme=G_demes[d];
        j=j1=j2=0;
        region=deme->now;
        ind=region-1000;
        village = G_villages[ind];
        b=betas[ind];
        malp=village->post;
        mpst=malp;
        level=deme->immunity_level;
        moi=deme->moi;
        cumm=deme->cummulative_exposures;
        age=deme->age;
        deme->clintoday=0;
//        /** village here now **/
//        village=G_villages[ind];
//        village->herenow.push_back(d);

        /** setprobdeath **/
        k1 = age/G_mu2;
        deathprob = (0.5*(exp(-age*G_sig)))+(2.5*G_sig*exp((1/G_mu1*age)/(1/1.8)));
        deme->death = deathprob*G_mu;
//        if (age>75){
//            deme->death = deathprob*G_mu/0.5;
//        }
//        else{deme->death = deathprob*G_mu/0.99;}

        /** birth death **/
        mu=ran2();
        if (mu<= deme->death){

            if(deme->drug == ARTESUNATE){
                G_groups->SUM_ART = G_groups->SUM_ART - 1;
                G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART - 1;
            }
            else if(deme->drug == PIPERAQUINE){
                G_groups->SUM_PIP = G_groups->SUM_PIP - 1;
                G_villdata[ind]->SUM_PIP = G_villdata[ind]->SUM_PIP - 1;
            }
            else if(deme->drug == PIPARTE){
                G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE - 1;
                G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE - 1;
            }
            else if(deme->drug == AL){
                G_groups->SUM_AL = G_groups->SUM_AL - 1;
                G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL - 1;
            }
            else if(deme->drug == LUMEFANTRINE){
                G_groups->SUM_LUM = G_groups->SUM_LUM - 1;
                G_villdata[ind]->SUM_LUM = G_villdata[ind]->SUM_LUM - 1;
            }

            deme->transmission = SUSCEPTIBLE;
            deme->resistance = RESISTRO;
            deme->drug = NODRUG;
            deme->prim = NO;
            deme->strategy = NOSTRAT;
            deme->age = 0;
            deme->immunity_level = 0;
            deme->immunity = NONIMMUNE;
            deme->moi = 0;
            deme->mda_days = 0;
            deme->mda_times = 0;
            deme->cummulative_exposures = 0;
            deme->clinical = NOINFECTION;
            deme->infectious_qeue.clear();
            deme->infections.clear();

            randomEvent = ran2();
            if (randomEvent <= 0.482){ //in accordance with 1998 census
                deme->gender = MALE;
            }
            else{
                deme->gender = FEMALE;
            }

            randomEvent2 = ran2();
            if (randomEvent2 <= G_itnprop){
               if ( deme->itn==NOITN){
                  G_groups->SUM_ITN = G_groups->SUM_ITN + 1;
                  G_groups->SUM_NOITN = G_groups->SUM_NOITN - 1;
               }
               deme->itn = ITN;
            }
            else{
                 if ( deme->itn==ITN){
                    G_groups->SUM_ITN = G_groups->SUM_ITN - 1;
                    G_groups->SUM_NOITN = G_groups->SUM_NOITN + 1;
                 }
                 deme->itn = NOITN;
            }
        }

        /**  aging **/
        if (timestep%365==0){
            deme->age = deme->age + 1;
        }

        /** drug loss **/
        if(deme->drug == ARTESUNATE ){
            randomEvent=ran2();
            if (randomEvent<=G_xa0){
                deme->drug = NODRUG;
                G_groups->SUM_ART = G_groups->SUM_ART - 1;
                G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART - 1;
                if (deme->infectious_qeue.empty()==false ||  deme->infections.empty()==false){
                   G_groups->SUM_FAIL = G_groups->SUM_FAIL + 1;
                }
            }
        }
        else if(deme->drug == PIPARTE ){ // need to reset deme->prim ? 
            randomEvent=ran2();
            randomEvent2=ran2(); // what is the relation between the two events?
            if (randomEvent<=G_xai){
                deme->drug = PIPERAQUINE;
                G_groups->SUM_PIP = G_groups->SUM_PIP + 1;
                G_villdata[ind]->SUM_PIP = G_villdata[ind]->SUM_PIP + 1;
                G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE - 1;
                G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE - 1;
            }
            else if (randomEvent2<=G_xab){
                deme->drug = ARTESUNATE;
                G_groups->SUM_ART = G_groups->SUM_ART + 1;
                G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE - 1;
                G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART + 1;
                G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE - 1;
            }
        }
        else if(deme->drug == PIPERAQUINE ){
            randomEvent=ran2();
            if (randomEvent<=G_xab){
                deme->drug = NODRUG;
                G_groups->SUM_PIP = G_groups->SUM_PIP - 1;
                G_villdata[ind]->SUM_PIP = G_villdata[ind]->SUM_PIP - 1;
                if (deme->infectious_qeue.empty()==false ||  deme->infections.empty()==false){
                    G_groups->SUM_FAIL = G_groups->SUM_FAIL + 1;
                }
            }
        }
        else if(deme->drug == AL ){
            randomEvent=ran2();
            randomEvent2=ran2();
            if (randomEvent<=G_xala){
                deme->drug = LUMEFANTRINE;
                G_groups->SUM_LUM = G_groups->SUM_LUM + 1;
                G_villdata[ind]->SUM_LUM = G_villdata[ind]->SUM_LUM + 1;
                G_groups->SUM_AL = G_groups->SUM_AL - 1;
                G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL - 1;
            }
            else if (randomEvent2<=G_xall){
                deme->drug = ARTESUNATE;
                G_groups->SUM_ART = G_groups->SUM_ART + 1;
                G_groups->SUM_AL = G_groups->SUM_AL - 1;
                G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART + 1;
                G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL - 1;
            }
        }
        else if(deme->drug == LUMEFANTRINE ){
            randomEvent=ran2();
            if (randomEvent<=G_xlum){
                deme->drug = NODRUG;
                G_groups->SUM_LUM = G_groups->SUM_LUM - 1;
                G_villdata[ind]->SUM_LUM = G_villdata[ind]->SUM_LUM - 1;
                if (deme->infectious_qeue.empty()==false ||  deme->infections.empty()==false){
                    G_groups->SUM_FAIL = G_groups->SUM_FAIL + 1;
                }
            }
        }
        if (deme->prim == YES ){
            randomEvent=ran2();
            if (randomEvent<=G_xprim){
                deme->prim = NO;
                G_groups->SUM_PRIM = G_groups->SUM_PRIM - 1;
                G_villdata[ind]->SUM_PRIM = G_villdata[ind]->SUM_PRIM - 1;
            }
        }


        /**  set susceptibility **/
//        susceptibility = (1-G_r*exp(-G_k*age));
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
           deme->infectiousness =  deme->infectiousness*G_phic;
           }
           else if (deme->clinical == ASYMPTOMATIC){
           deme->susceptibility =  deme->susceptibility;
           deme->infectiousness =  deme->infectiousness*G_phia;
           }
        }
        // prophylaxis from artesunate uptake
        if (deme->drug==ARTESUNATE || deme->drug==PIPARTE){
           deme->susceptibility =  deme->susceptibility* G_artprophylaxis;
        }

        /**  set clinical outcome probability **/
//        deme->probclinicallb=(0.4*exp(-(cumm-1)*0.1)+exp((-0.5*cumm)))/sqrt(level);
        deme->probclinicallb=(0.1*exp(-(cumm-2)*0.1)+exp((-0.9*cumm)))/sqrt(level);
        //deme->probclinicallb=exp(-0.1*(pow(cumm,1.65)));
        deme->probclinicalbb=exp(-0.15*(moi-1))*deme->probclinicallb;

        /**  set treatment rates **/ // Same values every step? Village property
        if (mpst == 1){
            treatclinical[d] = (1.0/G_timecomp)*G_fullcourse*G_covab; //clinical == asymptomatic? 
            //cout<<"postratec\t"<<treatclinical[d]<<endl;
            treatasympt[d] = (1.0/G_timecomp)*G_fullcourse*G_covab;
            }
        else{
            treatclinical[d] = (1.0/G_timecomp)*G_fullcourse*G_covab*G_nomp;
            //cout<<"treatclin: "<<treatclinical[d]<<endl;
            treatasympt[d] = (1.0/G_timecomp)*G_fullcourse*G_covab*G_nomp;
            //cout<<"treatssym: "<<treatasympt[d]<<endl;
        }

        /**  infection to bloodstage **/
        if (deme->blood_qeue.empty()==false){
            vector<int>::iterator it1 = deme->blood_qeue.begin();
            for ( ; it1 != deme->blood_qeue.end(); ) {
                j1+=1;
                t1=ran2();
                aux1 = deme->blood_qeue.at(j1-1);
                if (t1<G_gamma){
                    deme->transmission=BLOODSTAGE;
                    it1 = deme->blood_qeue.erase(it1);
                    deme->infectious_qeue.push_back(aux1);
                    deme->moi = deme->moi+1;
                    j1-=1;
                } else {
                    ++it1;
                }
            }
        }

        // Does this mean infections can go from blood_q to infectious_q
        // and then to infections (G_gamma * G_sigma) in one time step?

        // clinical != NOINFECTION iff infectious?

        /**  infection to infectious stage **/ 
        if (deme->infectious_qeue.empty()==false){
            vector<int>::iterator it2 = deme->infectious_qeue.begin();
            for ( ; it2 != deme->infectious_qeue.end(); ) {
                j2+=1;
                t2=ran2();
                aux2 = deme->infectious_qeue.at(j2-1);
                if (t2<G_sigma){
                    deme->transmission=INFECTIOUS;
                    deme->infections.push_back(aux2);
                    if (deme->moi>1){
                        pci=ran2();
                        if (pci<= deme->probclinicalbb){
                            deme->clinical = CLINICAL;
                            deme->clintoday=1;
                        }
                        else{deme->clinical = ASYMPTOMATIC;
                        deme->clintoday=0;}
                    }
                    else if (deme->moi == 1){
                        pci=ran2();
                        if (pci<= deme->probclinicallb){
                            deme->clinical = CLINICAL;
                            deme->clintoday=1;
                        }
                        else{deme->clinical = ASYMPTOMATIC;
                        deme->clintoday=0;}
                    }
                    it2 = deme->infectious_qeue.erase(it2);
                    j2-=1;
                } else {
                    ++it2;
                }
            }
        }


        /** recovery **/
        if (G_demes[d]->infections.empty()==false){ // what about blood_q
            vector<int>::iterator it = G_demes[d]->infections.begin();
            for ( ; it != G_demes[d]->infections.end(); ) {
                j+=1;
                delta=ran2();
                if (delta<=G_delta){
                    it = G_demes[d]->infections.erase(it);
                    G_demes[d]->moi-=1;
                    j-=1;
                    if (G_demes[d]->infectious_qeue.empty()==true && G_demes[d]->infections.empty()==true){
                        G_demes[d]->moi=0;
                        rec+=1;
                        G_demes[d]->transmission=SUSCEPTIBLE;
                        G_demes[d]->immunity_level= G_demes[d]->immunity_level+1;
                        G_demes[d]->immunity_days=1;
                        G_demes[d]->immunity=IMMUNE;
                        G_demes[d]->resistance=RESISTRO;
                        G_demes[d]->clinical = NOINFECTION;
                    }
                }
                else{++it;}
            }
        } //close recovery

        /** clinical resolution **/
        if (G_demes[d]->clinical== CLINICAL ){
            mellow=ran2();
            if (mellow<=G_mellow){
              G_demes[d]->clinical=ASYMPTOMATIC;
              G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
              G_groups->SUM_ASYM = G_groups->SUM_ASYM + 1;
            }
        }

        /** immunity loss **/
        if (deme->immunity==IMMUNE && deme->immunity_days>40 && deme->immunity_level>1 ){
            alpha=ran2();
            if (alpha<=G_alpha){
                deme->immunity_level=deme->immunity_level-1;
                if (deme->immunity_level==0){
                    deme->immunity==NONIMMUNE;
                    G_groups->SUM_IMM = G_groups->SUM_IMM - 1;
                    G_groups->SUM_NIMM = G_groups->SUM_NIMM + 1;
                    deme->immunity_days=0;
                }
            }
         }

        /** immunity count **/
        if (G_demes[d]->immunity==IMMUNE){
            G_demes[d]->immunity_days=G_demes[d]->immunity_days+1;
        }


        /** within host dominance **/
        int jr=0;
        int totres=0;
        if (deme->infections.empty()==false){
            vector<int>::iterator itt = deme->infections.begin();
            for ( ; itt != deme->infections.end(); ) {
                jr+=1;
                auxt = deme->infections.at(jr-1);
                totres=totres+(auxt-41);
                ++itt;
            }
//            int rr = int(ran2()*deme->infections.size());
//            deme->dominant=deme->infections.at(rr);

            //cout<<"within: "<<totres<<endl;
            if (totres==0){
                deme->dominant=RESISTRO;
            }
            else{
                if (deme->drug==ARTESUNATE){
                    deme->dominant=RESISTRA;
                }
                else if (deme->drug==PIPARTE){
                    deme->dominant=RESISTRA;
                }
                else if (deme->drug==NODRUG || deme->drug==PIPERAQUINE){
                    // cout << "tid:" << tid << ",dem:" << d  << ",size:" << deme->infections.size();
                    int rr = int(ran2()*deme->infections.size());
                    // cout << ",rr:" << rr << endl << std::flush;
                    // try{
                        deme->dominant=deme->infections.at(rr);
                    // } catch (const std::out_of_range& e) {
                    //     cout << "tid:" << tid << ",dem:" << d << endl << e.what() << endl << std::flush;
                    //     exit(0);
                    // }
                }
            }
        }


        /** drug effect **/
        // int j=0;
        j=0;
        int jj=0;

        if (deme->drug == NODRUG){
            continue;
        }
        else{
            if (deme->infectious_qeue.empty()==false){
                //cout<<"qeue"<<deme->infectious_qeue.size()<<endl;
                vector<int>::iterator it3 = deme->infectious_qeue.begin();
                for ( ; it3 != deme->infectious_qeue.end(); ) {
                    j+=1;
                    aux2 = deme->infectious_qeue.at(j-1);
                    if(deme->drug == ARTESUNATE){
                        if(aux2==RESISTRA){
                            t1=ran2();
                            if (t1<G_cBrada){
                               //cout<<"CLEAR 1"<<endl;
                               it3=deme->infectious_qeue.erase(it3);
                               j-=1;
                               if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                               //continue;
                            }
                            else{++it3;}//continue;}}
                        }
                        else if(aux2==RESISTRB){
                            t1=ran2();
                            if (t1<G_cBrbda){
                                //cout<<"CLEAR 2"<<endl;
                                it3=deme->infectious_qeue.erase(it3);
                                j-=1;
                                if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                //continue;
                            }
                            else{++it3;}//continue;}
                        }
                        else if(aux2==RESISTRO){
                            t1=ran2();
                            if (t1<G_cBroda){
                                //cout<<"CLEAR 3"<<endl;
                                it3=deme->infectious_qeue.erase(it3);
                                j-=1;
                                if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                //continue;
                            }
                            else{++it3;}//continue;}
                        }
                    }
                    else if (deme->drug == PIPERAQUINE){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradb){
                                    //cout<<"CLEAR 4"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdb){
                                    //cout<<"CLEAR 5"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodb){
                                    //cout<<"CLEAR 6"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                    }
                    else if (deme->drug == PIPARTE){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradab){
                                    //cout<<"CLEAR 7"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdab){
                                    //cout<<"CLEAR 8"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodab){
                                    //cout<<"CLEAR 9"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                    }
                    else if (deme->drug == LUMEFANTRINE){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradlf){
                                    //cout<<"CLEAR 10"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdlf){
                                    //cout<<"CLEAR 11"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodlf){
                                    //cout<<"CLEAR 12"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                    }
                    else if (deme->drug == AL){
                            if(aux2==RESISTRA){
                                t1=ran2();
                                if (t1<G_cBradal){
                                    //cout<<"CLEAR 13"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRB){
                                t1=ran2();
                                if (t1<G_cBrbdal){
                                    //cout<<"CLEAR 14"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it3;}//continue;}
                            }
                            else if(aux2==RESISTRO){
                                t1=ran2();
                                if (t1<G_cBrodal){
                                    //cout<<"CLEAR 15"<<endl;
                                    it3=deme->infectious_qeue.erase(it3);
                                    j-=1;
                                    //continue;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         //cout<<"RECOVERED 1"<<endl;
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                }
                                else{++it3;}//continue;}
                            }
                    }
                    else{++it3;}
                }

            } //close bloostage

            if (deme->infections.empty()==false){
                //cout<<"infections"<<deme->infections.size()<<endl;
                vector<int>::iterator it4 = deme->infections.begin();
                for ( ; it4 != deme->infections.end(); ) {
                    jj+=1;
                    aux3 = deme->infections.at(jj-1);
                    //cout<<"what value"<<aux3<<endl;
                    //cout<<"what drug"<<deme->drug<<endl;
                    if(deme->drug == ARTESUNATE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clrada){
                                    //cout<<"BLEAR 1"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbda){
                                    //cout<<"BLEAR 2"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clroda){
                                    //cout<<"BLEAR 3"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                    }
                    else if (deme->drug == PIPERAQUINE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradb){
                                    //cout<<"BLEAR 4"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdb){
                                    //cout<<"BLEAR 5"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodb){
                                    //cout<<"BLEAR 6"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                    }
                    else if (deme->drug == PIPARTE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradab){
                                    //cout<<"BLEAR 7"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdab){
                                    //cout<<"BLEAR 8"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodab){
                                    //cout<<"BLEAR 9"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                    }
                    else if (deme->drug == LUMEFANTRINE){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradlf){
                                    //cout<<"BLEAR 10"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdlf){
                                    //cout<<"BLEAR 11"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodlf){
                                    //cout<<"BLEAR 12"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                    }
                    else if (deme->drug == AL){
                            if(aux3==RESISTRA){
                                t1=ran2();
                                if (t1<G_clradal){
                                    //cout<<"BLEAR 13"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRB){
                                t1=ran2();
                                if (t1<G_clrbdal){
                                    //cout<<"BLEAR 14"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            else if(aux3==RESISTRO){
                                t1=ran2();
                                if (t1<G_clrodal){
                                    //cout<<"BLEAR 15"<<endl;
                                    it4=deme->infections.erase(it4);
                                    jj-=1;
                                    if (deme->infectious_qeue.empty()==true && deme->infections.empty()==true){
                                         cc=cc+1;
                                         G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                                         G_groups->SUM_S = G_groups->SUM_S + 1;
                                         G_groups->SUM_SUCC = G_groups->SUM_SUCC + 1;

                                         deme->transmission=SUSCEPTIBLE;
                                         deme->immunity=IMMUNE;
                                         deme->immunity_days=1;
                                         deme->immunity_level=deme->immunity_level+1;
                                         deme->moi=0;
                                         deme->clinical=NOINFECTION;
                                    }
                                    //continue;
                                }
                                else{++it4;}//continue;}
                            }
                            ++it4;
                    }
                    else{++it4;}
                } //iterate across all infections

            } // close infectious
        } // close drug effect
    }//close loop for individuals
    cout<<"natural recovery:"<<rec<<endl;
    cout<<"drug effect recovery:"<<cc<<endl;
} //close natural history


void clinicalresolution(){
    Deme* deme;
    int d;
    float mellow;

    for (d=0;d<G_N;d++){
        deme=G_demes[d];
        if (G_demes[d]->clinical== CLINICAL ){
            mellow=ran2();
            if (mellow<=G_mellow){
              G_demes[d]->clinical=ASYMPTOMATIC;
              G_groups->SUM_CLIN = G_groups->SUM_CLIN - 1;
              G_groups->SUM_ASYM = G_groups->SUM_ASYM + 1;
            }
        }
    }
}


void birthdeath(){
     Deme* deme;
     int d;
     float mu;
     float randomEvent;
     float randomEvent2;

     for (d=0;d<G_N;d++){
         deme=G_demes[d];
         mu=ran2();
         if (mu<= G_demes[d]->death){
         int region=G_demes[d]->home;
         int ind=region-1000;
         //cout<<"hey1"<<endl;
               if (G_demes[d]->transmission==INFECTIOUS){
                   if(G_demes[d]->resistance == RESISTRO ){
                        G_groups->SUM_IRO = G_groups->SUM_IRO - 1;
                        G_villdata[ind]->SUM_IRO = G_villdata[ind]->SUM_IRO - 1;
                   }

                   else if(G_demes[d]->resistance == RESISTRA ){
                        G_groups->SUM_IRA = G_groups->SUM_IRA - 1;
                        G_villdata[ind]->SUM_IRA = G_villdata[ind]->SUM_IRA - 1;
                   }

                   else if(G_demes[d]->resistance == RESISTRB ){
                        G_groups->SUM_IRB = G_groups->SUM_IRB - 1;
                        G_villdata[ind]->SUM_IRB = G_villdata[ind]->SUM_IRB - 1;

                   }
               }

               else if (G_demes[d]->transmission==BLOODSTAGE){

                   if(G_demes[d]->resistance == RESISTRO ){
                        G_groups->SUM_BRO = G_groups->SUM_BRO - 1;
                        G_villdata[ind]->SUM_BRO = G_villdata[ind]->SUM_BRO - 1;
                   }

                   else if(G_demes[d]->resistance == RESISTRA ){
                        G_groups->SUM_BRA = G_groups->SUM_BRA - 1;
                        G_villdata[ind]->SUM_BRA = G_villdata[ind]->SUM_BRA - 1;
                   }

                   else if(G_demes[d]->resistance == RESISTRB ){
                        G_groups->SUM_BRB = G_groups->SUM_BRB - 1;
                        G_villdata[ind]->SUM_BRB = G_villdata[ind]->SUM_BRB - 1;
                   }
              }

              else if (G_demes[d]->transmission==LIVERSTAGE){
                   if(G_demes[d]->resistance == RESISTRO ){
                        G_groups->SUM_LRO = G_groups->SUM_LRO - 1;
                        G_villdata[ind]->SUM_LRO = G_villdata[ind]->SUM_LRO - 1;
                   }

                   else if(G_demes[d]->resistance == RESISTRA ){
                        G_groups->SUM_LRA = G_groups->SUM_LRA - 1;
                        G_villdata[ind]->SUM_LRA = G_villdata[ind]->SUM_LRA - 1;
                   }

                   else if(G_demes[d]->resistance == RESISTRB ){
                        G_groups->SUM_LRB = G_groups->SUM_LRB - 1;
                        G_villdata[ind]->SUM_LRB = G_villdata[ind]->SUM_LRB - 1;
                   }
              }

              else if (G_demes[d]->transmission==SUSCEPTIBLE){
                    G_groups->SUM_S = G_groups->SUM_S - 1;
                    G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S - 1;
              }

              else { cout << "birthdeath: couldnt find transmission state" << endl; int c = getchar();}

              if(G_demes[d]->drug == ARTESUNATE){
                  G_groups->SUM_ART = G_groups->SUM_ART - 1;
                  G_villdata[ind]->SUM_ART = G_villdata[ind]->SUM_ART - 1;
              }
              else if(G_demes[d]->drug == PIPERAQUINE){
                  G_groups->SUM_PIP = G_groups->SUM_PIP - 1;
                  G_villdata[ind]->SUM_PIP = G_villdata[ind]->SUM_PIP - 1;
              }
              else if(G_demes[d]->drug == PIPARTE){
                  G_groups->SUM_PIPARTE = G_groups->SUM_PIPARTE - 1;
                  G_villdata[ind]->SUM_PIPARTE = G_villdata[ind]->SUM_PIPARTE - 1;
              }
              else if(G_demes[d]->drug == AL){
                  G_groups->SUM_AL = G_groups->SUM_AL - 1;
                  G_villdata[ind]->SUM_AL = G_villdata[ind]->SUM_AL - 1;
              }
              else if(G_demes[d]->drug == LUMEFANTRINE){
                  G_groups->SUM_LUM = G_groups->SUM_LUM - 1;
                  G_villdata[ind]->SUM_LUM = G_villdata[ind]->SUM_LUM - 1;
              }
              if(G_demes[d]->prim == YES){
                  G_groups->SUM_PRIM= G_groups->SUM_PRIM - 1;
                  G_villdata[ind]->SUM_PRIM = G_villdata[ind]->SUM_PRIM - 1;
              }

             G_groups->SUM_S = G_groups->SUM_S + 1;
             G_villdata[ind]->SUM_S = G_villdata[ind]->SUM_S + 1;
             G_demes[d]->transmission = SUSCEPTIBLE;
             G_demes[d]->resistance = RESISTRO;
             G_demes[d]->drug = NODRUG;
             G_demes[d]->prim = NO;
             G_demes[d]->strategy = NOSTRAT;
             G_demes[d]->age = 0;
             G_demes[d]->immunity_level = 0;
             G_demes[d]->moi = 0;
             G_demes[d]->cummulative_exposures = 0;

      //}
             randomEvent = ran2();
             if (randomEvent <= 0.482){ //in accordance with 1998 census
                 G_demes[d]->gender = MALE;
             }
             else{
                 G_demes[d]->gender = FEMALE;
             }


             randomEvent2 = ran2();
             if (randomEvent2 <= G_itnprop){
                if ( G_demes[d]->itn==NOITN){
                   G_groups->SUM_ITN = G_groups->SUM_ITN + 1;
                   G_groups->SUM_NOITN = G_groups->SUM_NOITN - 1;
                }
                G_demes[d]->itn = ITN;
             }
             else{
                  if ( G_demes[d]->itn==ITN){
                     G_groups->SUM_ITN = G_groups->SUM_ITN - 1;
                     G_groups->SUM_NOITN = G_groups->SUM_NOITN + 1;
                  }
                 G_demes[d]->itn = NOITN;
             }
          } // close death
      }   // close individuals' loop
}  // close function

void printResistance(int timestep){
     Deme* deme;
     vill* village;
     int RO =0; int RA=0; int RB=0;
     int ALLR0;
     int ALLRA;
     int ALLRB;
     int FAIL = 0;
     int CLEAR = 0;

     int I=0;
     int B=0;
     int L=0;
     int S=0;
     int CLIN = 0;
     int ASYM = 0;
     int NOI = 0;
     int postclin = 0;
     int nonpostclin = 0;
     int postasym = 0;
     int nonpostasym= 0;
     int popclin = 0;
     int popnonclin = 0;

     for (int d=0;d<G_xregion;d++){
        G_villdata[d]->SUM_IRO=0;
        G_villdata[d]->SUM_IRA=0;
        G_villdata[d]->SUM_IRB=0;
        G_villdata[d]->SUM_IRB=0;
        G_villdata[d]->SUM_ART=0;
        G_villdata[d]->SUM_PIP=0;
        G_villdata[d]->SUM_PIPARTE=0;
        G_villdata[d]->SUM_AL=0;
        G_villdata[d]->SUM_LUM=0;
        if (mp[d]>0){
            popclin=popclin+G_villdata[d]->popsize;
        }
        else{popnonclin=popnonclin+G_villdata[d]->popsize;}
     }

     for (int i=0; i<G_N; i++){
        deme=G_demes[i];
        int region=deme->now;
        int ind=region-1000;
        village=G_villages[ind];
        int j=0;
        int temp=0;

        if (deme->infections.empty()==false){
            I+=1;
            vector<int>::iterator it = deme->infections.begin();
            for ( ; it != deme->infections.end(); ) {
                j+=1;
                int aux2 = deme->infections.at(j-1);
                if (aux2>temp){
                    temp=aux2;
                }
                ++it;
            }
            if (temp==RESISTRO){
                RO+=1;
                G_villdata[ind]->SUM_IRO+=1;
            }
            else if (temp==RESISTRA){
                RA+=1;
                G_villdata[ind]->SUM_IRA+=1;
            }
            else if (temp==RESISTRB){
                RB+=1;
                G_villdata[ind]->SUM_IRB+=1;
            }

            if (deme->clinical==CLINICAL){
                CLIN+=1;
                if (mp[ind]>0){
                    postclin+=1;
                }
                else if (mp[ind]==0){
                    nonpostclin+=1;
                }
            }
            else if (deme->clinical==ASYMPTOMATIC){
                ASYM+=1;
                if (mp[ind]>0){
                    postasym+=1;
                }
                else if (mp[ind]==0){
                    nonpostasym+=1;
                }
            }
            else if (deme->clinical==NOINFECTION){
                NOI+=1;
            }
        }

        else if (deme->transmission==LIVERSTAGE){
            L+=1;
        }
        else if (deme->transmission==SUSCEPTIBLE){
            S+=1;
        }
        else if (deme->infectious_qeue.empty()==false && deme->infections.empty()==true){
            B+=1;
        }

        //else {cout<<"wrong transmission state:"<< deme->transmission <<endl;int gh = getchar();}
    }
    //cout<<"clins"<<CLIN<<endl;

    FAIL = G_groups->printFAILS();
    CLEAR = G_groups->printCLEARS();

    int SART = 0;
    int SPIP = 0;
    int SPIPARTE = 0;
    int SLUM = 0;
    int SAL = 0;
    int SSUM = 0;
    int SNO = 0;
    int MDAS =0;
    int TODAY=0;
//    SART=G_groups->printTreatA();
//    SPIP=G_groups->printTreatB();
//    SPIPARTE=G_groups->printTreatC();
//    SAL=G_groups->printTreatD();
//    SLUM=G_groups->printTreatE();

    int AGE;
    int SUMAGE = 0;
    int SUMINFAGE=0;
    float PREVAGE=0.0;
    int tot=0;
    int clins=0;


    for (int d=0; d<G_N; d++){
        deme=G_demes[d];
        //cout<<"AGE: " << deme->age<<endl;
        //see the age of a specific deme
        if (deme->age>=2 && deme->age <=10){
            SUMAGE=SUMAGE+1;
            if (deme->infections.empty()==false){
                SUMINFAGE=SUMINFAGE+1;
            }
        }
        if (deme->mda_times>0){
            MDAS+=1;
        }
        if (deme->drug==ARTESUNATE){
            SART+=1;
            tot+=1;
        }
        else if (deme->drug==PIPARTE){
            SPIPARTE+=1;
            tot+=1;
        }
        else if (deme->drug==PIPERAQUINE){
            SPIP+=1;
            tot+=1;
        }
        else if (deme->drug==LUMEFANTRINE){
            SLUM+=1;
            tot+=1;
        }
        else if (deme->drug==AL){
            SAL+=1;
            tot+=1;
        }
        else if (deme->drug==NODRUG){
            SNO+=1;
            tot+=1;
        }
        if (deme->clintoday==1){
            clins+=1;
        }
    }
    PREVAGE=float(SUMINFAGE/SUMAGE);
    SSUM=SART+SPIP+SPIPARTE+SLUM+SAL;


    //cout << "SUMINFAGE: "<<SUMINFAGE<< "\tSUMAGE: "<<SUMAGE<<endl;
    cout << "S:" << S << "\tL:" << L << "\tB:" << B <<  "\tI:" << I << "\tSUM:" << S+L+B+I<<endl;
    cout << "RO:" << RO << "\tRA:" << RA << "\tRB:" << RB << "\tSUM:" << RO+RA+RB<<endl;
    cout << "CLIN:" << CLIN <<"\tASYM: "<< ASYM << "\tSUMI: "<< CLIN+ASYM << "\tSUMALL: "<< CLIN+ASYM+NOI <<endl;
    cout << "FAIL:" << FAIL <<"\tCLEAR: "<< CLEAR <<"\tpopCLin: "<< popclin <<"\tpopnonCLin: "<< popnonclin << endl;
    cout << "ART:" << SART << "\tPIP:" << SPIP << "\tPIPARTE:" << SPIPARTE << "\tAL:" << SAL << "\tLUM:" << SLUM << "\tSUMDRUG:" << SSUM <<"\tCheck total: "<<tot<< endl;
    fout_resistance << RO << "\t" << RA << "\t" << RB << "\t" << RO+RA+RB<<"\t" << CLIN<<"\t" << ASYM<<"\t" << SUMAGE<<"\t" << SUMINFAGE<<"\t"<<  SART<<" \t"<<SPIP << "\t"<<SPIPARTE<<"\t"<<SAL <<"\t"<<SLUM <<"\t"<<SSUM << "\t"<< clins << "\t"<<FAIL<<"\t"<<CLEAR << "\t"<<postclin<<"\t"<<nonpostclin<<"\t"<<postasym<<"\t"<<nonpostasym <<endl;
    //fout_resistance<< " \t " << ALLR0 <<" \t "<< ALLRA << " \t "<< ALLRB << "\t"<< ALLR0+ALLRA+ALLRB<< endl;
}

void printTaus(){
     cout << "tau:" << G_tau << "\ttau1:" << G_tau1 << "\ttauab:" << G_tauab<< endl;
     fout_taus << G_tau << "\t" << G_tau1 << "\t" << G_tauab <<  endl;
}

void printTreat(){
     int SART = 0;
     int SPIP = 0;
     int SPIPARTE = 0;
     int SSUM = 0;
     int SAL = 0;
     int SLUM = 0;

     SART=G_groups->printTreatA();
     SPIP=G_groups->printTreatB();
     SPIPARTE=G_groups->printTreatC();
     SAL=G_groups->printTreatD();
     SLUM=G_groups->printTreatE();
     SSUM = SART+SPIP+SPIPARTE+SAL+SLUM;

     //fout_resistance << " ART: " << SART << " PIP: " << SPIP << " PIPARTE: " << SPIPARTE <<endl;
     cout << "ART:" << SART << "\tPIP:" << SPIP << "\tPIPARTE:" << SPIPARTE << "\tAL:" << SAL << "\tLUM:" << SLUM << "\tSUM:" << SSUM << endl;
}


void printbednet(){
     int SITN=0;
     int SNOITN =0;

     SITN = G_groups->printsumitn();
     SNOITN = G_groups->printsumnoitn();
     cout << " ITN: " << SITN << " NOITN: "<< SNOITN << " SUMITN: " << SITN+SNOITN <<endl;
}


void printage(){
     int AGE = 0;
     float SUSC = 0.0;
     float DEATH = 0.0;
     int DRUG = 13;
     float LAT;
     float LONG;

     Deme * deme;
     deme = G_demes[1100];

     AGE = deme->age;              //see the age of a specific deme
     SUSC = deme->susceptibility; //see susceptibility of a specific deme
     DEATH = deme->death;         //see dying prob of a specific deme
     DRUG = deme->drug;
     LAT = deme->latitude;
     LONG = deme->longitude;

     cout<< "Age:"<< AGE << "\tSUSCEPTIBILITY:" << SUSC << " \tDEATH:" << DEATH << endl;
     cout<< "\tDRUG:" << DRUG << "\tLONG:" << LONG << "\tLAT:"<< LAT<<endl;
}

void printinfo(int timestep){
     float LAT;
     float LONG;
     int AGE;
     int Tstate;
     int d;
     int RESIS;

     using namespace std;
   //  std::ofstream fout_info;
//     std::ostringstream strfile_info;
     ofstream fout_info;
     ostringstream strfile_info;

     strfile_info << "D:/Users/ricaguas/Desktop/Individual_Based/infoasym" ;
     strfile_info << "-time-"<<timestep ;
     strfile_info << ".txt";
     fout_info.open(strfile_info.str().c_str());

     Deme * deme;
     for (d=0;d<G_N;d++){

         deme = G_demes[d];

         //if (deme->transmission==INFECTIOUS){
         AGE = deme->age;
         LAT = deme->latitude;
         LONG = deme->longitude;
         Tstate=deme->transmission;
         RESIS=deme->resistance;


         fout_info  <<  LONG << "\t" << LAT <<  "\t" << Tstate << endl;
         //}
     }
     fout_info.close();
}

void printinfo2(int timestep){
     float LAT;
     float LONG;
     int Tstate;
     int RESIS;
     int d;
     int AGE;

     using namespace std;
   //  std::ofstream fout_info;
//     std::ostringstream strfile_info;
     ofstream fout_info;
     ostringstream strfile_info;

     strfile_info << "D:/Users/ricaguas/Desktop/Individual_Based/info" ;
     strfile_info << "-time-"<<timestep ;
     strfile_info << ".txt";
     fout_info.open(strfile_info.str().c_str());

     Deme * deme;
     for (d=0;d<G_N;d++){
         float G_trtasym=treatasympt[d]*0.1; // 0.0001?
         float G_trtclin=treatclinical[d];

         deme = G_demes[d];

         AGE = deme->age;
         LAT = deme->latitude;
         LONG = deme->longitude;
         Tstate=deme->transmission;
         RESIS=deme->resistance;


         fout_info <<   LONG << "\t" << LAT  << "\t" << Tstate <<"\t" << RESIS << "\t" << G_trtasym  << "\t" << G_trtclin<< endl;
     }
     fout_info.close();
}

void updatevilldata(int timestep){

//     if (timestep % 30 == 0){
//     Deme* deme;
//     float LAT;
//     float LONG;
//     int POP;
//     int RO;
//     int RA;
//     int RB;
//     int SART = 0;
//     int SPIP = 0;
//     int SPIPARTE = 0;
//     int SSUM = 0;
//     int SAL = 0;
//     int SLUM = 0;
//
//     for (int d=0;d<G_xregion;d++){
//        G_villdata[d]->SUM_IRO=0;
//        G_villdata[d]->SUM_IRA=0;
//        G_villdata[d]->SUM_IRB=0;
//        G_villdata[d]->SUM_IRB=0;
//        G_villdata[d]->SUM_ART=0;
//        G_villdata[d]->SUM_PIP=0;
//        G_villdata[d]->SUM_PIPARTE=0;
//        G_villdata[d]->SUM_AL=0;
//        G_villdata[d]->SUM_LUM=0;
//     }
//
     using namespace std;
     ofstream fout_info;
     ostringstream strfile_info;
     strfile_info << "outputs/info" ;
     strfile_info << "-trigg-"<<G_trigger<< "-delta-"<< G_delta;
     strfile_info << "-G_phia-"<<G_phia ;
     strfile_info << "-mobility-"<<G_mobility_static ;
     strfile_info << ".txt";
     fout_info.open(strfile_info.str().c_str());
//
//     for (int i=0; i<G_N; i++){
//        deme=G_demes[i];
//        int region=deme->now;
//        int ind=region-1000;
//
//        int j=0;
//        int temp=0;
//
//        if (deme->infections.empty()==false){
//            //cout<<"state"<<deme->transmission<<endl;
//            vector<int>::iterator it = deme->infections.begin();
//            for ( ; it != deme->infections.end(); ) {
//                j+=1;
//                int aux2 = deme->infections.at(j-1);
//                if (aux2>temp){
//                    temp=aux2;
//                }
//                ++it;
//            }
//            if (temp==RESISTRO){
//                G_villdata[ind]->SUM_IRO+=1;
//            }
//            else if (temp==RESISTRA){
//                G_villdata[ind]->SUM_IRA+=1;
//            }
//            else if (temp==RESISTRB){
//                G_villdata[ind]->SUM_IRB+=1;
//            }
//        }
//        if (deme->drug==ARTESUNATE){
//            G_villdata[ind]->SUM_ART+=1;
//        }
//        else if (deme->drug==PIPARTE){
//            G_villdata[ind]->SUM_PIPARTE+=1;
//        }
//        else if (deme->drug==PIPERAQUINE){
//            G_villdata[ind]->SUM_PIP+=1;
//        }
//        else if (deme->drug==LUMEFANTRINE){
//            G_villdata[ind]->SUM_LUM+=1;
//        }
//        else if (deme->drug==AL){
//            G_villdata[ind]->SUM_AL+=1;
//        }
//    }

     int d;
     vill* village;
     for (d=0;d<G_xregion;d++){
         village=G_villages[d];
         int days =village->mda_start;
         int open = village->mptime;

         fout_info<<days<< "\t" << open << endl;
     }

     fout_info.close();
}

void printmosqs(){
     int MEXP=0;
     int MINF =0;

     mosq* mosquito;

     for (int d=0;d<G_xregion;d++){
        mosquito=G_mosquitoes[d];
        MEXP=MEXP+mosquito->infected;
        MINF=MINF+mosquito->infectious;
     }
     cout << " MEXP: " << MEXP << " MINF: "<< MINF << endl;
}

void printagestruct(int timestep){
    Deme* deme;
    int i;

    if (timestep == 18250){
        cout<<"print age"<<endl;
        for (i=0;i<G_N;i++){
            deme = G_demes[i];
            int age=deme->age;
            int moi=deme->moi;
            int clin=deme->clinical;
            int level=deme->immunity_level;
            int cumm = deme->cummulative_exposures;
            //cout<<age<<"\t"<<moi<<"\t"<<clin<<"\t"<<level<<"\t"<<cumm<<endl;
            fout_age<<age<<"\t"<<moi<<"\t"<<clin<<"\t"<<level<<"\t"<<cumm<<endl;
        }
        fout_age.close();
    }
}


void printregion(){
    Deme* deme;
    vill* village;
    inf=vector<int>(G_xregion,0);
    tot=vector<int>(G_xregion,0);

    for ( int i=0;i<G_xregion;i++){
        village=G_villages[i];
        int villsize = village->herenow.size();
        tot[i]=villsize;
        for (int j=0; j<villsize; j++){
            int indexp = village->herenow.at(j);
            deme = G_demes[indexp];
            if (deme->clinical==CLINICAL){
               inf[i]+=1;
            }
        }
        float prop = float(float(inf[i])/float(tot[i]));
        //cout<< "contas: "<<prop<<endl;
        fout_regions<<prop<<"\t";
    }
    fout_regions<<endl;
}
