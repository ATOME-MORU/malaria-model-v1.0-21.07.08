
#define SEGS 2

#include <vector>
#include <list>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>

using namespace std;

void space_vector( ifstream& ifs, vector<float>& v ){
 string s;
 getline( ifs, s );
 istringstream iss( s );
 copy( istream_iterator<float>( iss ), istream_iterator<float>(), back_inserter( v ) );
}
void space_vector2( ifstream& ifs2, vector<int>& listvill ){
 string s;
 getline( ifs2, s );
 istringstream iss( s );
 copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( listvill ) );
}
void migrate_vector( ifstream& ifs3, vector<float>& v2 ){
 string s;
 getline( ifs3, s );
 istringstream iss( s );
 copy( istream_iterator<float>( iss ), istream_iterator<float>(), back_inserter( v2 ) );
}
void static_vector( ifstream& ifs4, vector<float>& v22 ){
 string s;
 getline( ifs4, s );
 istringstream iss( s );
 copy( istream_iterator<float>( iss ), istream_iterator<float>(), back_inserter( v22 ) );
}
void season_vector( ifstream& ifseason, vector<float>& vseas ){
 string s;
 getline(ifseason, s );
 istringstream iss( s );
 copy( istream_iterator<float>( iss ), istream_iterator<float>(), back_inserter( vseas ) );
}

void age_vector( ifstream& ifage, vector<int>& vage ){
 string s;
 getline(ifage, s );
 istringstream iss( s );
 copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( vage ) );
}


void data_pres( ifstream& ifdatap, vector<int>& vdatap ){
 string s;
 getline(ifdatap, s );
 istringstream iss( s );
 copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( vdatap ) );
}

void data_values( ifstream& ifdata, vector<int>& vdata ){
 string s;
 getline(ifdata, s );
 istringstream iss( s );
 copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( vdata ) );
}

void data_tests( ifstream& iftests, vector<int>& vtests ){
 string s;
 getline(iftests, s );
 istringstream iss( s );
 copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( vtests) );
}

void data_pop( ifstream& ifdatapop, vector<float>& vpop ){
 string s;
 getline(ifdatapop, s );
 istringstream iss( s );
 copy( istream_iterator<float>( iss ), istream_iterator<float>(), back_inserter( vpop ) );
}


void immune_vector( ifstream& ifimmune, vector<int>& vimmune ){
 string s;
 getline(ifimmune, s );
 istringstream iss( s );
 copy( istream_iterator<int>( iss ), istream_iterator<int>(), back_inserter( vimmune ) );
}

typedef struct vs {
        vector<float> v;
}VectorSpace;

typedef struct vm {
        vector<float> v2;
}VectorSpace2;

typedef struct vst {
        vector<float> v22;
}VectorSpace3;

typedef struct vss {
        vector<float> vseas;
}VectorSpace5;

typedef struct vage {
        vector<int> v2a;
}VectorSpaceage;

typedef struct vdatap {
        vector<int> vdtp;
}VectorSpacedatap;

typedef struct vdata {
        vector<int> vdt;
}VectorSpacedata;

typedef struct vtests {
        vector<int> vtst;
}VectorSpacetests;

typedef struct vpop {
        vector<float> vp;
}VectorSpacepop;

typedef struct vimmune {
        vector<int> vimm;
}VectorSpaceimmune;

typedef struct listv_{
        vector<int> listvill;
        int listorder;
}listv;

typedef struct list_neighs {
        vector<int> neighbours;
}Neighbourhood;


typedef struct datav {
        int SUM_ART;
        int SUM_PIP;
        int SUM_PIPARTE;
        int SUM_LUM;
        int SUM_AL;
        int SUM_PRIM;
        int SUM_S ;
        int	SUM_LRO	;
        int	SUM_LRA	;
        int	SUM_LRB	;
        int	SUM_IRO	;
        int	SUM_IRA	;
        int	SUM_IRB	;
        int	SUM_BRO	;
        int	SUM_BRA	;
        int	SUM_BRB	;
        int SUM_INF;

        void init(){
         SUM_S      =   0;
         SUM_LRA	=	0;
         SUM_LRB	=	0;
         SUM_LRO	=	0;
         SUM_BRA	=	0;
         SUM_BRB	=	0;
         SUM_BRO	=	0;
         SUM_IRA	=	0;
         SUM_IRB	=	0;
         SUM_IRO	=	0;
         SUM_ART = 0;
         SUM_PIP = 0;
         SUM_PIPARTE = 0;
         SUM_AL = 0;
         SUM_LUM = 0;
         SUM_PRIM= 0;
         SUM_INF = 0;
         }

         int getS(){
             return
             SUM_S;
             }
         int getI(){
             return
             SUM_IRA	+
             SUM_IRB	+
             SUM_IRO	;
             }
         int getL(){
             return
             SUM_LRA	+
             SUM_LRB	+
             SUM_LRO	;
             }
         int getB(){
                return
                SUM_BRA	+
                SUM_BRB	+
                SUM_BRO	;
                }

         int printTreatA(){
              return
               SUM_ART;
               }

        int printTreatB(){
              return
              SUM_PIP;
              }

        int printTreatC(){
              return
              SUM_PIPARTE;
              }

        int printTreatD(){
              return
              SUM_AL;
              }

        int printTreatE(){
              return
              SUM_LUM;
              }

        int printALL(){
              return
              SUM_S     +
              SUM_BRA	+
              SUM_BRB	+
              SUM_BRO	+
              SUM_LRA	+
              SUM_LRB	+
              SUM_LRO	+
              SUM_IRA	+
              SUM_IRB	+
              SUM_IRO	;
              }

        float lat;
        float lon;
        int popsize;
}villdata;



typedef struct state_ {
        int	SUM_S	;
        int	SUM_LRO	;
        int	SUM_LRA	;
        int	SUM_LRB	;
        int	SUM_BRO	;
        int	SUM_BRA	;
        int	SUM_BRB	;
        int	SUM_IRO	;
        int	SUM_IRA	;
        int	SUM_IRB	;

        int SUM_ALLR0 ;
        int SUM_ALLRA ;
        int SUM_ALLRB ;

        int SUM_ART;
        int SUM_PIP;
        int SUM_PIPARTE;
        int SUM_PRIM;
        int SUM_AL;
        int SUM_LUM;
        int SUM_FAIL;
        int SUM_SUCC;

        int SUM_ITN;
        int SUM_NOITN;

        int SUM_PREG;
        int SUM_CLIN;
        int SUM_ASYM;
        int SUM_NOINF;
        int SUM_IMM;
        int SUM_NIMM;

        void init(){
                    SUM_BRA	=	0;
                    SUM_BRB	=	0;
                    SUM_BRO	=	0;

                    SUM_IRA	=	0;
                    SUM_IRB	=	0;
                    SUM_IRO	=	0;

                    SUM_LRA	=	0;
                    SUM_LRB	=	0;
                    SUM_LRO	=	0;

                    SUM_S = 0;

                    SUM_ALLR0 = 0;
                    SUM_ALLRA = 0;
                    SUM_ALLRB = 0;

                    SUM_ART = 0;
                    SUM_PIP = 0;
                    SUM_PIPARTE = 0;
                    SUM_AL = 0;
                    SUM_PRIM = 0;
                    SUM_LUM = 0;

                    SUM_FAIL = 0;
                    SUM_SUCC = 0;

                    SUM_ITN = 0;
                    SUM_NOITN = 0;

                    SUM_PREG = 0;
                    SUM_CLIN = 0;
                    SUM_ASYM = 0;
                    SUM_NOINF = 0;
                    SUM_IMM = 0;
                    SUM_NIMM = 0;
             }

	     int getB(){
             return
                SUM_BRA	+
                SUM_BRB	+
                SUM_BRO	;
            }

            int getI(){
                return
                SUM_IRA	+
                SUM_IRB	+
                SUM_IRO	;
                }

                int getL(){
                    return
                    SUM_LRA	+
                    SUM_LRB	+
                    SUM_LRO	;
                }

                int getS(){
                    return
                    SUM_S
                  	;
                    }
                int getC(){
                    return
                    SUM_CLIN
                  	;
                    }
                int getASYM(){
                    return
                    SUM_ASYM
                  	;
                    }
                int getNI(){
                    return
                    SUM_NOINF
                  	;
                    }
                int getIM(){
                    return
                    SUM_IMM
                  	;
                    }
                int getNIM(){
                    return
                    SUM_NIMM
                  	;
                    }

        void printTransmission(){
             cout << "S:" << this->getS() << "\tL:" << this->getL() << "\tB:" << this->getB() << "\tI:" << this->getI() << "\tSUM:" << getI()+getS()+getL()+getB() << endl;
        }
        void printClinical(){
             cout << "Clin:" << this->getC()<< "\tAsym:" << this->getASYM()<< "\tnotinf:" << this->getNI() << "\tSUM:" << getC()+getASYM()+getNI() << endl;
        }
        void printImmune(){
             cout << "Immune:" << this->getIM()<< "\tnonImmune:" << this->getNIM() << "\tSUM:" << getIM()+getNIM() << endl;
        }

        void printSum(){
             int sum =

                        SUM_S	+
                        SUM_LRO	+
                        SUM_LRA	+
                        SUM_LRB	+
                        SUM_BRO	+
                        SUM_BRA	+
                        SUM_BRB	+
                        SUM_IRO	+
                        SUM_IRA	+
                        SUM_IRB	;

             cout << "SUM:" << sum << endl;
        }


        int printsumitn(){
            return
            SUM_ITN;
            }

        int printsumnoitn(){
            return
            SUM_NOITN;
            }

        int printIRes0(){
            return
           	SUM_IRO;
        //cout<<" res0 "<< Res0;
        }

        int printIResA(){
            return
         	SUM_IRA;
      	  // cout<<" resA "<< ResA;
          }

        int printIResB(){
            return
         	SUM_IRB;
          //cout<<" resB "<< ResB;
          }

        int printTreatA(){
              return
               SUM_ART;
        }

        int printTreatB(){
              return
              SUM_PIP;
              }

        int printTreatC(){
              return
              SUM_PIPARTE;
              }

        int printTreatD(){
              return
              SUM_AL;
              }

        int printTreatE(){
              return
              SUM_LUM;
              }

        int printFAILS(){
              return
              SUM_FAIL;
              }
        int printCLEARS(){
              return
              SUM_SUCC;
              }
} StateDeme;


typedef struct Deme_ {
        int drug;
        int prim;
        int strategy;
        int resistance;
        int transmission;
        int age;
        int gender;
        int itn;
        float susceptibility;
        float infectiousness;
        float probclinicallb;
        float probclinicalbb;
        float death;
        int pregnancy;
        float latitude;
        float longitude;
        int pregnancy_days;
        int parity;
        list<int> nearest;
        int lstsize;
        int region;
        int clinical;
        int immunity;
        int immunity_days;
        int mda_days;
        int mda_times;
        int home;
        int oldhome;
        int now;
        int mobile;
        int visitdays;
        int immunity_level;
        int cummulative_exposures;
        int moi;
        int super;
        int dominant;
        vector<int> blood_qeue;
        vector<int> infectious_qeue;
        vector<int> infections;
        int clintoday;
        int asymtoday;
		void init(){}
}Deme;

typedef struct Village_{
        int npeople;
        int mda;
        int mda_days;
        int mda_times;
        int counter;
        int delay;
        int mobile;
        float latitude;
        float longitude;
        int post;
        int mda_start;
        int mptime;
        int postfail;
        int clinical;
        int mdastart;
        std::vector<int> nearest;
        std::vector<int> seasonalnet;
        std::vector<int> migrantnet;
        std::vector<int> staticnet;
        std::vector<int> herenow;
        std::vector<float> migrantflows;
        std::vector<float> staticflows;
        std::vector<int> incidence;
        std::vector<float> incidencecorr;
        std::vector<float> incidencecorrmda;
        std::vector<int> datapres;
        std::vector<int> data;
        std::vector<float> datacorr;
        std::vector<int> tests;
        std::vector<float> poptest;
}vill;


typedef struct Mosquitoes_ {
        int infected;
        int infectious;
        vector<int>resistance;
        vector<int>resistancei;
}mosq;

