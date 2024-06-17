#ifndef STGRAIN_H
#define STGRAIN_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

typedef double position;
typedef const position cpos;
typedef const double cdouble;
typedef std::string str;

struct StAtom;
class StGrain;
typedef vector<StAtom> vatoms;

enum EFTYPE{nxyz,txyz};

struct StFileNameType{

string fileName;
EFTYPE type;
    StFileNameType(){ }
    StFileNameType(string fn__,EFTYPE type__=EFTYPE::nxyz)
        :fileName(fn__),type(type__){ }

};

//-----------------------------------------------------------------------------
struct StAtom{
position x,y,z;
size_t id;
size_t atype;  ///  position of an atom type in array
double r2;
vector<size_t> nID; //neighbor idh
size_t nOfn;
bool fcc=false;

    StAtom(cpos &x__, cpos &y__, cpos &z__, const size_t &atype__)
        {x=x__;y=y__;z=z__; atype=atype__; nOfn=0;}

    StAtom(cpos &x__, cpos &y__, cpos &z__, const size_t &atype__,const size_t id__)
        {x=x__;y=y__;z=z__; atype=atype__; nOfn=0;id=id__;nID.reserve(20);}

    void set_r2(){r2=x*x+y*y+z*z;}

    void FCC_NN_Angle_Analysis(const StGrain &grain,cpos &toleranceA);  // analysis of angles and neighbors for FCC lattice

};

struct StAtomType{
std::string name,charge;

    StAtomType() { name ="?"; charge="?";}
    StAtomType(const string name__):name (name__) { }

    bool operator == (const StAtomType &a){
        return a.name==name;
    }

};

struct StBoxPlane{
double A,B,C,D;

static double maxX,minX;
static double maxY,minY;
static double maxZ,minZ;
};

//-----------------------------------------------------------------------------
class StGrain
{
public:

    vatoms atoms;
    vector<StAtomType> atomTypes;
    vector<StBoxPlane> boxWalls;
    int threads;
    bool AAEnabled=false; //angle analysis on/off


    StGrain();
    bool openFile(const string & fileName);

private:
    bool openXYZFile(const string &fileName);

    int findAtomName(const string &aname__);
};


//-----------------------------------------------------------------------------

double atomDistR2(const StAtom &a, const StAtom &b);
bool CNA(StGrain &grain,cpos &distance, cpos &tolerance, cpos &toleranceA);

enum EPNF{pos,neg,fcc,nfcc};
bool saveAtoms(string &fileName,const StGrain &grain, const size_t nOfB,EFTYPE ftype, EPNF pnf);

#endif // STGRAIN_H
