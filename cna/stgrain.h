/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * stgrain.h
* Copyright (C) 2019 Kazimierz Skrobas <kskrobas@unipress.waw.pl>
 *
 * "cna" is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cna is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef STGRAIN_H
#define STGRAIN_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "inparams.h"
using namespace std;



struct StAtom;
class StGrain;
typedef vector<StAtom> vatoms;

ostream & operator <<(ostream &o,const StAtom &a);





//-----------------------------------------------------------------------------
struct StAtom{
position x,y,z;
size_t id;
size_t atype;  ///  position of an atom type in array
double r2;
vector<size_t> nID; //neighbor idh
size_t nOfn;  //NumberOfNeighbors
bool fcc=false;
bool zb=false;


enum EREGTYPE{OUT,BULK,MARGIN};
size_t rtype;  /// position type: in/out/boundary of box

    StAtom(cpos &x__, cpos &y__, cpos &z__, const size_t &atype__)
        {x=x__;y=y__;z=z__; atype=atype__; nOfn=0;}

    StAtom(cpos &x__, cpos &y__, cpos &z__, const size_t &atype__,const size_t id__)
        {x=x__;y=y__;z=z__; atype=atype__; nOfn=0;id=id__;nID.reserve(24);}

    void set_r2(){r2=x*x+y*y+z*z;}

    void FCC_NN_Angle_Analysis(const StGrain &grain,cpos &toleranceA);  // analysis of angles and neighbors for FCC lattice
    void ZB_NN_Angle_Analysis(const StGrain &grain, cpos &tolerancA);

    friend ostream & operator <<(ostream &o,const StAtom &a);
    void fullInfo();


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
    StInParams *inparams;

    vatoms atoms;    
    vector<size_t> count_OBM; /// counter of outside, bulk and margin atoms
    vector<StAtomType> atomTypes;        
    vector<StBoxPlane> boxWalls;

    StGrain();
    bool openFile(const string & fileName);
    int findAtomName(const string &aname__) const ;

private:
    bool openXYZFile(const string &fileName);
    bool openLMPFile(const string &fileName);


};


//-----------------------------------------------------------------------------

double atomDistR2(const StAtom &a, const StAtom &b);
bool CNA(StGrain &grain,cpos &distance, cpos &tolerance, cpos &toleranceA);


bool saveAtoms(string &fileName, const StGrain &grain, const size_t nOfB,
               EFTYPE ftype, EPNF pnf, const StBox &box,const vector<string> &ignore);
bool saveAtomsAndNeighbors(string &fileName, const StGrain &grain,
                           const size_t nOfN,EFTYPE ftype);

#endif // STGRAIN_H
