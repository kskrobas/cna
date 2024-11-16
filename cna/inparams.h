/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * inparams.h
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
#ifndef INPARAMS_H
#define INPARAMS_H

#include <algorithm>
#include <string>
#include <vector>

using namespace std;

typedef double position;
typedef const position cpos;
typedef const double cdouble;
typedef std::string str;


//-----------------------------------------------------------------------------
// trim from start
static inline std::string &ltrim(std::string &s)
{
auto isNotSpace=[](char c){return !std::isspace(c);};
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), isNotSpace));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s)
{
auto isNotSpace=[](char c){return !std::isspace(c);};
        s.erase(std::find_if(s.rbegin(), s.rend(), isNotSpace).base(), s.end());
        return s;
}




// trim from both ends
static inline std::string &trim(std::string &s)
{
        return ltrim(rtrim(s));
}

//-----------------------------------------------------------------------------
template<typename T>
vector<T>  split(const T & str, const T & delimiters)
{
vector<T> v;
typename T::size_type start = 0;
auto pos = str.find_first_of(delimiters, start);

        while(pos != T::npos) {
            if(pos != start) // ignore empty tokens
                v.emplace_back(str, start, pos - start);
            start = pos + 1;
            pos = str.find_first_of(delimiters, start);
        }
        if(start < str.length()) // ignore trailing delimiter
            v.emplace_back(str, start, str.length() - start); // add what's left of the string
return v;
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

enum EFTYPE{nxyz,txyz};
enum EPNF{pos,neg,fcc,nfcc,zb,nzb};
enum EAAN{ignore,neighb};  // save atoms and neighborhood (awn)

struct StFileNameType{

string fileName;
EFTYPE f_type;
EPNF   l_type;

EAAN   aan=EAAN::ignore;
size_t nOfnb;

    StFileNameType(){ }
    StFileNameType(string fn__,EFTYPE type__=EFTYPE::nxyz)
        :fileName(fn__),f_type(type__){ }
    StFileNameType(string fn__,EFTYPE type__,EPNF ltype__)
        :fileName(fn__),f_type(type__),l_type(ltype__){ }

    StFileNameType(EFTYPE type__,EPNF ltype__)
        :f_type(type__),l_type(ltype__){ }

    StFileNameType(EAAN aan__,const size_t nOfnb__)
        :aan(aan__),nOfnb(nOfnb__) { }

    void operator() (EFTYPE type__,EPNF ltype__)
            { f_type=type__; l_type=ltype__; }

};


//-----------------------------------------------------------------------------


struct StOutFileNames{
vector<StFileNameType> posVerifiedAtoms;
vector<StFileNameType> lattVerifiedAtoms;
vector<StFileNameType> negVerifiedAtoms;
vector<StFileNameType> lattNegVerifiedAtoms;
vector<StFileNameType> selAtomsAndNeighbs;

    bool empty(){ return posVerifiedAtoms.empty() &&
                        negVerifiedAtoms.empty() &&
                        lattVerifiedAtoms.empty() &&
                        lattNegVerifiedAtoms.empty() &&
                        selAtomsAndNeighbs.empty(); }
};
//------------------------------------------------------------------

//-----------------------------------------------------------------------------
struct StBox{
enum ETYPE{IGN,CUB,CYL,SPH} btype=IGN;  //default ignore box boundaries
    union{
        struct{position xlo,xhi,ylo,yhi,zlo,zhi;};
        position bounds[6];
    };

    bool isPointInside(const position &x,const position &y,const position &z ) const ;

};


//-----------------------------------------------------------------------------

struct StInParams{

string  *inFileName, *outStatFileName;
StOutFileNames *outFileNames;
bool *printStat;
string *scmd,*sval;
position *tol,*dst,*tolA;
size_t *nb;
StBox *box;
vector<string> ignoreKeyValue;
string selectedRegion;
//enum EIGNORECNA{CNAOFF,CYL} ignorecna=CNAOFF;
enum ANGLEDISTR{ADOFF,FCC,ZB} adistr=ADOFF;
int threads=1;



};


class Cnocna{


public:




};


//-----------------------------------------------------------------------------

bool parseInputScript(const string &fileName,StInParams &sparams,bool &verb);


#endif // INPARAMS_H
