#ifndef INPARAMS_H
#define INPARAMS_H

#include "stgrain.h"

#include <algorithm>




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


struct StOutFileNames{
vector<StFileNameType> posVerifiedAtoms;
vector<StFileNameType> lattVerifiedAtoms;
vector<StFileNameType> negVerifiedAtoms;
vector<StFileNameType> lattNegVerifiedAtoms;
    bool empty(){ return posVerifiedAtoms.empty() && negVerifiedAtoms.empty() && lattVerifiedAtoms.empty() && lattNegVerifiedAtoms.empty(); }
};
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

struct StInParams{
StGrain *grain;
string  *inFileName, *outStatFileName;
StOutFileNames *outFileNames;
bool *printStat;
string *scmd,*sval;
position *tol,*dst,*tolA;
size_t *nb;
StBox *box;

};

//-----------------------------------------------------------------------------

bool parseInputScript(const string &fileName,StInParams &sparams,bool &verb);


#endif // INPARAMS_H
