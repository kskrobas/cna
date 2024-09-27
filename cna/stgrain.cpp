#include "stgrain.h"
#include "colormsg.h"
#include "affinemat.h"
#include "cprogress.h"
#include <omp.h>
#include <functional>
#include <sstream>

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
inline position sqr(cpos &x){return x*x;}

std::ostream &operator <<(std::ostream &o, const StAtom &v )
{
return    o<<v.x<<"    "<<v.y<<"    "<<v.z;
}

//-----------------------------------------------------------------------------
int StGrain::findAtomName(const string &aname__) const
{
int apos=0;

            for(auto &atype :atomTypes){
            const std::string & aname(atype.name);
                if( aname==aname__ )
                return apos;

                apos++;
            }

return -1;
}
//---------------------------------------------------------
StGrain::StGrain()
{
        threads=1;
}
//---------------------------------------------------------
bool StGrain::openFile(const string &fileName)
{
        if(fileName.find(".xyz")!=string::npos){
        return openXYZFile(fileName);
        }
        if(fileName.find(".lmp")!=string::npos){
        return openLMPFile(fileName);
        }

        errMsg(" unrecognized file format/extenion: "+fileName);

return false;
}
//---------------------------------------------------------
bool StGrain::openXYZFile(const string &fileName)
{
fstream fin(fileName,ios::in);

            if(!fin){
                errMsg(" file doesn't exist: "+fileName);
            return false;
            }

CProgress progress;
int row=0,arows=-1;
position cx,cy,cz;  // center position

            cx=cy=cz=0;


            try{
                    fin.exceptions(ifstream::failbit | ifstream::badbit | ifstream::eofbit);

                    fin>>arows;
                    while(fin.get()!='\n' ) ;


                    if(arows<1){
                        cerr<<" rows <1"<<endl;
                    return false;
                    }


                    while(fin.get()!='\n' ) ;  //ignore a comment row

                    atoms.clear();
                    atoms.reserve(arows);

            position x,y,z;
            string aname;
            int atype;
            std::function<StAtom(position &, position &, position &, int &, int &) > createAtom;

                    if(AAEnabled || ZBAAEnabled){
                     createAtom=[](position &x, position &y, position &z , int &at, int &id)
                                    { return StAtom(x,y,z,at,id);} ;
                    }else{
                     createAtom=[](position &x, position &y, position &z , int &at, int &id)
                                       { (void) id; return StAtom(x,y,z,at);} ;
                    }


                    if(arows>=1e6){
                        progress.title=" progress reading: ";
                        progress.start(arows);
                    }


                    fin.exceptions(ifstream::failbit | ifstream::badbit);


                    for(row=0 ;row<arows;row++,progress++){
                        fin>>aname>>x>>y>>z;
                        cx+=x; cy+=y; cz+=z;

                        if( (atype=findAtomName(aname))<0){
                            atomTypes.push_back(StAtomType(aname));
                            atype=atomTypes.size()-1;
                        }

                        atoms.emplace_back(createAtom(x,y,z,atype,row));
                    }

                    atoms.shrink_to_fit();
                    fin.close();
                    progress.stop();
                    cout<<"\r";
                    cout.flush();

            }
            catch(std::ifstream::failure &e){
                cerr<<" exceptions during file procesing, row: "<<row<<", e.what(): "<<e.what()<<endl;
                fin.close();
                progress.stop();
                return false;
            }

cpos iN=1.0/arows;

            cx*=iN; cy*=iN; cz*=iN;

            omp_set_num_threads(threads);
            #pragma omp parallel for
            for(int i=0;i<arows;i++){
                atoms[i].x-=cx;
                atoms[i].y-=cy;
                atoms[i].z-=cz;
                atoms[i].set_r2();
            }
            return true;
}
//-----------------------------------------------------------------------------
bool StGrain::openLMPFile(const string &fileName)
{
fstream fin(fileName,ios::in);

        if(!fin){
            errMsg(" file doesn't exist: "+fileName);
        return false;
        }

CProgress progress;
int nAtoms,row=2,atypes,id=-1;
position cx,cy,cz;  // center position

            cx=cy=cz=0;


            try{
            std::string fline;

                fin.exceptions(ifstream::failbit | ifstream::badbit | ifstream::eofbit);
                //ignore 1st, 2nd lines
                std::getline(fin,fline);
                std::getline(fin,fline);
                ////////////////////////

                do{
                    row++;
                    std::getline(fin,fline);
                    if(fline.empty()) continue;

                vector<string> toks{split<string>(fline," \t")};

                    if(toks[0]=="Atoms")
                        break;

                    if( toks.size()<2 || toks.size()>3)
                        continue;


                    if(toks[1]=="atoms"){
                        nAtoms=std::stoi(toks[0]);
                    continue;
                    }

                    if(toks[1]=="atom"){
                        atypes=std::stoi(toks[0]);
                        atomTypes.resize(atypes);

                        for(int i=0;i<atypes;i++)
                            atomTypes[i].name=std::to_string(i+1);

                    continue;
                    }

                } while(true);


                //ignore the line
                std::getline(fin,fline);

            int atype;
            position x,y,z;
            std::function<StAtom(position &, position &, position &, int &, int &) > createAtom;

                    if(AAEnabled || ZBAAEnabled){
                     createAtom=[](position &x, position &y, position &z , int &at, int &id)
                                    { return StAtom(x,y,z,at,id);} ;
                    }else{
                     createAtom=[](position &x, position &y, position &z , int &at, int &id)
                                       { (void) id; return StAtom(x,y,z,at);} ;
                    }

                atoms.reserve(nAtoms);
                row++;

                for(int atom=0;atom<nAtoms;atom++,row++){
                    fin>>id; // eof error detection


                    std::getline(fin,fline);
                vector<string> toks{split<string>(fline," \t")};

                    atype=std::stoi(toks[0])-1;
                    x=std::stod(toks[1]);
                    y=std::stod(toks[2]);
                    z=std::stod(toks[3]);

                    cx+=x; cy+=y; cz+=z;

                    atoms.emplace_back(createAtom(x,y,z,atype,atom));
                }

                atoms.shrink_to_fit();
            }
            catch(std::ifstream::failure &e){
            std::stringstream sstream;
                sstream<<" exceptions during file procesing, row: "<<row<<", "<<e.what();
                errMsg(sstream.str());
                fin.close();
                progress.stop();
            return false;
            }
            catch(const std::invalid_argument &ia){
            std::stringstream sstream;
                sstream<<" invalid argument exception, row:  "<<row<<", "<<ia.what();
                errMsg(sstream.str());
                fin.close();
                progress.stop();
            return false;
            }



cpos iN=1.0/nAtoms;

            cx*=iN; cy*=iN; cz*=iN;

            omp_set_num_threads(threads);
            #pragma omp parallel for
            for(int i=0;i<nAtoms;i++){
                atoms[i].x-=cx;
                atoms[i].y-=cy;
                atoms[i].z-=cz;
                atoms[i].set_r2();
            }

return true;
}





//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

double atomDistR2(const StAtom &a, const StAtom &b)
{
return sqr(a.x-b.x)+sqr(a.y-b.y)+sqr(a.z-b.z);
}
//-----------------------------------------------------------------------------

bool CNA(StGrain &grain, cpos &distance, cpos &tolerance, cpos &toleranceA)
{
//cpos dist2=sqr(distance);
const size_t N=grain.atoms.size();
const size_t Nh=N/2;

cpos dist2Min=sqr(distance-tolerance);
cpos dist2Max=sqr(distance+tolerance);
CProgress progress;

                 progress.title=" progress CNA: ";
                 progress.start(N);

std::function<void (StAtom &a, StAtom &b)> updateNOfN;

                if(grain.AAEnabled || grain.ZBAAEnabled){
                    updateNOfN=[](StAtom &a, StAtom &b)
                                { a.nOfn++; b.nOfn++; a.nID.push_back(b.id); b.nID.push_back(a.id);};
                }
                else{
                    updateNOfN=[](StAtom &a, StAtom &b)
                                { a.nOfn++; b.nOfn++;};
                }


                omp_set_num_threads(grain.threads);
                #pragma omp parallel
                {
                size_t i,j;
                position adist2;                

                    #pragma omp for
                    for(i=0;i<Nh;i++){
                    auto &atomA=grain.atoms[i];

                        for(j=i+1;j<N;j++){
                        auto &atomB=grain.atoms[j];

                            adist2=atomDistR2(atomA,atomB);

                            if(dist2Min<adist2 && adist2<dist2Max)
                                #pragma omp critical
                                updateNOfN(atomA,atomB);
                        }
                        #pragma omp critical
                        progress++;
                    }


                    #pragma omp for
                    for(i=N-2;i>=Nh;i--){
                    auto &atomA=grain.atoms[i];

                        for(j=i+1;j<N;j++){
                        auto &atomB=grain.atoms[j];

                            adist2=atomDistR2(atomA,atomB);

                            if(dist2Min<adist2 && adist2<dist2Max)
                                #pragma omp critical
                                updateNOfN(atomA,atomB);
                        }

                        #pragma omp critical
                        progress++;
                    }
                }

                progress.stop();
                cout<<"\r";
                cout.flush();

                if(grain.AAEnabled){
                        #pragma omp parallel for
                        for(size_t i=0; i<N; i++)
                            grain.atoms[i].FCC_NN_Angle_Analysis(grain,toleranceA);
                }
                else{
                    if(grain.ZBAAEnabled){
                        #pragma omp parallel for
                        for(size_t i=0;i<N;i++)
                            grain.atoms[i].ZB_NN_Angle_Analysis(grain,toleranceA);
                    }
                }


return true;
}

//-----------------------------------------------------------------------------


class CSaveOptions{
public:
    virtual ~CSaveOptions() {  }
    virtual bool check(const StAtom & ) { return true; }
};
//.........................................................
class COptionsInterface: public CSaveOptions{
protected:
        CSaveOptions *cso;
public:

        COptionsInterface(CSaveOptions *cso__):cso(cso__) {  }
};
//.........................................................
class COpt_AtomTypeIgnore: public COptionsInterface{

public:
        size_t type;
        COpt_AtomTypeIgnore(CSaveOptions *cso):COptionsInterface(cso) { }
        COpt_AtomTypeIgnore(CSaveOptions *cso,const size_t type__):COptionsInterface(cso),type(type__) { }

        bool check(const StAtom & atom)
        { return atom.atype!=type && cso->check(atom); }

};
//.........................................................
class COpt_NumberOfNeighborsIgnore: public COptionsInterface{
public:
        size_t nOfn;
        COpt_NumberOfNeighborsIgnore(CSaveOptions *cso): COptionsInterface(cso) {  }
        COpt_NumberOfNeighborsIgnore(CSaveOptions *cso,const size_t nOfn__): COptionsInterface(cso),nOfn(nOfn__) {  }

        bool check(const StAtom &atom)
        { return atom.nOfn!=nOfn && cso->check(atom); }
};
//.........................................................
class COpt_OutsideBoxIgnore: public COptionsInterface{
public:
        StBox box;
        COpt_OutsideBoxIgnore(CSaveOptions *cso): COptionsInterface(cso) {  }
        COpt_OutsideBoxIgnore(CSaveOptions *cso,const StBox box__):COptionsInterface(cso),box(box__){ }

        bool check(const StAtom &atom)
        {return box.isAtomInside(atom) && cso->check(atom); }
};
//.........................................................


//-----------------------------------------------------------------------------
bool saveAtoms(string &fileName, const StGrain &grain, const size_t nOfB,
               EFTYPE ftype, EPNF pnf, const StBox &box, const vector<string> &ignore)
{
fstream fout(fileName,ios::out);

            if(!fout){
                errMsg("file not saved "+fileName);
            return false;
            }


auto &atoms=grain.atoms;
const size_t numOfatoms=atoms.size();
size_t nOfrows=0;
std::function<bool(const StAtom &)> np_condition;
std::function<string(const StAtom &)> atomNT;
std::function<string(const StAtom &)> typeAcc;
CProgress progress;

            fout<<"          "<<endl;

            if(pnf==EPNF::pos)
                    {   np_condition=[&nOfB](const StAtom &a){ return a.nOfn==nOfB;} ;
                        fout<<"--- positively verified atoms ---"<<endl; }
            else{
                    if(pnf==EPNF::neg)
                    {   np_condition=[&nOfB](const StAtom &a){ return a.nOfn!=nOfB;} ;
                        fout<<"--- negatively verified atoms ---"<<endl; }
                    else{
                        if(pnf==EPNF::fcc)
                            {   np_condition=[](const StAtom &a){ return a.fcc;} ;
                                fout<<"--- FCC verified atoms ---"<<endl; }
                        else
                            if(pnf==EPNF::nfcc)
                            {   np_condition=[](const StAtom &a){ return !a.fcc;} ;
                                fout<<"--- non FCC verified atoms ---"<<endl; }
                            else
                                if(pnf==EPNF::zb)
                                {   np_condition=[](const StAtom &a){ return a.zb;} ;
                                    fout<<"--- ZB verified atoms ---"<<endl; }
                                else
                                {   np_condition=[](const StAtom &a){ return !a.zb;} ;
                                    fout<<"--- non ZB verified atoms ---"<<endl; }
                        }
            }
            switch (ftype){
            case EFTYPE::nxyz:  atomNT=[&grain](const StAtom &a) { return grain.atomTypes[a.atype].name; } ; break;
            case EFTYPE::txyz:  atomNT=[](const StAtom &a) { return std::to_string(a.nOfn); } ; break;
            }

CSaveOptions *cso=new CSaveOptions;
vector<CSaveOptions *> ptr_cso;

            ptr_cso.reserve(ignore.size()+1);
            ptr_cso.push_back(cso);

            if(!ignore.empty()){            
                    for(auto & iopt: ignore){
                    const vector<string> toks{split<string>(iopt," \t")};

                            if (toks[0]=="atype" ){                            
                                if(const int atype=grain.findAtomName(toks[1]); atype>=0)
                                    cso=new COpt_AtomTypeIgnore(ptr_cso.back(),atype);
                                else{
                                    warnMsg("nosave atype "+toks[1]+" unknown");
                                    continue;
                                }
                            }

                            if (toks[0]=="nb")
                                cso=new COpt_NumberOfNeighborsIgnore(ptr_cso.back(),std::stoi(toks[1]));

                            ptr_cso.push_back(cso);
                    }//for
            }//if


            if(box.btype!=StBox::IGN){
                cso=new COpt_OutsideBoxIgnore(ptr_cso.back(),box);
                ptr_cso.push_back(cso);
            }


            if(numOfatoms>1e6){
                progress.title=(" progress writing: ");
                progress.start(numOfatoms);
            }

            /////////////////////////////////////////////////
            /// save atoms
            ///
            for(size_t i=0;i<numOfatoms;i++,progress++){

                if( !ptr_cso.back()->check(atoms[i]))
                    continue;;

                if(np_condition(atoms[i])){
                    fout<<atomNT(atoms[i])<<"    "
                        <<atoms[i].x<<"    "<<atoms[i].y<<"    "<<atoms[i].z
                        <<"    "<<atoms[i].nOfn<<endl;
                    nOfrows++;
                }
            }
            /////////////////////////////////////////////////////


            fout.seekg(0);
            fout<<nOfrows;
            fout.close();

            progress.stop();
            cout<<"\r";

            for(auto &a: ptr_cso)
                delete a;

return true;
}
//-----------------------------------------------------------------------------
void StAtom::FCC_NN_Angle_Analysis(const StGrain &grain, cpos &toleranceA)
{
        if(nOfn!=12) return;

const auto &atoms=grain.atoms;
auto & refAtom=atoms[nID[0]];
StVector ref(refAtom.x-x,refAtom.y-y,refAtom.z-z);
StVector b;
double vc;
size_t nOf_90=0;
size_t nOf_60=0;
size_t nOf_180=0;

cpos tol_60min=0.5-toleranceA;
cpos tol_60max=0.5+toleranceA;
cpos tol_90=toleranceA;
cpos tol_180min=1-toleranceA;
//cpos tol_180max=0=/.5+toleranceA;


        for(size_t i=1;i<12;i++){
            b.x=atoms[nID[i]].x-x;
            b.y=atoms[nID[i]].y-y;
            b.z=atoms[nID[i]].z-z;

            vc=std::fabs(cosa(ref,b));

            if(tol_60min<vc && vc<tol_60max) nOf_60++;
            else{
                if(vc<tol_90) nOf_90++;
                else
                    if(tol_180min<vc &&vc<=1) nOf_180++;
            }

        }

        fcc= (nOf_180==1  && nOf_90==2 && nOf_60==8);
}
//-----------------------------------------------------------------------------
void StAtom::ZB_NN_Angle_Analysis(const StGrain &grain, cpos &toleranceA)
{
            if(nOfn!=4) return;

const auto &atoms=grain.atoms;
auto & refAtom=atoms[nID[0]];
StVector ref(refAtom.x-x,refAtom.y-y,refAtom.z-z);
StVector b;
double vc;

size_t nOf_109=0;

cpos tol_109min=1.0/3.0-toleranceA;
cpos tol_109max=1.0/3.0+toleranceA;


            for(size_t i=1;i<4;i++){
                b.x=atoms[nID[i]].x-x;
                b.y=atoms[nID[i]].y-y;
                b.z=atoms[nID[i]].z-z;
                vc=std::fabs(cosa(ref,b));

                if(tol_109min<vc && vc<tol_109max)
                    nOf_109++;

            }

            zb=(nOf_109 == 3);
}

//-----------------------------------------------------------------------------
bool StBox::isAtomInside(const StAtom &atom) const
{
bool testX=( (atom.x>xlo) & (atom.x<xhi) );
bool testY=( (atom.y>ylo) & (atom.y<yhi) );
bool testZ=( (atom.z>zlo) & (atom.z<zhi) );

return testX && testY && testZ;
}


/*
if(box.btype==StBox::IGN){ //
    for(size_t i=0;i<numOfatoms;i++,progress++){

        if( !ptr_cso.back()->check(atoms[i]))
            continue;;

        if(np_condition(atoms[i])){
            fout<<atomNT(atoms[i])<<"    "
                <<atoms[i].x<<"    "<<atoms[i].y<<"    "<<atoms[i].z
                <<"    "<<atoms[i].nOfn<<"    "<<atoms[i].fcc<<endl;
            nOfrows++;
        }
    }
}
else{
    for(size_t i=0;i<numOfatoms;i++,progress++){

        if( !ptr_cso.back()->check(atoms[i]))
            continue;;

        if(np_condition(atoms[i]) )
            if(box.isAtomInside(atoms[i])){
                fout<<atomNT(atoms[i])<<"    "
                    <<atoms[i].x<<"    "<<atoms[i].y<<"    "<<atoms[i].z
                    <<"    "<<atoms[i].nOfn<<endl;
                nOfrows++;
            }
        }
    }
    */

