#include "stgrain.h"
#include "colormsg.h"
#include "affinemat.h"
#include "cprogress.h"
#include <omp.h>
#include <functional>


//-----------------------------------------------------------------------------
inline position sqr(cpos &x){return x*x;}

std::ostream &operator <<(std::ostream &o, const StAtom &v )
{
return    o<<v.x<<"    "<<v.y<<"    "<<v.z;
}

//-----------------------------------------------------------------------------
int StGrain::findAtomName(const string &aname__)
{
int apos=0;

            for(auto &atype :atomTypes){
            std::string & aname(atype.name);
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

return false;
}
//---------------------------------------------------------
bool StGrain::openXYZFile(const string &fileName)
{
fstream fin(fileName,ios::in);

        if(!fin){
            errMsg(" file doesn't exist");
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
bool saveAtoms(string &fileName, const StGrain &grain, const size_t nOfB, EFTYPE ftype, EPNF pnf)
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

            if(numOfatoms>1e6){
                progress.title=(" progress writing: ");
                progress.start(numOfatoms);
            }

            for(size_t i=0;i<numOfatoms;i++,progress++){
                if(np_condition(atoms[i])){
                    fout<<atomNT(atoms[i])<<"    "
                        <<atoms[i].x<<"    "<<atoms[i].y<<"    "<<atoms[i].z
                        <<"    "<<atoms[i].nOfn<<"    "<<atoms[i].fcc<<endl;
                    nOfrows++;
                }
            }


            fout.seekg(0);
            fout<<nOfrows;
            fout.close();

            progress.stop();
            cout<<"\r";

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












