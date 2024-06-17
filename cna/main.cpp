#include <iostream>
#include <algorithm>
#include <iomanip>
#include"stgrain.h"
#include"colormsg.h"
#include"help.h"

using namespace std;

#define vline "━━"

//-----------------------------------------------------------------------------


struct StOutFileNames{
vector<StFileNameType> posVerifiedAtoms;
vector<StFileNameType> fccVerifiedAtoms;
vector<StFileNameType> negVerifiedAtoms;
vector<StFileNameType> fccNegVerifiedAtoms;
    bool empty(){ return posVerifiedAtoms.empty() && negVerifiedAtoms.empty() && fccVerifiedAtoms.empty() && fccNegVerifiedAtoms.empty(); }
};
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
int main(int argc, char *argv[])
{
        if(argc<2){
            errMsg(" not enough parameters");
            help();
        return 1;
        }

StGrain grain;

string inFileName,outStatFileName;
StOutFileNames outFileNames;
bool verb=false,printStat=false;
string scmd,sval;
position tol=-1,dst=-1,tolA=-1;
size_t nb=0;

        try{

            for(int i=1;i<argc;){
                scmd=std::string(argv[i++]);

                //------------------------
                if(scmd=="-h"){
                    help();
                return 0;
                }

                if(scmd=="-aa"){
                    grain.AAEnabled=true;
                continue;
                }

                if(scmd=="-ps"){
                    printStat=true;
                continue;
                }

                if (scmd=="-v"){
                    verb=true;
                continue;
                }
                //------------------------
                if(i>=argc){
                    errMsg(" unrecognized argument "+scmd);
                return 1;
                }

                sval=std::string(argv[i++]);

                if(scmd=="-d"){
                    dst=std::stod(sval);
                continue;
                }

                if(scmd=="-i"){
                    inFileName=std::move(sval);
                continue;
                }

                if(scmd=="-nb"){
                    nb=(size_t) std::stoi(sval);
                continue;
                }

                if(scmd=="-ofcc"){
                    outFileNames.fccVerifiedAtoms.emplace_back(StFileNameType(sval));
                continue;
                }

                if(scmd=="-ofcct"){
                    outFileNames.fccVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz));
                continue;
                }

                if(scmd=="-onfcc"){
                    outFileNames.fccNegVerifiedAtoms.emplace_back(StFileNameType(sval));
                continue;
                }

                if(scmd=="-onfcct"){
                    outFileNames.fccNegVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz));
                continue;
                }

                if(scmd=="-on"){
                    outFileNames.negVerifiedAtoms.emplace_back(StFileNameType(sval));
                continue;
                }

                if(scmd=="-ont"){
                    outFileNames.negVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz));
                continue;
                }

                if(scmd=="-op"){
                    outFileNames.posVerifiedAtoms.emplace_back(StFileNameType(sval));
                continue;
                }

                if(scmd=="-opt"){
                    outFileNames.posVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz));
                continue;
                }

                if(scmd=="-psf"){
                    outStatFileName=sval;
                continue;
                }

                if(scmd=="-tol"){
                    tol=std::stod(sval);
                continue;
                }

                if(scmd=="-tolA"){
                    tolA=std::stod(sval);
                continue;
                }

                if(scmd=="-th"){
                   grain.threads=std::stoi(sval);
                continue;
                }
                //------------------------

                errMsg(" unrecognized argument "+scmd);
                throw 1;
            }

            if(dst<0)               {errMsg(" distance not given "); throw 2;}
            if(inFileName.empty())  {errMsg(" input file not given "); throw 2;}
            if(outFileNames.empty()){errMsg(" output file not given ");throw 2;}
            if(tol<0)               {warnMsg(" tolerance not given, assumed default value (0.1)"); tol=0.125;}
            if(grain.threads<1)     {errMsg(" number of threads musn't be less than 1"); throw 2;}

            if(grain.AAEnabled && tolA <0) {errMsg(" angle tolerance (-tolA) not given "); throw 2;}

        }
        catch(const int &err){
        return 2;
        }
        catch (const std::invalid_argument& ia){
            errMsg("invalid argument:  "+scmd);
        return 2;
        }
        catch(...){
        return -1;
        }
        //---------------------------------------------------------

        if(verb){
            infoMsg("parameters OK");
        }

        if(verb){ infoMsg("file reading :"+inFileName);}
        if(!grain.openFile(inFileName)) return 3;

        if(verb){ infoMsg("CNA analysis");}

        //---------------------------------------------------------

        CNA(grain,dst,tol,tolA);

        //---------------------------------------------------------

        if(printStat || ! outStatFileName.empty()){
        //auto cmp=[](const StAtom &a, const StAtom &b){return a.nOfn<b.nOfn;};
        vector<size_t> total_nOfn;
        const size_t N=grain.atoms.size();
        size_t nOffcc=0;
        constexpr int nN=20;

                total_nOfn.resize(nN);
                for(auto & v: total_nOfn) v=0;

                for(size_t i=0;i<N;i++){
                    total_nOfn.at(grain.atoms[i].nOfn)++;
                    nOffcc+=(size_t) grain.atoms[i].fcc;
                }

                if(printStat){
                    infoMsg(" statistic data :");

                    for(size_t i=0;i<nN;i++){
                        if(i%5==0)
                        { cout<<endl<<" "; for(size_t i=0;i<45;i++) cout<<vline; cout<<"\n┃"; }

                        cout<<" "<<setw(3)<<i<<" : "<<setw(8)<<total_nOfn[i]<<"  ┃";
                    }

                    cout<<endl; for(size_t i=0;i<45;i++) cout<<vline; cout<<"\n";

                    cout<<"\n total number of atoms: "<<N;

                    if(grain.AAEnabled)
                        cout<<"\n fcc atoms : "<<nOffcc<<endl;

                    cout<<endl;
                }

                if(!outStatFileName.empty()){
                fstream fout(outStatFileName,ios::out);

                        if(!fout){
                            errMsg(" couldn't save statistic to file "+outStatFileName);
                        }
                        else{
                        std::time_t datetime = std::time(nullptr);
                            fout<<"#statistic data "<<endl;
                            fout<<"#date: "<<std::asctime(std::localtime(&datetime));
                            fout<<"#input: "<<inFileName<<endl;
                            fout<<"#size: 20"<<endl;
                            fout<<"#fcc atoms: "<<nOffcc<<endl;
                            fout<<"#n.type   abundance"<<endl;

                            for(size_t i=0;i<nN;i++)
                                fout<<" "<<setw(3)<<i<<"    "<<setw(8)<<total_nOfn[i]<<endl;
                        }

                        fout.close();
                }

        }



        //-----------------------------------------------------------------------------------------

        if(verb && !outFileNames.negVerifiedAtoms.empty()){ infoMsg("save nv. outcomes");}

        for(auto & fn: outFileNames.negVerifiedAtoms)
            saveAtoms(fn.fileName,grain,nb,fn.type,EPNF::neg);


        if(verb && !outFileNames.posVerifiedAtoms.empty() ){ infoMsg("save pv. outcomes");}

        for(auto & fn: outFileNames.posVerifiedAtoms)
            saveAtoms(fn.fileName,grain,nb,fn.type,EPNF::pos);


        if(verb && !outFileNames.fccVerifiedAtoms.empty()){infoMsg("save fcc outcomes");}

        for(auto &fn:outFileNames.fccVerifiedAtoms)
            saveAtoms(fn.fileName,grain,nb,fn.type,EPNF::fcc);


        if(verb && !outFileNames.fccNegVerifiedAtoms.empty()){infoMsg("save non fcc outcomes");}

        for(auto &fn:outFileNames.fccNegVerifiedAtoms)
            saveAtoms(fn.fileName,grain,nb,fn.type,EPNF::nfcc);

        //-----------------------------------------------------------------------------------------
        cout<<endl;

return 0;
}
