/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cpp
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
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <cmath>
#include "stgrain.h"
#include "colormsg.h"
#include "help.h"
#include "inparams.h"

using namespace std;

#define vline "━━"


//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
        if(argc<2){
            errMsg(" not enough parameters");
            help();
        return 1;
        }


StInParams inParams;
StGrain grain;
StBox box;
string inFileName,outStatFileName;
StOutFileNames outFileNames;
bool verb=false,printStat=false;
string scmd,sval;
position tol=-1,dst=-1,tolA=-1;
size_t nb=0;


        inParams.inFileName=&inFileName;
        inParams.outStatFileName=&outStatFileName;
        inParams.outFileNames=&outFileNames;
        inParams.printStat=&printStat;
        inParams.tol=&tol;
        inParams.dst=&dst;
        inParams.tolA=&tolA;
        inParams.nb=&nb;
        inParams.box=&box;

        grain.inparams=&inParams;


        try{

            for(int i=1;i<argc;){
                scmd=std::string(argv[i++]);

                //------------------------
                if(scmd=="-inp"){
                    sval=std::string(argv[i++]);

                    if(!parseInputScript(sval,inParams,verb))
                        return 1;

                continue;
                }

                //------------------------



                if(scmd=="-h"){
                    help();
                return 0;
                }

                if(scmd=="-aafcc"){
                    inParams.adistr=StInParams::FCC;
                continue;
                }

                if(scmd=="-aazb"){
                    inParams.adistr=StInParams::ZB;
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
                    outFileNames.lattVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::nxyz,EPNF::fcc));
                continue;
                }

                if(scmd=="-ofcct"){
                    outFileNames.lattVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz,EPNF::fcc));
                continue;
                }

                if(scmd=="-ozb"){
                    outFileNames.lattVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::nxyz,EPNF::zb));
                continue;
                }

                if(scmd=="-ozbt"){
                    outFileNames.lattVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz,EPNF::zb));
                continue;
                }

                if(scmd=="-onfcc"){
                    outFileNames.lattNegVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::nxyz,EPNF::nfcc));
                continue;
                }

                if(scmd=="-onfcct"){
                    outFileNames.lattNegVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz,EPNF::nfcc));
                continue;
                }

                if(scmd=="-onzb"){
                    outFileNames.lattNegVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::nxyz,EPNF::nzb));
                continue;
                }

                if(scmd=="-onzbt"){
                    outFileNames.lattNegVerifiedAtoms.emplace_back(StFileNameType(sval,EFTYPE::txyz,EPNF::nzb));
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
                   inParams.threads=std::stoi(sval);
                continue;
                }
                //------------------------

                errMsg(" unrecognized argument "+scmd);
                throw 1;
            }

            if(dst<0)               {errMsg(" distance not given "); throw 2;}
            if(inFileName.empty())  {errMsg(" input file not given "); throw 2;}
            if(outFileNames.empty()){errMsg(" output file not given ");throw 2;}
            if(tol<0)               {warnMsg(" tolerance not given, assumed default value (0.125)"); tol=0.125;}
            if(inParams.threads<1)     {errMsg(" number of threads musn't be less than 1"); throw 2;}

            if(inParams.adistr==StInParams::FCC && tolA <0) {errMsg(" angle tolerance (-tolA) not given "); throw 2;}

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

        if(verb) { infoMsg("show stat");}
        if(printStat || ! outStatFileName.empty()){        
        vector<size_t> total_nOfn;
        const size_t N=grain.atoms.size();
        size_t nOffcc=0,nOfzb=0,negAtoms=0;
        size_t max_nOfn=1;
        constexpr size_t nOfcols=5;
        
                omp_set_num_threads(inParams.threads);
                #pragma omp parallel for reduction(max:max_nOfn)
                for(size_t i=0;i<N;i++){
                    if(grain.atoms[i].nOfn>max_nOfn)
                        max_nOfn=grain.atoms[i].nOfn;
                }                
                
                max_nOfn=nOfcols*std::ceil( (float) (max_nOfn+1)/nOfcols);
                
                
        const int nN(max_nOfn);

                total_nOfn.resize(nN);
                for(auto & v: total_nOfn) v=0;

                for(size_t i=0;i<N;i++){
                    total_nOfn.at(grain.atoms[i].nOfn)++;
                    nOffcc+=(size_t) grain.atoms[i].fcc;
                    nOfzb+=(size_t) grain.atoms[i].zb;
                }

                //-----------------------------------
                if(printStat){
                    infoMsg("statistic data :");

                    for(size_t i=0;i<nN;i++){
                        if(i%nOfcols==0){
                            if(i==0)
                                { cout<<endl<<"┏"; for(size_t i=0;i<42;i++) cout<<vline; cout<<"┓\n┃"; }
                            else
                            {     cout<<endl<<"┣"; for(size_t i=0;i<42;i++) cout<<vline; cout<<"┫\n┃"; }
                        }


                        cout<<" "<<setw(3)<<i<<" : "<<setw(8);
                        if(i==nb)
                            cout<<cgreen<<setw(8)<<total_nOfn[i]<<cdef" ┃";
                        else
                            cout<<total_nOfn[i]<<" ┃";
                    }

                    cout<<endl<<"┗"; for(size_t i=0;i<42;i++) cout<<vline; cout<<"┛\n";
                //------------------------------------

                for(int i=0;i<nN;i++)
                    if(i!=nb)  negAtoms+=total_nOfn[i];

                cout<<"\n total number of atoms: "<<N;
                cout<<"\n number of neg. ver. atoms: "<<negAtoms;

                    if(grain.inparams->adistr==StInParams::FCC)
                        cout<<"\n fcc atoms : "<<nOffcc<<endl;
                    if(grain.inparams->adistr==StInParams::ZB)
                        cout<<"\n zb atoms : "<<nOfzb<<endl;

                    if(!grain.inparams->selectedRegion.empty()){
                        cout<<"\n number of atoms by the region";
                        cout<<"\n margin  : "<<grain.count_OBM[2];
                        cout<<"\n bulk  : "<<grain.count_OBM[1];
                        cout<<"\n outside : "<<grain.count_OBM[0];
                    }

                cout<<endl;
                }

                if(!outStatFileName.empty()){
                fstream fout(outStatFileName,ios::out);

                        if(!fout){
                            errMsg(" couldn't save statistic to file "+outStatFileName);
                        }
                        else{
                        std::time_t datetime = std::time(nullptr);
                            fout<<"#ver: 0"<<endl;
                            fout<<"#title: statistic data "<<endl;
                            fout<<"#date: "<<std::asctime(std::localtime(&datetime));
                            fout<<"#input: "<<inFileName<<endl;
                            fout<<"#totNumAtoms: "<<N<<endl;
                            fout<<"#size: 20"<<endl;
                            fout<<"#nv. atoms: "<<negAtoms<<endl;
                            fout<<"#fcc atoms: "<<nOffcc<<endl;
                            fout<<"#zbb atoms: "<<nOfzb<<endl;
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
            saveAtoms(fn.fileName,grain,nb,fn.f_type,EPNF::neg,box,inParams.ignoreKeyValue);


        if(verb && !outFileNames.posVerifiedAtoms.empty() ){ infoMsg("save pv. outcomes");}

        for(auto & fn: outFileNames.posVerifiedAtoms)
            saveAtoms(fn.fileName,grain,nb,fn.f_type,EPNF::pos,box,inParams.ignoreKeyValue);




        for(auto &fn:outFileNames.lattVerifiedAtoms){
            if(verb){
            std::string ltype;
                switch (fn.l_type){
                case EPNF::fcc : ltype="fcc"; break;
                case EPNF::zb  : ltype="zb";break;
                default        : ltype="unk"; continue;
                }
                infoMsg("save "+ltype+" atoms");
            }
            saveAtoms(fn.fileName,grain,nb,fn.f_type,fn.l_type,box,inParams.ignoreKeyValue);
        }


        for(auto &fn:outFileNames.lattNegVerifiedAtoms){
            if(verb){
            std::string ltype;
                switch (fn.l_type){
                case EPNF::nfcc : ltype="non fcc"; break;
                case EPNF::nzb  : ltype="non zb";break;
                default         : ltype="unk";   continue;
                }
                infoMsg("save "+ltype+" atoms");
            }
            saveAtoms(fn.fileName,grain,nb,fn.f_type,fn.l_type,box,inParams.ignoreKeyValue);
        }


        for(auto &fn:outFileNames.selAtomsAndNeighbs){
            saveAtomsAndNeighbors(fn.fileName,grain,fn.nOfnb,fn.f_type);

        }



        //-----------------------------------------------------------------------------------------
        cout<<endl;

return 0;
}
