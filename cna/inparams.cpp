
#include "inparams.h"
#include "colormsg.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <regex>


using namespace std;
typedef vector<string> vstring;


#define RE_NUMBER_0  "[[:s:]]+[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define RE_NUMBER_1  "[[:s:]]*[+-]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define PRE_NUMBER_0 "[[:s:]]+[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define PRE_NUMBER_1 "[[:s:]]*[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?"
#define PRE_NUMBER_2 "[[:s:]]*[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?[[:s:]]*"
#define PROB_NUMBER "[[:s:]]+(0|1[.]?|0[.][0-9]*)"
#define UINT_NUMBER  "[[:s:]]+[0-9]+"
#define VAR          "\\$\\{\\w+\\}"

const std::string sRE_NUMBER(RE_NUMBER_0);    //real number
const std::string sRE_NUMBER_1(RE_NUMBER_1);    //real number
const std::string sPRE_NUMBER(PRE_NUMBER_0); //positive, real number
const std::string sPRE_NUMBER_1(PRE_NUMBER_1); //positive, real number
const std::string sPRE_NUMBER_2(PRE_NUMBER_2); //positive, real number
const std::string sPROB_NUMBER(PROB_NUMBER);  //positive number for probabilty purposes from (0...1) range
const std::string sUINT_NUMBER(UINT_NUMBER);
const std::string sVAR(VAR);


bool parseInputScript(const string &fileName, StInParams &sparams, bool &verb)
{
fstream fin(fileName,ios::in);

        if(!fin){
            errMsg("couldn't open "+fileName);
        return false;
        }


string fline;
size_t flineNr=0;

        try{
            while( !fin.eof() ){
                std::getline(fin,fline);
                trim(fline);
                flineNr++;
                if(fline.empty() || fline[0]=='#') continue;

                if(verb) cout<<setw(4)<<" "<<flineNr<<": "<<fline<<endl;


                // replace tabs by spaces
                std::replace(std::begin(fline),std::end(fline),'\t',' ');


                if(regex_match(fline,std::regex("avedist"+sPRE_NUMBER))){
                vstring toks{split<string>(fline," ")};

                        *sparams.dst=std::stod(toks[1]);
                continue;
                }


                if(regex_match(fline,std::regex("box[[:s:]]+cuboid("+sRE_NUMBER_1+"){6}"))){
                vstring toks{split<string>(fline," ")};
                size_t i=0,j=2;

                        for(;i<6;i++,j++)
                            sparams.box->bounds[i]=std::stod(toks[j]);

                        sparams.box->btype=StBox::CUB;
                continue;
                }


                if(regex_match(fline,std::regex("ifile[[:s:]]+[[:print:]]+"))){
                vstring toks{split<string>(fline," ")};

                        *sparams.inFileName=toks[1];
                continue;
                }

                if(regex_match(fline,std::regex("nosave[[:s:]]+atype[[:s:]]+[[:print:]]+"))){
                const std::string keyValue{fline.substr(fline.find("atype"))};

                        sparams.ignoreKeyValue.emplace_back(keyValue);
                continue;
                }


                if(regex_match(fline,std::regex("nosave[[:s:]]+nb[[:s:]]+[0-9]+"))){
                const std::string keyValue{fline.substr(fline.find("nb"))};

                        sparams.ignoreKeyValue.emplace_back(keyValue);
                continue;
                }

                if(regex_match(fline,std::regex("mode[[:s:]]+(fcc|zb)"))){
                vstring toks{split<string>(fline," ")};

                        sparams.grain->AAEnabled=(toks[1]=="fcc");
                continue;
                }

                if(regex_match(fline,std::regex("nbnum"+sUINT_NUMBER))){
                vstring toks{split<string>(fline," ")};

                        *sparams.nb=std::stoi(toks[1]);
                continue;
                }

                if(regex_match(fline,std::regex("print[[:s:]]+stat"))){
                        *sparams.printStat=true;
                continue;
                }

                if(regex_match(fline,std::regex("save[[:s:]]+n(t|fcc|fcct|zb|zbt)?[[:s:]]+[[:print:]]+"))){
                vstring toks{split<string>(fline," ")};
                const string _mode_{toks[1]};
                StFileNameType fnt;

                        fnt.fileName=toks[2];

                        if(_mode_=="n"){  // all, negatively verified atoms
                            fnt(EFTYPE::nxyz,EPNF::neg);
                            sparams.outFileNames->negVerifiedAtoms.emplace_back(fnt);
                        }
                        else{
                            if(_mode_=="nfcc")
                                fnt(EFTYPE::nxyz,EPNF::nfcc);
                            else{
                                if(_mode_=="nfcct")
                                    fnt(EFTYPE::txyz,EPNF::nfcc);
                            }
                            sparams.outFileNames->lattNegVerifiedAtoms.emplace_back(fnt);
                        }

                continue;
                }


                if(regex_match(fline,std::regex("save[[:s:]]+p(t|fcc|fcct|zb|zbt)?[[:s:]]+[[:print:]]+"))){
                vstring toks{split<string>(fline," ")};
                const string _mode_{toks[1]};
                StFileNameType fnt;

                        fnt.fileName=toks[2];

                        if(_mode_=="p"){  // all, negatively verified atoms
                            fnt(EFTYPE::nxyz,EPNF::pos);
                            sparams.outFileNames->posVerifiedAtoms.emplace_back(fnt);
                        }
                        else{
                            if(_mode_=="pfcc")
                                fnt(EFTYPE::nxyz,EPNF::fcc);
                            else{
                                if(_mode_=="pfcct")
                                    fnt(EFTYPE::txyz,EPNF::fcc);
                            }
                            sparams.outFileNames->lattVerifiedAtoms.emplace_back(fnt);
                        }

                continue;
                }



                if(regex_match(fline,std::regex("threads"+sUINT_NUMBER))){
                vstring toks{split<string>(fline," ")};

                        sparams.grain->threads=std::stoi(toks[1]);

                continue;
                }


                if(regex_match(fline,std::regex("toldist"+sPRE_NUMBER))){
                vstring toks{split<string>(fline," ")};

                        *sparams.tol=std::stod(toks[1]);
                continue;
                }


                if(regex_match(fline,std::regex("tolang"+sPRE_NUMBER))){
                vstring toks{split<string>(fline," ")};

                        *sparams.tolA=std::stod(toks[1]);
                continue;
                }



                errMsg("  unknown last command or its argument(s) ");
                throw 1;


            }
        }
        catch(int i){
                fin.close();
                errMsg(" line 206, exception, code: "+std::to_string(i));
                return false;
        }
        catch(...){
                fin.close();
                errMsg(" line 211, exception, code: unknown ");
                return false;
        }

return true;
}

