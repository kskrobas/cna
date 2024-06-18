#include "help.h"

#include<iostream>
using namespace std;


void help()
{

    cerr<<"usage:  cna -i <inputFile> [options]"<<endl;
    cerr<<"    OPTIONS:"<<endl;
    cerr<<" -aafcc - enable fcc angle analysis\n";
    cerr<<" -aazb  - enable zb angle analysis\n";
    cerr<<" -v   - verbose\n";
    cerr<<" -ps  - print statistic\n";
    cerr<<" -psf - print statistic to file\n";
    cerr<<" -d    <real number> - average distance to neighbor \n";
    cerr<<" -tol  <real number> - tolerance of distance \n";
    cerr<<" -tolA <real number> - tolerance of angle \n";
    cerr<<" -nb  <int number>  - number of neighbors\n";
    cerr<<" -op  <fileName>    - output file with positively verified atoms and format nxyz\n";
    cerr<<" -opt <fileName>    - output file with positively verified atoms and format txyz\n";
    cerr<<" -ofcc  <fileName>    - output file with fcc verified atoms and format nxyz\n";
    cerr<<" -ofcct <fileName>    - output file with fcc verified atoms and format txyz\n";
    cerr<<" -ozb  <fileName>    - output file with zbb verified atoms and format nxyz\n";
    cerr<<" -ozbt <fileName>    - output file with zbb verified atoms and format txyz\n";

    cerr<<" -on  <fileName>    - output file with negatively verified atoms and format nxyz\n";
    cerr<<" -ont <fileName>    - output file with negatively verified atoms and format txyz\n";
    cerr<<" -onfcc  <fileName>    - output file with non fcc verified atoms and format nxyz\n";
    cerr<<" -onfcct <fileName>    - output file with non fcc verified atoms and format txyz\n";
    cerr<<" -onzb  <fileName>    - output file with non zbb verified atoms and format nxyz\n";
    cerr<<" -onzbt <fileName>    - output file with non zbb verified atoms and format txyz\n";

    cerr<<" -th  <int number>  - number of OMP threads\n";

    cerr<<"\n example usage:  cna -i atoms37044.xyz -d 2.5 -tol 0.5 -th 8 -opt ppp.xyz -ont nnn.xyz -v -nb 12 -aa \n";

    cerr<<" author/email: Kazimierz.Skrobas@ncbj.gov.pl\n";
    cerr<<" date: "<<__DATE__<<endl;
    cerr<<endl;
}
