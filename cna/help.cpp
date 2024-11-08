/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * help.cpp
* Copyright (C) 2019 Kazimierz Skrobas <kskrobas@unipress.waw.pl>
 *
 * "cna" is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * npcl is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "help.h"

#include<iostream>
using namespace std;


void help()
{

    cerr<<"example usages:  \n cna -i <inputFile> [options]"<<endl;
    cerr<<"\n cna -inp parameters.cna [options]\n\n";
    cerr<<"    command line arguments:"<<endl;
    cerr<<" -aafcc - enable fcc angle analysis\n";
    cerr<<" -aazb  - enable zb angle analysis\n";
    cerr<<" -v   - verbose\n";
    cerr<<" -ps  - print statistic\n";
    cerr<<" -psf - print statistic to file\n";
    cerr<<" -d    <real number> - average distance to neighbor \n";
    cerr<<" -i <fileName> - input file (*.xyz or *.lmp format) with atomic positions\n";
    cerr<<" -inp <scriptName> - alternative for command line options, here arguments are read from script\n";
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

    cerr<<" input script commands: \n";
    cerr<<" infile <fileName> - input file name (equivalent of -i) \n";
    cerr<<" box cube xlo xhi ylo yhi zlo zhi - boundaries of box; if enabled, atoms from box are saved only\n";
    cerr<<" threads <number> - number of threads\n";
    cerr<<" print stat - print statistic\n";
    cerr<<" save (n|nt|nfcc|nfcct|nzb|nzbt) <fileName> - save negatively verified atoms\n";
    cerr<<" save (p|pt|pfcc|pfcct|pzb|pzbt) <fileName> - save positively verified atoms\n";
    cerr<<" toldist <number> - tolerance for  distance between atomic positions\n";
    cerr<<" tolangle <number> - tolerance for angle between atomic positions\n";
    cerr<<" nbh - number of neighbors\n";
    cerr<<" dist <number> - distance to the 1st neighbor\n";

    cerr<<"WARNING: using both command line and script instructions may cause unpredictable behavior\n";
    cerr<<" author/email: Kazimierz.Skrobas@ncbj.gov.pl\n";
    cerr<<" date: "<<__DATE__<<endl;
    cerr<<endl;
}
