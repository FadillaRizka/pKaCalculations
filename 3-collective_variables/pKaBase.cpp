/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "pKaBase.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/Matrix.h"
#include "tools/Vector.h"
#include "ActionRegister.h"
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>


using namespace std;

namespace PLMD{
namespace colvar{

void pKaBase::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First Acid/Base group");
  keys.add("atoms","GROUPB","List of Hydrogen Atoms");
  keys.add("atoms", "GROUPC", "Selected oxygen atoms (e.g., solute)");
  keys.add("compulsory","LAMBDA","1","The lambda parameter of the sum_exp function; 0 implies 1");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory", "D_1","0.0", "Charge shift for GROUPC atoms");
  keys.addOutputComponent("sd","default","Acid-base distance order parameter");
  keys.addOutputComponent("tc","default","Total charge order parameter");
}

pKaBase::pKaBase(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista,gc_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);
  parseAtomList("GROUPC",gc_lista);

  list_a = ga_lista;
  list_b = gb_lista;
  list_c = gc_lista;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parse("D_0",d0);
  parse("D_1", d1);
  parse("LAMBDA",lambda);

  bool dopair=false;
  parseFlag("PAIR",dopair);

  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
   if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }
  
  addComponentWithDerivatives("sd"); componentIsNotPeriodic("sd");
  addComponentWithDerivatives("tc"); componentIsNotPeriodic("tc");

  if (gb_lista.size() > 0) {
    if (doneigh)
      nl = Tools::make_unique<NeighborList>(ga_lista, gb_lista, serial, dopair, pbc, getPbc(), comm, nl_cut, nl_st);
    else
      nl = Tools::make_unique<NeighborList>(ga_lista, gb_lista, serial, dopair, pbc, getPbc(), comm);
  } else {
    if (doneigh)
      nl = Tools::make_unique<NeighborList>(ga_lista, serial, pbc, getPbc(), comm, nl_cut, nl_st);
    else
      nl = Tools::make_unique<NeighborList>(ga_lista, serial, pbc, getPbc(), comm);
  }

  std::vector<AtomNumber> atoms;
  atoms.insert(atoms.end(),list_a.begin(),list_a.end());
  atoms.insert(atoms.end(),list_b.begin(),list_b.end());
  requestAtoms(atoms);
  
  for (unsigned i = 0; i < atoms.size(); ++i)
	  idxMap[ atoms[i] ] = i;

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh){
   log.printf("  using neighbor lists with\n");
   log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }

}

pKaBase::~pKaBase(){
  //  delete nl;
}

void pKaBase::prepare(){
  if(nl->getStride()>0){
    if(firsttime || (getStep()%nl->getStride()==0)){
      invalidateList=true;
      firsttime=false;
    }else{
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

void pKaBase::calculate()
{
 double IonDistance=0.0;
 double TotalCharge=0.0;
 double alpha=0.0001;
 double rcut=3.0; 

 unsigned len_acids = list_a.size();
 unsigned len_acids_hyd = len_acids + list_b.size();
 unsigned Nc   = list_c.size();        
 unsigned NaC  = len_acids - Nc; 

 vector<bool> isInCidx(len_acids, false); //mapping group C atoms
 for (auto cNum : list_c)
     isInCidx[idxMap[cNum]] = true;

 vector<unsigned> nonC; //declare nonC vector
 nonC.reserve(len_acids - Nc); //reserve vector nonC (member of A but not C)
 for (unsigned ia = 0; ia < len_acids; ia++) {
     if (!isInCidx[ia])
         nonC.push_back(ia);
 }

 vector<double> sum_exp(len_acids_hyd);
 fill(sum_exp.begin(),sum_exp.end(),0.);

 Tensor virial;
 Tensor virial_dist; 
 std::vector<Vector> deriv_dist(len_acids_hyd);
 std::vector<Vector> deriv_tc(len_acids_hyd);
 Vector zeros;
 zeros.zero();
 fill(deriv_dist.begin(), deriv_dist.end(), zeros);
 fill(deriv_tc.begin(), deriv_tc.end(), zeros);

 if(nl->getStride()>0 && invalidateList){
   nl->update(getPositions());
 }

unsigned nt=OpenMP::getNumThreads();

const unsigned nn=nl->size();

if(nt==0)nt=1;

 std::vector<Vector> omp_deriv_tc(len_acids_hyd);
 fill(omp_deriv_tc.begin(), omp_deriv_tc.end(), zeros);
 Tensor omp_virial;

 vector<vector<Vector>> dist(len_acids, vector<Vector>(len_acids_hyd));
 Matrix<double> distmod(len_acids,len_acids_hyd);

// 1. precompute distances
#pragma omp parallel for
for(unsigned int i=0;i<len_acids;i++) {   
  for(unsigned int j=i+1;j<len_acids;j++) {   //O-O distance
     if(pbc){
        dist[i][j]=pbcDistance(getPosition(i),getPosition(j));
     } else {
        dist[i][j]=delta(getPosition(i),getPosition(j));
     }
     dist[j][i] = -dist[i][j];
     distmod[i][j] = dist[i][j].modulo();
     distmod[j][i] = distmod[i][j];
  }
  for(unsigned int j=len_acids;j<len_acids_hyd;j++) {   //O-H distance
     if(pbc){
        dist[i][j]=pbcDistance(getPosition(i),getPosition(j));
     } else {
        dist[i][j]=delta(getPosition(i),getPosition(j));
     }
     distmod[i][j] = dist[i][j].modulo();
  }
  distmod[i][i]=0.;
}

// 2. calculates sum_exp for charge calcultions
for(unsigned int j=len_acids;j<len_acids_hyd;j++) {   
   double sum_tmp = 0.0;
   #pragma omp parallel for reduction(+:sum_tmp)
   for(unsigned int i=0;i<len_acids;i++) {   
      sum_tmp += exp(lambda * distmod[i][j]);
   }
   sum_exp[j] = sum_tmp;
}

// 3. assignment of d0 and d1 for charge calculations
std::vector<double> d(list_a.size(), d0);
for (std::size_t i = 0; i < list_a.size(); ++i) { //for each atom in group A
    if (std::find(list_c.begin(), list_c.end(), list_a[i]) != list_c.end()) { // if also in group C
        d[i] = d1;
    }
}

// 4. charge calculations
Matrix<double> c(len_acids_hyd,len_acids_hyd);
vector<double> coord(len_acids);
vector<double> charge(len_acids);
fill(coord.begin(),coord.end(),0.);
for(unsigned int i=0;i<len_acids;i++) {   
   double sum_tmp = 0.0;
   #pragma omp parallel for reduction(+:sum_tmp)
   for(unsigned int j=len_acids;j<len_acids_hyd;j++) {   

      c[i][j] = exp( lambda * distmod[i][j] ) / sum_exp[j];
      sum_tmp += c[i][j];

   }
   coord[i] = sum_tmp;
   charge[i] = coord[i] - d[i];
}

// 5. derivation of charge
std::vector<vector<Vector>> ompdfunc_delta(len_acids, vector<Vector>(len_acids_hyd));
std::vector<vector<Vector>> dfunc_delta(len_acids, vector<Vector>(len_acids_hyd));

#pragma omp parallel for 
for(unsigned i=0;i<len_acids;i++) {
  for(unsigned j=0;j<len_acids_hyd;j++) {
    dfunc_delta[i][j]=zeros;
    ompdfunc_delta[i][j]=zeros;
  }
}

double dfunc_coord;
#pragma omp parallel for private(dfunc_coord)
for(unsigned int i=0;i<len_acids;i++) {
  for(unsigned int j=len_acids;j<len_acids_hyd;j++){
    if(distmod[i][j]<rcut) {
      for(unsigned int k=0;k<len_acids;k++) {
        dfunc_coord = -lambda *  c[k][j] * c[i][j];
        ompdfunc_delta[i][k] += dfunc_coord * dist[k][j]/distmod[k][j];
        ompdfunc_delta[i][j] -= dfunc_coord * dist[k][j]/distmod[k][j];
      }
      dfunc_coord = lambda *  c[i][j];
      ompdfunc_delta[i][i] += dfunc_coord * dist[i][j]/distmod[i][j];
      ompdfunc_delta[i][j] -= dfunc_coord * dist[i][j]/distmod[i][j];
    }
  }
}

#pragma omp critical
for(unsigned i=0;i<len_acids;i++) {
  for(unsigned j=0;j<len_acids_hyd;j++) {
    dfunc_delta[i][j]=ompdfunc_delta[i][j]; // gather all results
  }
}


// 6. Calculation of TC (probably not needed) and its derivative
std::vector<double> dfunc_theta(len_acids);

#pragma omp parallel for reduction(+:TotalCharge)
 for(unsigned int i=0;i<len_acids;i++) {
     TotalCharge    += sqrt(pow(charge[i],2)+alpha);
     dfunc_theta[i]  = charge[i]/sqrt(charge[i]*charge[i]+alpha);
 }
TotalCharge -= len_acids * sqrt(alpha);

#pragma omp parallel for
 for(unsigned int m=0;m<len_acids_hyd;m++) {
   for(unsigned int i=0;i<len_acids;i++) {
     omp_deriv_tc[m] -= dfunc_theta[i] * dfunc_delta[i][m];
  }
}

// 7. IonDistance calculations
#pragma omp parallel for reduction(+:IonDistance)
for(auto ia : nonC) {
    for(auto cNum : list_c) {
        unsigned k = idxMap[cNum];
        IonDistance -= distmod[ia][k] * charge[ia] * charge[k];
    }
}
IonDistance *= 1.0 / static_cast<double>(Nc);

// 8. Derivative of IonDistance
std::vector<Vector> omp_deriv_dist(len_acids_hyd);
fill(omp_deriv_dist.begin(), omp_deriv_dist.end(), zeros);
#pragma omp parallel for
for (unsigned m = 0; m < len_acids_hyd; ++m) {
    for (unsigned i = 0; i < NaC; ++i) {
        for (unsigned k = 0; k < Nc; ++k) {
            unsigned idx_k = idxMap[list_c[k]];
                // 1st term
            if (m == i)
               omp_deriv_dist[m] += charge[i] * charge[idx_k] * dist[i][idx_k] / distmod[i][idx_k];
            else if (m == idx_k)
               omp_deriv_dist[m] -= charge[i] * charge[idx_k] * dist[i][idx_k] / distmod[i][idx_k];
                // 2nd and 3rd terms
            omp_deriv_dist[m] += distmod[i][idx_k] * (charge[idx_k] * dfunc_delta[i][m] + charge[i] * dfunc_delta[idx_k][m]);
        }
    }
    omp_deriv_dist[m] *= 1.0 / static_cast<double>(Nc);
}

#pragma omp critical
for(unsigned i=0;i<len_acids_hyd;i++) deriv_dist[i]+=omp_deriv_dist[i];
#pragma omp critical
for(unsigned i=0;i<len_acids_hyd;i++) deriv_tc[i]+=omp_deriv_tc[i];

 Value* vsd = getPntrToComponent("sd");
 Value* vtc = getPntrToComponent("tc");

 for(unsigned i=0;i<deriv_dist.size();++i) setAtomsDerivatives(vsd,i,deriv_dist[i]);
 vsd->set(IonDistance);

 for(unsigned i=0;i<deriv_tc.size();++i) setAtomsDerivatives(vtc,i,deriv_tc[i]);
 vtc->set(TotalCharge);

 }
}
}
