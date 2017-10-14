// #include <Rcpp.h>
// #include <iostream>
// #include <fstream>
// #include <eigen3/Eigen/Dense>
#include "parent.h"
#include "function_sup.h"
#include <random>
#include <stdint.h>
#include <RcppEigen.h>

// using namespace Eigen;
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List ibm(const int nbInd, const int length, const double prob_mut,const int loci,
               const Eigen::MatrixXf & CovarMatrix,
               bool pursue, const Eigen::MatrixXf ParentGenetValue,
               bool intratrait,const Eigen::MatrixXf matrixEpistasisIntra,
               bool intertraits,const Eigen::MatrixXf matrixEpistasisInter){





const int damDesign =4;
const int offspringDesign=4;
const int sireDesign=nbInd;
const int traits=2;
const int alleles=2;
const int lociXalleles=loci*alleles ;
const int offspringTotal=sireDesign*damDesign*offspringDesign;
const int maleMeanSlice=damDesign*offspringDesign;
const int femaleMeanNumber=sireDesign*damDesign;
const int numberVarCovar=((traits*traits)-traits)/2+traits;

//////////////////////////////////////////////////////////
// Initialisation of Rmat used to generate multivariate normal //
//////////////////////////////////////////////////////////

SelfAdjointEigenSolver<MatrixXf> eigensolver(CovarMatrix);
Eigen::MatrixXf eigenValuesMat=  eigensolver.eigenvalues().array().sqrt().replicate<1,traits>();
Eigen::MatrixXf RightPartRmat=eigensolver.eigenvectors().transpose().array()*eigenValuesMat.array();
Eigen::MatrixXf allPartRmat = eigensolver.eigenvectors()* RightPartRmat;
Eigen::MatrixXf Rmat = allPartRmat.transpose();

// // // // // //////////////////////////////////////////////////////////////
// // // Initialisation of Epistatic matrix (in a purpose of testing function //
// // // // // // // // // // // // // // // // // // // // // // // // // //

// MatrixXf matrixEpistasis =matrixEpis(loci);

MatrixXf matEpistasisIntra =matrixEpistasisIntra;
MatrixXf matEpistasisInter =matrixEpistasisInter;


//////////////////////////////////////////////////////////////////////////////////////
//////////////// Create the base population of parents and offsprings ///////////////
////////////////  OR get the genetic value from another run /////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

std::vector<parent> p;
for(int i = 0; i < nbInd; ++i)  p.push_back(parent(lociXalleles, traits));

if(pursue){
  int  storeRowPass=0;
  Eigen::MatrixXf passGenet(lociXalleles,traits);
  for (int pass=0;pass<nbInd;++pass){
          passGenet=ParentGenetValue.block(storeRowPass,0,lociXalleles,traits);
          p[pass].setGenetMatrix(passGenet);
          storeRowPass+=lociXalleles;
  }
}


std::vector<parent> o;
for(int i = 0; i < offspringTotal; ++i) o.push_back(parent(lociXalleles, traits));


////////////////////////////////////////////////////////////////
////////////////// Initialise G-matrix collection //////////////
////////////////////////////////////////////////////////////////

Eigen::MatrixXf Gmat(length,traits*traits);
Eigen::MatrixXf Gmat2(length,numberVarCovar);
Eigen::MatrixXf OffspringsValueMeanEpi(length,traits);
Eigen::MatrixXf OffspringsValueMeanNoEpi(length,traits);


std::random_device rd;
std::mt19937 generator(rd());
std::uniform_int_distribution<int> uni_distribution(0, sireDesign-1 );

std::vector<int> femaleRepro;
for (int k= 0; k < nbInd; ++k) femaleRepro.push_back(k);
std::vector<int> randomoffspring;
for (int k= 0; k < offspringTotal; ++k) randomoffspring.push_back(k);


///////////START OF CREATING POPULATION GENERATION

for(int i = 0; i < length; ++i) {

Eigen::MatrixXf offspringGenetValue(offspringTotal,traits);
Eigen::MatrixXf offspringGenetValueKept(nbInd,traits);
Eigen::MatrixXf offspringGenetValue4sum(nbInd,traits);

std::vector<int> maleReproducers;
std::vector<int> femaleReproducers;
std::vector<int> femaleForMale;
femaleForMale.reserve(damDesign*offspringDesign);
std::vector<int> femaleForReproduction;
femaleForReproduction.reserve(offspringTotal);
std::vector<int> maleForReproduction;
maleForReproduction.reserve(offspringTotal);


///// Find the male & female reproductors

for (int f=0;f<sireDesign;++f) {
    maleReproducers.push_back(uni_distribution(generator));
    femaleForMale=RandSamplWithoutReplac(femaleRepro,damDesign,maleReproducers[f],offspringDesign);
    femaleForReproduction.insert(std::end(femaleForReproduction), std::begin(femaleForMale), std::end(femaleForMale));
    maleForReproduction.insert(maleForReproduction.end(), damDesign*offspringDesign, maleReproducers[f]);
}

///reproduction recombination

for (int re=0;re<offspringTotal;++re){
    o[re].reproduction_recombinaison(traits,loci,p[maleForReproduction[re]],p[femaleForReproduction[re]]);
}

///// mutation , get genetic value, swap indivi

if(!intratrait && !intertraits){
  for (int mu=0;mu<offspringTotal;++mu){
    o[mu].mutation(prob_mut,traits,lociXalleles,Rmat);
    offspringGenetValue.row(mu)=o[mu].getGeneticValue();
  }
}else if(intratrait && !intertraits){
  for (int mu=0;mu<offspringTotal;++mu){
    o[mu].mutation(prob_mut,traits,lociXalleles,Rmat);
    offspringGenetValue.row(mu)=o[mu].getGeneticValue()+o[mu].getGeneticValueEpistasisIntratrait(loci,traits,matEpistasisIntra);
  }
}else if(intratrait && intertraits){
  for (int mu=0;mu<offspringTotal;++mu){
    o[mu].mutation(prob_mut,traits,lociXalleles,Rmat);
    offspringGenetValue.row(mu)=o[mu].getGeneticValue()+o[mu].getGeneticValueEpistasisIntratraitIntertraits(loci,traits,matEpistasisIntra,matEpistasisInter);
  }
}


///// get random offsprings
int co=0;
for (int next=0;next<nbInd;++next){
    offspringGenetValueKept.row(next)=offspringGenetValue.row(co);
    offspringGenetValue4sum.row(next)=o[co].getGeneticValue();
    p[next].swapIndividual(o[co]);
    co+=maleMeanSlice;
}

MatrixXf centered = offspringGenetValue4sum.rowwise() - offspringGenetValue4sum.colwise().mean();
MatrixXf cov = (centered.adjoint() * centered) / double(offspringGenetValue4sum.rows() - 1);
Map<RowVectorXf> VectorG(cov.data(), cov.size());
Gmat.row(i)=VectorG;
Gmat2.row(i)=GmatrixBreedDesign(maleMeanSlice,femaleMeanNumber,numberVarCovar,offspringTotal,damDesign,offspringDesign,sireDesign,offspringGenetValue);

OffspringsValueMeanEpi.row(i)=offspringGenetValueKept.colwise().mean();
OffspringsValueMeanNoEpi.row(i)=offspringGenetValue4sum.colwise().mean();
}

///////// STORE GENETIQUE VALUE  OF ALL PARENTS In ONE MATRIX for the last iteration ///////////
int  storeRow=0;
Eigen::MatrixXf parent_genet(lociXalleles*nbInd,traits);
  for (int con=0;con<nbInd;++con){
    parent_genet.block(storeRow,0,lociXalleles,traits)=p[con].getGenet();
    storeRow+=lociXalleles;
}

return Rcpp::List::create(
                          Rcpp::Named("epimatInter")=matrixEpistasisInter,
                          Rcpp::Named("OffspringsValueMeanEpi")=OffspringsValueMeanEpi,
                          Rcpp::Named("OffspringsValueMeanNoEpi")=OffspringsValueMeanNoEpi,
                          Rcpp::Named("Gmat_bd")=Gmat2,
                          Rcpp::Named("Gmat_sum")=Gmat,
                          Rcpp::Named("Parents_geneticValues")=parent_genet
                          );

}
