/*
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2017  stephane <email>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "parent.h"
#include <iostream>
#include <random>
#include <ctime>
#include "function_sup.h"
// #include <eigen3/Eigen/Dense>
#include <RcppEigen.h>
#include <fstream>



using namespace std;

parent::parent(int lociXalleles, int traits)
{
    m_genet=MatrixXf::Zero(lociXalleles,traits);
    geneticValue=VectorXf(traits);
}

parent::~parent(){}


void parent::reproduction_recombinaison(int numberTraits, int loci, parent &sire, parent &dam){
std::random_device rd;
std::mt19937 generator(rd());
std::bernoulli_distribution randomAllele(0.5);

Eigen::MatrixXf m_genetSire=sire.getGenet();
Eigen::MatrixXf m_genetDam=dam.getGenet();

   for (int i=0;i<loci;++i){
            if (randomAllele(generator)){m_genetSire.row(i)=m_genetSire.row(i+loci);}
            if (randomAllele(generator)){m_genetDam.row(i)=m_genetDam.row(i+loci);}
    }
   m_genet.block(0,0,loci,numberTraits)=m_genetSire.block(0,0,loci,numberTraits);
   m_genet.block(loci,0,loci,numberTraits)=m_genetDam.block(0,0,loci,numberTraits);

}


void parent::mutation(double probaMutation,int numberTraits, int lociXalleles, Eigen::MatrixXf& Rmat){
std::random_device rd;
std::mt19937 generator(rd());
std::bernoulli_distribution distributionMutations(probaMutation);
std::normal_distribution<double> normsample(0,1);

Eigen::MatrixXf Rgener=Eigen::MatrixXf::Zero(lociXalleles,numberTraits);

int countMut=0;
for (int i=0;i<lociXalleles;++i){
        if(distributionMutations(generator)){
            ++countMut;
            for (int j=0;j<numberTraits;++j){
                Rgener(i,j)= normsample(generator);
            }
        }
}

if(countMut!=0){
    Eigen::MatrixXf addMutations=Rgener*Rmat;
    m_genet+=addMutations;
}

}


void parent::setGenetMatrix(Eigen::MatrixXf& genetVal){
  m_genet=genetVal;
}

void parent::swapIndividual(parent &offspring)
{
    m_genet=offspring.getGenet();
    geneticValue=offspring.getGeneticValue();
}


Eigen::VectorXf parent::getGeneticValueEpistasisIntratrait(int loci, int numberTraits, Eigen::MatrixXf& Epis)
{
  Eigen::MatrixXf geneticValueLoci=m_genet.block(0,0,loci,numberTraits)+m_genet.block(loci,0,loci,numberTraits);
  Eigen::VectorXf geneticValueEpi(numberTraits);
  Eigen::VectorXf phenoValueEpi(numberTraits);
  for (int i=0;i<numberTraits;i++){
  Eigen::MatrixXf lociMat=geneticValueLoci.col(i)*geneticValueLoci.col(i).transpose();
  Eigen::MatrixXf geneticMatrixEpiTrait=lociMat.array()*Epis.array();
  geneticValueEpi(i)=geneticMatrixEpiTrait.sum();
  }
  phenoValueEpi=geneticValueEpi;
  return phenoValueEpi;
}



Eigen::VectorXf parent::getGeneticValueEpistasisIntratraitIntertraits(int loci, int numberTraits, Eigen::MatrixXf& Episintra,Eigen::MatrixXf& Episinter)
{
  Eigen::MatrixXf geneticValueLoci=m_genet.block(0,0,loci,numberTraits)+m_genet.block(loci,0,loci,numberTraits);
  Eigen::VectorXf geneticValueEpiIntra(numberTraits);
  Eigen::VectorXf geneticValueEpiInter(numberTraits);
  Eigen::VectorXf phenoValueEpi(numberTraits);
  for (int i=0;i<numberTraits;i++){
    Eigen::MatrixXf lociMatintra=geneticValueLoci.col(i)*geneticValueLoci.col(i).transpose();
    Eigen::MatrixXf geneticMatrixEpiTraitintra=lociMatintra.array()*Episintra.array();
    geneticValueEpiIntra(i)=geneticMatrixEpiTraitintra.sum();

    for (int j=i+2;j<numberTraits;j++){

      Eigen::MatrixXf lociMatinter=geneticValueLoci.col(i)*geneticValueLoci.col(j).transpose();
      Eigen::MatrixXf geneticMatrixEpiTraitinter=lociMatinter.array()*Episinter.array();
      geneticValueEpiInter(i)=geneticMatrixEpiTraitinter.sum();

    }
  }
  phenoValueEpi= geneticValueEpiIntra+geneticValueEpiInter;
  return phenoValueEpi;
}


Eigen::MatrixXf& parent::getGenet()
{
    return m_genet;
}

Eigen::VectorXf& parent::getGeneticValue()
{
    geneticValue=m_genet.colwise().sum();
    return geneticValue;
}

void parent::afficherEtat()
{
    cout << "m_genet : \n" <<  m_genet << endl;
    cout << "genval" << geneticValue << endl;
}



