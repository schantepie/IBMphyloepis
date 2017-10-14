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

#include "function_sup.h"
// #include <eigen3/Eigen/Dense>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <RcppEigen.h>
#include <random>

using namespace Eigen;
using namespace std;



////////// Fast pseudo random generator xorshift128+
//
//random_device rd;
//
///* The state must be seeded so that it is not everywhere zero. */
//uint64_t s[2] = { (uint64_t(rd()) << 32) ^ (rd()),(uint64_t(rd()) << 32) ^ (rd()) };
//
//uint64_t xorshift128plus(void) {
//	uint64_t x = s[0];
//	uint64_t const y = s[1];
//	s[0] = y;
//	x ^= x << 23; // a
//	s[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
//	return s[1] + y;
//}



////////// Random female without selfing
//std::pair<std::vector<int>,std::vector<int>> RandSamplWithoutReplac(std::vector<int>& vData, int numberDraws, int omit, int numberReplication){
std::vector<int> RandSamplWithoutReplac(std::vector<int>& vData, int numberDraws, int omit, int numberReplication){
std::random_device rd;
std::mt19937 generator(rd());
std::vector<int> vData_tmp;
vData_tmp=vData;
vData_tmp.erase(vData_tmp.begin()+omit);
std::shuffle(vData_tmp.begin(),vData_tmp.end(),generator);
vData_tmp.resize(numberDraws);
std::vector<int> sample;
sample.reserve(numberDraws*numberReplication);

if (numberReplication>1){
        for (int tot= 0; tot < numberReplication; ++tot) sample.insert(sample.end(),numberReplication,vData_tmp[tot]);
//        return std::make_pair(sample, vData_tmp);
        return sample;

}else{
        sample=vData_tmp;
//        return std::make_pair(sample, sample);
        return sample;

}
}



Eigen::VectorXf GmatrixBreedDesign(int maleMeanSlice,int femaleMeanNumber, int numberVarCovar, int offspringTotal, int damDesign, int offspringDesign, int nbInd, Eigen::MatrixXf& offspringGeneticValue){
Eigen::MatrixXf meanMale(nbInd,numberVarCovar);
Eigen::MatrixXf meanFemaleMinusMeanMale(femaleMeanNumber,numberVarCovar);
Eigen::MatrixXf offspringCorrelation=offspringGeneticValue.rowwise().sum(); //only work for 2 traits
Eigen::MatrixXf offspringGeneticValueAndCorrel(offspringTotal,numberVarCovar);
offspringGeneticValueAndCorrel<<offspringGeneticValue,offspringCorrelation;
int increMale=0;
int increFemale=0;
int comptFemale=0;
    for (int i=0; i<nbInd; ++i){//nbInd or nbMale according design
        meanMale.row(i)=offspringGeneticValueAndCorrel.block(increMale,0,maleMeanSlice,numberVarCovar).colwise().mean();
        increMale+=maleMeanSlice;
         for (int j=0; j<damDesign; ++j){//nbInd or nbMale according design
            meanFemaleMinusMeanMale.row(comptFemale)=offspringGeneticValueAndCorrel.block(increFemale,0,offspringDesign,numberVarCovar).colwise().mean()-meanMale.row(i);
            increFemale+=offspringDesign;
            ++comptFemale;
        }

    }
Eigen::MatrixXf Error=meanMale.rowwise()-offspringGeneticValueAndCorrel.colwise().mean();
Eigen::VectorXf MSs=(Error.array().pow(2).colwise().sum()*damDesign*offspringDesign)/(nbInd-1);
Eigen::VectorXf MSd=(meanFemaleMinusMeanMale.array().pow(2).colwise().sum()*offspringDesign)/(nbInd*(damDesign-1));
Eigen::VectorXf Va=4*(MSs-MSd)/(damDesign*offspringDesign);
Va[2]=(Va[2]-Va[1]-Va[0])/2;//only work for 2 traits
return Va;
}





Eigen::MatrixXf matrixEpis(int loci){
std::random_device rd;
std::mt19937 generator(rd());
std::normal_distribution<double> normsample(0,0.1);
Eigen::MatrixXf matEpi(loci,loci);
        for (int i=0; i<loci; ++i){
                for (int j=0; j<loci; ++j){
                        matEpi(i,j)=normsample(generator);
                }
        }
        matEpi.diagonal()=VectorXf::Zero(loci);
return matEpi;
}

