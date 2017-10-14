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

#ifndef FUNCTION_SUP_H
#define FUNCTION_SUP_H

#include <random>
#include <stdint.h>
#include <RcppEigen.h>
// #include <eigen3/Eigen/Dense>
#include <vector>


using namespace std;

std::vector<int> RandSamplWithoutReplac(std::vector<int>& vData, int numberDraws, int omit, int numberReplication);

Eigen::VectorXf GmatrixBreedDesign(int maleMeanSlice,int femaleMeanNumber, int numberVarCovar, int offspringTotal, int damDesign, int offspringDesign, int sireDesign, Eigen::MatrixXf& offspringGeneticValue);

Eigen::MatrixXf matrixEpis(int loci);


//uint64_t xorshift128plus(void);

#endif // FUNCTION_SUP_H
