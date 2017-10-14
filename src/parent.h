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

#ifndef PARENT_H
#define PARENT_H
#include <iostream>
#include <vector>
// #include <eigen3/Eigen/Dense>
#include <RcppEigen.h>
#include <random>
using namespace Eigen;
using namespace std;

class parent
{
    public:

     parent(int lociXalleles, int traits);
     ~parent();

    void mutation(double probaMutation, int numberTraits, int lociXalleles, Eigen::MatrixXf& Rmat);

    void reproduction_recombinaison(int numberTraits, int loci, parent &sire, parent &dam);

    void setGenetMatrix(Eigen::MatrixXf& genetVal);

    void swapIndividual(parent &offspring);

    void afficherEtat();

    Eigen::VectorXf getGeneticValueEpistasisIntratrait (int loci, int numberTraits, Eigen::MatrixXf& Epis);

    Eigen::VectorXf getGeneticValueEpistasisIntratraitIntertraits (int loci, int numberTraits, Eigen::MatrixXf& Episintra,Eigen::MatrixXf& Episinter);

    Eigen::MatrixXf& getGenet();

    Eigen::VectorXf& getGeneticValue();

    private:

    Eigen::MatrixXf m_genet;
    Eigen::VectorXf geneticValue;
};

#endif // PARENT_H
