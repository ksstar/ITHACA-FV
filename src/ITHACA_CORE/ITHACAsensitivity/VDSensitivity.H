/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

License
    This file is part of ITHACA-FV

    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef VDSensitivity_H
#define VDSensitivity_H
#include <iostream> 
#include <vector>
#include <string>
#include "Eigen/Eigen" //poi controlla la location
#include "FofM.H"
#include "ITHACAstream.H"
#include "ITHACAutilities.H"
#include "ITHACAsampling.H"


class VDSensitivity
{
	public:
		//Constructors
		VDSensitivity();
		explicit VDSensitivity(int Npara,int Np, std::vector<std::string>& elenco);
		~VDSensitivity() {};

		//Members
		/// Matrix to store the range used in training/offline stage
		Eigen::MatrixXd trainingRange;
		/// Matrices used for sensitivity indeces calculation
		Eigen::MatrixXd Mat1;
		Eigen::MatrixXd Mat2;
		Eigen::MatrixXd Nmat;
		Eigen::MatrixXd Nminus;
		/// Figure of merit object
		FofM M;
		/// Number of parameters used for the analysis
		int No_parameters;
		/// Number of points used for Monte Carlo sampling
		int Npoints;
		/// List of parameters' names
		std::vector<std::string> paramlist;
		/// Vectors to store the output of the model
		Eigen::VectorXd fM1;
		Eigen::VectorXd fM2;
		Eigen::VectorXd	fNmat;
		Eigen::VectorXd fNminus;
		/// Mean value of the output
		double Ey;
		/// Variance of the output
		double Vy;
		/// first order effect index
		double Sj;
		/// total effect index
		double STj;


		//Methods

		///----------------------------------------------------
		/// method to build M1 and M2 matrices, it internally calls samplingMC
		/// pdflist is the list with the names of the distributions of the parameters
		/// plist is a vector ps requested in samplingMC

		void buildSamplingSets(std::vector<std::string>& pdflist, Eigen::MatrixXd plist);

		///----------------------------------------------------
		/// method to export M1 and M2 matrices

		void exportSamplingSets();

		///----------------------------------------------------
		/// method to build Nmat e Nminus used for sensitivity index calculation
		/// param is an int to choose among the parameters adopted for the analysis, it must be less than No_parameters

		void buildNmats(int& param);

		///----------------------------------------------------
		/// method to export Nmat and Nminus matrices

		void exportNmats(int& param);

		///----------------------------------------------------
		/// method to load in fM1,fM2, fNmat and fNminus the M.modelOutput values. M.modelOutput must be already computed i.e. M.MObuilt=true. Each vector mentioned must be loaded separately in place of f
		/// startat is an int to choose from which element of M.modelOutput start at, Npoints values are picked.

		Eigen::VectorXd load_foutput(int& startat);

		///----------------------------------------------------
		/// method to compute the first order effect index of a parameter
		/// parameter is the name of the parameter
		/// x to choose M1 or M2 to compute the variance and the mean value of the output, it can be only 1 or 2
		/// in case the folder corresponding to parameter already exist, or x!=1 & x!=2, -1 is returned

		double getSindex(std::string parameter,int x);

		///----------------------------------------------------
		/// method to compute the total effect index of a parameter
		/// parameter is the name of the parameter
		/// x to choose M1 or M2 to compute the variance and the mean value of the output, it can be only 1 or 2
		/// in case the folder corresponding to parameter already exist, or x!=1 & x!=2, -1 is returned

		double getSTindex(std::string parameter,int x);


	private:

		/// method used in the VDSensitivity object's constructor to set the size of all matrices and vectors of the class
		void setAll();
		
};

#endif

