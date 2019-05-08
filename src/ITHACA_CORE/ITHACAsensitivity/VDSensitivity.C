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

#include "VDSensitivity.H"


using std::vector; using std::string;

// Constructors
VDSensitivity::VDSensitivity() {};

VDSensitivity::VDSensitivity(int Npara,int Np, std::vector<std::string>& elenco)
{
	No_parameters=Npara;
	paramlist=elenco;
	Npoints=Np;	
	setAll();
}




void VDSensitivity::buildSamplingSets(std::vector<std::string>& pdflist,Eigen::MatrixXd plist)
{
	for (int i=0; i<No_parameters; i++)
	{	
		Mat1.col(i)=ITHACAsampling::samplingMC(pdflist[i],trainingRange(i,0),trainingRange(i,1), plist(i,0), plist(i,1), Npoints);
		Mat2.col(i)=ITHACAsampling::samplingMC(pdflist[i],trainingRange(i,0),trainingRange(i,1), plist(i,0), plist(i,1), Npoints);
        }

}
void VDSensitivity::exportSamplingSets()
{
		
	ITHACAstream::exportMatrix(Mat1,"Mat1","eigen","./ITHACAoutput/Sensitivity/SamplingMC/");
	ITHACAstream::exportMatrix(Mat2,"Mat2","eigen","./ITHACAoutput/Sensitivity/SamplingMC/");

}

void VDSensitivity::buildNmats(int& param)
{
	if (param>No_parameters)
	{
		std::cout<<"Number of parameters is: "<<No_parameters<<", you asked for param. #: "<<param<<", ERROR, program aborted"<<std::endl;
		exit(0);
	}
	else
	{
		Nmat=Mat2;
		Nmat.col(param-1)=Mat1.col(param-1);
		Nminus=Mat1; 
		Nminus.col(param-1)=Mat2.col(param-1);
	}
}

void VDSensitivity::exportNmats(int& param)
{
	word name={"Nmat_"}; 
	name.append(std::to_string(param));
	ITHACAstream::exportMatrix(Nmat,name,"eigen","./ITHACAoutput/Sensitivity/SamplingMC/");
	name={"Nminus_"};
	name.append(std::to_string(param));
	ITHACAstream::exportMatrix(Nminus,name,"eigen","./ITHACAoutput/Sensitivity/SamplingMC/");
}

double VDSensitivity::getSindex(std::string parameter,int x)
{
	word name; 
	name.append(M.name);	
	word folder="./ITHACAoutput/Sensitivity/VDI/";
	folder.append(name); 
	folder.append("/");
	folder.append("FirstOrderEffect/");
	Eigen::VectorXd values;
	values.resize(Npoints);
	values.setZero();
	double Uj=0;
	if (x==1){values=fM1;}
	else if (x==2){values=fM2;}
	else 
	{	
		std::cout<<"Option not available, x can be only 1 or 2,returned -1"<<std::endl;
		return -1;
	}

	if (ITHACAutilities::check_folder(folder)==false)
	{
		Ey=values.sum()/Npoints;
		Vy=0; 
		for (int i=0; i<Npoints; i++)
		{
			Vy+=values(i)*values(i)-(Ey*Ey);
			Uj+=fM1(i)*fNmat(i);
		}
		Vy=Vy/(Npoints-1);
		Uj=Uj/(Npoints-1);
		Sj=(Uj-(Ey*Ey))/Vy;
		Eigen::MatrixXd tmp;
		tmp.resize(1,1);
		tmp(0,0)=Sj;
		word index={"S_"};
		index.append(parameter);
		ITHACAstream::exportMatrix(tmp,index,"eigen",folder);
		return Sj;
	}
	else
	{
		std::cout<<"Variance decomposition index of the couple "<<M.name<<" "<<parameter<<" was already computed, returned -1"<<std::endl;
		return -1;
	} 
}

double VDSensitivity::getSTindex(std::string parameter,int x)
{
	word name; 
	name.append(M.name);	
	word folder="./ITHACAoutput/Sensitivity/VDI/";
	folder.append(name); 
	folder.append("/");
	folder.append("TotalEffect/");
	Eigen::VectorXd values;
	values.resize(Npoints);
	values.setZero();
	double U_j=0;
	if (x==1){values=fM1;}
	else if (x==2){values=fM2;}
	else 
	{	
		std::cout<<"Option not available, x can be only 1 or 2,returned -1"<<std::endl;
		return -1;
	}

	if (ITHACAutilities::check_folder(folder)==false)
	{
		Ey=values.sum()/Npoints;
		Vy=0; 
		for (int i=0; i<Npoints; i++)
		{
			Vy+=values(i)*values(i)-Ey*Ey;
			U_j+=fM1(i)*fNminus(i);
		}
		Vy=Vy/(Npoints-1);
		U_j=U_j/(Npoints-1);
		STj=1-(U_j-Ey*Ey)/Vy;
		Eigen::MatrixXd tmp;
		tmp.resize(1,1);
		tmp(0,0)=STj;
		word index={"ST_"};
		index.append(parameter);
		ITHACAstream::exportMatrix(tmp,index,"eigen",folder);
		return STj;
	}
	else
	{
		std::cout<<"Variance decomposition index of the couple "<<M.name<<" "<<parameter<<" was already computed, returned -1"<<std::endl;
		return -1;
	} 
}

void VDSensitivity::setAll()
{

	Mat1.resize(Npoints,No_parameters);
	Mat1.setZero();
	Mat2.resize(Npoints,No_parameters);
	Mat2.setZero();
	Nmat.resize(Npoints,No_parameters);
	Nmat.setZero();
	Nminus.resize(Npoints,No_parameters);
	Nminus.setZero();
	fM1.resize(Npoints);
	fM1.setZero();
	fM2.resize(Npoints); 
	fM2.setZero();
	fNmat.resize(Npoints);	
	fNmat.setZero();
	fNminus.resize(Npoints);
	fNminus.setZero();
	trainingRange.resize(No_parameters,2);
	trainingRange.setZero();
	
}

Eigen::VectorXd VDSensitivity::load_foutput( int& startat)
{
	Eigen::VectorXd f;
	f.resize(Npoints);	
	if (M.MObuilt==true)
	{
		f=M.modelOutput.segment(startat,Npoints);
	}
	else{std::cout<<"Model output not computed or assigned, program aborted"<<std::endl; exit(0);}

	return f;

}

