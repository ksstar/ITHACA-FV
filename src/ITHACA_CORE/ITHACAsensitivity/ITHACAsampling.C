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

#include "ITHACAsampling.H"

std::vector<std::string> ITHACAsampling::distributions={"UNIFORM","NORMAL","POISSON","EXPONENTIAL"};

Eigen::VectorXd ITHACAsampling::samplingMC(std::string pdftype,double& lowerE, double& upperE, double& distpara1,double& distpara2, int& Npoints)
{       
	std::random_device rd;
	std::mt19937 generator(rd());
		
	//to make it non-case sensitive
        for (int i=0; i< pdftype.size(); i++)
	{
		pdftype[i]=toupper(pdftype[i]);
	}

        int pos=-1;
	bool pdf_found=false;
	Eigen::VectorXd samplingVector;
        samplingVector.resize(Npoints);
	std::uniform_real_distribution<> dist0(distpara1,distpara2);
	std::normal_distribution<> dist1(distpara1,distpara2);
	std::poisson_distribution<> dist2(distpara1);
	std::exponential_distribution<> dist3(distpara1);

	for (int j=0; j<distributions.size(); j++)
	{
		if(pdftype==distributions[j])
		{
			pdf_found=true;
                        pos=j;
		}
	}
	double random;
        int p_counter=0;
	if (pdf_found==true)
	{
		while (p_counter<samplingVector.size())
			{
			  switch(pos)
			  {
		            case 0:
				 random=dist0(generator);
			    break;
			   
			    case 1:
				 random=dist1(generator);
			    break;
			   
			    case 2:
				 random=dist2(generator);
			    break;
			   
			    case 3:
				 random=dist3(generator);
			    break;
			  }
			  if (random>=lowerE && random<=upperE) 
			  {
			    samplingVector(p_counter)=random;
			    p_counter++;
			  }
			}
	}
	else
	{
		std::cout<<"pdf '"<<pdftype<<"' not implemented, programm aborted"<<std::endl;
		exit(0);
	}

	return samplingVector;
}

