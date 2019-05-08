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
#include "Tm.H"

Tm::Tm() {}

Tm::Tm(int argc, char *argv[], int Nsampled)
{
	Npoints=Nsampled;
	modelOutput.resize(Npoints);
	modelOutput.setZero();	
	name={"Tm"};	
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"	
	#include "createT.H"
}

void Tm::buildMO(std::string dir)
{
	Time& runTime = _runTime();
	fvMesh& mesh = _mesh();
	volScalarField& T=_T();
	std::string folder=dir;

	if (ITHACAutilities::check_folder(folder)==true)
	{
		for (int j=0; j<Npoints; j++)
		{
			folder.append(std::to_string(j));
			folder.append("/");
			ITHACAstream::read_fields(ptrfield,T,folder); 
			auto k=ptrfield.last();
			modelOutput(j)=k.weightedAverage(mesh.V()).value();
			folder=dir;	
			ptrfield.clear();			
		}
		MObuilt=true;
	} 
		

	else
	{
		std::cout<<"Outputs of the model are not computed yet, programm aborted"<<std::endl;
		exit(0);		
	}
}

