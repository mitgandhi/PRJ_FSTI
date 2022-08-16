#include "CPistonGap.h"
#include "../../caspar_input\input.h"
#include "logger.h"
#pragma once

extern class input myinput;
extern int hd_levels;


bool CPistonGap::RefineMesh(int gridref,bool forcebalance){
	//cout << "Refine Mesh" << endl;
	bool refineflag = false;
	//Variable declaration
	int GridN,GridS,GridE,GridW;//indices of neighboring points
	int pointcheck;				//check if a point already exists
	int counter;
	int tracepoint1,tracepoint2;				//for tracing diagonal neighbors
	double X,nX,Y,nY,tol,coeff,cutoffheight,cutoffpress,toltemp;		//variables for finding duplicate points
	
	//include all old grid points (important)
	for(int i = 0;i<HDlevels[gridref].GridH.size();i++){

		//Don't define any neighbors at this time. At the new grid level, none of the existing points will be neighbors to each other.
		HDlevels[gridref+1].GridE.push_back(0);
		HDlevels[gridref+1].GridN.push_back(0);
		HDlevels[gridref+1].GridS.push_back(0);
		HDlevels[gridref+1].GridW.push_back(0);

		
		//Transfer the grid height
		HDlevels[gridref+1].GridH.push_back(HDlevels[gridref].GridH[i]);

		//Transfer the fluid pressure
		HDlevels[gridref+1].GridP.push_back(HDlevels[gridref].GridP[i]);
		
		//Transfer x coordinate
		HDlevels[gridref+1].GridX.push_back(HDlevels[gridref].GridX[i]);

		//Transfer y coordinate
		HDlevels[gridref+1].GridY.push_back(HDlevels[gridref].GridY[i]);

		//Transfer density
		HDlevels[gridref+1].GridD.push_back(HDlevels[gridref].GridD[i]);

		//Transfer expansion
		HDlevels[gridref+1].GridA.push_back(HDlevels[gridref].GridA[i]);

		//Transfer squeeze
		HDlevels[gridref+1].GridQ.push_back(HDlevels[gridref].GridQ[i]);

		//Transfer viscosity
		HDlevels[gridref+1].GridV.push_back(HDlevels[gridref].GridV[i]);

	};

	//Interpolate new points and define connectivity
	for(int i = 0;i<HDlevels[gridref].GridH.size();i++){
		//cout << gridref << "\t" << i << "\t" << HDlevels[gridref].GridH.size() << endl;
		bool refine = false;

		

		//Check the neighbors of point i in the current grid
		GridN = HDlevels[gridref].GridN[i];
		GridS = HDlevels[gridref].GridS[i];
		GridE = HDlevels[gridref].GridE[i];
		GridW = HDlevels[gridref].GridW[i];

		double maxnormalizedtilt = 0.3;
		//Get the location of point i
		X = HDlevels[gridref+1].GridX[i];
		Y = HDlevels[gridref+1].GridY[i];
		//(HDlevels[gridref].GridP[i]>cutoffpress)||(HDlevels[gridref].GridH[i]+def[gridref][i]<cutoffheight))
		//Current logic is to refine every grid point that has at least one neighbor.

		if(GridN)
			if(fabs(HDlevels[gridref].GridH[GridN-1]-HDlevels[gridref].GridH[i])/(HDlevels[gridref].GridH[GridN-1]+HDlevels[gridref].GridH[i])>maxnormalizedtilt)
				refine = true;
		if(GridS)
			if(fabs(HDlevels[gridref].GridH[GridS-1]-HDlevels[gridref].GridH[i])/(HDlevels[gridref].GridH[GridS-1]+HDlevels[gridref].GridH[i])>maxnormalizedtilt)
				refine = true;
		if(GridE)
			if(fabs(HDlevels[gridref].GridH[GridE-1]-HDlevels[gridref].GridH[i])/(HDlevels[gridref].GridH[GridE-1]+HDlevels[gridref].GridH[i])>maxnormalizedtilt)
				refine = true;
		if(GridW)
			if(fabs(HDlevels[gridref].GridH[GridW-1]-HDlevels[gridref].GridH[i])/(HDlevels[gridref].GridH[GridW-1]+HDlevels[gridref].GridH[i])>maxnormalizedtilt)
				refine = true;

		//refine = true;
		
		//cout << gridref << "\t" << hd_levels << endl;
		if(gridref >= hd_levels)
			refine = false;

		if(gridref == 0)
			//refine = true;
			if(i%M == 0)
				refine = true;
			else if(i%M == (M-1))
				refine = true;

		/*if(forcebalance){
			refine = refinememory[refinecounter];
			refinecounter++;
		}
		else
			refinememory.push_back(refine);*/

		

		if(refine){
			refineflag = true;


			//N
			if(GridN){

				//coordinates of the new point
				nX = 0.5*(HDlevels[gridref].GridX[GridN-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridN-1]+Y);

				//set tolerance for duplicate checking
				tol = 0.1*sqrt(pow(X-nX,2.0)+pow(Y-nY,2.0));

				//check for a duplicate point (already defined)
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);

				//if a point has already been defined in this location, update it's records
				if(pointcheck){

					//point i is to its south
					HDlevels[gridref+1].GridS[pointcheck] = i+1;

					//point GridN is to its north
					HDlevels[gridref+1].GridN[pointcheck] = GridN;
				}

				//if no point has already been defined in this location, define one by interpolating between point i and its Northern neighbor.
				else{

					//This point is now going to be the Northern neighbor of point i
					pointcheck = HDlevels[gridref+1].GridX.size();

					//Interpolate the cylinder bore surface temperature to the new grid point
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[GridN-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[GridN-1]));
					
					//Interpolate the original film thickness to the new grid point
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[GridN-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[GridN-1]));
					
					//Interpolate the fluid pressure to the new grid point
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[GridN-1]));
					
					//Interpolate the dh/dt to the new grid point
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[GridN-1]));
					
					//Assign the new x coordinate to the point
					HDlevels[gridref+1].GridX.push_back(nX);
					
					//Assign the new y coordinate to the point
					HDlevels[gridref+1].GridY.push_back(nY);
					
					//Define the new point's Northern Neighbor
					HDlevels[gridref+1].GridN.push_back(GridN);
					
					//Define the new point's Southern Neighbor
					HDlevels[gridref+1].GridS.push_back(i+1);

					//Try to find the new point's Eastern Neighbor
					//Reset the trace indices
					tracepoint1 = 0;
					tracepoint2 = 0;

					//First look to the east of point i for tracepoint1
					if(HDlevels[gridref+1].GridE.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridE[i];

					//if tracepoint1 exists, look to its north for tracepoint2
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						tracepoint2 = HDlevels[gridref+1].GridN[tracepoint1-1];

					//add tracepoint2 to the connectivity of the new point.
					HDlevels[gridref+1].GridE.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridW.size())){
						HDlevels[gridref+1].GridW[tracepoint2-1] = pointcheck+1;
					}
					//cout << "W1 " << i << "\t" << pointcheck << endl;
					
					//Try to find the new point's Western Neighbor
					//Reset the trace indices
					tracepoint1 = 0;
					tracepoint2 = 0;

					//First look to the West of point i for tracepoint1
					if(HDlevels[gridref+1].GridW.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridW[i];
					
					//Then trace North to tracepoint2
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						tracepoint2 = HDlevels[gridref+1].GridN[tracepoint1-1];

					//Add tracepoint2 to the connectivity of the new point.
					HDlevels[gridref+1].GridW.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridE.size())){
						HDlevels[gridref+1].GridE[tracepoint2-1] = pointcheck+1;
					}
				}

				//Update point i to indicate it's new Northern Neighbor
				HDlevels[gridref+1].GridN[i] = pointcheck+1;

				//Update point GridN to indicate it's new Southern Neighbor
				HDlevels[gridref+1].GridS[GridN-1] = pointcheck+1;
				//cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}
			else if((gridref == 0)&&(i%M == (M-1))){

				//coordinates of the new point
				nX = X;
				nY = 0.5*(3.0*Y-HDlevels[gridref].GridY[GridS-1]);

				//set tolerance for duplicate checking
				tol = 0.1*sqrt(pow(X-nX,2.0)+pow(Y-nY,2.0));

				//check for a duplicate point (already defined)
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);

				//if a point has already been defined in this location, update it's records
				if(pointcheck){

					//point i is to its south
					HDlevels[gridref+1].GridS[pointcheck] = i+1;
					
				}

				//if no point has already been defined in this location, define one by interpolating between point i and its Northern neighbor.
				else{

					//This point is now going to be the Northern neighbor of point i
					pointcheck = HDlevels[gridref+1].GridX.size();

					//Interpolate the cylinder bore surface temperature to the new grid point
					HDlevels[gridref+1].GridA.push_back(0.5*(3.0*HDlevels[gridref].GridA[i]-HDlevels[gridref].GridA[GridS-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(3.0*HDlevels[gridref].GridQ[i]-HDlevels[gridref].GridQ[GridS-1]));
					
					//Interpolate the original film thickness to the new grid point
					HDlevels[gridref+1].GridH.push_back(0.5*(3.0*HDlevels[gridref].GridH[i]-HDlevels[gridref].GridH[GridS-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(3.0*HDlevels[gridref].GridV[i]-HDlevels[gridref].GridV[GridS-1]));
					
					//Interpolate the fluid pressure to the new grid point
					HDlevels[gridref+1].GridP.push_back(operatingpistongap.pCase);
					
					//Interpolate the dh/dt to the new grid point
					HDlevels[gridref+1].GridD.push_back(0.5*(3.0*HDlevels[gridref].GridD[i]-HDlevels[gridref].GridD[GridS-1]));
					
					//Assign the new x coordinate to the point
					HDlevels[gridref+1].GridX.push_back(nX);
					
					//Assign the new y coordinate to the point
					HDlevels[gridref+1].GridY.push_back(nY);
					
					//Define the new point's Northern Neighbor
					HDlevels[gridref+1].GridN.push_back(0);
					
					//Define the new point's Southern Neighbor
					HDlevels[gridref+1].GridS.push_back(i+1);
					
					//Try to find the new point's Eastern Neighbor
					//Reset the trace indices
					tracepoint1 = 0;
					tracepoint2 = 0;

					//First look to the east of point i for tracepoint1
					if(HDlevels[gridref+1].GridE.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridE[i];

					//if tracepoint1 exists, look to its north for tracepoint2
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						tracepoint2 = HDlevels[gridref+1].GridN[tracepoint1-1];

					//add tracepoint2 to the connectivity of the new point.
					HDlevels[gridref+1].GridE.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridW.size()))
						HDlevels[gridref+1].GridW[tracepoint2-1] = pointcheck+1;
					
					//Try to find the new point's Western Neighbor
					//Reset the trace indices
					tracepoint1 = 0;
					tracepoint2 = 0;

					//First look to the West of point i for tracepoint1
					if(HDlevels[gridref+1].GridW.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridW[i];
					
					//Then trace North to tracepoint2
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						tracepoint2 = HDlevels[gridref+1].GridN[tracepoint1-1];

					//Add tracepoint2 to the connectivity of the new point.
					HDlevels[gridref+1].GridW.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridE.size()))
						HDlevels[gridref+1].GridE[tracepoint2-1] = pointcheck+1;
				}

				//Update point i to indicate it's new Northern Neighbor
				HDlevels[gridref+1].GridN[i] = pointcheck+1;
				boundarynodes.push_back(pointcheck);
				//Log << "New Point " << pointcheck << endl;
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}

			//E (similar to above)
			if(GridE){
				nX = 0.5*(HDlevels[gridref].GridX[GridE-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridE-1]+Y);
				if(nX < X)
					nX += PI*rK;
				if(nX < 0)
					nX += 2.0*PI*rK;
				nX = fmod(nX,2.0*PI*rK);
				toltemp = abs(X-nX);
				toltemp = min(toltemp,abs(X-nX-2.0*PI*rK));
				toltemp = min(toltemp,abs(X-nX+2.0*PI*rK));
				tol = 0.1*sqrt(pow(toltemp,2.0)+pow(Y-nY,2.0));
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);
				if(pointcheck){
					HDlevels[gridref+1].GridW[pointcheck] = i+1;
					HDlevels[gridref+1].GridE[pointcheck] = GridE;
				}
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[GridE-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[GridE-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[GridE-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[GridE-1]));
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[GridE-1]));
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[GridE-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridN.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridN[i];
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridE.size()))
						tracepoint2 = HDlevels[gridref+1].GridE[tracepoint1-1];
					HDlevels[gridref+1].GridN.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridS.size()))
						HDlevels[gridref+1].GridS[tracepoint2-1] = pointcheck+1;
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridS.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridS[i];
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridE.size()))
						tracepoint2 = HDlevels[gridref+1].GridE[tracepoint1-1];
					//cout << gridref << "\t" << i << "\t" << HDlevels[gridref].GridH.size() << "\t" << HDlevels[gridref+1].GridH.size() << "\t" << tracepoint1 << "\t" << tracepoint2 << endl;
					HDlevels[gridref+1].GridS.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridN.size()))
						HDlevels[gridref+1].GridN[tracepoint2-1] = pointcheck+1;
					HDlevels[gridref+1].GridE.push_back(GridE);
					HDlevels[gridref+1].GridW.push_back(i+1);
				}
				HDlevels[gridref+1].GridE[i] = pointcheck+1;
				HDlevels[gridref+1].GridW[GridE-1] = pointcheck+1;
				//cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}

			//S (similar to above)
			if(GridS){
				nX = 0.5*(HDlevels[gridref].GridX[GridS-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridS-1]+Y);
				tol = 0.1*sqrt(pow(X-nX,2.0)+pow(Y-nY,2.0));
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);
				if(pointcheck){
					HDlevels[gridref+1].GridN[pointcheck] = i+1;
					HDlevels[gridref+1].GridS[pointcheck] = GridS;
				}
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[GridS-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[GridS-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[GridS-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[GridS-1]));
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[GridS-1]));
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[GridS-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					HDlevels[gridref+1].GridN.push_back(i+1);
					HDlevels[gridref+1].GridS.push_back(GridS);
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridE.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridE[i];
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						tracepoint2 = HDlevels[gridref+1].GridS[tracepoint1-1];
					HDlevels[gridref+1].GridE.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridW.size()))
						HDlevels[gridref+1].GridW[tracepoint2-1] = pointcheck + 1;
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridW.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridW[i];
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						tracepoint2 = HDlevels[gridref+1].GridS[tracepoint1-1];
					HDlevels[gridref+1].GridW.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridE.size()))
						HDlevels[gridref+1].GridE[tracepoint2-1] = pointcheck+1;
				}
				HDlevels[gridref+1].GridS[i] = pointcheck+1;
				HDlevels[gridref+1].GridN[GridS-1] = pointcheck+1;
				//cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}
			else if((gridref == 0)&&(i%M == 0)){
				nX = X;
				nY = 0.5*(3.0*Y-HDlevels[gridref].GridY[GridN-1]);
				tol = 0.1*sqrt(pow(X-nX,2.0)+pow(Y-nY,2.0));
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);
				if(pointcheck){
					HDlevels[gridref+1].GridN[pointcheck] = i+1;
				}
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(3.0*HDlevels[gridref].GridA[i]-HDlevels[gridref].GridA[GridN-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(3.0*HDlevels[gridref].GridQ[i]-HDlevels[gridref].GridQ[GridN-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(3.0*HDlevels[gridref].GridH[i]-HDlevels[gridref].GridH[GridN-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(3.0*HDlevels[gridref].GridV[i]-HDlevels[gridref].GridV[GridN-1]));
					HDlevels[gridref+1].GridP.push_back(operatingpistongap.pDC);
					HDlevels[gridref+1].GridD.push_back(0.5*(3.0*HDlevels[gridref].GridD[i]-HDlevels[gridref].GridD[GridN-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					HDlevels[gridref+1].GridN.push_back(i+1);
					HDlevels[gridref+1].GridS.push_back(0);
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridE.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridE[i];
					if((tracepoint1 >0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						tracepoint2 = HDlevels[gridref+1].GridS[tracepoint1-1];
					HDlevels[gridref+1].GridE.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridW.size()))
						HDlevels[gridref+1].GridW[tracepoint2-1] = pointcheck + 1;
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridW.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridW[i];
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						tracepoint2 = HDlevels[gridref+1].GridS[tracepoint1-1];
					HDlevels[gridref+1].GridW.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridE.size()))
						HDlevels[gridref+1].GridE[tracepoint2-1] = pointcheck+1;
				}
				HDlevels[gridref+1].GridS[i] = pointcheck+1;
				boundarynodes.push_back(pointcheck);
				//Log << "New Point " << pointcheck << endl;
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}

			//W (similar to above)
			if(GridW){
				nX = 0.5*(HDlevels[gridref].GridX[GridW-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridW-1]+Y);
				
				if(nX > X){
					nX -= PI*rK;
				}
				
				if(nX < 0)
					nX += 2.0*PI*rK;
				nX = fmod(nX,2.0*PI*rK);
				//if(nX <= 0)
					//nX += 2.0*PI*rK;
				toltemp = abs(X-nX);
				toltemp = min(toltemp,abs(X-nX-2.0*PI*rK));
				toltemp = min(toltemp,abs(X-nX+2.0*PI*rK));
				tol = 0.1*sqrt(pow(toltemp,2.0)+pow(Y-nY,2.0));
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);
				//cout << i << "\t" << pointcheck << endl;
				if(pointcheck){
					HDlevels[gridref+1].GridE[pointcheck] = i+1;
					HDlevels[gridref+1].GridW[pointcheck] = GridW;
				}
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[GridW-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[GridW-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[GridW-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[GridW-1]));
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[GridW-1]));
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[GridW-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridN.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridN[i];
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridW.size()))
						tracepoint2 = HDlevels[gridref+1].GridW[tracepoint1-1];
					HDlevels[gridref+1].GridN.push_back(tracepoint2);
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridS.size()))
						HDlevels[gridref+1].GridS[tracepoint2-1] = pointcheck+1;
					tracepoint1 = 0;
					tracepoint2 = 0;
					if(HDlevels[gridref+1].GridS.size()>i)
						tracepoint1 = HDlevels[gridref+1].GridS[i];
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridW.size()))
						tracepoint2 = HDlevels[gridref+1].GridW[tracepoint1-1];
					HDlevels[gridref+1].GridS.push_back(tracepoint2);
					//cout << tracepoint1 << "\t" << tracepoint2 << "\t" << HDlevels[gridref+1].GridS.size() << "\t" << i << endl;
					if((tracepoint2 > 0)&&(tracepoint2 <= HDlevels[gridref+1].GridN.size()))
						HDlevels[gridref+1].GridN[tracepoint2-1] = pointcheck+1;
					HDlevels[gridref+1].GridE.push_back(i+1);
					HDlevels[gridref+1].GridW.push_back(GridW);
				}
				HDlevels[gridref+1].GridW[i] = pointcheck+1;
				HDlevels[gridref+1].GridE[GridW-1] = pointcheck+1;
				//cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}

			//NE
			if(GridN&&GridE){

				//calculate the coordinates for the new point
				nX = 0.5*(HDlevels[gridref].GridX[GridE-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridN-1]+Y);
				if(nX < X)
					nX += PI*rK;
				if(nX < 0)
					nX += 2.0*PI*rK;
				nX = fmod(nX,2.0*PI*rK);
				//calculate the tolerance for duplicate checking
				toltemp = abs(X-nX);
				toltemp = min(toltemp,abs(X-nX-2.0*PI*rK));
				toltemp = min(toltemp,abs(X-nX+2.0*PI*rK));
				tol = 0.1*sqrt(pow(toltemp,2.0)+pow(Y-nY,2.0));

				//check for a point already defined in this location
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);

				//if a point already exists in this location, update the connectivity.
				if(pointcheck){

					//Define the North and East neighbors of the new point.
					HDlevels[gridref+1].GridS[pointcheck] = HDlevels[gridref+1].GridE[i];
					HDlevels[gridref+1].GridW[pointcheck] = HDlevels[gridref+1].GridN[i];
				}

				//if no point exists yet, define one based on the four surrounding points (if possible)
				else if(HDlevels[gridref].GridN[GridE-1]){

					//Find the index of the new point
					pointcheck = HDlevels[gridref+1].GridX.size();

					//Interpolate the cylinder surface temperature
					HDlevels[gridref+1].GridA.push_back(0.25*(HDlevels[gridref].GridA[GridN-1]+HDlevels[gridref].GridA[GridE-1]+HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[HDlevels[gridref].GridN[GridE-1]-1]));
					HDlevels[gridref+1].GridQ.push_back(0.25*(HDlevels[gridref].GridQ[GridN-1]+HDlevels[gridref].GridQ[GridE-1]+HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[HDlevels[gridref].GridN[GridE-1]-1]));
					
					//interpolate the film thickness
					HDlevels[gridref+1].GridH.push_back(0.25*(HDlevels[gridref].GridH[GridN-1]+HDlevels[gridref].GridH[GridE-1]+HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[HDlevels[gridref].GridN[GridE-1]-1]));
					HDlevels[gridref+1].GridV.push_back(0.25*(HDlevels[gridref].GridV[GridN-1]+HDlevels[gridref].GridV[GridE-1]+HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[HDlevels[gridref].GridN[GridE-1]-1]));
					
					//interpolate the pressure
					HDlevels[gridref+1].GridP.push_back(0.25*(HDlevels[gridref].GridP[GridN-1]+HDlevels[gridref].GridP[GridE-1]+HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[HDlevels[gridref].GridN[GridE-1]-1]));
					
					//interpolate the dh/dt
					HDlevels[gridref+1].GridD.push_back(0.25*(HDlevels[gridref].GridD[GridN-1]+HDlevels[gridref].GridD[GridE-1]+HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[HDlevels[gridref].GridN[GridE-1]-1]));

					//set the x coordinate
					HDlevels[gridref+1].GridX.push_back(nX);
					
					//set the Y coordinate
					HDlevels[gridref+1].GridY.push_back(nY);
					
					//Look for a neighbor to the North
					tracepoint1 = HDlevels[gridref+1].GridE[GridN-1];
					HDlevels[gridref+1].GridN.push_back(tracepoint1);
					
					//If a neighbor has been found, update it's connectivity as well.
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						HDlevels[gridref+1].GridS[tracepoint1-1] = pointcheck+1;
					
					//Define the neighbor to the South
					HDlevels[gridref+1].GridS.push_back(HDlevels[gridref+1].GridE[i]);
					
					//Look for a neighbor to the East
					tracepoint1 = HDlevels[gridref+1].GridN[GridE-1];
					HDlevels[gridref+1].GridE.push_back(tracepoint1);
					
					//If a neighbor has been found, update it's connectivity as well.
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridW.size()))
						HDlevels[gridref+1].GridW[tracepoint1-1] = pointcheck+1;
					
					//Define the neighbor to the North
					HDlevels[gridref+1].GridW.push_back(HDlevels[gridref+1].GridN[i]);
				}

				//If there aren't four surrounding points, then define the new point just from the two diagonals, similar to above
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[GridN-1]+HDlevels[gridref].GridA[GridE-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[GridN-1]+HDlevels[gridref].GridQ[GridE-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[GridN-1]+HDlevels[gridref].GridH[GridE-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[GridN-1]+HDlevels[gridref].GridV[GridE-1]));
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[GridN-1]+HDlevels[gridref].GridP[GridE-1]));
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[GridN-1]+HDlevels[gridref].GridD[GridE-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					tracepoint1 = HDlevels[gridref+1].GridE[GridN-1];
					HDlevels[gridref+1].GridN.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						HDlevels[gridref+1].GridS[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridS.push_back(HDlevels[gridref+1].GridE[i]);
					tracepoint1 = HDlevels[gridref+1].GridN[GridE-1];
					HDlevels[gridref+1].GridE.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridW.size()))
						HDlevels[gridref+1].GridW[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridW.push_back(HDlevels[gridref+1].GridN[i]);
				}

				//Update the connectivity of the points to the West and South
				HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridE[i]-1]=pointcheck+1;
				HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridN[i]-1]=pointcheck+1;
				//cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}

			//SE (similar to above)
			if(GridS&&GridE){
				nX = 0.5*(HDlevels[gridref].GridX[GridE-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridS-1]+Y);
				if(nX < X)
					nX += PI*rK;
				if(nX < 0)
					nX += 2.0*PI*rK;
				nX = fmod(nX,2.0*PI*rK);
				toltemp = abs(X-nX);
				toltemp = min(toltemp,abs(X-nX-2.0*PI*rK));
				toltemp = min(toltemp,abs(X-nX+2.0*PI*rK));
				tol = 0.1*sqrt(pow(toltemp,2.0)+pow(Y-nY,2.0));
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);
				if(pointcheck){
					HDlevels[gridref+1].GridN[pointcheck] = HDlevels[gridref+1].GridE[i];
					HDlevels[gridref+1].GridW[pointcheck] = HDlevels[gridref+1].GridS[i];
				}
				else if(HDlevels[gridref].GridS[GridE-1]){
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.25*(HDlevels[gridref].GridA[GridS-1]+HDlevels[gridref].GridA[GridE-1]+HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[HDlevels[gridref].GridS[GridE-1]-1]));
					HDlevels[gridref+1].GridQ.push_back(0.25*(HDlevels[gridref].GridQ[GridS-1]+HDlevels[gridref].GridQ[GridE-1]+HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[HDlevels[gridref].GridS[GridE-1]-1]));
					HDlevels[gridref+1].GridH.push_back(0.25*(HDlevels[gridref].GridH[GridS-1]+HDlevels[gridref].GridH[GridE-1]+HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[HDlevels[gridref].GridS[GridE-1]-1]));
					HDlevels[gridref+1].GridV.push_back(0.25*(HDlevels[gridref].GridV[GridS-1]+HDlevels[gridref].GridV[GridE-1]+HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[HDlevels[gridref].GridS[GridE-1]-1]));
					HDlevels[gridref+1].GridP.push_back(0.25*(HDlevels[gridref].GridP[GridS-1]+HDlevels[gridref].GridP[GridE-1]+HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[HDlevels[gridref].GridS[GridE-1]-1]));
					HDlevels[gridref+1].GridD.push_back(0.25*(HDlevels[gridref].GridD[GridS-1]+HDlevels[gridref].GridD[GridE-1]+HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[HDlevels[gridref].GridS[GridE-1]-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					HDlevels[gridref+1].GridN.push_back(HDlevels[gridref+1].GridE[i]);
					tracepoint1 = HDlevels[gridref+1].GridE[GridS-1];
					HDlevels[gridref+1].GridS.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						HDlevels[gridref+1].GridN[tracepoint1-1] = pointcheck+1;
					tracepoint1 = HDlevels[gridref+1].GridS[GridE-1];
					HDlevels[gridref+1].GridE.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridW.size()))
						HDlevels[gridref+1].GridW[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridW.push_back(HDlevels[gridref+1].GridS[i]);
				}
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[GridS-1]+HDlevels[gridref].GridA[GridE-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[GridS-1]+HDlevels[gridref].GridQ[GridE-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[GridS-1]+HDlevels[gridref].GridH[GridE-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[GridS-1]+HDlevels[gridref].GridV[GridE-1]));
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[GridS-1]+HDlevels[gridref].GridP[GridE-1]));
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[GridS-1]+HDlevels[gridref].GridD[GridE-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					HDlevels[gridref+1].GridN.push_back(HDlevels[gridref+1].GridE[i]);
					tracepoint1 = HDlevels[gridref+1].GridE[GridS-1];
					HDlevels[gridref+1].GridS.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						HDlevels[gridref+1].GridN[tracepoint1-1] = pointcheck+1;
					tracepoint1 = HDlevels[gridref+1].GridS[GridE-1];
					HDlevels[gridref+1].GridE.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridW.size()))
						HDlevels[gridref+1].GridW[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridW.push_back(HDlevels[gridref+1].GridS[i]);
				}
				HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridE[i]-1]=pointcheck+1;
				HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridS[i]-1]=pointcheck+1;
				///cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}

			//SW (similar to above)
			if(GridS&&GridW){
				nX = 0.5*(HDlevels[gridref].GridX[GridW-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridS-1]+Y);
				if(nX > X)
					nX -= PI*rK;
				if(nX < 0)
					nX += 2.0*PI*rK;
				nX = fmod(nX,2.0*PI*rK);
				//if(nX < 0)
					//nX += 2.0*PI*rK;
				
				//if(nX <= 0)
				//	nX += 2.0*PI*rK;
				toltemp = min(abs(X-nX),abs(X-nX-2.0*PI*rK));
				toltemp = min(toltemp,abs(X-nX+2.0*PI*rK));
				tol = 0.1*sqrt(pow(toltemp,2.0)+pow(Y-nY,2.0));
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);
				if(pointcheck){
					HDlevels[gridref+1].GridN[pointcheck] = HDlevels[gridref+1].GridW[i];
					HDlevels[gridref+1].GridE[pointcheck] = HDlevels[gridref+1].GridS[i];
				}
				else if(HDlevels[gridref].GridS[GridW-1]){
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.25*(HDlevels[gridref].GridA[GridS-1]+HDlevels[gridref].GridA[GridW-1]+HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[HDlevels[gridref].GridS[GridW-1]-1]));
					HDlevels[gridref+1].GridQ.push_back(0.25*(HDlevels[gridref].GridQ[GridS-1]+HDlevels[gridref].GridQ[GridW-1]+HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[HDlevels[gridref].GridS[GridW-1]-1]));
					HDlevels[gridref+1].GridH.push_back(0.25*(HDlevels[gridref].GridH[GridS-1]+HDlevels[gridref].GridH[GridW-1]+HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[HDlevels[gridref].GridS[GridW-1]-1]));
					HDlevels[gridref+1].GridV.push_back(0.25*(HDlevels[gridref].GridV[GridS-1]+HDlevels[gridref].GridV[GridW-1]+HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[HDlevels[gridref].GridS[GridW-1]-1]));
					HDlevels[gridref+1].GridP.push_back(0.25*(HDlevels[gridref].GridP[GridS-1]+HDlevels[gridref].GridP[GridW-1]+HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[HDlevels[gridref].GridS[GridW-1]-1]));
					HDlevels[gridref+1].GridD.push_back(0.25*(HDlevels[gridref].GridD[GridS-1]+HDlevels[gridref].GridD[GridW-1]+HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[HDlevels[gridref].GridS[GridW-1]-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					HDlevels[gridref+1].GridN.push_back(HDlevels[gridref+1].GridW[i]);
					tracepoint1 = HDlevels[gridref+1].GridW[GridS-1];
					HDlevels[gridref+1].GridS.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						HDlevels[gridref+1].GridN[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridE.push_back(HDlevels[gridref+1].GridS[i]);
					tracepoint1 = HDlevels[gridref+1].GridS[GridW-1];
					HDlevels[gridref+1].GridW.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridE.size()))
						HDlevels[gridref+1].GridE[tracepoint1-1] = pointcheck+1;
				}
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[GridS-1]+HDlevels[gridref].GridA[GridW-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[GridS-1]+HDlevels[gridref].GridQ[GridW-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[GridS-1]+HDlevels[gridref].GridH[GridW-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[GridS-1]+HDlevels[gridref].GridV[GridW-1]));
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[GridS-1]+HDlevels[gridref].GridP[GridW-1]));
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[GridS-1]+HDlevels[gridref].GridD[GridW-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					HDlevels[gridref+1].GridN.push_back(HDlevels[gridref+1].GridW[i]);
					tracepoint1 = HDlevels[gridref+1].GridW[GridS-1];
					HDlevels[gridref+1].GridS.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridN.size()))
						HDlevels[gridref+1].GridN[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridE.push_back(HDlevels[gridref+1].GridS[i]);
					tracepoint1 = HDlevels[gridref+1].GridS[GridW-1];
					HDlevels[gridref+1].GridW.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridE.size()))
						HDlevels[gridref+1].GridE[tracepoint1-1] = pointcheck+1;
				}
				HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridW[i]-1]=pointcheck+1;
				HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridS[i]-1]=pointcheck+1;
				//cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}

			//NW (similar to above)
			if(GridN&&GridW){
				nX = 0.5*(HDlevels[gridref].GridX[GridW-1]+X);
				nY = 0.5*(HDlevels[gridref].GridY[GridN-1]+Y);
				if(nX > X)
					nX -= PI*rK;
				
				if(nX < 0)
					nX += 2.0*PI*rK;
				nX = fmod(nX,2.0*PI*rK);
				
				//if(nX <= 0)
					//nX += 2.0*PI*rK;
				toltemp = min(abs(X-nX),abs(X-nX-2.0*PI*rK));
				toltemp = min(toltemp,abs(X-nX+2.0*PI*rK));
				tol = 0.1*sqrt(pow(toltemp,2.0)+pow(Y-nY,2.0));
				pointcheck = CheckPoint(nX,nY,tol,gridref+1);
				if(pointcheck){
					HDlevels[gridref+1].GridS[pointcheck] = HDlevels[gridref+1].GridW[i];
					HDlevels[gridref+1].GridE[pointcheck] = HDlevels[gridref+1].GridN[i];
				}
				else if(HDlevels[gridref].GridN[GridW-1]){
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.25*(HDlevels[gridref].GridA[GridN-1]+HDlevels[gridref].GridA[GridW-1]+HDlevels[gridref].GridA[i]+HDlevels[gridref].GridA[HDlevels[gridref].GridN[GridW-1]-1]));
					HDlevels[gridref+1].GridQ.push_back(0.25*(HDlevels[gridref].GridQ[GridN-1]+HDlevels[gridref].GridQ[GridW-1]+HDlevels[gridref].GridQ[i]+HDlevels[gridref].GridQ[HDlevels[gridref].GridN[GridW-1]-1]));
					HDlevels[gridref+1].GridH.push_back(0.25*(HDlevels[gridref].GridH[GridN-1]+HDlevels[gridref].GridH[GridW-1]+HDlevels[gridref].GridH[i]+HDlevels[gridref].GridH[HDlevels[gridref].GridN[GridW-1]-1]));
					HDlevels[gridref+1].GridV.push_back(0.25*(HDlevels[gridref].GridV[GridN-1]+HDlevels[gridref].GridV[GridW-1]+HDlevels[gridref].GridV[i]+HDlevels[gridref].GridV[HDlevels[gridref].GridN[GridW-1]-1]));
					HDlevels[gridref+1].GridP.push_back(0.25*(HDlevels[gridref].GridP[GridN-1]+HDlevels[gridref].GridP[GridW-1]+HDlevels[gridref].GridP[i]+HDlevels[gridref].GridP[HDlevels[gridref].GridN[GridW-1]-1]));
					HDlevels[gridref+1].GridD.push_back(0.25*(HDlevels[gridref].GridD[GridN-1]+HDlevels[gridref].GridD[GridW-1]+HDlevels[gridref].GridD[i]+HDlevels[gridref].GridD[HDlevels[gridref].GridN[GridW-1]-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					tracepoint1 = HDlevels[gridref+1].GridW[GridN-1];
					HDlevels[gridref+1].GridN.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						HDlevels[gridref+1].GridS[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridS.push_back(HDlevels[gridref+1].GridW[i]);
					HDlevels[gridref+1].GridE.push_back(HDlevels[gridref+1].GridN[i]);
					tracepoint1 = HDlevels[gridref+1].GridN[GridW-1];
					HDlevels[gridref+1].GridW.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridE.size()))
						HDlevels[gridref+1].GridE[tracepoint1-1] = pointcheck+1;
				}
				else{
					pointcheck = HDlevels[gridref+1].GridX.size();
					HDlevels[gridref+1].GridA.push_back(0.5*(HDlevels[gridref].GridA[GridN-1]+HDlevels[gridref].GridA[GridW-1]));
					HDlevels[gridref+1].GridQ.push_back(0.5*(HDlevels[gridref].GridQ[GridN-1]+HDlevels[gridref].GridQ[GridW-1]));
					HDlevels[gridref+1].GridH.push_back(0.5*(HDlevels[gridref].GridH[GridN-1]+HDlevels[gridref].GridH[GridW-1]));
					HDlevels[gridref+1].GridV.push_back(0.5*(HDlevels[gridref].GridV[GridN-1]+HDlevels[gridref].GridV[GridW-1]));
					HDlevels[gridref+1].GridP.push_back(0.5*(HDlevels[gridref].GridP[GridN-1]+HDlevels[gridref].GridP[GridW-1]));
					HDlevels[gridref+1].GridD.push_back(0.5*(HDlevels[gridref].GridD[GridN-1]+HDlevels[gridref].GridD[GridW-1]));
					HDlevels[gridref+1].GridX.push_back(nX);
					HDlevels[gridref+1].GridY.push_back(nY);
					tracepoint1 = HDlevels[gridref+1].GridW[GridN-1];
					HDlevels[gridref+1].GridN.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridS.size()))
						HDlevels[gridref+1].GridS[tracepoint1-1] = pointcheck+1;
					HDlevels[gridref+1].GridS.push_back(HDlevels[gridref+1].GridW[i]);
					HDlevels[gridref+1].GridE.push_back(HDlevels[gridref+1].GridN[i]);
					tracepoint1 = HDlevels[gridref+1].GridN[GridW-1];
					HDlevels[gridref+1].GridW.push_back(tracepoint1);
					if((tracepoint1 > 0)&&(tracepoint1 <= HDlevels[gridref+1].GridE.size()))
						HDlevels[gridref+1].GridE[tracepoint1-1] = pointcheck+1;
				}
				HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridW[i]-1]=pointcheck+1;
				HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridN[i]-1]=pointcheck+1;
				//cout << "New Point " << HDlevels[gridref+1].GridB.size() << ":\n";
				//cout << GridN << "\t" << GridS << "\t" << GridE << "\t" << GridW << "\n";
				//cout << HDlevels[gridref+1].GridN[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridS[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridE[HDlevels[gridref+1].GridB.size()-1] << "\t" << HDlevels[gridref+1].GridW[HDlevels[gridref+1].GridB.size()-1] << "\n";
			}
			
		};
	}

	

	return refineflag;
};
int CPistonGap::CheckPoint(double nX,double nY, double tol, int gridref){

	//Variable declaration
	vector<double> dist;	//distance from each grid point to the coordinates nX,nY
	double mindist;			//minimum distance found
	int toreturn;			//return value

	if(nX < 0.0)
		nX += 2.0*PI*rK;

	nX = fmod(nX,2.0*PI*rK);

	//variable initialization
	mindist = 1e11;											//set minimum distance somewhere far away.
	dist.resize(HDlevels[gridref].GridX.size());	//resize the distance vector
	
	//find the distance from each point in the grid to the new coordinates
	for(int i = 0;i<HDlevels[gridref].GridX.size();i++){
		dist[i] = sqrt(pow(HDlevels[gridref].GridX[i]-nX,2.0)+pow(HDlevels[gridref].GridY[i]-nY,2.0));

		if(dist[i] <= mindist)
			mindist = dist[i];
	}
	
	//if a point has been found within the specified tolerance, output it's index
	if(mindist < tol){

		//loop through the points, looking for the one that is closest to the new coordinates.
		for(int i = 0;i<dist.size();i++){

			//once it's been found, return it's index.
			if(mindist == dist[i]){
				toreturn = i;
				break;
			}
		}
	}
	
	//If a point hasn't been found within the tolerance, return 0
	else{
		toreturn = 0;
	}

	//return the final value.
	//cout << "c" << toreturn << endl;
	return(toreturn);
};
void CPistonGap::SolveReynoldsHD(int gridref,int maxlev,double alphapress,double dt){

	//Indices for neighboring cells
	int GridN,GridS,GridE,GridW;
	//double alphapress = AlphaP;

	//Loop through every grid point in this level and solve the Reynolds equation
	for(int i = 0;i<HDlevels[gridref].GridP.size();i++){
		//Get the indices for the neighboring points
		GridN = HDlevels[gridref].GridN[i];
		GridS = HDlevels[gridref].GridS[i];
		GridE = HDlevels[gridref].GridE[i];
		GridW = HDlevels[gridref].GridW[i];
		
		// This loop should only solve for the pressure of a point if it is interior to the region. Do not update pressure on boundaries!
		if(((GridN && GridS) && (GridE && GridW))&&((GridN!=GridS)&&(GridE!=GridW))){
			double pn,ps,pe,pw,pp;			//pressure (North, South, East, West, and Point) [bar]
			double hn,hs,he,hw,hp;			//film thickness [m]
			double mun,mus,mue,muw,mup;		//viscosity [Pa-s]
			double dx,dy;					//grid spacing [m]
			double e;						//residual error from Reynolds equation [???]
			double dedp;					//derivative of error with respect to pressure [???/Pa]
			double rhon,rhos,rhoe,rhow,rhop;//fluid density [Kg/m^3]
			
			//Calculate the grid spacing based on the neighboring points.
			dy = (HDlevels[gridref].GridY[GridN-1]-HDlevels[gridref].GridY[GridS-1])/2.0;
			dx = (HDlevels[gridref].GridX[GridE-1]-HDlevels[gridref].GridX[GridW-1])/2.0;
			if(dx < 0.0)
				dx += PI*rK;
			dx = fmod(dx,2.0*PI*rK);

			if(dx*dy<=0.0){
				Log << "Error at: " << i << "\n";
			}
			else
			{
				//If reshere is > gridref, that means the pressure value in the levels structure is actually a residual, 
				//and we need to use updatepress instead.
				if(reshere[i] > gridref){
					pn = updatepress[GridN-1];
					ps = updatepress[GridS-1];
					pe = updatepress[GridE-1];
					pw = updatepress[GridW-1];
					pp = updatepress[i];
				}

				//Otherwise we can just pull the pressure value out of the levels structure and use that.
				else{
					pp = HDlevels[maxlev].GridP[i];
					pn = HDlevels[maxlev].GridP[GridN-1];
					ps = HDlevels[maxlev].GridP[GridS-1];
					pe = HDlevels[maxlev].GridP[GridE-1];
					pw = HDlevels[maxlev].GridP[GridW-1];
				}

				//If this level is the finest refinement level for this point, update the pressure.
				if(reshere[i] == gridref){
					HDlevels[maxlev].GridP[i] += alphapress * updatepress[i];
					pp = HDlevels[maxlev].GridP[i];
				}
			
				//get the film thickness, and ensure that it's above the minimum allowed.
				hp = HDlevels[gridref].GridH[i];
				//hp = max(hp,minh);

				//get density from the oil model
				rhop = HDlevels[gridref].GridD[i];
			
				//get viscosity from the oil model
				mup = HDlevels[gridref].GridV[i];
			
				//Gather information about the point to the North
				hn = HDlevels[gridref].GridH[GridN-1];											//height		//saturate height to minimum allowable
				hn = (hn+hp)/2.0;
				mun = HDlevels[gridref].GridV[GridN-1];							//viscosity
				mun = (mun+mup)/2.0;
				rhon = HDlevels[gridref].GridD[GridN-1];	//density
				rhon = (rhon+rhop)/2.0;
			
				//Repeat of above for the point to the South
				hs = HDlevels[gridref].GridH[GridS-1];
				hs = (hs+hp)/2.0;
				mus = HDlevels[gridref].GridV[GridS-1];
				mus = (mus+mup)/2.0;
				rhos = HDlevels[gridref].GridD[GridS-1];
				rhos = (rhos+rhop)/2.0;
			
				//Repeat of above for the point to the East
				he = HDlevels[gridref].GridH[GridE-1];
				he = (he+hp)/2.0;
				mue = HDlevels[gridref].GridV[GridE-1];
				mue = (mue+mup)/2.0;
				rhoe = HDlevels[gridref].GridD[GridE-1];
				rhoe = (rhoe+rhop)/2.0;
			
				//Repeat of above for the point to the West
				hw = HDlevels[gridref].GridH[GridW-1];
				hw = (hw+hp)/2.0;
				muw = HDlevels[gridref].GridV[GridW-1];
				muw = (muw+mup)/2.0;
				rhow = HDlevels[gridref].GridD[GridW-1];
				rhow = (rhow+rhop)/2.0;
			
				//If the residual already exists here, we use that to replace the source term, but solve for a pressure correction in the normal way.
				if(reshere[i] > gridref){
				
					/*e = (1.0/12.0)*(rhon-rhos)/(2.0*dy)*(rhop*pow(hp,3.0)/mup)*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(pow(hn,3.0)-pow(hs,3.0))/(2.0*dy)*(rhop/mup)*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(1.0/mun-1.0/mus)/(2.0*dy)*(rhop*pow(hp,3.0))*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(ps+pn-2.0*pp)/pow(dy,2.0)
						+(1.0/12.0)*(rhoe-rhow)/(2.0*dx)*(rhop*pow(hp,3.0)/mup)*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(pow(he,3.0)-pow(hw,3.0))/(2.0*dx)*(rhop/mup)*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(1.0/mue-1.0/muw)/(2.0*dx)*(rhop*pow(hp,3.0))*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(pe+pw-2.0*pp)/pow(dx,2.0)
						+HDlevels[gridref].GridP[i];
					dedp = (1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(-2.0)/pow(dy,2.0)	//change of error with change in pressure
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(-2.0)/pow(dx,2.0);*/

					e =  (((rhon*(pn-pp)*pow(hn,3.0)/(12.0*dy*mun))-(rhos*(pp-ps)*pow(hs,3.0)/(12.0*dy*mus)))/(1.0*dy))	//error with Reynolds source terms
						+(((rhoe*(pe-pp)*pow(he,3.0)/(12.0*dx*mue))-(rhow*(pp-pw)*pow(hw,3.0)/(12.0*dx*muw)))/(1.0*dx))
						+HDlevels[gridref].GridP[i];
					dedp = (((-1.0*rhon*pow(hn,3.0)/(12.0*dy*mun))-(rhos*pow(hs,3.0)/(12.0*dy*mus)))/(1.0*dy))	//change of error with change in pressure
						+(((-1.0*rhoe*pow(he,3.0)/(12.0*dx*mue))-(rhow*pow(hw,3.0)/(12.0*dx*muw)))/(1.0*dx));
				
					updatepress[i] -= alphapress * e / dedp;									//updating the pressure
					//updatepress[i] = updatepress[i]<1e4?1e4:updatepress[i];		//limit the pressure to positive values
					//updatepress[i] = updatepress[i]>1e9?1e9:updatepress[i];
					if(updatepress[i]!=updatepress[i]){
						cout << "1 " << e << "\t" << dedp << "\t" << updatepress[i] << "\n";
						cout << rhon << "\t" << hn << "\t" << mun << "\t" << pn << "\t" << pp << "\t" << dy << "\t" << rhos << "\t" << hs << "\t" << mus << "\t" << ps << "\t" << rhoe << "\t" << he << "\t" << mue << "\t" << pe << "\t" << dx << "\t" << rhow << "\t" << hw << "\t" << muw << "\t" << HDlevels[gridref].GridP[i] << "\n";
					
						system("PAUSE");
					}
				}
				else{

					/*e = (1.0/12.0)*(rhon-rhos)/(2.0*dy)*(rhop*pow(hp,3.0)/mup)*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(pow(hn,3.0)-pow(hs,3.0))/(2.0*dy)*(rhop/mup)*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(1.0/mun-1.0/mus)/(2.0*dy)*(rhop*pow(hp,3.0))*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(ps+pn-2.0*pp)/pow(dy,2.0)
						+(1.0/12.0)*(rhoe-rhow)/(2.0*dx)*(rhop*pow(hp,3.0)/mup)*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(pow(he,3.0)-pow(hw,3.0))/(2.0*dx)*(rhop/mup)*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(1.0/mue-1.0/muw)/(2.0*dx)*(rhop*pow(hp,3.0))*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(pe+pw-2.0*pp)/pow(dx,2.0)
						+HDlevels[gridref].GridB[i];
					dedp = (1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(-2.0)/pow(dy,2.0)	//change of error with change in pressure
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(-2.0)/pow(dx,2.0);*/

					e =  (((rhon*(pn-pp)*pow(hn,3.0)/(12.0*dy*mun))-(rhos*(pp-ps)*pow(hs,3.0)/(12.0*dy*mus)))/(1.0*dy))	//error with Reynolds source terms
						+(((rhoe*(pe-pp)*pow(he,3.0)/(12.0*dx*mue))-(rhow*(pp-pw)*pow(hw,3.0)/(12.0*dx*muw)))/(1.0*dx))
						-0.25 * operatingpistongap.speedK * operatingpistongap.omega * geometrypistongap.dK * (rhoe*he - rhow*hw)/dx
						-0.5 * operatingpistongap.vK * (rhon*hn - rhos*hs)/dy
						+rhop * HDlevels[gridref].GridQ[i]
						+hp * (rhop - HDlevels[gridref].GridA[i])/dt;
					dedp = (((-1.0*rhon*pow(hn,3.0)/(12.0*dy*mun))-(rhos*pow(hs,3.0)/(12.0*dy*mus)))/(1.0*dy))	//change of error with change in pressure
						+(((-1.0*rhoe*pow(he,3.0)/(12.0*dx*mue))-(rhow*pow(hw,3.0)/(12.0*dx*muw)))/(1.0*dx));

					HDlevels[maxlev].GridP[i] -= alphapress * e/dedp;																//update the pressure
					HDlevels[maxlev].GridP[i] = HDlevels[maxlev].GridP[i]<1e4?1e4:HDlevels[maxlev].GridP[i];		//limit the pressure to positive values
					HDlevels[maxlev].GridP[i] = HDlevels[maxlev].GridP[i]>1e9?1e9:HDlevels[maxlev].GridP[i];	//limit the maximum pressure
					if(HDlevels[maxlev].GridP[i]!=HDlevels[maxlev].GridP[i]){
						cout << "2 " << e << "\t" << dedp << "\t" << HDlevels[maxlev].GridP[i] << "\n";
						//cout << i << "\t" << GridN << "\t" << GridS << "\t" << HDlevels[gridref].GridY[GridN-1] << "\t" << HDlevels[gridref].GridY[GridS-1] << "\t" << rhon << "\t" << hn << "\t" << mun << "\t" << pn << "\t" << pp << "\t" << dy << "\t" << rhos << "\t" << hs << "\t" << mus << "\t" << ps << "\t" << rhoe << "\t" << he << "\t" << mue << "\t" << pe << "\t" << dx << "\t" << rhow << "\t" << hw << "\t" << muw << "\t" << HDlevels[gridref].GridB[i] << "\n";
						system("PAUSE");
					}
				}
			}
		}
	}

	//interpolate the pressure update to the next finer level (if one exists).
	if(gridref!=maxlev){

		//it's only necessary to update the grid points that exist in finer grids, as the ones in coarser grids are in the same locations of those already solved.
		for(int i = HDlevels[gridref].GridP.size();i<HDlevels[gridref+1].GridP.size();i++){
				
			bool dothething = true;
			for(int z = 0;z<boundarynodes.size();z++){
				if(i==boundarynodes[z]){
					dothething = false;
					break;
				}
			}
			//don't update the gap boundaries!!!
			if(dothething){
				//if(HDlevels[gridref].GridC[i]==section){
				//variable declaration
				int GridN,GridS,GridE,GridW,GridNW,GridNE,GridSW,GridSE;	//indices for the surrounding points in all directions.
				int count;													//counter of how many surrounding cells there are.
				double avg = 0.0;											//variable for calculating the sum of pressures in surrounding cells
			
				//identify the surrounding grid points
				GridN = HDlevels[gridref+1].GridN[i];
				GridS = HDlevels[gridref+1].GridS[i];
				GridE = HDlevels[gridref+1].GridE[i];
				GridW = HDlevels[gridref+1].GridW[i];
				GridNE = 0;										//The grid points in non-cardinal directions have not yet been identified.
				GridNW = 0;
				GridSE = 0;
				GridSW = 0;


				//Use a logical search pattern to identify grid points in non-cardinal directions
				if(GridN){
					GridNE = HDlevels[gridref+1].GridE[GridN-1];
					GridNW = HDlevels[gridref+1].GridW[GridN-1];
				}
				if(GridE){
					if((GridNE==0)||(GridNE>HDlevels[gridref].GridP.size()))
						GridNE = HDlevels[gridref+1].GridN[GridE-1];
					GridSE = HDlevels[gridref+1].GridS[GridE-1];
				}
				if(GridW){
					if((GridNW==0)||(GridNW>HDlevels[gridref].GridP.size()))
						GridNW = HDlevels[gridref+1].GridN[GridW-1];
					GridSW = HDlevels[gridref+1].GridS[GridW-1];
				}
				if(GridS){
					if((GridSE==0)||(GridSE>HDlevels[gridref].GridP.size()))
						GridSE = HDlevels[gridref+1].GridE[GridS-1];
					if((GridSW==0)||(GridSW>HDlevels[gridref].GridP.size()))
						GridSW = HDlevels[gridref+1].GridW[GridS-1];
				}

			

				//Find the average of pressure in the surrounding cells
				count = 0;															//initialize the sum to 0
				if((GridN&&(GridN <= HDlevels[gridref].GridP.size()))){		//if there is a cell to the north, and it's pressure has been solved above...
					count++;														//increment the counter
					avg += updatepress[GridN-1];									//add its pressure to the sum.
				}
				if((GridNE&&(GridNE <= HDlevels[gridref].GridP.size()))){	//similar to above for the North East
					count++;
					avg += updatepress[GridNE-1];
				}
				if((GridE&&(GridE <= HDlevels[gridref].GridP.size()))){		//similar to above for the East
					count++;
					avg += updatepress[GridE-1];
				}
				if((GridSE&&(GridSE <= HDlevels[gridref].GridP.size()))){	//similar to above for the South East
					count++;
					avg += updatepress[GridSE-1];
				}
				if((GridS&&(GridS <= HDlevels[gridref].GridP.size()))){		//similar to above for the South
					count++;
					avg += updatepress[GridS-1];
				}
				if((GridSW&&(GridSW <= HDlevels[gridref].GridP.size()))){	//similar to above for the South West
					count++;
					avg += updatepress[GridSW-1];
				}
				if((GridW&&(GridW <= HDlevels[gridref].GridP.size()))){		//similar to above for the West
					count++;
					avg += updatepress[GridW-1];
				}
				if((GridNW&&(GridNW <= HDlevels[gridref].GridP.size()))){	//similar to above for the North West
					count++;
					avg += updatepress[GridNW-1];
				}

				//if there is more than one surrounding point, average the values to get the new value
				//Note: there must be at least two surrounding points because otherwise this point will not have been defined in the first place.
				if(count>1){
					updatepress[i] = avg/(1.0*count);
					if(updatepress[i]!=updatepress[i]){
						cout << "3 " << updatepress[i] << "\n";
					
						system("PAUSE");
					}
				}
			}
			/*else
			{
				Log << "Skipping boundary node " << i << endl;
			}*/
		}
	}
};
void CPistonGap::ReynoldsResidualHD(int gridref, int maxlev,double dt){
	
	//Indices for surrounding grid points.
	int GridN,GridS,GridE,GridW;

	//Solve for the residual of the points in the next coarser grid.
	for(int i = 0;i<HDlevels[gridref-1].GridP.size();i++){
		//Find the indices for the surrounding grid points in the current (fine) grid
		GridN = HDlevels[gridref].GridN[i];
		GridS = HDlevels[gridref].GridS[i];
		GridE = HDlevels[gridref].GridE[i];
		GridW = HDlevels[gridref].GridW[i];

		//This loop should only execute if point i is interior to the region. Do not update pressure on boundaries!
		if(((GridN && GridS) && (GridE && GridW))&&((GridN!=GridS)&&(GridE!=GridW))){

			//only if the point isn't solved for in a finer grid, calculate the residual value
			if(reshere[i]==0){			
				double pn,ps,pe,pw,pp;			//pressure [Pa]
				double hn,hs,he,hw,hp;			//Film thickness [m]
				double mun,mus,mue,muw,mup;		//Viscosity [Pa-s]
				double dx,dy;					//Grid spacing [m]
				double rhon,rhos,rhoe,rhow,rhop;//Density [kg/m^3]
				double e;						//Residual error from Reynolds equation. [???]
				
				//Calculate the grid spacing from the neighboring points.
				dy = (HDlevels[gridref].GridY[GridN-1]-HDlevels[gridref].GridY[GridS-1])/2.0;
				dx = (HDlevels[gridref].GridX[GridE-1]-HDlevels[gridref].GridX[GridW-1])/2.0;
				if(dx < 0.0)
					dx += PI*rK;
				dx = fmod(dx,2.0*PI*rK);
				//cout << i << "\t" << dx << "\t" << dy << "\n";

				if(dx*dy<=0.0){
					Log << "Problem at: " << i << "\n";
				}
				else
				{
					//Compile the properties at grid point i
					pp = HDlevels[maxlev].GridP[i];																	//Pressure
					hp = HDlevels[gridref].GridH[i];											//Film thickness
					mup = HDlevels[gridref].GridV[i];	//Get the viscosity from the oil model
					rhop = HDlevels[gridref].GridD[i];	//Get the density from the oil model
				
					//Gather Information about the point to the North.
					pn = HDlevels[maxlev].GridP[GridN-1];																			//presure
					hn = HDlevels[gridref].GridH[GridN-1];
					hn = (hn + hp)/2.0;//film thickness
					mun = HDlevels[gridref].GridV[GridN-1];	//viscosity
					mun = (mun + mup)/2.0;
					rhon = HDlevels[gridref].GridD[GridN-1];	//density
					rhon = (rhon+rhop)/2.0;
				
					//Similar to above for South
					ps = HDlevels[maxlev].GridP[GridS-1];
					hs = HDlevels[gridref].GridH[GridS-1];
					hs = (hs+hp)/2.0;
					mus = HDlevels[gridref].GridV[GridS-1];
					mus = (mus+mup)/2.0;
					rhos = HDlevels[gridref].GridD[GridS-1];
					rhos = (rhos+rhop)/2.0;
				
					//Similar to above for East
					pe = HDlevels[maxlev].GridP[GridE-1];
					he = HDlevels[gridref].GridH[GridE-1];
					he = (he+hp)/2.0;
					mue = HDlevels[gridref].GridV[GridE-1];
					mue = (mue+mup)/2.0;
					rhoe = HDlevels[gridref].GridD[GridE-1];
					rhoe = (rhoe+rhop)/2.0;
				
					//Similar to above for West
					pw = HDlevels[maxlev].GridP[GridW-1];
					hw = HDlevels[gridref].GridH[GridW-1];
					hw = (hw+hp)/2.0;
					muw = HDlevels[gridref].GridV[GridW-1];
					muw = (muw+mup)/2.0;
					rhow = HDlevels[gridref].GridD[GridW-1];
					rhow = (rhow+rhop)/2.0;
				
					//Calculating the residual error from the Reynolds equation.
					e =  (((rhon*(pn-pp)*pow(hn,3.0)/(12.0*dy*mun))-(rhos*(pp-ps)*pow(hs,3.0)/(12.0*dy*mus)))/(1.0*dy))	//error with Reynolds source terms
						+(((rhoe*(pe-pp)*pow(he,3.0)/(12.0*dx*mue))-(rhow*(pp-pw)*pow(hw,3.0)/(12.0*dx*muw)))/(1.0*dx))
						-0.25 * operatingpistongap.speedK * operatingpistongap.omega * geometrypistongap.dK * (rhoe*he - rhow*hw)/dx
						-0.5 * operatingpistongap.vK * (rhon*hn - rhos*hs)/dy
						+rhop * HDlevels[gridref].GridQ[i]
						+hp * (rhop - HDlevels[gridref].GridA[i])/dt;
					/*e = (1.0/12.0)*(rhon-rhos)/(2.0*dy)*(rhop*pow(hp,3.0)/mup)*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(pow(hn,3.0)-pow(hs,3.0))/(2.0*dy)*(rhop/mup)*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(1.0/mun-1.0/mus)/(2.0*dy)*(rhop*pow(hp,3.0))*(pn-ps)/(2.0*dy)
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(ps+pn-2.0*pp)/pow(dy,2.0)
						+(1.0/12.0)*(rhoe-rhow)/(2.0*dx)*(rhop*pow(hp,3.0)/mup)*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(pow(he,3.0)-pow(hw,3.0))/(2.0*dx)*(rhop/mup)*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(1.0/mue-1.0/muw)/(2.0*dx)*(rhop*pow(hp,3.0))*(pe-pw)/(2.0*dx)
						+(1.0/12.0)*(rhop*pow(hp,3.0)/mup)*(pe+pw-2.0*pp)/pow(dx,2.0)
						+HDlevels[gridref].GridB[i];*/

					//store the residual error in the pressure vector at the next coarser level.
					HDlevels[gridref-1].GridP[i] = e;	

					//mark this level as actually containing pressure rather than residual.
					reshere[i] = gridref;		
					if(e!=e){
						cout << "4 " << e << "\n";
						//cout << rhon << "\t" << hn << "\t" << mun << "\t" << pn << "\t" << pp << "\t" << dy << "\t" << rhos << "\t" << hs << "\t" << mus << "\t" << ps << "\t" << rhoe << "\t" << he << "\t" << mue << "\t" << pe << "\t" << dx << "\t" << rhow << "\t" << hw << "\t" << muw << "\t" << HDlevels[gridref].GridB[i] << "\n";
					
						system("PAUSE");
					}
				}
			}

			//If this point already has a residual calculated, just pass it on.
			else
				HDlevels[gridref-1].GridP[i] = HDlevels[gridref].GridP[i];
		
		}
	}
};
double CPistonGap::CalcRes(vector<double> pold,vector<double> pnew){
	
	//define and initialize a variable for the residual
	double res; res = 0.0;

	//the residual is defined by the point with the largest difference between pold and pnew.
	for(int i = 0;i<pold.size();i++){
		res = max(res,abs(pold[i]-pnew[i]));
	};

	//return the residual value.
	return(res);
};
void CPistonGap::CalcPressureHD(double dt,bool forcebalance){
	pold = p;
	
	//Define Base Grid
	HDlevels.resize(0);
	HDlevels.resize(1);
	//Initialize Vectors
	HDlevels[0].GridD.resize(M*N);
	HDlevels[0].GridH.resize(M*N);
	HDlevels[0].GridP.resize(M*N);
	HDlevels[0].GridV.resize(M*N);
	HDlevels[0].GridX.resize(M*N);
	HDlevels[0].GridY.resize(M*N);
	HDlevels[0].GridN.resize(M*N);
	HDlevels[0].GridS.resize(M*N);
	HDlevels[0].GridE.resize(M*N);
	HDlevels[0].GridW.resize(M*N);
	HDlevels[0].GridQ.resize(M*N);
	HDlevels[0].GridA.resize(M*N);
	//HDlevels[0].GridC.resize(M*N);
	for(int i = 0;i<M*N;i++){
		//Fluid and height values
		HDlevels[0].GridQ[i] = (dht_total(i));
		HDlevels[0].GridA[i] = (density_expansion(i));
		HDlevels[0].GridD[i] = (rho2d(i));
		HDlevels[0].GridH[i] = (h(i));
		HDlevels[0].GridP[i] = (p(i));
		HDlevels[0].GridV[i] = (mu(i));
		HDlevels[0].GridX[i] = (rK*phi(i));
		HDlevels[0].GridY[i] = (dy * (double)(i % M));
		//maxgridy = dy * (double) (M-1);
		//mingridy = 0.0;
		//Connectivity
		if(i%M==(M-1))
			HDlevels[0].GridN[i] = (0);
		else
			HDlevels[0].GridN[i] = (i+2);
		if(i%M==0)
			HDlevels[0].GridS[i] = (0);
		else
			HDlevels[0].GridS[i] = (i);
		if(i>=(N-1)*M)
			HDlevels[0].GridE[i] = (i-(N-1)*M+1);
		else
			HDlevels[0].GridE[i] = (i+M+1);
		if(i<M)
			HDlevels[0].GridW[i] = (i+(N-1)*M+1);
		else
			HDlevels[0].GridW[i] = (i-M+1);
		//HDlevels[0].GridC[i] = 1;
	}
	//Set neighbors to zero when in groove to prevent solving or refining.
	for(int i = 0;i<N*M;i++){
		for(int j = 0;j<numgvlizhi;j++){
			if((i%M < nlimitlizhi(j))&&(i%M > slimitlizhi(j))){
				HDlevels[0].GridN[i] = 0;
				HDlevels[0].GridS[i] = 0;
				HDlevels[0].GridE[i] = 0;
				HDlevels[0].GridW[i] = 0;
			}
		}
	}

	//if force balance loop, use previous grid refinement
	if(!forcebalance){
		refinememory.resize(0);
	}
	refinecounter = 0;
	//refine mesh
	boundarynodes.resize(0);
	bool refined = true;
	int gridref = 0;
	do{
		HDlevels.resize(gridref+2);
		refined = RefineMesh(gridref,forcebalance);
		gridref++;
		//if(refined>5)
			//refined = false;
		//Log << "Refining\t";
	}while(refined);
	HDlevels.resize(gridref);
	//cout << gridref << " Refinements\n";
		
	gridref--;
	levelrequest = levelrequest<gridref? gridref:levelrequest;
	//Log << "Done\n";

	/*fout.open("./output/piston/matlab/GridCheck.txt",ios::app);
	for(int j = 0;j<=gridref;j++)
		for( int i = 0;i<HDlevels[j].GridN.size();i++)
			fout << j << "\t" << i << "\t" << HDlevels[j].GridN[i] << "\t" << HDlevels[j].GridS[i] << "\t" << HDlevels[j].GridE[i] << "\t" << HDlevels[j].GridW[i] << "\t" << HDlevels[j].GridX[i] << "\t" << HDlevels[j].GridY[i] << "\n";
	fout.close();
	fout.clear();
	system("PAUSE");*/

	//calculate pressure
	//Resize variables
	reshere.resize(HDlevels[gridref].GridP.size());
	updatepress.resize(HDlevels[gridref].GridP.size());

	bool keepgoing = true;
	double intres;
	vector<double> oldP;
	int j = 0;
		
	oldP.resize(HDlevels[gridref].GridP.size());
	oldP = HDlevels[gridref].GridP;
	//Make five passes through the multigrid method.
	do{
	//for(int trash = 0;trash<5;trash++){
		//Reset the book-keeping
		for(int i = 0;i<HDlevels[gridref].GridP.size();i++){
			reshere[i] = 0;
			updatepress[i] = 0.0;
		}
		
		//Calculate the residual value at each level, going from the finest level to the coarsest level.
		for(int i = gridref;i>0;i--)
			ReynoldsResidualHD(i,gridref,dt);
		
		//Solve the updated pressure value at each level going from the coarsest level to the finest.
		for(int i = 0;i<=gridref;i++)
			SolveReynoldsHD(i,gridref,0.7,dt);

		j++;
		intres = CalcRes(oldP,HDlevels[gridref].GridP);
		//Log << intres << "\n";
		if(intres<(Rmin_R*1e6))
			keepgoing = false;
		if(j>100000)
			keepgoing = false;

		oldP = HDlevels[gridref].GridP;
	//}
	}while(keepgoing);
		

	//interpolate back pressure
	int GridN,GridS,GridE,GridW,GridNW,GridNE,GridSW,GridSE;
	double avg;
	for(int j = gridref;j>0;j--){
		for(int i = 0;i<HDlevels[j-1].GridP.size();i++){
			GridN = 0;
			GridS = 0;
			GridE = 0;
			GridW = 0;
			GridNW = 0;
			GridNE = 0;
			GridSE = 0;
			GridSW = 0;
			GridN = HDlevels[j].GridN[i];
			GridS = HDlevels[j].GridS[i];
			GridE = HDlevels[j].GridE[i];
			GridW = HDlevels[j].GridW[i];
			
			if(GridN)
				GridNE = HDlevels[j].GridE[GridN-1];
			else if(GridE)
				GridNE = HDlevels[j].GridN[GridE-1];
			if(GridN)
				GridNW = HDlevels[j].GridW[GridN-1];
			else if(GridW)
				GridNW = HDlevels[j].GridN[GridW-1];
			if(GridS)
				GridSE = HDlevels[j].GridE[GridS-1];
			else if(GridE)
				GridSE = HDlevels[j].GridS[GridE-1];
			if(GridS)
				GridSW = HDlevels[j].GridW[GridS-1];
			else if(GridW)
				GridSW = HDlevels[j].GridS[GridW-1];
			
			avg = 1.0;

			HDlevels[j-1].GridP[i] = 0.0;

			if(GridN){
				avg -= 0.125;
				HDlevels[j-1].GridP[i] += 0.125 * HDlevels[j].GridP[GridN-1];
			}
			if(GridS){
				avg -= 0.125;
				HDlevels[j-1].GridP[i] += 0.125 * HDlevels[j].GridP[GridS-1];
			}
			if(GridE){
				avg -= 0.125;
				HDlevels[j-1].GridP[i] += 0.125 * HDlevels[j].GridP[GridE-1];
			}
			if(GridW){
				avg -= 0.125;
				HDlevels[j-1].GridP[i] += 0.125 * HDlevels[j].GridP[GridW-1];
			}
			if(GridNW){
				avg -= 0.0625;
				HDlevels[j-1].GridP[i] += 0.0625 * HDlevels[j].GridP[GridNW-1];
			}
			if(GridNE){
				avg -= 0.0625;
				HDlevels[j-1].GridP[i] += 0.0625 * HDlevels[j].GridP[GridNE-1];
			}
			if(GridSW){
				avg -= 0.0625;
				HDlevels[j-1].GridP[i] += 0.0625 * HDlevels[j].GridP[GridSW-1];
			}
			if(GridSE){
				avg -= 0.0625;
				HDlevels[j-1].GridP[i] += 0.0625 * HDlevels[j].GridP[GridSE-1];
			}
			HDlevels[j-1].GridP[i] += avg * HDlevels[j].GridP[i];
			//cout << avg << endl;
		}
	}
	for(int i = 0;i<N*M;i++){
		//if(HDlevels[gridref].GridP[i] == HDlevels[gridref].GridP[i])
			p(i) = HDlevels[0].GridP[i];
	}

};