#include "CSlipperGap.h"
#include "../../caspar_slipper_fem/src/CasparSlipperFEM.h"
#include <numeric>

using namespace CasparSlipperFEM;

void CSlipperGap::CalcSlipperThermal(const double alpha)
{
	GapLog << endl << "Slipper Thermal Analysis..." << endl;

	//fluid2slipper - no longer needed because the interpolation happens at the end of each timestep
	//Fluid2Slipper_Thermal();

	//Create the average gap flux values
	for(int i=0; i<t_slipper.faceFlux.size(); i++)
	{
		if(t_slipper.faceFluxcnt > 0)
		{
			//add in the flux relaxation
			t_slipper.faceFlux[i] = t_slipper.faceFlux_old[i] + alpha*(t_slipper.faceFlux[i]/t_slipper.faceFluxcnt - t_slipper.faceFlux_old[i]);
		} else {
			//there hasn't been any flux applied, so use a zero flux condition
			t_slipper.faceFlux[i] = 0.0;
		}
	}
	
	//can update the old flux values now
	t_slipper.faceFlux_old = t_slipper.faceFlux;

	//perform the analysis
	thermoelastic slipper_te;
	slipper_te.option_file = t_slipper.option_file;
	slipper_te.caspar_gap_input = gapinput;
	slipper_te.vtk_file = "./output/slipper/vtk/slipper_thermal." + ::n2s(operatingslippergap.phi_deg) + ".vtk";
	
	//resize the communication vectors
	slipper_te.heat_flux.resize(t_slipper.facecnt);

	//set the faceFlux from t_slipper to the analysis heat_flux
	for(int i=0; i<t_slipper.facecnt; i++)
	{
		slipper_te.heat_flux[i] = t_slipper.faceFlux[i];
	}

	//If this is the first thermoelastic solution we will accecpt the full temperature solution and not underrelax
	if(accumulate(t_slipper.nodeTemperature_old_FullMesh.begin(), t_slipper.nodeTemperature_old_FullMesh.end(), static_cast<double>(0)) == 0)
	{
		//the heat flux is being applied to the gap faceset
		slipper_te.run_analysis("gap");
	} else {
		//the heat flux is being applied to the gap faceset
		slipper_te.run_analysis("gap", gapinput->options_slipper.numeric.AlphaTEHD, t_slipper.nodeTemperature_old_FullMesh);
	}

	//update the old full nodal temperature
	t_slipper.nodeTemperature_old_FullMesh = slipper_te.node_temp;

	//copy the gap temp / deform from the analysis back to t_slipper
	for(int i=0; i<t_slipper.nodecnt; i++)
	{
		const int Gnid = t_slipper.nodes[i].Gid;
		t_slipper.nodeTemperature[i] = slipper_te.node_temp[Gnid];
		
		//get the z deformation
		t_slipper.nodeDeformation[i] = slipper_te.node_deform[Gnid][2];
	}

	//interpolate the temp & deformation back to the fluid grid (slipper2fluid)
	Slipper2Fluid_Thermal(alpha);

	//Zero the gap flux vector and counter
	for(int i=0; i<t_slipper.faceFlux.size(); i++)
	{
		t_slipper.faceFlux[i] = 0;
	}
	t_slipper.faceFluxcnt = 0;

	GapLog << endl << "Slipper Thermal Analysis Finished!" << endl;
}
void CSlipperGap::CalcSwashplateThermal(const double alpha)
{
	GapLog << endl << "Swashplate Thermal Analysis..." << endl;

	//fluid2slipper - no longer needed because the interpolation happens at the end of each timestep
	//Fluid2Slipper_Thermal();

	//Create the average gap flux values
	for(int i=0; i<t_swashplate.faceFlux.size(); i++)
	{
		if(t_swashplate.faceFluxcnt > 0)
		{
			//add in the flux relaxation
			double newflux = double(gapinput->operating_conditions.npistons)*t_swashplate.faceFlux[i]/t_swashplate.faceFluxcnt;
			t_swashplate.faceFlux[i] = t_swashplate.faceFlux_old[i] + alpha*(newflux - t_swashplate.faceFlux_old[i]);
		} else {
			//there hasn't been any flux applied, so use a zero flux condition
			t_swashplate.faceFlux[i] = 0.0;
		}
	}
	
	//can update the old flux values now
	t_swashplate.faceFlux_old = t_swashplate.faceFlux;

	//perform the analysis
	thermoelastic swashplate_te;
	swashplate_te.caspar_gap_input = gapinput;
	swashplate_te.option_file = t_swashplate.option_file;
	swashplate_te.vtk_file = "./output/slipper/vtk/swashplate_thermal." + ::n2s(operatingslippergap.phi_deg) + ".vtk";
	
	//resize the communication vectors
	swashplate_te.heat_flux.resize(t_swashplate.facecnt);

	//set the faceFlux from t_swashplate to the analysis heat_flux
	for(int i=0; i<t_swashplate.facecnt; i++)
	{
		swashplate_te.heat_flux[i] = t_swashplate.faceFlux[i];
	}

	//If this is the first thermoelastic solution we will accecpt the full temperature solution and not underrelax
	if(accumulate(t_swashplate.nodeTemperature_old_FullMesh.begin(), t_swashplate.nodeTemperature_old_FullMesh.end(), static_cast<double>(0)) == 0)
	{
		//the heat flux is being applied to the gap faceset
		swashplate_te.run_analysis("gap");
	} else {
		//the heat flux is being applied to the gap faceset
		swashplate_te.run_analysis("gap", gapinput->options_slipper.numeric.AlphaTEHD, t_swashplate.nodeTemperature_old_FullMesh);
	}

	//update the old full nodal temperature
	t_swashplate.nodeTemperature_old_FullMesh = swashplate_te.node_temp;

	//copy the gap temp / deform from the analysis back to t_slipper
	for(int i=0; i<t_swashplate.nodecnt; i++)
	{
		const int Gnid = t_swashplate.nodes[i].Gid;
		t_swashplate.nodeTemperature[i] = swashplate_te.node_temp[Gnid];
		
		//get the z deformation
		t_swashplate.nodeDeformation[i] = swashplate_te.node_deform[Gnid][2];
	}

	//interpolate the temp & deformation back to the fluid grid (swashplate2fluid)
	Swashplate2Fluid_Thermal(alpha);

	//Zero the gap flux vector and counter
	for(int i=0; i<t_swashplate.faceFlux.size(); i++)
	{
		t_swashplate.faceFlux[i] = 0;
	}
	t_swashplate.faceFluxcnt = 0;

	GapLog << endl << "Swashplate Thermal Analysis Finished!" << endl;
}
void CSlipperGap::Fluid2Slipper_Thermal(void)
{
	ANNpoint qPt = annAllocPt(2);
	//Interpolation points
	int k = 1;
	ANNidxArray	  nnIdx = new ANNidx[k];
	ANNdistArray  dists = new ANNdist[k];

	for(int f=0; f<t_slipper.facecnt; f++)
	{

		qPt[0] = t_slipper.faces[f].x;
		qPt[1] = t_slipper.faces[f].y;

		Fluid.KDslip.kdtree->annkSearch
		(
			qPt,
			k,
			nnIdx,
			dists,
			Fluid.KDslip.eps
		);	

		//The closest fluid volume centroid
		int c1 = nnIdx[0];
		int m1 = nnIdx[0]/Fluid.N;
		int n1 = nnIdx[0]%Fluid.N;

		//We need to find two other points to use
		int m2, n2;
		int m3, n3;

		double x = qPt[0];
		double y = qPt[1];
		double x1 = Fluid.LRx(m1,n1);
		double y1 = Fluid.LRy(m1,n1);

		//to polar!
		double r = sqrt(x*x+y*y);
		double t = atan2(y,x);
		double r1 = sqrt(x1*x1+y1*y1);
		double t1 = atan2(y1,x1);

		if(r>r1)
		{
			m2 = m1+1;
			//Special checks to use the 2nd inner point if at a boundary
			m2 = m2<0 ? 1 : m2;
			m2 = m2>=Fluid.M ? Fluid.M-2 : m2;

			n2 = n1;
		} else {
			m2 = m1-1;
			//Special checks to use the 2nd inner point if at a boundary
			m2 = m2<0 ? 1 : m2;
			m2 = m2>=Fluid.M ? Fluid.M-2 : m2;

			n2 = n1;
		}
			
		double x2 = Fluid.LRx(m2,n2);
		double y2 = Fluid.LRy(m2,n2);

		if(t>t1)
		{
			m3 = m1;
			n3 = n1+1;

			//Special trick if n is at a boundary since n is "wraped"
			n3 = (n3+Fluid.N)%Fluid.N;
		} else {
			m3 = m1;			
			n3 = n1-1;

			//Special trick if n is at a boundary since n is "wraped"
			n3 = (n3+Fluid.N)%Fluid.N;
		}

		double x3 = Fluid.LRx(m3,n3);
		double y3 = Fluid.LRy(m3,n3);

		//Barycentric Interpolation
		//http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)

		const double detT = (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
		const double lam1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
		const double lam2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
		const double lam3 = 1.0 - lam1 - lam2;

		if(detT == 0)
		{
			//Two of the points are the same
			//This really should not happen with the current scheme

			//If it does an error needs to be thrown as interpolation will fail with NaN
		}

		if(lam1 < 0 || lam2 < 0 || lam3 < 0)
		{
			//The query point is not bounded by the interpolation points
				
			const double lamS = fabs(lam1) + fabs(lam2) + fabs(lam3);
			
			//This is only a big deal if lamS is significantly greater than 1.
			//A warning should be thrown if it is greater than 2-3??

			if(lamS > 3.0)
			{
				//cout << "WARNING: Forcing F2S f: " << f << " with lamS: " << lamS << endl;
				GapLog.message("WARNING: There is a mismatch between the fluid and slipper thermal mesh. Forcing F2ST interpolation at f=" + ::n2s(f) + ", lamS=" + ::n2s(lamS));
			}
		}

		//interp
		t_slipper.faceFlux [f] +=	lam1*Fluid.Qflux(m1,n1) +
											lam2*Fluid.Qflux(m2,n2) +
											lam3*Fluid.Qflux(m3,n3);
	}

	t_slipper.faceFluxcnt += 1;

	annDeallocPt(qPt);
	delete [] nnIdx;
	delete [] dists;
}
void CSlipperGap::Slipper2Fluid_Thermal(const double alpha)
{
	//how many pts should be tried for interpolation boundedness
	int maxBoundTry = 10;
	ANNpoint qPt = annAllocPt(2);

	for(int m=0; m<Fluid.M; m++)
	{
		for(int n=0; n<Fluid.N; n++)
		{

			qPt[0] = Fluid.LRx(m,n);
			qPt[1] = Fluid.LRy(m,n);

			//Used to find the "best" pts in case of extrapolation
			vector<double> lamS;

			//Interpolation points
			for(int k = 1; k<=maxBoundTry+1; k++)
			{
				bool force = false;
				if(k == maxBoundTry+1)
				{
					//The closest `maxBoundTry` points didn't bound the query pt
					//so just force the min(sum(abs(lam))) one and extrapolate

					force = true;

					//Use whatever k value had the lowest lamS
						
					double minlamS = lamS[0];
					k = 1;
					for(unsigned int i=1; i<lamS.size(); i++)
					{
						if(lamS[i] < minlamS)
						{
							k = i+1;
							minlamS = lamS[i];
						}
					}
						
				}

				ANNidxArray	  nnIdx = new ANNidx[k];
				ANNdistArray  dists = new ANNdist[k];

				t_slipper.KDfaces.kdtree->annkSearch
				(
					qPt,
					k,
					nnIdx,
					dists,
					t_slipper.KDfaces.eps
				);	


				//Barycentric Interpolation
				//http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)

				const double x = qPt[0];
				const double y = qPt[1];
				const double x1 = t_slipper.nodes[t_slipper.faces[nnIdx[k-1]].nodes[0]].x;
				const double y1 = t_slipper.nodes[t_slipper.faces[nnIdx[k-1]].nodes[0]].y;
				const double x2 = t_slipper.nodes[t_slipper.faces[nnIdx[k-1]].nodes[1]].x;
				const double y2 = t_slipper.nodes[t_slipper.faces[nnIdx[k-1]].nodes[1]].y;
				const double x3 = t_slipper.nodes[t_slipper.faces[nnIdx[k-1]].nodes[2]].x;
				const double y3 = t_slipper.nodes[t_slipper.faces[nnIdx[k-1]].nodes[2]].y;

				const double detT = (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
				const double lam1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
				const double lam2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
				const double lam3 = 1.0 - lam1 - lam2;

				if(force)
				{
					const double lamS = fabs(lam1) + fabs(lam2) + fabs(lam3);
						
					//This is only a big deal if lamS is significantly greater than 1.
					//A warning should be thrown if it is greater than 2-3??
						
					if(lamS > 3.0)
					{
						//Only throw an error message for non-boundary volumes
						if(Fluid.boundary(m,n) == -1)
						{
							//cout << "WARNING: Forcing S2F m: " << m << " n: " << n << " with lamS: " << lamS << endl;
							GapLog.message("WARNING: There is a mismatch between the fluid and slipper thermal mesh. Forcing ST2F interpolation at m=" + ::n2s(m) + ", n=" + ::n2s(n) + ", lamS=" + ::n2s(lamS));
						}
					}
				}

				if( (lam1 < 0 || lam2 < 0 || lam3 < 0 || detT == 0) && !force )
				{
					lamS.push_back(fabs(lam1) + fabs(lam2) + fabs(lam3));
						
					delete [] nnIdx;
					delete [] dists;
					
					continue;
				}

				//Deformation
				double newdeform =
								lam1*t_slipper.nodeDeformation[t_slipper.faces[nnIdx[k-1]].nodes[0]] + 
								lam2*t_slipper.nodeDeformation[t_slipper.faces[nnIdx[k-1]].nodes[1]] + 
								lam3*t_slipper.nodeDeformation[t_slipper.faces[nnIdx[k-1]].nodes[2]];

				//Temperature
				double newtemp =
								lam1*t_slipper.nodeTemperature[t_slipper.faces[nnIdx[k-1]].nodes[0]] + 
								lam2*t_slipper.nodeTemperature[t_slipper.faces[nnIdx[k-1]].nodes[1]] + 
								lam3*t_slipper.nodeTemperature[t_slipper.faces[nnIdx[k-1]].nodes[2]];

				//damp the new deformation
				t_slipper.deform(m,n) = newdeform;
				t_slipper.temp(m,n) = newtemp;

				//

				//sucuessful so
				delete [] nnIdx;
				delete [] dists;
				break;
			}

		}
	}

	//clean up the ann pt
	annDeallocPt(qPt);

	//calculate the pdeform (shifted deform)
	t_slipper.pdeform = t_slipper.deform - min(t_slipper.deform);
}
void CSlipperGap::Fluid2Swashplate_Thermal(void)
{
	ANNpoint qPt = annAllocPt(2);

	vector<double> tmpflux(t_swashplate.facecnt, 0);
	vector<double> tmpfluxcnt(t_swashplate.facecnt, 0);
	
	for(int m=0; m<Fluid.M; m++)
	{
		for(int n=0; n<Fluid.N; n++)
		{

			qPt[0] = Fluid.Gx(m,n);
			qPt[1] = Fluid.Gy(m,n);

			int k=1;
			ANNidxArray	  nnIdx = new ANNidx[k];
			ANNdistArray  dists = new ANNdist[k];

			t_swashplate.KDfaces.kdtree->annkSearch
			(
				qPt,
				k,
				nnIdx,
				dists,
				t_swashplate.KDfaces.eps
			);	

			tmpflux[nnIdx[k-1]] += Fluid.Qflux(m,n);
			tmpfluxcnt[nnIdx[k-1]] += 1;
			
			delete [] nnIdx;
			delete [] dists;
		}
	}

	//take the 'interpolation' flux average and assign it to each face in the thermal swashplate object
	for(int i=0; i<t_swashplate.facecnt; i++)
	{
		if(tmpfluxcnt[i] > 0)
		{
			t_swashplate.faceFlux[i] += (tmpflux[i] / tmpfluxcnt[i]);
		}
	}

	t_swashplate.faceFluxcnt += 1;
	
	//clean up the ann pt
	annDeallocPt(qPt);
}
void CSlipperGap::Swashplate2Fluid_Thermal(const double alpha)
{
	//how many pts should be tried for interpolation boundedness
	int maxBoundTry = 10;
	ANNpoint qPt = annAllocPt(2);

	for(int m=0; m<Fluid.M; m++)
	{
		for(int n=0; n<Fluid.N; n++)
		{

			qPt[0] = Fluid.Gx(m,n);
			qPt[1] = Fluid.Gy(m,n);

			//Used to find the "best" pts in case of extrapolation
			vector<double> lamS;

			//Interpolation points
			for(int k = 1; k<=maxBoundTry+1; k++)
			{
				bool force = false;
				if(k == maxBoundTry+1)
				{
					//The closest `maxBoundTry` points didn't bound the query pt
					//so just force the min(sum(abs(lam))) one and extrapolate

					force = true;

					//Use whatever k value had the lowest lamS
						
					double minlamS = lamS[0];
					k = 1;
					for(unsigned int i=1; i<lamS.size(); i++)
					{
						if(lamS[i] < minlamS)
						{
							k = i+1;
							minlamS = lamS[i];
						}
					}
						
				}

				ANNidxArray	  nnIdx = new ANNidx[k];
				ANNdistArray  dists = new ANNdist[k];

				t_swashplate.KDfaces.kdtree->annkSearch
				(
					qPt,
					k,
					nnIdx,
					dists,
					t_swashplate.KDfaces.eps
				);	


				//Barycentric Interpolation
				//http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)

				const double x = qPt[0];
				const double y = qPt[1];
				const double x1 = t_swashplate.nodes[t_swashplate.faces[nnIdx[k-1]].nodes[0]].x;
				const double y1 = t_swashplate.nodes[t_swashplate.faces[nnIdx[k-1]].nodes[0]].y;
				const double x2 = t_swashplate.nodes[t_swashplate.faces[nnIdx[k-1]].nodes[1]].x;
				const double y2 = t_swashplate.nodes[t_swashplate.faces[nnIdx[k-1]].nodes[1]].y;
				const double x3 = t_swashplate.nodes[t_swashplate.faces[nnIdx[k-1]].nodes[2]].x;
				const double y3 = t_swashplate.nodes[t_swashplate.faces[nnIdx[k-1]].nodes[2]].y;

				const double detT = (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
				const double lam1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
				const double lam2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
				const double lam3 = 1.0 - lam1 - lam2;

				if(force)
				{
					const double lamS = fabs(lam1) + fabs(lam2) + fabs(lam3);
						
					//This is only a big deal if lamS is significantly greater than 1.
					//A warning should be thrown if it is greater than 2-3??
						
					if(lamS > 3.0)
					{
						//Only throw an error message for non-boundary volumes
						if(Fluid.boundary(m,n) == -1)
						{
							//cout << "WARNING: Forcing S2F m: " << m << " n: " << n << " with lamS: " << lamS << endl;
							GapLog.message("WARNING: There is a mismatch between the fluid and swashplate thermal mesh. Forcing ST2F interpolation at m=" + ::n2s(m) + ", n=" + ::n2s(n) + ", lamS=" + ::n2s(lamS));
						}
					}
				}

				if( (lam1 < 0 || lam2 < 0 || lam3 < 0 || detT == 0) && !force )
				{
					lamS.push_back(fabs(lam1) + fabs(lam2) + fabs(lam3));
						
					delete [] nnIdx;
					delete [] dists;
					
					continue;
				}

				//Deformation
				double newdeform =
								lam1*t_swashplate.nodeDeformation[t_swashplate.faces[nnIdx[k-1]].nodes[0]] + 
								lam2*t_swashplate.nodeDeformation[t_swashplate.faces[nnIdx[k-1]].nodes[1]] + 
								lam3*t_swashplate.nodeDeformation[t_swashplate.faces[nnIdx[k-1]].nodes[2]];

				//Temperature
				double newtemp =
								lam1*t_swashplate.nodeTemperature[t_swashplate.faces[nnIdx[k-1]].nodes[0]] + 
								lam2*t_swashplate.nodeTemperature[t_swashplate.faces[nnIdx[k-1]].nodes[1]] + 
								lam3*t_swashplate.nodeTemperature[t_swashplate.faces[nnIdx[k-1]].nodes[2]];


				//damp the new deformation
				t_swashplate.deform(m,n) = newdeform;
				t_swashplate.temp(m,n) = newtemp;

				//

				//sucuessful so
				delete [] nnIdx;
				delete [] dists;
				break;
			}

		}
	}

	//clean up the ann pt
	annDeallocPt(qPt);

	//calculate the pdeform (shifted deform)
	//not active for now.. if so remember that the 'shift' needs to be constant over a revolution
	//t_swashplate.pdeform = t_swashplate.deform - min(t_swashplate.deform);
	t_swashplate.pdeform = t_swashplate.deform;
}
void CSlipperGap::sTSolid::load(const string Option_File)
{
	//set the option file for this thermal solid object
	option_file = Option_File;

	*GapLog << endl << "Loading Thermal Inputs from " << option_file << " ..." << endl;
	thermoelastic slipper_te;
	slipper_te.option_file = option_file;
	
	//used to form our internal faceset
	vector<int> nid;
	vector<vector<double> > nodexyz;
	vector<vector<int> > setfaces;

	fem * body = slipper_te.open();
	slipper_te.get_faceset(body, "gap", nid, nodexyz, setfaces);
	unsigned int total_body_nodecnt = slipper_te.get_nodecnt(body);
	slipper_te.close(body);

	//populate the nodes
	for(int i=0; i<nid.size(); i++)
	{
		nodes.push_back(node(nodexyz[i][0], nodexyz[i][1], nid[i]));
	}

	//populate the faces
	for(int i=0; i<setfaces.size(); i++)
	{
		//find the centroid
		double cx = 0;
		double cy = 0;
		for(int j=0; j<3; j++)
		{
			const int id = setfaces[i][j];
			cx += nodes[id].x / 3.0;
			cy += nodes[id].y / 3.0;
		}

		faces.push_back(face(setfaces[i][0], setfaces[i][1], setfaces[i][2], cx, cy));
	}

	facecnt = (int) faces.size();
	nodecnt = (int) nodes.size();

	nodeDeformation.resize(nodecnt, 0);
	nodeTemperature.resize(nodecnt, 0);
	faceFlux.resize(facecnt, 0);
	faceFlux_old.resize(facecnt, 0);
	faceFluxcnt = 0;

	nodeTemperature_old_FullMesh.resize(total_body_nodecnt, 0);

	//ReSize the faces KD tree
	KDfaces.points = annAllocPts(facecnt,KDfaces.dim);

	//Load the points
	for(int f=0; f<facecnt; f++)
	{
		KDfaces.points[f][0] = faces[f].x;
		KDfaces.points[f][1] = faces[f].y;
	}

	//Build the tree
	KDfaces.kdtree = new ANNkd_tree(KDfaces.points, facecnt, KDfaces.dim);

	*GapLog << "Loaded " << facecnt << " thermal gap faces and " << nodecnt << " thermal gap nodes." << endl << endl;



}
