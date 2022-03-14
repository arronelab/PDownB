#include "ktlMoleculeRandom.h"
#include "hydrationShellRandom.h"
#include "experimentalData.h"
#include <string.h> 


/***********************************************
  argv[1] scattering data file
  argv[2] sequence file
  argv[3] initial prediction coords  file (can be empty)
  argv[4] paired distances file (can be empty)
  argv[5] fixed sections file (again can be empty)
  argv[6] crystal symmetry file again can be empty
  argv[7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it. Will say none if not.
  argv[8] request to apply hydrophobic covering BETWEEN monomers will be a list of pairs to try to hydropobically pair. Will say none if not.
  argv[9] kmin
  argv[10] kmax
  argv[11] prediction file 
 **********************************************/

bool checkTransition(double &chiSqVal,double &chiSqCurr,double &uniformProb,int index,int &maxSteps){
  /*double tempFrac;
  if(index < maxSteps/2){
    tempFrac= 1.0-double(index)/double(maxSteps/2);
  }else{
    tempFrac=0.0000000001;
  }
  double annealProb;
  double tempVal;
  if(chiSqVal<chiSqCurr){
    annealProb=1.0;
  }else{
    annealProb = std::exp(-(chiSqVal-chiSqCurr)/tempFrac);
  }
  //std::cout<<index<<" "<<maxSteps<<" "<<chiSqVal<<" "<<chiSqCurr<<" "<<annealProb<<" "<<uniformProb<<"\n";
  if(annealProb>uniformProb){
    return true;
  }else{
    return false;
  }*/
  if(chiSqVal<chiSqCurr){
    return true;
  }else{
    return false;
    }
}

int main( int argc, const char* argv[] )
{  
 /*************************************

  set up model parameters 
 
  *************************************/
  
  double lmin=4.0; // closest distance two non adjactent local (same secondary unit) moelcules can get
  double rmin=3.7;double rmax=3.9; // max and min calpha-Dists
  double closestApproachDist=3.9; // closest distance two non adjactent non local moelcules (different secondary unit) can get

  /*************************************
  
   determine initial model 

   *************************************/

  ktlMolecule mol;
  if(strcmp(argv[3],"no_initial_prediction")==0){
    // read in from sequernce and sc structure pred
    mol.readInSequence(argv[2],rmin,rmax,lmin);
    // generate random start (with no overlap)
    mol.getRandomMolecule();
    mol.writeMoleculeToFile("testMol.dat");
    // identify hydrophobic residues
    mol.getHydrophobicResidues();
  }else{
    // read in from sequernce and sc structure pred
    mol.readInSequence(argv[2],rmin,rmax,lmin);
    // read in coordinates
    mol.readInCoordinates(argv[3]);
    //here would be the filling in missing section routine
    // identify hydrophobic residues
    mol.writeMoleculeToFile("testMol.dat");
    mol.getHydrophobicResidues();
  }
  bool doAll = false;
  /*********************************
   
   check if we are restricting sections which are to be altered
   
   ********************************************/
  std::vector<int> fixedSecList;
   if(strcmp(argv[5],"none")==0){
    // here we make no restrictions
    doAll =true;
  }else{
    std::ifstream fixedSecFile;
    fixedSecFile.open(argv[5]);
    std::string line;int index;
    if(fixedSecFile.is_open()){
      while(!fixedSecFile.eof()){
	std::getline(fixedSecFile,line);
	std::stringstream ss(line);
	ss>>index;
	fixedSecList.push_back(index);
      }
    }else{
    }
    fixedSecFile.close();
  }
  
  /*********************************************
   
   check if we need to apply hydrophobic covering
   
   ********************************************/

  //read in the internal hydrophobic list
  std::vector<int> internalHydrophicChecklist;
  if(strcmp(argv[7],"none")==-1){
    for(const char* it=argv[7];*it;++it){
      internalHydrophicChecklist.push_back(*it-'0');
    }
  }
  // enforce any requested hyrdophobic covering
  // std::cout<<" any hydro aminos ? "<<internalHydrophicChecklist.size()<<"\n";
  for(int i=0;i<internalHydrophicChecklist.size();i++){
    // start random genrator
    std::random_device rdev{};
    std::default_random_engine generator1{rdev()};
    std::uniform_real_distribution<double> distributionR(0.0,1.0);
    // select max number of steps
    int noHydrationCoverSteps=100;
    // get initial hydration value first choose no of neighbours to check for nearest distances
    int noNeighbours =10;
    double  currGlobalRad = mol.getGlobalRadiusOfCurvatureWithinSec(internalHydrophicChecklist[i],noNeighbours);
    double globalRadFit = 0.01*(currGlobalRad-7.5)*(currGlobalRad-7.5);
    std::vector<double> ovelps = mol.checkOverlapWithRad(closestApproachDist,internalHydrophicChecklist[i]); 
    double distSum=0.0;
    for(int l=0;l<ovelps.size();l++){
      distSum = distSum +  std::abs(closestApproachDist-ovelps[l]);
    }   
    globalRadFit =globalRadFit + 0.01*distSum;
    int k=0;
    // now  
    while(k<noHydrationCoverSteps && globalRadFit>0.0001){
      for(int j=0;j<mol.getSubsecSize(internalHydrophicChecklist[i]);j++){
	ktlMolecule molCopy = mol;
	molCopy.changeMoleculeSingleMulti(j,internalHydrophicChecklist[i]);
	double newGlobalRad = molCopy.getGlobalRadiusOfCurvatureWithinSec(internalHydrophicChecklist[i],noNeighbours); 
        std::vector<double> ovelps = molCopy.checkOverlapWithRad(closestApproachDist,internalHydrophicChecklist[i]);   
	double newGlobalRadFit = 0.01*(newGlobalRad-7.5)*(newGlobalRad-7.5);
        double distSum=0.0;
        for(int l=0;l<ovelps.size();l++){
	  distSum = distSum + std::abs(closestApproachDist-ovelps[l]);
	}   
        newGlobalRadFit =newGlobalRadFit + 0.01*distSum;
        double uProb = distributionR(generator1);
	if(checkTransition(newGlobalRadFit,globalRadFit,uProb,k,noHydrationCoverSteps)){
	  currGlobalRad = newGlobalRad;
	  globalRadFit=  newGlobalRadFit;       
	  mol=molCopy;
	}
      }
      k++;
    }
  }
  mol.writeMoleculeToFile(argv[11]);

  /******************************************

     read in the scattering and set up the scattering model

   ******************************************/
  experimentalData ed(argv[1]);
  
  /*
   ed.generatePR();
   */

  /*************************************************

    generate the hydration layer model
   
   *************************************************/
  // the hydration shell parameters
  double Rin= 6.0;
  double Rout=7.0;
  double RShell = 5.5;
  int ntrivs=6;
  double helixRatio=1.0;
  int solventsPerLink =1;
  
  hydrationShellMinimal hydrationShell(mol,Rin,Rout,RShell,ntrivs,helixRatio,solventsPerLink,closestApproachDist,rmin,rmax,lmin);

  // generate the hydration shell
  
  hydrationShell.tubeParamList();
  hydrationShell.constructInitialState();
  hydrationShell.allOverlap();

  hydrationShell.writeHydrationShellToFile("testLyzScat.dat");
  
  
  /********************************************

     load in additional information

   *******************************************/

  // contact predictions or similar distance constraints
  mol.loadContactPredictions(argv[4]);
  
  // To do: crystal symmetry
  

  /*******************************************
    establish the initial fit
   ******************************************/

  // Generate the initial distance distribution, grab all distances
  
  // check internal overlap of sections and fill overlap list
  
  std::vector<double> molDists;
  std::vector<double> overlapDists= mol.checkOverlapWithRad(closestApproachDist);

  molDists = mol.getDistSet();
  int molSize = mol.getNoAminos();
  std::sort(molDists.begin(),molDists.end());
   
  // get solvent molecules
  
  std::vector<std::vector<point> > solpts = hydrationShell.returnFlatSolList();
  
  // calculate withing solvent list;
  
  std::vector<double> solDists = mol.getExactDistSet(solpts);
  std::sort(solDists.begin(),solDists.end());
  int noSol=0;
  for(int i=0;i<solpts.size();i++){
    noSol = noSol + solpts[i].size();
  }
  // now calculate the molecule/solvent distances
  
  double maxDist = solDists[solDists.size()-1];
  
  std::vector<double> solMolDists = mol.solMolDists(solpts);
  std::sort(solMolDists.begin(),solMolDists.end());
  // set the phases of the scattering model
  
  double  kmin = std::atof(argv[9]);
  double  kmax = std::atof(argv[10]);
  
  int noDistBins = int(1.1*std::ceil((kmax-kmin)*maxDist/3.14159265359));
  ed.setPhases(maxDist,kmin,kmax);
  double scatterFit = ed.fitToScattering(molDists,solDists,solMolDists,molSize,noSol);
  //check for overlaps
  double distSumCurr=0.0;
  for(int l=0;l<overlapDists.size();l++){
    distSumCurr = distSumCurr + std::exp(std::abs(closestApproachDist-overlapDists[l]))-1.0;
  }
  if(overlapDists.size()>0){
    distSumCurr =0.1*(1.0/overlapDists.size())*distSumCurr;
  }
  double distSum=0.0;
  scatterFit =scatterFit + distSumCurr;
  scatterFit =scatterFit + mol.getLennardJonesContact();
  /******************************************************************
  
                      Fit to scattering data
  
  ******************************************************************/ 
  int k=0;
  int noSections = mol.noChains();
  int noScatterFitSteps=1;
  std::cout<<"No. chains "<<noSections<<"\n";
  while(k<noScatterFitSteps && scatterFit>0.0001){
    // start random genrator
    std::random_device rdev{};
    std::default_random_engine generator1{rdev()};
    std::uniform_real_distribution<double> distributionR(0.0,1.0);
    // select max number of steps
    // check for overlap
    // now
    int netIndex=0;
    for(int i=1;i<=noSections;i++){
      if(i>1){
	netIndex=netIndex+mol.getSubsecSize(i-1);
      }
      for(int j=0;j<mol.getSubsecSize(i)-1;j++){
	int totalIndex = netIndex+j;
	if((doAll==true) || (std::find(fixedSecList.begin(),fixedSecList.end(),totalIndex)!=fixedSecList.end())){
	  ktlMolecule molCopy = mol;
	  molCopy.changeMoleculeSingleMulti(totalIndex,i);
	  bool cacaDist= molCopy.checkCalphas(i);
	  if(cacaDist==false){
	    // get molecule distances
	    std::vector<double> molDistsCopy;
	    std::vector<double> ovelaps = molCopy.checkOverlapWithRad(closestApproachDist);   
	    molDistsCopy = molCopy.getDistSet();
	    std::sort(molDistsCopy.begin(),molDistsCopy.end()); 
	    distSum=0.0;
	    for(int l=0;l<ovelaps.size();l++){
	      distSum = distSum + std::exp(abs(closestApproachDist-ovelaps[l]))-1.0;
	    }
	    if(ovelaps.size()>0){
	      distSum = 0.1*distSum*(1.0/ovelaps.size());
	    }       
	    mol = molCopy;
	    molDists = molDistsCopy;
	    distSumCurr = distSum;
	  }
	}
      }
      k++;
    }
  }

  mol.writeMoleculeToFile(argv[11]);
}
      
