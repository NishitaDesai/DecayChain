#include "Pythia8/Pythia.h"

using namespace Pythia8;

//===========================================================

class State{

public:
  State():psize(0), bratio(1.0), next(0), previous(0), dbr(0.0), csize(1) {}; 
  //~State();
  
  int particles[25];
  int psize;
  float bratio;
  State* next;
  State* previous;
  float dbr;
  int csize;
};

//===========================================================

class DecayChain {

public:
  DecayChain(): DEBUG(false), PDG(pythia.particleData), tol(1.0e-4), noSMdecay(true) {};

  bool init(int argc, char** argv);
  State* addParticle(int i);
  void displayParticle(State* Head);

  float getFrac(State* Head, int nLep = -1, int nJet = -1, int nTau = -1, bool MET = false);
  float get2Frac(State* Head1, State* Head2, int nLep = -1, int nJet = -1, int nTau = -1, bool MET = false);
  float bosonDecay(int nW = 0, int nZ = 0, int nlep = -1, int njet = -1, int ntau = -1, bool MET = false);  
  
  bool isJ(int i);
  bool isL(int i);
  bool isI(int i);

  float tol;
  bool DEBUG;
  bool noSMdecay;
  
protected:
  Pythia pythia;
  ParticleData& PDG;
  ofstream clog;

  
  
};

//-----------------------------------------------------------

bool DecayChain::isJ(int i){
  if (abs(i) < 6) return true;
  return false;
}

//-----------------------------------------------------------

bool DecayChain::isL(int i){
  if(abs(i) == 11 || abs(i) == 13) return true;
  return false;
}

//-----------------------------------------------------------

bool DecayChain::isI(int i){
  int ii = abs(i);
  if (ii == 12 || ii == 14 || ii == 16) return true;
  return false;
}

//-----------------------------------------------------------

bool DecayChain::init(int argc, char** argv){

  PDG = pythia.particleData;
  clog.open("DChain.log",ios::out);
  
  if(argc > 1) {
    stringstream ss;
    ss << "SLHA:file = " << argv[1] << endl;
    pythia.readString(ss.str());
    pythia.readString("SUSY:gg2gluinogluino = on");
    pythia.readString("SLHA:minMassSM = 10.0");
    pythia.readString("SLHA:allowUserOverride = true");
    pythia.readString("ResonanceWidths:minWidth = 1.0e-30");
    pythia.readString("Init:showChangedSettings = off      ! list changed settings");
    pythia.readString("Init:showChangedParticleData = off ! list changed particle data");

  }
  else {
    clog << "Please give an SLHA file" <<endl;
    return false;
  }

  // Initialize.
  pythia.init();

  bool DEBUG = false;
  return true;
  
}

//-----------------------------------------------------------

State* DecayChain::addParticle(int iIn){

  // TODO: populate initial list with final state particles
  State* Head = new State();
  Head->particles[0] = iIn;
  Head->psize = 1;
  int stateSize = 1;

  State* current = Head;
  bool checkDecay = true;
  float discardedBR = 0.0;

  if (DEBUG) cout << "Starting particle " << iIn << endl;
  while(checkDecay) {

    int temp[25];
    int tsize = 0;
    int stateSizeOld = stateSize;
    bool headChanged = false;
    
    for (int is = 0; is < current->psize; is++){
      ParticleDataEntry* part = PDG.particleDataEntryPtr(current->particles[is]);
      int pid = abs(part->id());
      int pidmax = 23;
      if(noSMdecay) pidmax = 25;
      if(!part->mayDecay() ||  pid < 6 || (pid > 10 && pid < pidmax)){
	temp[tsize] = current->particles[is];
	tsize++;
      
	// If full list contains only stable particles
	if (tsize == current->psize) {
	  // Node complete
	  if(stateSizeOld == stateSize){
	    if(!current->next) {
	      // Case one: full chain done; exit
	      checkDecay = false;
	      break;
	    }
	    else {
	      // Finished current node only

	      if (DEBUG) {
		clog << "Stable node : ";
		for (int it = 0; it < tsize; it++)
		  clog << temp[it] << "  ";
		clog << endl;
		clog << "  Current: " << current << " next: " << current->next << endl;
	      }
	      
	      current = current->next;
	      
	      break;
	    }
	  }
	}
      }
      else {
	
	if(DEBUG) clog << "Performing particle " << part->name() << "  " << pid << endl;

	// Repace the current node with list of nodes corresponding to
	// the decay table
	
	State *oldstate = current;
	State* top = current;
	State* last = current;

	if(DEBUG) clog << "   Oldstate: " << oldstate
		       << " Top: " << top
		       << " Head: " << Head
		       << endl;

	int iChan = 0;
	// Add new nodes for a single particle
	for(int im = 0; im < part->sizeChannels(); im++){

	  if(part->channel(im).bRatio() * top->bratio < tol) {
	    discardedBR += part->channel(im).bRatio() * top->bratio;
	    continue;
	  }
	  iChan++;
	  
	  State* newstate = new State();
	  if (DEBUG) clog << "    Generated newstate: " << newstate << endl;
	  
	  
	  int ip = 0;
	  for (; ip < tsize; ip++)
	    newstate->particles[ip] = temp[ip];
	  
	  DecayChannel& ch = part->channel(im); 
	  if (DEBUG)  clog << "    Adding bratio: " << im << "  " << ch.bRatio() << "  " ;
	  for (int j = 0; j < ch.multiplicity(); j++){
	    newstate->particles[ip] = ch.product(j);
	    ip++;
	    if (DEBUG) clog << ch.product(j) << "  ";
	  }
	  if(DEBUG) clog << endl;
	  
	  // Copy rest of the particle list
	  for (int j = tsize + 1; j < top->psize; j++){
	    newstate->particles[ip] = top->particles[j];
	    ip++;
	    
	  }
	  newstate->psize = ip;
	  newstate->bratio = top->bratio * ch.bRatio();

	  // Reset pointers
	  if (iChan == 1) {
	    newstate->previous = top->previous;

	    // Change head for first particle
	    if (top == Head) {
	      Head = newstate;
	      headChanged = true;
	      if(DEBUG) clog << "     Head pointing to newstate" << endl;
	    }
	    else {
	      if(DEBUG) {
		clog << "     Changing next of " << top->previous << " to " << newstate <<endl;
		clog << "     top = " << top << " Head= " << Head <<endl;
	      }
	      (top->previous)->next = newstate;
	    }
	  }
	  else {
	    newstate->previous = oldstate;
	    oldstate->next = newstate;
	    stateSize++;
	  }

	  // Finish up before loop ends
	  oldstate = newstate;
	  last = newstate;
	}
	if(DEBUG) clog << "  Last: " << last
		       << " Top: " << top
		       << " Head: " << Head
		       << endl;

	if(DEBUG) clog << "  ----- Channel Loop Ends ----" <<endl;

	
	// This finishes adding of new nodes for a single particle
	if(last == top) {
	  // No decay mode found
	  // Delete node completely
	  if(top->previous) (top->previous)->next = top->next;
	  if(top->next) (top->next)->previous = top->previous;
	  if(top == Head) Head = top->next;
	  stateSize--;
	}
	else {
	    last->next = top->next;
	    if(top->next) (top->next)->previous = last;
	}

	if(top) delete(top);
	headChanged = true;
	break;
      }
    }


    if(stateSizeOld !=stateSize || headChanged) {

      current = Head;

      if (DEBUG) {
	clog <<"Total states: "<< stateSize <<endl;
        State* cur = Head;
	int ssize = 0;
	while(true){
	  ssize++;
	  clog << "  " << ssize << "  ";
	  for (int i = 0; i< cur->psize;i++)
	    clog << cur->particles[i] << "  " ;
	  clog << cur->bratio;
	  if (DEBUG) clog << "  " << cur << "  " << cur->previous << "  " << cur->next;
	  clog << endl;
	  if(!cur->next) {
	    clog << "  Error after " << ssize << endl;
	    break;
	  }
	  cur = cur->next;
	}
      }
    }

  }

  Head->dbr = discardedBR;
  Head->csize = stateSize;
  return Head;
}

//-----------------------------------------------------------

void DecayChain::displayParticle(State* Head){
  
  State* cur = Head;
  int ssize = 0;
  float brsum = 0.0;
  while(true){
    ssize++;
    cout << "  " << ssize << "  ";
    for (int i = 0; i< cur->psize;i++)
      cout << cur->particles[i] << "  " ;
    cout << cur->bratio;
    brsum += cur->bratio;
    if (DEBUG) cout << "  " << cur << "  " << cur->previous << "  " << cur->next;
    cout << endl;
    if(!cur->next) break;
    cur = cur->next;
  }
  cout << "BRsum = " << brsum << endl;
  cout << "Discarded BR = " << Head->dbr <<endl;
  
}

//-----------------------------------------------------------
float DecayChain::bosonDecay(int nW, int nZ, int nlep, int njet, int ntau, bool MET){

  cout << " nW+nZ=" << nW+nZ << "njet = " << njet <<endl;
  // Add W and Z decay tables, make combinations
  float tempbr = 0.0;
  // lep, jets, tau, MET
  float Wbr[4] = {0.213, 0.6, 0.107, 0.321};
  float Zbr[4] = {0.04, 0.6, 0.34, 0.02};

  // You can only get jets in multiples of two
  if(njet%2 != 0) return 0.0;
  int nhad = njet/2;
  if (nhad > nW + nZ) return 0.0;
  
  if (nhad == nW && nZ == 0) return Wbr[1];
  if (nhad == nZ && nW == 0) return Zbr[1];

  int bosarr[20], bSize;
  for (int i = 0; i < nW ;i++){
    bosarr[bSize] = Wbr[1];
    bSize++;
  }
  for (int i = 0; i < nZ ;i++){
    bosarr[bSize] = Zbr[1];
    bSize++;
  }
  
  if (nlep < 0) {

    // If exactly as many jets as from bosons
    if (bSize == nhad) {
      tempbr = 1.0;
      for (int ib = 0; ib < bSize; ib++)
	tempbr *= bosarr[ib];
      return tempbr;
    }

    // Else pick combination
    // Write code for getting all combinations (nbos, nhad)
    // Implementation of Knuth's combinatoric Algorithm L: n = bSize, t = nhad   

    int c[25] = {0};
    
    // Initialise
    for(int i = 1; i < nhad+1; i++){
      c[i] = i-1;
    }
    c[nhad+1] = bSize;
    c[nhad+2] = 0;

    // Visit
    int iComb;
    while (true){
      float tempbr2 = 1.0;
      for (int j = nhad; j > 0; j--)
	tempbr2 *= bosarr[c[j]] ;

      tempbr += tempbr2;

      iComb = 1;
      
      while (c[iComb] + 1 == c[iComb+1]) {
	c[iComb] = iComb-1;
	iComb++;
      }
  
      if(iComb > nhad) break;
      c[iComb] = c[iComb] + 1;
    }

  }
  else if (njet < 0) {
    // Pure leptonic states
    
  }
  else {
    // Both leptons and jets specified
    if (nlep == 0){
      // TODO: Z can go into MET
      
      if(bSize > nhad && nZ == 0) return 0.0;
      float tempbr2 = 0.0;
      for (int i = 0; i < bSize; i++){
	tempbr2 *= bosarr[i];
      }
      return tempbr2;
    }
  }

  return 0.0;
}

//-----------------------------------------------------------
float DecayChain::getFrac(State* Head, int nLep, int nJet, int nTau, bool MET){

  if(!Head){
    clog << "DecayChain::getFrac : No head node found. Exiting. " << endl;
    return 0.0;
  }
  
  if (nLep < 0  && nJet < 0 && nTau < 0 && !MET) {
    clog << "DecayChain::getFrac : No final states given. Exiting. " << endl;
    return 0.0;
  }  
  
  State* cur = Head;
  float brsum = 0.0;
  
  while(true){
    int njet = 0, nlep = 0, ninv = 0, ntau = 0;
    int nW = 0, nZ = 0, nphot = 0;
    for (int i = 0; i< cur->psize;i++) {
      int pid = abs(cur->particles[i]);
      if(pid < 6 || pid == 21) njet++;
      else if(pid == 11 || pid == 13) nlep++;
      else if(pid == 15) ntau++;
      else if(pid == 12 || pid == 14 || pid == 16) ninv++;
      else if(pid > 1000000) ninv++;
      else if(pid == 22) nphot++;
      else if(pid == 24) nW++;
      else if(pid == 25) nZ++;
      else cout << " Not asigned:" << pid << endl;
      }
    bool accept = true;
    float tempbr = 0.0;

    // If extra jets or leptons than required
    if((nlep > nLep && nLep > -1) ||
       (njet > nJet && nJet > -1) ||
       (ntau > nTau && nTau > -1) )  accept = false;

    if (accept) {
      if (nW + nZ == 0) {
	if (MET && ninv == 0) accept = false;
	if (nLep > -1 && nLep != nlep) accept = false;
	if (nJet > -1 && nJet != njet) accept = false;
	if (nTau > -1 && nTau != ntau) accept = false;
      }
      else  {
	// Get BR for getting excess jets/leptons/met from gauge bosons
	tempbr = bosonDecay(nW, nZ, nLep-nlep, nJet-njet, (MET && ninv == 0));
	if (tempbr > 0) accept = true;
	cout << "Looking for " << nJet-njet << "jets from " << nW + nZ << "bosons, br = " << tempbr << endl; 
      }
    }

    if (accept){
      if(tempbr > tol) brsum *= cur->bratio * tempbr;
      else brsum += cur->bratio;
    }
    
    if(!cur->next) break;
    cur = cur->next;
  }

  clog << "getFrac::Calculating head = "<< Head << " nLep = " << nLep << ", nJet = " << nJet <<" brsum = " << brsum <<endl; 
  
  return brsum;
}

//-----------------------------------------------------------

float DecayChain::get2Frac(State* Head1, State* Head2, int nLep, int nJet, int nTau, bool MET){

  if(!Head1 || !Head2){
    clog << "DecayChain::get2Frac : No head node found. Exiting. " << endl;
    return 0.0;
  }
  if (nLep < 0 && nJet < 0 && !MET) {
    clog << "DecayChain::get2Frac : No final states given. Exiting. " << endl;
    return 0.0;
  }

  float brfrac = 0;

  if (nLep > -1 && nJet < 1){
    // Only Leptons
    for (int il = 0; il <= nLep; il++)
      brfrac += getFrac(Head1, il, nJet, nTau, MET) * getFrac(Head2, nLep-il, nJet, nTau, MET); 
  }
  else if (nLep < 1 && nJet > -1) {
    // Only jets exclusive
    for (int ij = 0; ij <= nJet; ij++){
      if (DEBUG) {

      }
      clog << "get2Frac:: Calculating combination: "<< ij << ", " << nJet - ij << endl;
      double frac1 = getFrac(Head1, nLep, ij, -1, MET);
      double frac2 = getFrac(Head2, nLep, nJet - ij, -1, MET);
      clog << "get2Frac:: frac1 = " <<  frac1 << "  frac2 = " << frac2 << endl;
	
      brfrac += getFrac(Head1, nLep, ij, nTau, MET) * getFrac(Head2, nLep, nJet - ij, nTau, MET); 
    }
  }
  else if (nLep + nJet == -2){
    // Only MET
    brfrac += getFrac(Head1, nLep, nJet, nTau, MET) * getFrac(Head2, nLep, nJet, nTau, MET);     
  }
  else {
    for (int il = 0; il <= nLep; il++){
      for (int ij = 0; ij <= nJet; ij++){
	brfrac += getFrac(Head1, il, ij, nTau, MET) * getFrac(Head2, nLep-il, nJet-ij, nTau, MET); 
      }
    }
  }

  return brfrac;
  
}

//===========================================================

int main(int argc, char** argv) {

  DecayChain dchain;
  if(!dchain.init(argc, argv)){
    cout << "Could not init " << endl;
    return -1;
  }

  dchain.tol = 1.0e-5;
  dchain.noSMdecay = false;

  State* part1 = dchain.addParticle(1000021);
  State* part2 = part1; // dchain.addParticle(1000021);
  
  //dchain.displayParticle(part1);

  cout.precision(6);
  cout << "Discarded BR = " << part1->dbr << endl;
  cout << "tol = " << dchain.tol << " ssize = " << part1->csize << endl;
  cout << endl;
  
  // cout << "0l + 2j + MET: " << dchain.getFrac(part1, 0, 2, -1, true) << endl;
  // cout << "0l + 3j + MET: " << dchain.getFrac(part1, 0, 3, -1, true) << endl;
  // cout << "0l + 4j + MET: " << dchain.getFrac(part1, 0, 4, -1, true) << endl;
  // cout << "1l + 2j + MET: " << dchain.getFrac(part1, 1, 2, -1, true) << endl;
  // cout << "1l + 3j + MET: " << dchain.getFrac(part1, 1, 3, -1, true) << endl;
  // cout << "1l + 4j + MET: " << dchain.getFrac(part1, 1, 4, -1, true) << endl;

  // for (int il = 0; il < 5; il++) {
  //   cout << il << "l: " << dchain.getFrac(part1, il, -1, -1, false) << endl;
  // }
  // cout << endl;

  // for (int ij = 0; ij < 10; ij++) {
  //   cout << ij << "j: " << dchain.getFrac(part1, -1, ij, -1, false) << endl;
  // }
  // cout << endl;

  // for (int il = 0; il < 5; il++) {
  //   cout << il << "tau: " << dchain.getFrac(part1, -1, -1, il, false) << endl;
  // }
  // cout << endl;

  // for (int ij = 0; ij < 9; ij++) {
  //   cout << ij << "j (2 part): " << dchain.get2Frac(part1, part2, -1, ij, -1, true) << endl;
  // }
  // cout << endl;


  cout << "2j: " << dchain.getFrac(part1,-1,2,-1,false) << endl;
  cout << "2j + MET: " << dchain.getFrac(part1,-1,2,-1,true) << endl;
  cout << "4j: " << dchain.get2Frac(part1,part2,-1, 4,-1,false) << endl;
  cout << "4j + MET: " << dchain.get2Frac(part1,part2,-1, 4,-1,true) << endl;  
}

/******** TODO ************:

1. Add method to delete particle

**************************/
