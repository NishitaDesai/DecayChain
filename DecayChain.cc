#include "Pythia8/Pythia.h"

using namespace Pythia8;

class State{

public:
  State():psize(0), bratio(1.0), next(0), previous(0) {}; 
  //~State();
  
  int particles[25];
  int psize;
  float bratio ;
  State* next;
  State* previous;
};


class DecayChain {

public:
  DecayChain(): DEBUG(false), PDG(pythia.particleData) {};

  bool init(int argc, char** argv);
  State* addParticle(int i);
  void displayParticle(State* Head);
  
  bool isJ(int i);
  bool isL(int i);
  bool isI(int i);

protected:
  Pythia pythia;
  ParticleData& PDG;
  ofstream clog;

private:
  bool DEBUG;
  
  
};

bool DecayChain::isJ(int i){
  if (abs(i) < 6) return true;
  return false;
}

bool DecayChain::isL(int i){
  if(abs(i) == 11 || abs(i) == 13) return true;
  return false;
}

bool DecayChain::isI(int i){
  int ii = abs(i);
  if (ii == 12 || ii == 14 || ii == 16) return true;
  return false;
}

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


State* DecayChain::addParticle(int iIn){

  // TODO: populate initial list with final state particles
  State* Head = new State();
  Head->particles[0] = iIn;
  Head->psize = 1;
  int stateSize = 1;

  State* current = Head;
  bool checkDecay = true;

  if (DEBUG) cout << "Starting particle " << iIn << endl;
  while(checkDecay) {

    int temp[25];
    int tsize = 0;
    int stateSizeOld = stateSize;
    bool headChanged = false;
    
    for (int is = 0; is < current->psize; is++){
      ParticleDataEntry* part = PDG.particleDataEntryPtr(current->particles[is]);
      int pid = abs(part->id());
      if(!part->mayDecay() ||  pid < 6 || (pid > 10 && pid < 23)){
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

	  if(part->channel(im).bRatio() * top->bratio < 1.0e-4) continue;
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

  return Head;
}

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

  
}


int main(int argc, char** argv) {

  DecayChain dchain;
  if(!dchain.init(argc, argv)){
    cout << "Could not init " << endl;
    return -1;
  }

  State* part1 = dchain.addParticle(6);
  dchain.displayParticle(part1);
  
}
