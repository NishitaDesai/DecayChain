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



int main(int argc, char** argv) {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;
  ParticleData& PDG = pythia.particleData;

  ofstream clog;
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
    return -1;
  }
  
  // Initialize.
  pythia.init();

  bool DEBUG = false;
  
  // TODO: populate initial list with final state particles
  int nStates = 1;
  State* Final1 = new State();
  Final1->particles[0] = 1000006;
  Final1->psize = 1;
  int stateSize = 1;

  State* current = Final1;
  bool checkDecay = true;
  
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
		       << " Final1: " << Final1
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
	    if (top == Final1) {
	      Final1 = newstate;
	      headChanged = true;
	      if(DEBUG) clog << "     Final1 pointing to newstate" << endl;
	    }
	    else {
	      if(DEBUG) {
		clog << "     Changing next of " << top->previous << " to " << newstate <<endl;
		clog << "     top = " << top << " Final1= " << Final1 <<endl;
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
		       << " Final1: " << Final1
		       << endl;

	if(DEBUG) clog << "  ----- Channel Loop Ends ----" <<endl;

	
	// This finishes adding of new nodes for a single particle
	if(last == top) {
	  // No decay mode found
	  // Delete node completely
	  if(top->previous) (top->previous)->next = top->next;
	  if(top->next) (top->next)->previous = top->previous;
	  if(top == Final1) Final1 = top->next;
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

      current = Final1;

      if (DEBUG) {
	clog <<"Total states: "<< stateSize <<endl;
        State* cur = Final1;
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

  clog <<"Total states: "<< stateSize <<endl;
  State* cur = Final1;
  int ssize = 0;
  float brsum = 0.0;
  while(true){
    ssize++;
    clog << "  " << ssize << "  ";
    for (int i = 0; i< cur->psize;i++)
      clog << cur->particles[i] << "  " ;
    clog << cur->bratio;
    brsum += cur->bratio;
    if (DEBUG) clog << "  " << cur << "  " << cur->previous << "  " << cur->next;
    clog << endl;
    if(!cur->next) {
      if (ssize != stateSize && brsum < 0.99) clog << "  Error after " << ssize << endl;
      break;
    }
    cur = cur->next;
  }
  clog << "BRsum = " << brsum << endl;


  return 0;
}
