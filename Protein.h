#ifndef PROTEIN_H
#define PROTEIN_H

#include "Acid.h"
#include <list>

using namespace std;

class Protein {
public:
	Protein() {}
	Protein(int id) {
	m_id = id;
	}
	~Protein() {
	}
    list<Acid>* getAcids() {
    	return m_acids;
	}	
	list<Acid>::iterator getIterator() {
		return m_iterator;
	}
	int getId() {
		return m_id;
	}
private:
	list<Acid>* m_acids = new list<Acid>;
	list<Acid>::iterator m_iterator;
	int m_id;
};

#endif
