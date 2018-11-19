/// CS361 amino acids project 1
/// collaborated with Joseph Rugh, Rapheal Sandor, and Derek Tiller
#include <iostream>
#include "Util.h"
#include "Protein.h"
#include <random>
#include <ctime>
#include <tuple>
#include <fstream>

typedef unsigned int uint; // Prevent warnings for, for loops who's comparison field is comparing against .size() from list

struct Locations
{
    list<tuple<double, double, double>> locations1;
    list<tuple<double, double, double>> locations2;
    list<tuple<double, double, double>> locations3;
};

int requestAcidAmount(Protein protein, default_random_engine * randomEngine, int lastId);

void display(Protein protein);

void insertSolution ( Protein p1, Protein p2,Protein p3, default_random_engine * randomEngine);

void linkSolutionNodes(list<Acid> & acids, default_random_engine * randomEngine, list<Acid>::iterator iterator);

void attachSolution (list<Acid> & aList, Protein &p, int id);

int createNodes(list<Acid> & acids, default_random_engine * randomEngine, int numAcids, int lastId);

void linkNodes(list<Acid> & acids, default_random_engine * randomEngine, list<Acid>::iterator iterator);

void translateNodes(list<Acid> * acids, list<Acid>::iterator iterator, int bottomVertexId);

void rotateNodes(list<Acid> * acids, list<Acid>::iterator iterator, tuple<int, int> aicdIds);

void rotateNodesByMat(list<Acid> * acids, list<Acid>::iterator iterator, double mat[3][3]);

void gnuPlotOutput(ofstream& p1plot, ofstream& p2plot, ofstream& p3plot, ofstream& comout, Locations * locations);

/*
 * Tuple containing 3 tuple pairs.
 * Each pair contains the id of the acid that sits at the bottom of the "V" shape for the key,
 * and the acid that makes up the free vector.
 *
 * From the first acid id the Acid that makes up the standing vector may be found.
 */
tuple<tuple<int, int>, tuple<int, int>, tuple<int, int>> findMatchingKey(Protein protein1, Protein protein2,
                                                                         Protein protein3, Protein newProteinOrder[3]);

void compareKeys(list<Acid>::iterator iterator, list<Acid> * acids,
                 double length1a, double length1b, double angle1,
                 int & ida, int & idb);

double getAngle(list<Acid>::iterator & acid, Acid* a1, Acid* a2);

Locations* findLocations(Protein protein1, Protein protein2, Protein protein3);

/*--------Constants----------*/

const double CHANCE_FOR_NONE = 0.2;  // 20% chance that a node will have no connections. Caps at 50% then goes in the reverse.
const int MAX_WINDOW         = 10;   // Window width for adding nodes is 10

int main()
{

    ofstream p1plot;
    ofstream p2plot;
    ofstream p3plot;
    ofstream comout;

	default_random_engine * randomEngine = new default_random_engine(time(NULL));
	Protein protein1(0), protein2(1), protein3(2);

	int lastId = 0;
	lastId = requestAcidAmount(protein1, randomEngine, lastId);
	lastId = requestAcidAmount(protein2, randomEngine, lastId);
	requestAcidAmount(protein3, randomEngine, lastId);

	insertSolution(protein1, protein2, protein3, randomEngine);

	Protein newProteinOrder[3];
	auto proteinIds = findMatchingKey(protein1, protein2, protein3, newProteinOrder);

	// If idb1 != -1 -> proceed to translate and rotation nodes
	int idb1 = get<1>(get<0>(proteinIds));
    if (idb1 != -1)
    {
        translateNodes(newProteinOrder[0].getAcids(), newProteinOrder[0].getIterator(), get<0>(get<0>(proteinIds)));
        translateNodes(newProteinOrder[1].getAcids(), newProteinOrder[1].getIterator(), get<0>(get<1>(proteinIds)));
        translateNodes(newProteinOrder[2].getAcids(), newProteinOrder[2].getIterator(), get<0>(get<2>(proteinIds)));
        rotateNodes(newProteinOrder[0].getAcids(), newProteinOrder[0].getIterator(), get<0>(proteinIds));
        rotateNodes(newProteinOrder[1].getAcids(), newProteinOrder[1].getIterator(), get<1>(proteinIds));
        rotateNodes(newProteinOrder[2].getAcids(), newProteinOrder[2].getIterator(), get<2>(proteinIds));

        auto locations = findLocations(newProteinOrder[0], newProteinOrder[1], newProteinOrder[2]);

        cout << "Matching Coordinates locations:" << endl;
        list<tuple<double, double, double>> locations1 = locations->locations1;
        list<tuple<double, double, double>>::iterator itr = locations1.begin();
        while (itr != locations1.end())
        {
            cout << "(" << get<0>(*itr) << ", " << get<1>(*itr) << ", " << get<2>(*itr) << ")" << endl;
            itr++;
        }




        gnuPlotOutput(p1plot, p2plot, p3plot, comout, locations);
        delete locations;
    }

	return 0;
}

int requestAcidAmount(Protein protein, default_random_engine * randomEngine, int lastId)
{
	int numAcids;
	cout << protein.getId() << ") How many acids:";
	cin >> numAcids;

	lastId = createNodes(*protein.getAcids(), randomEngine, numAcids, lastId);
	linkNodes(*protein.getAcids(), randomEngine, protein.getIterator());
	return lastId;
}

void display(Protein protein)
{
	cout << "---------------------------------\n";
	cout << "Protein -> " << protein.getId() << "\n";
	cout << "---------------------------------\n";

	list<Acid>::iterator iterator = protein.getIterator();
	list<Acid> acids = *protein.getAcids();
	iterator = acids.begin();
	while (iterator != acids.end())
	{
		iterator->display();
		iterator++;
	}
}

int createNodes(list<Acid> & acids, default_random_engine * randomEngine, int numAcids, int lastId)
{
	normal_distribution<double> locDistribution(5.0, 3.0);

	double ax = 0, ay = 0, az = 0;

	int currentId = lastId;
	for (int i = 0; i < numAcids; i++)
	{
		// Random x, y, and z coordinate generation
		double x = ax += locDistribution(*randomEngine);
		double y = ay += locDistribution(*randomEngine);
		double z = az += locDistribution(*randomEngine);

		Acid* acidPtr = new Acid(x, y, z);
		acidPtr->setId(currentId++);
		acids.push_back(*acidPtr);
	}

	return currentId;
}

void linkNodes(list<Acid> & acids, default_random_engine * randomEngine, list<Acid>::iterator iterator)
{
	int numAcids = acids.size();

	// Distribution used to provide a connection or no connection
	uniform_real_distribution<double> noNodeDistribution(0, 1);
	// Distribution used to determine how many connections there will be
	normal_distribution<double> connectionsDistributions(3.0, 1.0);

	iterator = acids.begin();
	for (int i = 0; i < numAcids - 1; i++)
	{
		if (!(noNodeDistribution(*randomEngine) > CHANCE_FOR_NONE))
		{
			iterator++;
			continue;
		}

		int numberOfConnections = (int)connectionsDistributions(*randomEngine);
		if (numberOfConnections < 1) numberOfConnections = 1;
		int currentConnections [numberOfConnections];
		for (int j = 0; j < numberOfConnections; j++) currentConnections[j] = -1; // Initialize to -1's
		for (int j = 0; j < numberOfConnections; j++)
		{
			int windowOffset = i + 1 + MAX_WINDOW;

			// Cut off the window view when we get close to the end of the chain
			if (windowOffset > numAcids) windowOffset = numAcids;

			// Distribution for which connection to connect to
			uniform_int_distribution<int> nodeDistribution(i + 1, windowOffset);
			int connectionId = nodeDistribution(*randomEngine);

            if(connectionId > numAcids - 1) connectionId = numAcids - 1;

			// Prevent duplicate connections
			bool preventDuplicate = false;
			for (int k = 0; k < numberOfConnections; k++)
			{
				// We already have that connection
				if (currentConnections [k] == connectionId  ) preventDuplicate = true;
			}

			if (!preventDuplicate)
			{
				currentConnections [j] = connectionId;
				Acid* acidToAdd = &getListElement(acids, connectionId);

				double dx = abs(iterator->getX() - acidToAdd->getX());
				double dy = abs(iterator->getY() - acidToAdd->getY());
				double dz = abs(iterator->getZ() - acidToAdd->getZ());

				iterator->addConnection(acidToAdd);
				iterator->addLength(sqrt(dx*dx + dy*dy + dz*dz));
			}
		}
		iterator++;
	}
}

void linkSolutionNodes(list<Acid> & acids, default_random_engine * randomEngine, list<Acid>::iterator iterator)
{
	int numAcids = acids.size();

	// Distribution used to provide a connection or no connection
	uniform_real_distribution<double> noNodeDistribution(0, 1);
	// Distribution used to determine how many connections there will be
	normal_distribution<double> connectionsDistributions(3.0, 1.0);

	iterator = acids.begin();
	for (int i = 0; i < numAcids - 1; i++)
	{
		if (!(noNodeDistribution(*randomEngine) > .10))
		{
			iterator++;
			continue;
		}

		int numberOfConnections = (int)connectionsDistributions(*randomEngine);
		if (numberOfConnections < 2 || numberOfConnections >= numAcids) numberOfConnections = 2;
		int currentConnections [numberOfConnections];
		for (int j = 0; j < numberOfConnections; j++) currentConnections[j] = -1; // Initialize to -1's
		for (int j = 0; j < numberOfConnections; j++)
		{
			int windowOffset = i + 2 + MAX_WINDOW;

			// Cut off the window view when we get close to the end of the chain
			if (windowOffset > numAcids) windowOffset = numAcids;

			// Distribution for which connection to connect to
			uniform_int_distribution<int> nodeDistribution(i+1 , windowOffset);
			int connectionId = nodeDistribution(*randomEngine);

            if(connectionId > numAcids - 1) connectionId = numAcids -1;

			// Prevent duplicate connections
			bool preventDuplicate = false;
			for (int k = 0; k < numberOfConnections; k++)
			{
				// We already have that connection
				if (currentConnections [k] == connectionId  ) preventDuplicate = true;
			}

			if (!preventDuplicate)
			{
				currentConnections [j] = connectionId;
				Acid* acidToAdd = &getListElement(acids, connectionId);

				double dx = abs(iterator->getX() - acidToAdd->getX());
				double dy = abs(iterator->getY() - acidToAdd->getY());
				double dz = abs(iterator->getZ() - acidToAdd->getZ());

				iterator->addConnection(acidToAdd);
				iterator->addLength(sqrt(dx*dx + dy*dy + dz*dz));
			}
		}
		iterator++;
	}
}

void translateNodes(list<Acid> * acids, list<Acid>::iterator iterator, int bottomVertexId)
{
	iterator = acids->begin();
	while (iterator != acids->end())
    {
        if (iterator->getId() == bottomVertexId) break;
        iterator++;
    }

	/*- Find the offset from the first acid -*/
	double dx = -iterator->getX();
	double dy = -iterator->getY();
	double dz = -iterator->getZ();

	iterator = acids->begin();
	while (iterator != acids->end())
	{
		iterator->ajustLocation(dx, dy, dz);
		iterator++;
	}
}

void rotateNodes(list<Acid> * acids, list<Acid>::iterator iterator, tuple<int, int> aicdIds)
{
	if (acids->size() < 3) return;

    auto getAcidFromPtrList = [](list<Acid>* acids, int id) -> Acid& {
        list<Acid>::iterator acidItr = acids->begin();
        while (acidItr != acids->end())
        {
            if (acidItr->getId() == id) return *acidItr;
            acidItr++;
        }
    };

    Acid bottomVertexAcid = getAcidFromPtrList(acids, get<0>(aicdIds));// Id of the acid at the bottom of the 'V' was used

    list<Acid*> tempAcids = bottomVertexAcid.getConnections();
    Acid* topOfStanding = getListElement(tempAcids, 0);              // Receiving the endpoint of the standing vector
    Acid* topOfFree = &getAcidFromPtrList(acids, get<1>(aicdIds));// Id of the free vector's acid endpoint was used

    //---------------------------------------------------------------------

	double x1 = topOfStanding->getX();
	double y1 = topOfStanding->getY();
	double z1 = topOfStanding->getZ();

	double angle1a = atan(y1/z1);

	double vec1[3] = { x1, y1, z1 };
	double mat1[3][3];
	double mat2[3][3];
	double mat3[3][3];

	zeroMat(mat1);

	// Canceling of the Y axis for vector 1
	mat1[0][0] = 1.0;
	mat1[1][1] = cos(angle1a);
	mat1[1][2] = -1.0*sin(angle1a);
	mat1[2][1] = sin(angle1a);
	mat1[2][2] = cos(angle1a);

	// Need to apply this rotation to all the acids
	double r1[3] = { 0.0, 0.0, 0.0 };
	multiplyMat2Vec(vec1, mat1, r1);

	// Storing to receive the second angle
	vec1[0] = r1[0];
	vec1[1] = r1[1];
	vec1[2] = r1[2];

	// Iterate over all the acids and apply the first rotation
    rotateNodesByMat(acids, iterator, mat1);

	zeroMat(mat2);
	double angle2 = atan(-vec1[2]/vec1[0]);

	// Canceling of the Z axis for vector 1
	mat2[0][0] = cos(angle2);
	mat2[0][2] = -1.0*sin(angle2);
	mat2[1][1] = 1.0;
	mat2[2][0] = sin(angle2);
	mat2[2][2] = cos(angle2);

	double r2[3];
	multiplyMat2Vec(vec1, mat2, r2);

	vec1[0] = r2[0];
	vec1[1] = r2[1];
	vec1[2] = r2[2];

	// Iterate over all the acids and apply the second rotation
	rotateNodesByMat(acids, iterator, mat2);

	double y2 = topOfFree->getY();
	double z2 = topOfFree->getZ();
	double angle1b = atan(y2/z2);

    // Canceling of the Y axis for vector 2
    zeroMat(mat3);

	mat3[0][0] = 1.0;
	mat3[1][1] = cos(angle1b);
	mat3[1][2] = -1.0*sin(angle1b);
	mat3[2][1] = sin(angle1b);
	mat3[2][2] = cos(angle1b);

	// Iterate over all the acids and apply the third rotation
	rotateNodesByMat(acids, iterator, mat3);
}

void rotateNodesByMat(list<Acid> * acids, list<Acid>::iterator iterator, double mat[3][3])
{
    iterator = acids->begin();
	while (iterator != acids->end())
    {
        double vect[3] = { iterator->getX(), iterator->getY(), iterator->getZ() };

        double result[3] = { 0.0, 0.0, 0.0 };
        multiplyMat2Vec(vect, mat, result);

        iterator->setLocation(result[0], result[1], result[2]);
        iterator++;
    }
}

tuple<tuple<int, int>, tuple<int, int>, tuple<int, int>> findMatchingKey(Protein protein1, Protein protein2,
                                                                         Protein protein3, Protein newProteinOrder[3])
{

	int ida1 = -1, ida2 = -1, ida3 = -1; // Standing vector ids
	int idb1 = -1, idb2 = -1, idb3 = -1; // Free vector ids

	// Rearranging the proteins so that the first protein in the array is the shortest

	auto size = [](Protein p) -> int { return p.getAcids()->size(); };
	int sid = size(protein1) > size(protein2) ?
								 size(protein2) > size(protein3) ? 2 : 1 :
								 size(protein1) > size(protein3) ? 2 : 0;

	Protein oldProteinOrder[3] = { protein1, protein2, protein3 };

	// Circulate and fill in remaining proteins
    newProteinOrder[0] = oldProteinOrder[sid];
    if (++sid > 2) sid = 0;
	newProteinOrder[1] = oldProteinOrder[sid];
    if (++sid > 2) sid = 0;
    newProteinOrder[2] = oldProteinOrder[sid];

	list<Acid>::iterator iterator1 = newProteinOrder[0].getIterator();
    list<Acid>::iterator iterator2 = newProteinOrder[1].getIterator();
    list<Acid>::iterator iterator3 = newProteinOrder[2].getIterator();

    list<Acid>* acids1 = newProteinOrder[0].getAcids();
    list<Acid>* acids2 = newProteinOrder[1].getAcids();
    list<Acid>* acids3 = newProteinOrder[2].getAcids();

    //----------------------------------------------------------------------------------

    iterator1 = acids1->begin();
    while (iterator1 != acids1->end())
    {
        // If less than 2 connections return
    	list<Acid*> connections1 = iterator1->getConnections();
    	if (connections1.size() < 2)
		{
			iterator1++;
			continue; // Only consider cases where there are 2 connections to make a key
    	}


    	// GOOD TILL THIS POINT


    	list<double> lengths = iterator1->getLengths();

    	// Standing vector data
    	// The vector that stands in place. This is the vector's length and the acid endpoint called a1.
    	double length1a = getListElement(lengths, 0);
    	Acid* a1 = getListElement(connections1, 0);

    	// For all lengths go through all other lengths of the free vector
    	for (uint i = 1; i < lengths.size(); i++)
    	{
    		// Length to compare
    		double length2a = getListElement(lengths, i); // This is the length of the free vector

    		// Angle to compare
    		Acid* a2 = getListElement(connections1, i); // This acid at the endpoint of the free vector called a2.
    		double angle1 = getAngle(iterator1, a1, a2); // Using the acid at the bottom of the 'V' shape called iterator1
                                                         // This is then used in conjunction with a1 and a2 (the end points)
                                                         // to find the angle between the vectors which is the angle for the
                                                         // first key

    		// If no matching keys are found aka id2 or id3 == 0 then continue to search
    		compareKeys(iterator2, acids2, length1a, length2a, angle1, ida2, idb2);
    		compareKeys(iterator3, acids3, length1a, length2a, angle1, ida3, idb3);

    		// Key match in all 3!
    		if (ida2 != -1 && ida3 != -1) {
    			ida1 = iterator1->getId();
    			idb1 = a2->getId(); // Acid ID for the free vector
    			return make_tuple(make_tuple(ida1, idb1), make_tuple(ida2, idb2), make_tuple(ida3, idb3));
			}
			ida2 = ida3 = -1; // Resetting since the keys are not correct
		}

    	iterator1++;
	}
	return make_tuple(make_tuple(ida1, idb1), make_tuple(ida2, idb2), make_tuple(ida3, idb3));
}

void compareKeys(list<Acid>::iterator iterator, list<Acid> * acids,
                 double length1a, double length2a, double angle1,
                 int & ida, int & idb)
{
    // We are trying to search for a second key by iterating through all the acids in the 2nd and 3rd protein
	iterator = acids->begin();
	while (iterator != acids->end())
	{
	    // If the connections in the 2nd and 3rd protein are less than two there is no key for that given acid
		list<Acid*> connections = iterator->getConnections();
		if (connections.size() < 2)
		{
			iterator++;
			continue; // Only consider cases where there are 2 connections to make a key
		}

		// The lengths of the 2nd and 3rd protein
		list<double> lengths = iterator->getLengths();

		// Standing vector data
		double length1b = getListElement(lengths, 0); // The standing vector
		Acid* b1 = getListElement(connections, 0);    // The standing vector's end point acid

		// Going through all the free vector lengths for protein 2 and 3
		for (uint i = 1; i < lengths.size(); i++)
		{
		    // grabbing the length of the free vector of the 2nd and 3rd protein
			double length2b = getListElement(lengths, i);

			// Check for matching lengths
			if (length1a == length1b && length2a == length2b)
			{
				// Check angles
				Acid* b2 = getListElement(connections, i); // The free vector's endpoint acid for 2nd and 3rd proteins
				double angle2 = getAngle(iterator, b1, b2); // angle for the 2nd and 3rd key

				// Matching keys!
				if (angle1 == angle2) {
					ida = iterator->getId(); // Acid ID at the bottom of the key
					idb = b2->getId(); // Acid ID for the free standing vector
					return;
				}
			}
		}

		iterator++;
	}
}

double getAngle(list<Acid>::iterator& acid, Acid* a1, Acid* a2)
{
	double startp[3] = { acid->getX(), acid->getY(), acid->getZ() };
	double endp1[3] = { a1->getX(), a1->getY(), a1->getZ() };
	double endp2[3] = { a2->getX(), a2->getY(), a2->getZ() };
	return angleBetweenVectors(startp, endp1, startp, endp2);
}

void attachSolution (list<Acid> & sPts, Protein & p, int id)
{
    list<Acid>::iterator proteinItr;
    list<Acid>::iterator solutionItr;

    list<Acid>* tempAcids = p.getAcids();
    proteinItr = tempAcids->begin();
    solutionItr = sPts.begin();
    while ((id != proteinItr->getId()) && proteinItr != tempAcids->end())
    {
        proteinItr++;

    }
    while (solutionItr != sPts.end())
    {
      // some local temp variables
      tempAcids->insert(proteinItr,* solutionItr);

      solutionItr++;
      //ptr2=ptr;
    }
}
// This starts the insertion function
void insertSolution(Protein p1, Protein p2, Protein p3, default_random_engine * randomEngine )
{

	list<Acid>::iterator iterator;  ///general purpose list iterator
    list<Acid>::iterator litr; ///solution iterator

    Acid * ptr;

	normal_distribution<double> ndist(2,2);

	list<Acid> solutionList;

    int numAcids;
	int selectedId;


    cout << "How many Acids in a solution?" << endl;
    cin >> numAcids;

    double x =0;
	double y =0;
	double z =0;
	while (numAcids <3)
	{
		cout << "Enter at least 3 acids" << endl;
		cin >> numAcids;
	}

	for(int i=0; i<numAcids ; i++)
    {
	   ptr=new Acid;
	   ptr->setId(100+i+1);
	   ptr->setLocation(x,y,z);
	   solutionList.push_back(*ptr);
	   x=x+ndist(*randomEngine);
       y=y+ndist(*randomEngine);
	   z=z+ndist(*randomEngine);
    }

	// Distribution used to provide a connection or no connection
	uniform_real_distribution<double> noNodeDistribution(0, 1);
	// Distribution used to determine how many connections there will be
	normal_distribution<double> connectionsDistributions(3.0, 1.0);

	// Linking the solutions
	iterator = solutionList.begin();
	for (int i = 0; i < numAcids - 2; i++)
	{
		if (!(noNodeDistribution(*randomEngine) > .10))
		{
			iterator++;
			continue;
		}

		int numberOfConnections = (int)connectionsDistributions(*randomEngine);
		if (numberOfConnections < 2 || numberOfConnections >= numAcids) numberOfConnections = 2;
		int currentConnections [numberOfConnections];

		for (int j = 0; j < numberOfConnections; j++) currentConnections[j] = -1; // Initialize to -1's

		for (int j = 0; j < numberOfConnections; j++)
		{
			int windowOffset = numAcids-2;

			// Distribution for which connection to connect to
			uniform_int_distribution<int> nodeDistribution(i+1 , windowOffset);
			int connectionId = nodeDistribution(*randomEngine);


			bool unableToAdd = false;
            if (connectionId > numAcids - 1 || connectionId < i+1) unableToAdd= true;
			// Prevent duplicate connections
			bool preventDuplicate = false;
			for (int k = 0; k < numberOfConnections; k++)
			{
				// We already have that connection
				if (currentConnections [k] == connectionId  ) preventDuplicate = true;
			}

			if (!preventDuplicate && !unableToAdd)
			{
				currentConnections [j] = connectionId;
				Acid* acidToAdd = &getListElement(solutionList, connectionId);

				double dx = abs(iterator->getX() - acidToAdd->getX());
				double dy = abs(iterator->getY() - acidToAdd->getY());
				double dz = abs(iterator->getZ() - acidToAdd->getZ());
				iterator->addConnection(acidToAdd);
				iterator->addLength(sqrt(dx*dx + dy*dy + dz*dz));
			}
		}
		iterator++;
	}


    list<Acid> * myAcids = p1.getAcids();
   	int idRange = myAcids->size();
    cout << "Select Id to attach to solution in the range of "  << idRange-idRange+1 << " to " << idRange-1
    << endl;
    cin >> selectedId;
    while (selectedId > idRange-1 || selectedId < idRange-idRange+1)
    {
    	cout << "Enter a valid location ";
   		cin >> selectedId;
   	}


    attachSolution(solutionList, p1, selectedId);

    myAcids = p2.getAcids();
   	int idRange2 = myAcids->size();
    cout << "Select Id to attach to solution in the range of "  << idRange+1 << " to "
    << idRange2+idRange-1 << endl;
    cin >> selectedId;
    while (selectedId > idRange2+idRange-1 || selectedId <  idRange+1)
    {
    	cout << "Enter a valid location ";
   		cin >> selectedId;
   	}

    attachSolution(solutionList,p2, selectedId);

   	myAcids = p3.getAcids();
   	int idRange3 = myAcids->size();
    cout << "Select Id to attach to solution in the range of "  << idRange2+idRange+1 << " to "
    << idRange3+idRange2+idRange-1 << endl;
    cin >> selectedId;
     while (selectedId > idRange+idRange2+idRange3-1 || selectedId <  idRange2+idRange+1)
    {
   			cout << "Enter a valid location: ";
   			cin >> selectedId;
   	}

    attachSolution(solutionList, p3, selectedId);

}

void gnuPlotOutput(ofstream& p1plot, ofstream& p2plot, ofstream& p3plot, ofstream& comout, Locations * locations)
{
    p1plot.open("aminoAcids1.txt"); //First chains file of points matched
    p2plot.open("aminoAcids2.txt"); //Second chains file of points matched
    p3plot.open("aminoAcids3.txt"); //Third chains file of points matched
    comout.open("gnuPlotCommands.txt"); //command file for gnuplot

    list<tuple<double, double, double>> locations1 = locations->locations1;
    auto p1itr = locations1.begin();
    while (p1itr != locations1.end())
        {
            p1plot << get<0>(*p1itr) << ", " << get<1>(*p1itr) << ", " << get<2>(*p1itr) << endl;
            p1itr++;
        }

    list<tuple<double, double, double>> locations2 = locations->locations2;
    auto p2itr = locations2.begin();
    while (p2itr != locations2.end())
        {
            p2plot << get<0>(*p2itr) << ", " << get<1>(*p2itr) << ", " << get<2>(*p2itr) << endl;
            p2itr++;
        }

    list<tuple<double, double, double>> locations3 = locations->locations3;
    auto p3itr = locations3.begin();
    while (p3itr != locations3.end())
        {
            p3plot << get<0>(*p3itr) << ", " << get<1>(*p3itr) << ", " << get<2>(*p3itr) << endl;
            p3itr++;
        }
    comout << "set term png" << endl;
    comout << "set output \"aminoAcids.png\"" << endl;
    comout << "set xrange [-25:25]" << endl;
    comout << "set yrange [-25:25]" << endl;
    comout << "set zrange [-25:25]" << endl;
    comout << "splot 'aminoAcids1.txt', 'aminoAcids2.txt', 'aminoAcids3.txt' with circles" << endl;

    p1plot.close();
    p2plot.close();
    p3plot.close();
    comout.close();
}

Locations* findLocations(Protein protein1, Protein protein2, Protein protein3)
{
    Locations* locations = new Locations;

    list<Acid>::iterator iterator1 = protein1.getIterator();
    list<Acid>* a1 = protein1.getAcids();
    list<Acid>* a2 = protein2.getAcids();

    list<Acid>* a3 = protein3.getAcids();
    iterator1 = a1->begin();

    list<Acid>::iterator iterator2;

    auto compareDistance = [](list<Acid>::iterator iterator2,
                              double x, double y, double z,
                              Protein otherProtein, list<Acid>* otherAcids) -> tuple<bool, double, double, double> {
        iterator2 = otherProtein.getIterator();
        iterator2 = otherAcids->begin();
        while (iterator2 != otherAcids->end())
        {
            double ox = iterator2->getX();
            double oy = iterator2->getY();
            double oz = iterator2->getZ();
            double dx = abs(x - ox);
            double dy = abs(y - oy);
            double dz = abs(z - oz);

            double distance = sqrt(dx*dx + dy+dy + dz*dz);

            // Found acid in same spacial region
            if (distance <= .15D) return make_tuple(true, ox, oy, oz);

            iterator2++;
        }
        return make_tuple(false, -1, -1, -1);
    };

    while (iterator1 != a1->end())
    {
        double x = iterator1->getX();
        double y = iterator1->getY();
        double z = iterator1->getZ();

        tuple<bool, double, double, double> found1 = compareDistance(iterator2, x, y, z, protein2, a2);
        tuple<bool, double, double, double> found2 = compareDistance(iterator2, x, y, z, protein3, a3);

        if (get<0>(found1) && get<0>(found2))
        {
            locations->locations1.push_back(make_tuple(x, y, z));
            locations->locations2.push_back(make_tuple(get<1>(found1), get<2>(found1), get<3>(found1)));
            locations->locations3.push_back(make_tuple(get<1>(found2), get<2>(found2), get<3>(found2)));
        }

        iterator1++;
    }

    return locations;
}
