#ifndef ACID_H
#define ACID_H

#include <list>

using namespace std;

class Acid
{
public:
    Acid();
    Acid(double x, double y, double z);
	~Acid();
	void display();

	void setId(int id);

	int getId();

	void addConnection(Acid * acid);
	void addLength(double length);

	void setLocation(double x, double y, double z);
	void ajustLocation(double x, double y, double z);

	double getX();
	double getY();
	double getZ();

	list<Acid*> getConnections();
	list<double> getLengths();

private:
	int m_id;
	double m_x;
	double m_y;
	double m_z;
	list<Acid*> m_connections;
	list<double> m_lengths;
};

#endif
