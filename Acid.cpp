#include "Acid.h"
#include <iostream>

Acid::Acid()
{
}

Acid::Acid(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

Acid::~Acid()
{
}

void Acid::display()
{
	cout << "Acid #" << m_id << endl;
	cout << "(" << m_x << ", " << m_y << ", " << m_z << ")\n";
	cout << "Connections:\n";
	if (m_connections.size() == 0) cout << "None\n";
	list<Acid*>::iterator iterator = m_connections.begin();
	while (iterator != m_connections.end())
	{
		Acid * p = *iterator;
		cout << p->getId() << endl;

		iterator++;
	}
	cout << endl;
}

void Acid::setId(int id)
{
	m_id = id;
}

int Acid::getId()
{
	return m_id;
}

void Acid::addConnection(Acid * acid)
{
	m_connections.push_back(acid);
}

void Acid::addLength(double length)
{
	m_lengths.push_back(length);
}

void Acid::setLocation(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

void Acid::ajustLocation(double x, double y, double z)
{
	m_x += x;
	m_y += y;
	m_z += z;
}

double Acid::getX()
{
	return m_x;
}

double Acid::getY()
{
	return m_y;
}

double Acid::getZ()
{
	return m_z;
}

list<Acid*> Acid::getConnections()
{
	return m_connections;
}

list<double> Acid::getLengths()
{
	return m_lengths;
}

