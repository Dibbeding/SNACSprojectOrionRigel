#include "rigel.h"
#include "node.h"

extern int DIMENSIONS;

void Node::Refresh()
{
	for(int i=0;i<DIMENSIONS;i++)
	{
		vec[i]=0;
	}
}
Node::Node () {
	vec = new double [DIMENSIONS];
	for(int i = 0;i<DIMENSIONS;i++)
	{
		vec[i]=0;
	}
}
Node::~Node () {
  delete [] vec;
}
