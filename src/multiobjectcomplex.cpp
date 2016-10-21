//--------------------------------------------------------------------------
// MultiObjectComplex: Collection of multiple ShapeObjects
//--------------------------------------------------------------------------

#include "multiobjectcomplex.h"

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------
// MultiObjectComplex
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------
MultiObjectComplex::MultiObjectComplex()
{

}

//----------------------------------------------------------------
// MultiObjectComplex
//----------------------------------------------------------------
// Inputs:
//   multiObject - MultiObjectComplex object to copy
//
// Outputs:
//----------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------
MultiObjectComplex::MultiObjectComplex(const MultiObjectComplex& multiObject)
{
}

//----------------------------------------------------------------
// MultiObjectComplex
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//----------------------------------------------------------------
// Destructor
//----------------------------------------------------------------
MultiObjectComplex::~MultiObjectComplex()
{
}

//----------------------------------------------------------------
// AddShape
//----------------------------------------------------------------
// Inputs:
//   shape - the ShapeObject to add
//
// Outputs:
//----------------------------------------------------------------
// Adds a shape to the end of the collection
//----------------------------------------------------------------
void MultiObjectComplex::AddShape(ShapeObject* shape)
{
    this->_shapes.push_back(shape);
}

//----------------------------------------------------------------
// SetShapeAt
//----------------------------------------------------------------
// Inputs:
//   index - the index in the collection to set
//   shape - the ShapeObject to set at the index
//
// Outputs:
//----------------------------------------------------------------
// Sets a shape at a given index
//----------------------------------------------------------------
void MultiObjectComplex::SetShapeAt(int index, ShapeObject *shape)
{
    this->_shapes[index] = shape;
}

//----------------------------------------------------------------
// GetNumberOfShapes
//----------------------------------------------------------------
// Inputs:
//
// Outputs:
//   return - the number of shapes in the multi-object complex
//----------------------------------------------------------------
// Returns the number of shapes in the multi-object complex
//----------------------------------------------------------------
int MultiObjectComplex::GetNumberOfShapes()
{
    return this->_shapes.size();
}

//----------------------------------------------------------------
// GetShapeAt
//----------------------------------------------------------------
// Inputs:
//   index - the index in the collection to return
//
// Outputs:
//----------------------------------------------------------------
// Returns the shape at the given index
//----------------------------------------------------------------
ShapeObject* MultiObjectComplex::GetShapeAt(int index)
{
    return this->_shapes[index];
}
