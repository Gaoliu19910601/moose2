/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef VESICLEVOLUME_H
#define VESICLEVOLUME_H

#include "ElementUserObject.h"

//Forward Declarations
class VesicleVolume;

template<>
InputParameters validParams<VesicleVolume>();

class VesicleVolume : public ElementUserObject
{
public:
  ElementUserObject(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual Real getValue();

protected:
  virtual Real computeQpIntegral();
  virtual Real computeIntegral();

  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  unsigned int _qp;

  Real _integral_value;
};

#endif
