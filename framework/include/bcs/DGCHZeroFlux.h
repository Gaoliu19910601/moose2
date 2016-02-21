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

#ifndef DGCHZEROFLUX_H
#define DGCHZEROFLUX_H

#include "IntegratedBC.h"

//Forward Declarations
class DGCHZeroFlux;

template<>
InputParameters validParams<DGCHZeroFlux>();

class DGCHZeroFlux : public IntegratedBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  DGCHZeroFlux(const InputParameters & parameters);

  virtual ~DGCHZeroFlux() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  VariableSecond & _second_u;

private:
  Real _alpha;
};

#endif //DGCHZEROFLUX_H
