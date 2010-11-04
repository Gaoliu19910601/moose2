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

#ifndef NEARESTNODEVALUEAUX_H
#define NEARESTNODEVALUEAUX_H

#include "AuxKernel.h"
#include "NearestNode.h"


//Forward Declarations
class NearestNodeValueAux;

template<>
InputParameters validParams<NearestNodeValueAux>();

/** 
 * Constant auxiliary value
 */
class NearestNodeValueAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  NearestNodeValueAux(const std::string & name, MooseSystem & moose_system, InputParameters parameters);

  virtual ~NearestNodeValueAux() {}

  virtual void setup();
  
protected:
  virtual Real computeValue();

  NearestNode _nearest_node;

  unsigned int _paired_variable;
};

#endif //NEARESTNODEVALUEAUX_H
