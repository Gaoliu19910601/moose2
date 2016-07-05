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

#include "VesicleVolume.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<VesicleVolume>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredCoupledVar("variable", "The name of the variable that this userobject applies to");
  return params;
}

VesicleVolume::VesicleVolume(const InputParameters & parameters) :
    ElementUserObject(parameters),
    MooseVariableInterface(this, false),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _qp(0),
    _integral_value(0)
{}

void
VesicleVolume::initialize()
{
  _integral_value = 0;
}

void
VesicleVolume::execute()
{
  _integral_value += computeIntegral();
}

Real
VesicleVolume::getValue()
{
  gatherSum(_integral_value);
  return _integral_value;
}

void
VesicleVolume::threadJoin(const UserObject & y)
{
  const VesicleVolume & pps = static_cast<const VesicleVolume &>(y);
  _integral_value += pps._integral_value;
}

Real
VesicleVolume::computeIntegral()
{
  Real sum = 0;

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
    sum += _JxW[_qp]*_coord[_qp]*computeQpIntegral();
  return sum;
}

Real
VesicleVolume::computeQpIntegral()
{
  Real value = _u[_qp];
  return value;
}
