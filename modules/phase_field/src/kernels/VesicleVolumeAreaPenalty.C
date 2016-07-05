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

#include "VesicleVolumeAreaPenalty.h"


template<>
InputParameters validParams<VesicleVolumeAreaPenalty>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("alpha_v", 1000, "The penalty parameter of vesicle volume.");
  params.addParam<Real>("alpha_a", 1000, "The penalty parameter of vesicle area."); 
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  return params;
}

VesicleVolumeAreaPenalty::VesicleVolumeAreaPenalty(const InputParameters & parameters) :
    Kernel(parameters),
    _alpha_v(getParam<Real>("alpha_v")),
    _alpha_a(getParam<Real>("alpha_a"))
    _epsilon(getParam<Real>("epsilon")), 
    _second_phi(secondPhiFace()),
    _second_test(secondTestFace()),
    _second_u(second())
{
}

VesicleVolumeAreaPenalty::~VesicleVolumeAreaPenalty()
{
}

Real
VesicleVolumeAreaPenalty::computeQpResidual()
{
  Real r = 0;

  r += _alpha_v * _test[_i][_qp] * (_volume - _volume_0);

  r += _alpha_a * _test[_i][_qp] * (_area - _area_0) * (-3.0/std::sqrt(2.0) * _epsilon * _second_u.tr());

  return r;
}

Real
VesicleVolumeAreaPenalty::computeQpJacobian()
{
  Real r = 0;

  r += _alpha_a * _test[_i][_qp] * (_area - _area_0) * (-3.0/std::sqrt(2.0) * _epsilon * _second_phi[_j].tr());

  return r;
}

