/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "VesicleIC.h"

template<>
InputParameters validParams<VesicleIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addParam<Point>("center", Point(0,0,0), "The center coordinate");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter."); 
  return params;
}

VesicleIC::VesicleIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    _center(parameters.get<Point>("center")),
    _epsilon(getParam<Real>("epsilon"))
{}

Real
VesicleIC::value(const Point & p)
{
  Real value = 0.0;

  value = pow((p - _center).norm_sq(),0.5) - 0.25;

  value = std::tanh(value/sqrt(2.0)/_epsilon);

  return value;
}

