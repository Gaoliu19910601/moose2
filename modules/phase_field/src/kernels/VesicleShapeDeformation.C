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

#include "VesicleShapeDeformation.h"


template<>
InputParameters validParams<VesicleShapeDeformation>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  params.addParam<Real>("flexural_rigidity", 100, "Flexural rigidity of vesicle.");
  return params;
}

VesicleShapeDeformation::VesicleShapeDeformation(const InputParameters & parameters) :
    Kernel(parameters),
    _epsilon(getParam<Real>("epsilon")),
    _flexural_rigidity(getParam<Real>("flexural_rigidity")),
    _second_phi(secondPhiFace()),
    _second_test(secondTestFace()),
    _second_u(second())
{
}

VesicleShapeDeformation::~VesicleShapeDeformation()
{
}

Real
VesicleShapeDeformation::computeQpResidual()
{
  Real r = 0;

  r += _epsilon * (_second_test[_i][_qp].tr() * _second_u[_qp].tr()); // eps*Lap_w*Lap_u

  r += 2.0/_epsilon * (_grad_test[_i][_qp] * _grad_u[_qp]) * (3.0 * _u[_qp] * _u[_qp] + 2.0 * _flexural_rigidity * _epsilon * _u[_qp] - 1.0); // 2/eps*(Grad_w*Grad_w)*(3*u^2 + 2*C*eps*u -1.0)

  r += 1.0/_epsilon * _test[_i][_qp] * (6.0 * _u[_qp] + 2.0 * _flexural_rigidity * _epsilon) * (_grad_u[_qp] * _grad_u[_qp]); // 1/eps*w*(6*u+2*C*eps)(Grad_u*Grad_u)

  r += pow(1.0/_epsilon, 3.0) * _test[_i][_qp] * (3.0 * _u[_qp] * _u[_qp] + 2.0 * _flexural_rigidity * _epsilon * _u[_qp]) * (_u[_qp] * _u[_qp] - 1.0) * (_u[_qp] + _flexural_rigidity * _epsilon); // 1/eps^3*w*(3*u^2+2*C*eps*u-1)(u^2-1)(u+C*eps)

  return r;
}

Real
VesicleShapeDeformation::computeQpJacobian()
{
  Real r = 0;
  // eps*Lap_w*Lap_u
  r += _epsilon * _second_test[_i][_qp].tr() * _second_phi[_j][_qp].tr();

  // 2/eps*(Grad_w*Grad_w)*(3*u^2 + 2*C*eps*u -1.0)
  r += 2.0/_epsilon * (_grad_test[_i][_qp] * _grad_phi[_j][_qp]) * (3.0 * _u[_qp] * _u[_qp] + 2.0 * _flexural_rigidity * _epsilon * _u[_qp] - 1.0);
  r += 2.0/_epsilon * (_grad_test[_i][_qp] * _grad_u[_qp]) * (6.0 * _phi[_j][_qp] * 2.0 * _flexural_rigidity * _epsilon * _phi[_j][_qp]);
  
  // 1/eps*w*(6*u+2*C*eps)(Grad_u*Grad_u)
  r += 1.0/_epsilon * _test[_i][_qp] * (6.0 * _phi[_j][_qp] + 2.0 * _flexural_rigidity * _epsilon) * (_grad_u[_qp] * _grad_u[_qp]);
  r += 1.0/_epsilon * _test[_i][_qp] * (6.0 * _u[_qp] + 2.0 * _flexural_rigidity * _epsilon) * (2.0 * _grad_phi[_j][_qp] * _grad_u[_qp]);

  // 1/eps^3*w*(3*u^2+2*C*eps*u-1)(u^2-1)(u+C*eps)
  r += pow(1.0/_epsilon, 3.0) * _test[_i][_qp] * (6.0 * _phi[_j][_qp] + 2.0 * _flexural_rigidity * _epsilon * _phi[_j][_qp]) * (_u[_qp] * _u[_qp] - 1.0) * (_u[_qp] + _flexural_rigidity * _epsilon);
  r += pow(1.0/_epsilon, 3.0) * _test[_i][_qp] * (3.0 * _u[_qp] * _u[_qp] + 2.0 * _flexural_rigidity * _epsilon * _u[_qp] - 1.0) * (2.0 * _u[_qp] * _phi[_j][_qp]) * (_u[_qp] + _flexural_rigidity * _epsilon);
  r += pow(1.0/_epsilon, 3.0) * _test[_i][_qp] * (3.0 * _u[_qp] * _u[_qp] + 2.0 * _flexural_rigidity * _epsilon * _u[_qp] - 1.0) * (_u[_qp] * _u[_qp] - 1.0) * _phi[_j][_qp];

  return r;
}

