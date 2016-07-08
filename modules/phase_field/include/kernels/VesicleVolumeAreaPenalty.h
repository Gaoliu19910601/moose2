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

#ifndef VESICLEVOLUMEAREAPENALTY_H
#define VESICLEVOLUMEAREAPENALTY_H

#include "Kernel.h"
#include "VesicleVolume.h"
#include "VesicleArea.h"

class VesicleVolumeAreaPenalty;

template<>
InputParameters validParams<VesicleVolumeAreaPenalty>();

class VesicleVolumeAreaPenalty : public Kernel
{
public:
  VesicleVolumeAreaPenalty(const InputParameters & parameters);

  virtual ~VesicleVolumeAreaPenalty();

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  Real _alpha_v;

  Real _alpha_a;

  Real _epsilon;

  Real _volume, _volume_0;
  Real _area, _area_0;

  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;

  const VesicleArea & _vesicle_area_uo;
  const VesicleVolume & _vesicle_volume_uo;

};

#endif /* VESICLEVOLUMEAREAPENALTY_H */
