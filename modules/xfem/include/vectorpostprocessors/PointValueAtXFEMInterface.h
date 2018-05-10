//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POINTVALUEATXFEMINTERFACE_H
#define POINTVALUEATXFEMINTERFACE_H

// MOOSE includes
#include "GeneralVectorPostprocessor.h"
#include "CoupleableMooseVariableDependencyIntermediateInterface.h"
#include "MooseVariableInterface.h"
#include "SamplerBase.h"
#include "ElementPairLocator.h"

// Forward Declarations
class PointValueAtXFEMInterface;

class XFEM;
class MovingLineSegmentCutSetUserObject;

template <>
InputParameters validParams<PointValueAtXFEMInterface>();

class PointValueAtXFEMInterface : public GeneralVectorPostprocessor,
                                  public CoupleableMooseVariableDependencyIntermediateInterface,
                                  public MooseVariableInterface<Real>,
                                  protected SamplerBase
{
public:
  PointValueAtXFEMInterface(const InputParameters & parameters);

  virtual ~PointValueAtXFEMInterface() {}

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

protected:
  /**
   * Find the local element that contains the point.  This will attempt to use a cached element to
   * speed things up.
   *
   * @param p The point in physical space
   * @return The Elem containing the point or NULL if this processor doesn't contain an element that
   * contains this point.
   */
  const Elem * getLocalElemContainingPoint(const Point & p);

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// The points to evaluate at
  std::vector<Point> _points;

  /// The ID to use for each point (yes, this is Real on purpose)
  std::vector<Real> _ids;

  /// Vector of values per point
  std::vector<std::vector<Real>> _point_values;

  /// Whether or not the Point was found on this processor (short because bool and char don't work with MPI wrappers)
  std::vector<short> _found_points;

  unsigned int _qp;

  std::unique_ptr<PointLocatorBase> _pl;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;

  const ElementPairLocator::ElementPairList * _elem_pairs;

  const MovingLineSegmentCutSetUserObject * _geo_cut;

  bool _positive_level_set;

  /// The variable number of the level set variable we are operating on
  const unsigned int _level_set_var_number;

  /// system reference
  const System & _system;

  /// the subproblem solution vector
  const NumericVector<Number> * _solution;
};

#endif