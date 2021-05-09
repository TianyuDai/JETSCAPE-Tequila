/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef IPGLASMAWRAPPER_H
#define IPGLASMAWRAPPER_H

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "JetScapeModuleBase.h"
#include "InitialState.h"
#include "JetScapeLogger.h"

using namespace Jetscape;

class IPGlasmaWrapper : public Jetscape::InitialState {
  // this is wrapper class to read external files that
  // stores initial number of binary collisions and corresponding
  // configurations
public:
  IPGlasmaWrapper();
  ~IPGlasmaWrapper();

  /** Reads the input parameters from the XML file under the tag  <IS>. Calls InitTask(); This explicit call of InitTask() can be used for actual initialization of modules such as @a Trento if attached as a @a polymorphic class. It also initializes the tasks within the current module.
      @sa Read about @a polymorphism in C++.
   */
  //void Init();

  /** Default Exec() function. It can be overridden by other tasks.
   */
  void Exec();

  /** Default Clear() function. It can be overridden by other tasks.
   */
  void Clear();

  void InitTask();

  /** Default Write() function. It can be overridden by other tasks.
      @param w A pointer to the JetScapeWriter class.
   */
  virtual void Write(weak_ptr<JetScapeWriter> w);

  /** Generated number of collision participants.
  */
  double GetNpart() { return npart; };

  /** Generated number of binary collisions.
  */
  double GetNcoll() { return ncoll; };

  /** Generated total entropy
  */
  double GetTotalEntropy() { return totalentropy; };

private:
  int dim_x_, dim_y_;

  double npart = -1;
  double ncoll = -1;
  double totalentropy = -1;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<IPGlasmaWrapper> reg;
};

#endif  // IPGLASMAWRAPPER_H
