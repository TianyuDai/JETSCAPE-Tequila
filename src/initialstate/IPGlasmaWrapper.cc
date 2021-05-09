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

#include "IPGlasmaWrapper.h"

// Register the module with the base class
RegisterJetScapeModule<IPGlasmaWrapper> IPGlasmaWrapper::reg("IPGlasma");

IPGlasmaWrapper::IPGlasmaWrapper() {
  SetId("IPGlasma");
  event_id_ = -1;
}

IPGlasmaWrapper::~IPGlasmaWrapper() {}

void IPGlasmaWrapper::InitTask() {}

void IPGlasmaWrapper::Exec() {
  Clear();
  Jetscape::JSINFO << "Run IPGlasma ...";
  try {
    event_id_++;
  } catch (std::exception &err) {
    Jetscape::JSWARN << err.what();
    std::exit(-1);
  }
}

void IPGlasmaWrapper::Clear() {
  Jetscape::JSINFO << "clear initial condition vectors";
}

void IPGlasmaWrapper::Write(weak_ptr<JetScapeWriter> w) {}
