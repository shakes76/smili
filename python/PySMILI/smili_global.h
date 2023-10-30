/*
  This file is part of SMILI.
*/

#pragma once

// Make "signals:", "slots:" visible as access specifiers
#define QT_ANNOTATE_ACCESS_SPECIFIER(a) __attribute__((annotate(#a)))

// Define PYTHON_BINDINGS
#define PYTHON_BINDINGS

//SMILI
#include "milxFile.h"
//milxQt
#include "milxQtWindow.h"
#include "milxQtRenderWindow.h"
#include "milxQtModel.h"
#include "milxQtImage.h"
#include "milxQtFile.h"

