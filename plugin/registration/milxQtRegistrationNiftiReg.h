/*=========================================================================
  The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO)
  ABN 41 687 119 230.
  All rights reserved.

  Licensed under the CSIRO BSD 3-Clause License
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License in the file LICENSE.md or at

  https://stash.csiro.au/projects/SMILI/repos/smili/browse/license.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
=========================================================================*/
#ifndef MILXQTRegistrationNiftiReg_H
#define MILXQTRegistrationNiftiReg_H
#include "milxQtMain.h"
#include <QObject>
#include <QtCore>
#include <QFuture>

#include "_reg_ReadWriteImage.h"
#include "_reg_f3d2.h"
#include "_reg_ReadWriteImage.h"
#include "_reg_resampling.h"
#include "_reg_globalTransformation.h"
#include "_reg_localTransformation.h"
#include "_reg_tools.h"
#include "_reg_thinPlateSpline.h"
#include "_reg_aladin.h"
#include <fstream>
#include <vector>
#include <iostream>
#include "milxQtRegistrationStructures.h"


#ifdef _USE_NR_DOUBLE
#define PrecisionTYPE double
#else
#define PrecisionTYPE float
#endif

#ifdef _USE_CUDA
#include "_reg_f3d_gpu.h"
#endif
#include "float.h"
#include <limits>

#ifdef _WINDOWS
#include <time.h>
#endif


/**
    \class milxQtRegistrationNifti
    \brief NiftiReg calls
    \author 
*/
class milxQtRegistrationNifti : public QObject
{
	Q_OBJECT

public:
    milxQtRegistrationNifti();
    virtual ~milxQtRegistrationNifti();
	int cpp2def(PARAMSCPP2DEF params);
	int f3d(PARAMSF3D params);
	int aladin(PARAMSALADIN params);
	
	void cpp2def_async(PARAMSCPP2DEF params);
	void f3d_async(PARAMSF3D params);
	void aladin_async(PARAMSALADIN params);

signals:
	void cpp2defFinished();
	void registrationFinished();


protected:
	QFuture<int> future;
};


#endif // MILXQTRegistrationPLUGIN_H
