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
    \brief Algorithms calls
*/
class milxQtRegistrationNifti : public QObject
{
    Q_OBJECT

public:
    /*!
        \fn milxQtRegistrationNifti::milxQtRegistrationNifti(QObject * parent)
        \brief The standard constructor
    */
    milxQtRegistrationNifti(QObject * parent);

    /*!
        \fn milxQtRegistrationNifti::cpp2def(QObject * parent)
        \brief Cpp to deformation field transformation using Nifti Reg Library
    */
    int cpp2def(ParamsF3D params);

    /*!
        \fn milxQtRegistrationNifti::f3d(ParamsF3D params)
        \brief Perform a f3d registration
    */
    int f3d(ParamsF3D params);

    /*!
        \fn milxQtRegistrationNifti::aladin(ParamsAladin params)
        \brief Perform a aladin registration
    */
    int aladin(ParamsAladin params);

    /*!
        \fn milxQtRegistrationNifti::average(QString outputName, QStringList filenames)
        \brief Perform the average of a list of images
    */
    int average(QString outputName, QStringList filenames);

    /*!
        \fn milxQtRegistrationNifti::cpp2def_async(ParamsF3D params)
        \brief Perform a transformation from cpp to deformation field asynchroniously
    */
    void cpp2def_async(ParamsF3D params);

    /*!
        \fn milxQtRegistrationNifti::cpp2def_async(ParamsF3D params)
        \brief Perform a f3d registration asynchroniously
    */
    void f3d_async(ParamsF3D params);

    /*!
        \fn milxQtRegistrationNifti::aladin_async(ParamsAladin params)
        \brief Perform a aladin registration asynchroniously
    */
    void aladin_async(ParamsAladin params);

    /*!
        \fn milxQtRegistrationNifti::average_async(ParamsF3D params)
        \brief Perform the average of a list of images asynchroniously
    */
    void average_async(QString outputName, QStringList filenames);


signals:
    /*!
        \fn milxQtRegistrationNifti::cpp2defCompleted()
        \brief The transformation from cpp to deformation field has been completed
    */
    void cpp2defCompleted();

    /*!
        \fn milxQtRegistrationNifti::registrationCompleted()
        \brief The registration has been completed
    */
    void registrationCompleted();

    /*!
        \fn milxQtRegistrationNifti::averageCompleted()
        \brief The average of the images has been completed
    */
    void averageCompleted();


protected:
    QFuture<int> future; //!< Used to start the functions asynchroniously
};


#endif // MILXQTRegistrationPLUGIN_H
