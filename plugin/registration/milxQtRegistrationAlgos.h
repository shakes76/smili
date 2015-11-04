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
#ifndef milxQtRegistrationAlgosAlgos_H
#define milxQtRegistrationAlgosAlgos_H
#include "milxQtMain.h"
#include <QObject>
#include <QtCore>
#include <QFuture>
#include "milxQtRegistrationParams.h"
#include "milxRegistration.h"

class milxQtRegistration;

#ifdef USE_NIFTI
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
#endif

#ifdef USE_ELASTIX
#include "itkUseMevisDicomTiff.h"
#include "elxElastixMain.h"
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include "itkObject.h"
#include "itkDataObject.h"
#include <itksys/SystemTools.hxx>
#include <itksys/SystemInformation.hxx>
#include "elxTimer.h"
#endif

/**
    \class milxQtRegistrationAlgos
    \brief Algorithms calls
*/
class milxQtRegistrationAlgos : public QObject
{
    Q_OBJECT

public:
    /*!
        \fn milxQtRegistrationAlgos::milxQtRegistrationAlgos(QObject * parent)
        \brief The standard constructor
    */
    milxQtRegistrationAlgos(QObject * parent);

    /*!
    \fn milxQtRegistrationAlgos::affine(milxQtRegistrationParams params)
    \brief Perform a ITK affine registration
    */
    int affine(milxQtRegistrationParams params);

    /*!
    \fn milxQtRegistrationAlgos::demon(milxQtRegistrationParams params)
    \brief Perform a ITK demon registration
    */
    int demon(milxQtRegistrationParams params);


    /*!
    \fn milxQtRegistrationAlgos::affine_async(milxQtRegistrationParams params)
    \brief Perform a ITK affine registration asynchroniously
    */
    void affine_async(milxQtRegistrationParams params);

    /*!
    \fn milxQtRegistrationAlgos::demon_async(milxQtRegistrationParams params)
    \brief Perform a ITK demon registration asynchroniously
    */
    void demon_async(milxQtRegistrationParams params);

#ifdef USE_NIFTI
    /*!
        \fn milxQtRegistrationAlgos::cpp2def(milxQtRegistrationParams params)
        \brief Cpp to deformation field transformation using Nifti Reg Library
    */
    int cpp2def(milxQtRegistrationParams params);

    /*!
        \fn milxQtRegistrationAlgos::f3d(milxQtRegistrationParams params)
        \brief Perform a f3d registration
    */
    int f3d(milxQtRegistrationParams params);

    /*!
        \fn milxQtRegistrationAlgos::aladin(milxQtRegistrationParams params)
        \brief Perform a aladin registration
    */
    int aladin(milxQtRegistrationParams params);

    /*!
        \fn milxQtRegistrationAlgos::average(QString outputName, QStringList filenames)
        \brief Perform the average of a list of images
    */
    int average(QString outputName, QStringList filenames);

    /*!
    	\fn milxQtRegistrationAlgos::similarities(milxQtRegistration * image)
    	\brief Perform similarities calculation
    */
    int similarities(milxQtRegistration * image);

    /*!
        \fn milxQtRegistrationAlgos::cpp2def_async(milxQtRegistrationParams params)
        \brief Perform a transformation from cpp to deformation field asynchroniously
    */
    void cpp2def_async(milxQtRegistrationParams params);

    /*!
        \fn milxQtRegistrationAlgos::cpp2def_async(milxQtRegistrationParams params)
        \brief Perform a f3d registration asynchroniously
    */
    void f3d_async(milxQtRegistrationParams params);

    /*!
        \fn milxQtRegistrationAlgos::aladin_async(milxQtRegistrationParams params)
        \brief Perform a aladin registration asynchroniously
    */
    void aladin_async(milxQtRegistrationParams params);

    /*!
        \fn milxQtRegistrationAlgos::average_async(QString outputName, QStringList filenames);
        \brief Perform the average of a list of images asynchroniously
    */
    void average_async(QString outputName, QStringList filenames);

    /*!
    	\fn milxQtRegistrationAlgos::similarities_async(milxQtRegistration * image)
    	\brief Perform similarities calculation asynchroniously
    */
    void similarities_async(milxQtRegistration * image);
#endif


#ifdef USE_ELASTIX
    /*!
    \fn milxQtRegistrationAlgos::elastix(milxQtRegistrationParams params)
    \brief Perform an elastix registration
    */
    int elastix(milxQtRegistrationParams params);

    /*!
    \fn milxQtRegistrationAlgos::elastix_async(milxQtRegistrationParams params)
    \brief Perform an elastix registration asynchroniously
    */
    void elastix_async(milxQtRegistrationParams params);

    /*!
    \fn milxQtRegistrationAlgos::elastixClean(milxQtRegistrationParams params)
    \brief Remove unused files created by elastix
    */
    void elastixClean(milxQtRegistrationParams params);
#endif

signals:
    /*!
    \fn milxQtRegistrationAlgos::registrationCompleted()
    \brief The registration has been completed
    */
    void registrationCompleted();

    /*!
    \fn milxQtRegistrationAlgos::error(QString functionName, QString errorMsg)
    \brief An error happened during the registration
    */
    void error(QString functionName, QString errorMsg);

#ifdef USE_NIFTI
    /*!
        \fn milxQtRegistrationAlgos::cpp2defCompleted()
        \brief The transformation from cpp to deformation field has been completed
    */
    void cpp2defCompleted();

    /*!
        \fn milxQtRegistrationAlgos::averageCompleted()
        \brief The average of the images has been completed
    */
    void averageCompleted();

    /*!
    \fn milxQtRegistrationAlgos::similaritiesComputed()
    \brief similarity of the image has been computed
    */
    void similaritiesComputed();

#endif

protected:
    QFuture<int> future; //!< Used to start the functions asynchroniously
};


#endif // MILXQTRegistrationPLUGIN_H
