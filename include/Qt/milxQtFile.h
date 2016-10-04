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
#ifndef MILXQTFILE
#define MILXQTFILE

#include <QFile>
#include <QPointer>
//VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataCollection.h>
#include <vtkImageData.h>
#include <vtkCommand.h>
#include <vtkEventQtSlotConnect.h>

#include "milxQtImage.h"
#include "milxQtModel.h"

#include "milxFile.h"

/*!
    \class milxQtFile
    \brief This class represents the MILX Qt File I/O object using VTK/ITK/Qt.
    \author Shekhar S. Chandra, 2013

    The class opens and saves various data into a number of different formats via the VTK, ITK and Qt libraries transparently.

    Supported Image Formats:
    Images (*.png *.jpeg *.jpg *.bmp *.tiff *.tif)
    Medical Images (*.nii *.ima *.dcm *.dicom *.mhd) with gz extensions
    Other Images (*.vti *.mrc *.rec)

    Supported Model Formats:
    Model or Polygonal Files (*.vtp *.vtk *.ply *.obj). OBJ file format has read-only support.

    Usage Examples:
    Open a model (Poly Data) file:
    \code
    QPointer<milxQtFile> reader = new milxQtFile;
    vtkSmartPointer<vtkPolyData> modelData;

    reader->openModel(filename, modelData);
    \endcode
*/
class MILXQT_EXPORT milxQtFile : public QFile
{
    Q_OBJECT

public:
    /*!
        \fn milxQtFile::milxQtFile(QObject *theParent = 0)
        \brief The standard constructor
    */
    milxQtFile(QObject *theParent = 0);
    /*!
        \fn milxQtFile::~milxQtFile()
        \brief The standard destructor
    */
    virtual ~milxQtFile();

    /*!
        \brief ITK interfacing member for exceptions.
    */
    virtual inline const char * GetNameOfClass() const
    {
      return "milxQtFile";
    }

public slots:
    /*!
        \fn milxQtFile::supportedImageFormats()
        \brief Returns a string of supported image formats from all libraries for in file dialogs in Qt.

        The string looks like "*.nii.gz *.img" etc.
    */
    QString supportedImageFormats();
    /*!
        \fn milxQtFile::is8BitFormat(const QString filename, bool &errorEncountered)
        \brief Returns if a image is an 8-bit image format or not.

        Also writes Boolean into errorEncountered if error was encountered during reading of header.
    */
    bool is8BitFormat(const QString filename, bool &errorEncountered);
    /*!
    \fn milxQtFile::is32BitFormat(const QString filename, bool &errorEncountered)
    \brief Returns if a image is an 32-bit image format or not.

    Also writes Boolean into errorEncountered if error was encountered during reading of header.
    */
    bool is32BitFormat(const QString filename, bool &errorEncountered);
    /*!
        \fn milxQtFile::isFieldFormat(const QString filename, bool &errorEncountered)
        \brief Returns if a image is a deformation field type image (a vector image) format or not.

        Also writes Boolean into errorEncountered if error was encountered during reading of header.
    */
    bool isFieldFormat(const QString filename, bool &errorEncountered);
    /*!
        \fn milxQtFile::openImage(const QString filename, vtkImageData* data)
        \brief Opens an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, HDR etc.

        You need to pre-allocate data, as the result is copied to data.

        Returns true if successful. ImageData is allocated within this member, so pass a NULL pointer.
        Image is also flipped since the ITK reader orientation is different to VTK images.
    */
    bool openImage(const QString filename, vtkImageData* data);
    /*!
        \fn milxQtFile::openImage(const QString filename, milxQtImage* data)
        \brief Opens an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, HDR etc. Overloaded for milxQtImage.

        Returns true if successful.
    */
    bool openImage(const QString filename, milxQtImage* data);
    /*!
        \fn milxQtFile::openImageSeries(milxQtImage* data, QString &directoryPath)
        \brief Opens an DICOM series.

        Returns true if successful.
    */
    bool openImageSeries(milxQtImage* data, QString directoryPath = "");

    /*!
        \fn milxQtFile::saveImage(const QString filename, vtkImageData* data)
        \brief Saves data as an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, HDR etc.

        Returns true if successful. ImageData is allocated within this member, so pass a NULL pointer.
        Image is also flipped since the ITK reader orientation is different to VTK images.
    */
    bool saveImage(const QString filename, vtkImageData* data);
    /*!
        \fn milxQtFile::saveImage(const QString filename, milxQtImage* data)
        \brief Saves data as an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, HDR etc. Overloaded for milxQtImage.

        Returns true if successful.
    */
    bool saveImage(const QString filename, milxQtImage* data);

    /*!
        \fn milxQtFile::openModel(const QString filename, vtkPolyData* data)
        \brief Opens a model file, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj).

        You need to pre-allocate data, as the result is copied to data.

        Returns true if successful. PolyData is allocated within this member, so pass a NULL pointer.
    */
    bool openModel(const QString filename, vtkPolyData* data);
    /*!
        \fn milxQtFile::openModel(const QString filename, milxQtModel* data)
        \brief Opens a model file, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj).

        Returns true if successful.
    */
    bool openModel(const QString filename, milxQtModel* data);
    /*!
        \fn milxQtFile::openModelCollection(vtkPolyDataCollection* collection, QStringList &filenames)
        \brief Opens a series of model files, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj) into a PolyData collection.

        Filenames are asked for via a file dialog. Assumes collection is allocated but not necessarily (re)sized. It is sized accordingly internally.

        Returns true if successful.
    */
    bool openModelCollection(vtkPolyDataCollection* collection, QStringList &filenames);
    inline bool openModelCollection(vtkPolyDataCollection* collection)
    {
        QStringList tmpList;
        return openModelCollection(collection, tmpList); //List is discarded
    }
    /*!
        \fn milxQtFile::saveModel(const QString filename, vtkPolyData* data, const bool binary = false)
        \brief Saves a model as a file, which can be a VTK XML or Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk) or a Polygonal File (*.ply).

        Returns true if successful.
    */
    bool saveModel(const QString filename, vtkPolyData* data, const bool binary = false);
    /*!
        \fn milxQtFile::saveModel(const QString filename, milxQtModel* data, const bool binary = false)
        \brief Saves a model as a file, which can be a VTK XML or Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk) or a Polygonal File (*.ply).

        Returns true if successful.
    */
    bool saveModel(const QString filename, milxQtModel* data, const bool binary = false);
    /*!
        \fn milxQtFile::saveScalarsOfModel(const QString filename, milxQtModel* data)
        \brief Saves the scalars of a model as a CSV or VTK Text file. The CSV format is output if the extension is *.csv or *.txt.

        Returns true if successful.
    */
    bool saveScalarsOfModel(const QString filename, milxQtModel* data);

    ///Convert an QStringList to a vector of STL strings.
    inline std::vector< std::string > convertQStringList(const QStringList filenames)
    {
        std::vector< std::string > names;
        foreach(QString name, filenames)
            names.push_back(name.toStdString());

        return names;
    }

    inline QString getPixelType()
    {
      return dataPixelType;
    }
    inline QString getComponentType()
    {
      return dataComponentType;
    }
    inline size_t getNumberOfDimensions()
    {
      return dataDimensions;
    }
    inline size_t getNumberOfComponents()
    {
      return dataComponents;
    }

    void linkProgressEventOf(vtkObject * obj);
    inline void updateQtEvents()
    {
        qApp->processEvents();
    }

protected:
    QString name; //!< Name of the last opened file
    QString dataPixelType; //!< pixel type of image read
    QString dataComponentType; //!< Component type of image read
    size_t dataComponents; //!< Components of vector of image read
    size_t dataDimensions; //!< Dimensions of image read

    //Observer for the Qt event loop
    itkEventQtObserver::Pointer observeProgress;

    vtkSmartPointer<vtkEventQtSlotConnect> Connector; //!< VTK Events to slots convertor

private:

};

#endif // MILXQTFILE
