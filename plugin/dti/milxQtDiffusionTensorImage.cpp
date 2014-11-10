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
#include "milxQtDiffusionTensorImage.h"

#include <vtkSphereSource.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkAppendPolyData.h>

#include "vtkDiffusionTensorGlyphFilter.h"

typedef vtkDiffusionTensorGlyphFilter::VectorImageType::SpacingType SpacingType;

void CalcGlyph(void *arg)
{
  vtkDiffusionTensorGlyphFilter *glyphFilter = (vtkDiffusionTensorGlyphFilter*) arg;

  if(!glyphFilter)
    {
    std::cerr << "vtkDiffusionTensorGlyphFilter is not valid!" << std::endl;
    }

  ///Associated Tensor data
  std::vector<vectorImageType::IndexType> indices = glyphFilter->GetTensorImageIndices();
  vtkDiffusionTensorGlyphFilter::VectorImageType::Pointer image = glyphFilter->GetTensorImage();
  const vectorImageType::PixelType normVector = glyphFilter->GetNormaliseVector();

  ///Get the SH index in image
  vectorImageType::IndexType index = indices[glyphFilter->GetPointId()];
  ///Get the SH itself
  const vectorImageType::PixelType harmonicsVector = image->GetPixel(index);

//  SpacingType spacing = glyphFilter->GetTensorImage()->GetSpacing();
  double radiusScaling = 2.5;

  double pointCoords[3];
  glyphFilter->GetPoint(pointCoords);

  typedef double projectionType;

//  std::cout << "Calling CalcGlyph for point "
//             << glyphFilter->GetPointId() << std::endl;
//  std::cout << "Point coords are: "
//            << pointCoords[0] << " "
//            << pointCoords[1] << " "
//            << pointCoords[2] << std::endl;

  vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
      sphereSource->SetThetaResolution( glyphFilter->GetResolution() );
      sphereSource->SetPhiResolution( glyphFilter->GetResolution() );
      sphereSource->SetRadius( 1.0 );
      sphereSource->Update();
  vtkSmartPointer<vtkPolyData> glyph = sphereSource->GetOutput();

  vtkSmartPointer<vtkFloatArray> projections = vtkSmartPointer<vtkFloatArray>::New();
      projections->SetNumberOfComponents(3);
      projections->SetNumberOfTuples(glyph->GetNumberOfPoints());
      projections->SetName("Fibre Projections");
      projections->FillComponent(0, 0.0);
      projections->FillComponent(1, 0.0);
      projections->FillComponent(2, 0.0);

  ///For every point in sphere mesh
  double m_x = pointCoords[0], m_y = pointCoords[1], m_z = pointCoords[2];
  for(vtkIdType i = 0; i < glyph->GetNumberOfPoints(); i ++)
  {
      double point_sphere[3];
      glyph->GetPoint(i, point_sphere);
      const double x_s = point_sphere[0];
      const double y_s = point_sphere[1];
      const double z_s = point_sphere[2];

      ///Get the SH itself
      std::vector<double> sHarmonicCoefficients;
      for(size_t j = 0; j < harmonicsVector.GetSize(); j ++)
          sHarmonicCoefficients.push_back(harmonicsVector[j]/normVector[j]);

      ///Compute spherical harmonic amplitude
      const int n_coefs = sHarmonicCoefficients.size();
      const int l_max   = vtkDiffusionTensorGlyphFilter::LforN( n_coefs );
      double amplitude = vtkDiffusionTensorGlyphFilter::computeAmplitude( sHarmonicCoefficients, x_s, y_s, z_s, l_max );
//      cout << "Amplitude: " << amplitude << ", ";

      if(amplitude < 0)
          amplitude = 0;

      ///use SH to shape sphere points accordingly
      double x_g = x_s * amplitude * radiusScaling;
      double y_g = y_s * amplitude * radiusScaling;
      double z_g = z_s * amplitude * radiusScaling;

      glyph->GetPoints()->SetPoint(i, x_g, y_g, z_g);

      ///Compute color based on direction
      for(size_t j = 0; j < 3; j ++)
      {
          projectionType axis[3] = {0.0, 0.0, 0.0};
          coordinate currentProjection(projections->GetTuple3(i)), position(glyph->GetPoint(i)); //coordinate is vnl_vector_fixed

          axis[j] = 1.0;
          projectionType projection = vtkMath::Dot(axis, position.data_block()); //project to axis being done,
          currentProjection[j] = projection; //projection in each direction

          projections->SetTuple3(i, static_cast<float>(currentProjection[0])
                                  , static_cast<float>(currentProjection[1])
                                  , static_cast<float>(currentProjection[2]));
      }

      ///use point in input to displace sphere points accordingly
      x_g += m_x;
      y_g += m_y;
      z_g += m_z;
      glyph->GetPoints()->SetPoint(i, x_g, y_g, z_g);

      qApp->processEvents(); //keep UI responsive
  }

  ///Colour based on each of the gradient values in that direction
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
      scalars->SetNumberOfComponents(3);
      scalars->SetNumberOfTuples(glyph->GetNumberOfPoints());
      scalars->SetName("Colours");
      scalars->FillComponent(0, 0);
      scalars->FillComponent(1, 0);
      scalars->FillComponent(2, 0);
  for(vtkIdType i = 0; i < glyph->GetNumberOfPoints(); i ++)
  {
      coordinate currentProjection(projections->GetTuple3(i));
      coordinate currentProjSquared = element_product(currentProjection, currentProjection);
      projectionType maxProjection = currentProjSquared.max_value();
      currentProjSquared /= maxProjection;

      unsigned char colourOfPoint[3] = {0, 0, 0};
      colourOfPoint[0] = static_cast<unsigned char>( currentProjSquared[0]*255.0 );
      colourOfPoint[1] = static_cast<unsigned char>( currentProjSquared[1]*255.0 );
      colourOfPoint[2] = static_cast<unsigned char>( currentProjSquared[2]*255.0 );

      scalars->SetTupleValue(i, colourOfPoint);

      qApp->processEvents(); //keep UI responsive
  }
  //  cout << endl;
  glyph->GetPointData()->SetScalars(scalars);
  glyph->Modified();

  vtkSmartPointer< vtkPolyDataNormals > normals = vtkSmartPointer< vtkPolyDataNormals >::New();
  #if VTK_MAJOR_VERSION <= 5
    normals->SetInput(glyph);
  #else
    normals->SetInputData(glyph);
  #endif
    normals->ComputePointNormalsOn();
    normals->SplittingOff();
    normals->ComputeCellNormalsOff();
    normals->Update();

#if VTK_MAJOR_VERSION <= 5
  glyphFilter->SetSource(normals->GetOutput());
#else
  glyphFilter->SetSourceConnection(normals->GetOutputPort());
#endif
}

milxQtDiffusionTensorImage::milxQtDiffusionTensorImage(QWidget *theParent) : milxQtImage(theParent)
{
    milxQtWindow::prefix = "DTI Img: ";

    createActions();

    createConnections();
}

milxQtDiffusionTensorImage::~milxQtDiffusionTensorImage()
{
    //dtor
}

void milxQtDiffusionTensorImage::diffusionGlyphs(size_t subsampleFactor, size_t resolution)
{
    if(ITK_VERSION_MAJOR < 4 || VTK_MAJOR_VERSION < 6)
    {
        printError("You are not using ITK 4 or above OR VTK 6 or above. Ignoring.");
        return;
    }

    printDebug("Computing Diffusion Glyphs");

    bool ok1 = false, ok2 = false;
    if(subsampleFactor == 0)
    {
        subsampleFactor = QInputDialog::getInt(this, tr("Please Provide the sub-sample factor of the data"),
                                              tr("Sub-sample Factor:"), 1, 1, 1000, 1, &ok1);
        resolution = QInputDialog::getInt(this, tr("Please Provide the resolution of the glyphs"),
                                              tr("Resolution:"), 16, 4, 256, 1, &ok2);
        if(!ok1 || !ok2)
            return;
    }

    int extent[6];
    milxQtImage::viewer->GetImageActor()->GetDisplayExtent(extent);

    if(flipped)
    {
        int actualExtent[6];
        milxQtImage::imageData->GetExtent(actualExtent);

        if(extent[3]-extent[2] == 0) //flip y extent
        {
            extent[2] = actualExtent[3]-extent[2];
            extent[3] = actualExtent[3]-extent[3];
        }
    }

    emit working(-1);
    ///Extract the viewed slice
    vectorImageType::Pointer sliceVector = milx::Image<vectorImageType>::ExtractSlice<vectorImageType>(milxQtImage::imageVector, extent);//!< assumes float vector image here

    //Which dimension is sliced?
    vectorImageType::SizeType imageSize = sliceVector->GetLargestPossibleRegion().GetSize(); ///assumes float image here
    size_t sliceDimension = 0;
    for(size_t j = 0; j < vectorImageType::ImageDimension; j ++)
    {
      if(imageSize[j] == 1)
        sliceDimension = j;
    }

    ///Check mag image if its a 2D slice
    vectorImageType::SizeType sliceSubsampleSizes;
        sliceSubsampleSizes.Fill(subsampleFactor);
        sliceSubsampleSizes[sliceDimension] = 1;

    ///Subsample image to reduce computation costs
    const size_t components = milxQtImage::imageVector->GetNumberOfComponentsPerPixel();
    vectorImageType::Pointer imgSubSampled = milx::Image<vectorImageType>::SubsampleImage(sliceVector, sliceSubsampleSizes);
    cout << "Slice Information: " << endl;
    milx::Image<vectorImageType>::Information(sliceVector);

    ///Create Polydata to hold vector field
    printDebug("Setting up seed points for glyphs using the current slice");
    vtkSmartPointer<vtkPoints> slicePoints = vtkSmartPointer<vtkPoints>::New();
    typedef itk::Point<double, 3> InputImagePointType;
    itk::ImageRegionConstIteratorWithIndex<vectorImageType> imageIterator(imgSubSampled, imgSubSampled->GetLargestPossibleRegion());

    ///Store indices
    std::vector<vectorImageType::IndexType> indices;

    ///Setup seed points
    InputImagePointType point;
    vectorImageType::PixelType maxVector(components);
    maxVector.Fill(0.0);
    while(!imageIterator.IsAtEnd())
    {
        double position[3];

        imgSubSampled->TransformIndexToPhysicalPoint(imageIterator.GetIndex(), point);

        position[0] = point[0];
        position[1] = point[1];
        position[2] = point[2];

        slicePoints->InsertNextPoint(position);
        indices.push_back(imageIterator.GetIndex());

        vectorImageType::PixelType pixelVector = imageIterator.Get();
        for(size_t j = 0; j < pixelVector.GetSize(); j ++)
        {
            if(maxVector[j] < pixelVector[j])
                maxVector[j] = pixelVector[j];
        }

        ++imageIterator;
    }
    cout << "Max Vector found: " << maxVector << endl;
    cout << "ID Type Size: " << sizeof(vtkIdType) << endl;
    cout << "Integer Type Size: " << sizeof(unsigned) << endl;

//    double magNormVector = maxVector.GetNorm();
//    double magNormVector = maxVector[0];
//    for(size_t j = 0; j < maxVector.GetSize(); j ++)
//        maxVector[j] /= magNormVector;
//    cout << "Max Vector normalised: " << maxVector << endl;

    vtkSmartPointer<vtkPolyData> slice = vtkSmartPointer<vtkPolyData>::New();
        slice->SetPoints(slicePoints);

    ///Although the following sphere is unused, its info is used to create the total output so don't remove
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetThetaResolution( resolution );
        sphereSource->SetPhiResolution( resolution );
        sphereSource->SetRadius( 1.0 );
        sphereSource->Update();
    vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
        scalars->SetNumberOfComponents(3);
        scalars->SetNumberOfTuples(sphereSource->GetOutput()->GetNumberOfPoints());
        scalars->SetName("Dummy Colours");
        scalars->FillComponent(0, 0);
        scalars->FillComponent(1, 0);
        scalars->FillComponent(2, 0);
    sphereSource->GetOutput()->GetPointData()->SetScalars(scalars); //dummy allocation to get uchar array type into output of glyph filter.
    vtkSmartPointer< vtkPolyDataNormals > normals = vtkSmartPointer< vtkPolyDataNormals >::New();
      #if VTK_MAJOR_VERSION <= 5
        normals->SetInput(sphereSource->GetOutput());
      #else
        normals->SetInputData(sphereSource->GetOutput());
      #endif
        normals->ComputePointNormalsOn();
        normals->SplittingOff();
        normals->ComputeCellNormalsOff();
        normals->Update();

    vtkSmartPointer<vtkDiffusionTensorGlyphFilter> diffusionTensorGlyphs = vtkSmartPointer<vtkDiffusionTensorGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        diffusionTensorGlyphs->SetInput(slice);
        diffusionTensorGlyphs->SetSource(normals->GetOutput());
    #else
        diffusionTensorGlyphs->SetInputData(slice);
        diffusionTensorGlyphs->SetSourceData(normals->GetOutput());
    #endif
        diffusionTensorGlyphs->SetTensorImage(sliceVector);
        diffusionTensorGlyphs->SetTensorImageIndices(indices);
        diffusionTensorGlyphs->SetNormaliseVector(maxVector);
        diffusionTensorGlyphs->SetResolution(resolution);
        diffusionTensorGlyphs->SetGlyphMethod(CalcGlyph, diffusionTensorGlyphs);
        linkProgressEventOf(diffusionTensorGlyphs);
        printDebug("Creating Glyphs");
        diffusionTensorGlyphs->SetColorModeToColorBySource();
        diffusionTensorGlyphs->Update();

    QPointer<milxQtModel> model = new milxQtModel;
        model->setName("Diffusion Glyphs");
        model->SetInput(diffusionTensorGlyphs->GetOutput());
        model->generateModel();
        model->ImmediateModeRenderingOn(); //for large datasets

    emit done(-1);
    emit resultAvailable(model);
}

void milxQtDiffusionTensorImage::diffusionGlyphs2(size_t subsampleFactor, size_t resolution)
{
    if(ITK_VERSION_MAJOR < 4)
    {
        printError("You are not using ITK 4 or above. Ignoring.");
        return;
    }

    printDebug("Computing Diffusion Glyphs via Brute Force");

    bool ok1 = false, ok2 = false;
    if(subsampleFactor == 0)
    {
        subsampleFactor = QInputDialog::getInt(this, tr("Please Provide the sub-sample factor of the data"),
                                              tr("Sub-sample Factor:"), 1, 1, 1000, 1, &ok1);
        resolution = QInputDialog::getInt(this, tr("Please Provide the resolution of the glyphs"),
                                              tr("Resolution:"), 16, 4, 256, 1, &ok2);
        if(!ok1 || !ok2)
            return;
    }

    int extent[6];
    milxQtImage::viewer->GetImageActor()->GetDisplayExtent(extent);

    if(flipped)
    {
        int actualExtent[6];
        milxQtImage::imageData->GetExtent(actualExtent);

        if(extent[3]-extent[2] == 0) //flip y extent
        {
            extent[2] = actualExtent[3]-extent[2];
            extent[3] = actualExtent[3]-extent[3];
        }
    }

    emit working(-1);
    ///Extract the viewed slice
    vectorImageType::Pointer sliceVector = milx::Image<vectorImageType>::ExtractSlice<vectorImageType>(milxQtImage::imageVector, extent);//!< assumes float vector image here

    //Which dimension is sliced?
    vectorImageType::SizeType imageSize = sliceVector->GetLargestPossibleRegion().GetSize(); ///assumes float image here
    size_t sliceDimension = 0;
    for(size_t j = 0; j < vectorImageType::ImageDimension; j ++)
    {
      if(imageSize[j] == 1)
        sliceDimension = j;
    }

    ///Check mag image if its a 2D slice
    vectorImageType::SizeType sliceSubsampleSizes;
        sliceSubsampleSizes.Fill(subsampleFactor);
        sliceSubsampleSizes[sliceDimension] = 1;

    ///Subsample image to reduce computation costs
    const size_t components = milxQtImage::imageVector->GetNumberOfComponentsPerPixel();
    vectorImageType::Pointer imgSubSampled = milx::Image<vectorImageType>::SubsampleImage(sliceVector, sliceSubsampleSizes);
    cout << "Slice Information: " << endl;
    milx::Image<vectorImageType>::Information(sliceVector);
    typedef float projectionType;

    ///Create Polydata to hold vector field
    printDebug("Setting up seed points for glyphs using the current slice");
    typedef itk::Point<double, 3> InputImagePointType;
    itk::ImageRegionConstIteratorWithIndex<vectorImageType> imageIterator(imgSubSampled, imgSubSampled->GetLargestPossibleRegion());

    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetThetaResolution( resolution );
        sphereSource->SetPhiResolution( resolution );
        sphereSource->SetRadius( 1.0 );
        sphereSource->Update();
        sphereSource->GetOutput()->GetPointData()->RemoveArray("Normals");

    ///Setup glyphs per voxel
    vtkSmartPointer<vtkAppendPolyData> appendedPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
    InputImagePointType point;
    const double radiusScaling = 2.5;
//    vectorImageType::PixelType maxVector(components);
//    maxVector.Fill(0.0);
    while(!imageIterator.IsAtEnd())
    {
        imgSubSampled->TransformIndexToPhysicalPoint(imageIterator.GetIndex(), point);

        const vectorImageType::PixelType harmonicsVector = imageIterator.Get();

//        vectorImageType::PixelType pixelVector = imageIterator.Get();
//        for(size_t j = 0; j < pixelVector.GetSize(); j ++)
//        {
//            if(maxVector[j] < pixelVector[j])
//                maxVector[j] = pixelVector[j];
//        }

        vtkSmartPointer<vtkPolyData> glyph = vtkSmartPointer<vtkPolyData>::New();
            glyph->DeepCopy(sphereSource->GetOutput());

        vtkSmartPointer<vtkFloatArray> projections = vtkSmartPointer<vtkFloatArray>::New();
            projections->SetNumberOfComponents(3);
            projections->SetNumberOfTuples(glyph->GetNumberOfPoints());
            projections->SetName("Fibre Projections");
            projections->FillComponent(0, 0.0);
            projections->FillComponent(1, 0.0);
            projections->FillComponent(2, 0.0);

        ///For every point in sphere mesh
        double m_x = point[0], m_y = point[1], m_z = point[2];
        for(vtkIdType i = 0; i < glyph->GetNumberOfPoints(); i ++)
        {
            double point_sphere[3];
            glyph->GetPoint(i, point_sphere);
            const double x_s = point_sphere[0];
            const double y_s = point_sphere[1];
            const double z_s = point_sphere[2];

            ///Get the SH itself
            std::vector<double> sHarmonicCoefficients;
            for(size_t j = 0; j < harmonicsVector.GetSize(); j ++)
                sHarmonicCoefficients.push_back(harmonicsVector[j]);

            ///Compute spherical harmonic amplitude
            const int n_coefs = sHarmonicCoefficients.size();
            const int l_max   = vtkDiffusionTensorGlyphFilter::LforN( n_coefs );
            double amplitude = vtkDiffusionTensorGlyphFilter::computeAmplitude( sHarmonicCoefficients, x_s, y_s, z_s, l_max );
      //      cout << "Amplitude: " << amplitude << ", ";

            if(amplitude < 0)
                amplitude = 0;

            ///use SH to shape sphere points accordingly
            double x_g = x_s * amplitude * radiusScaling;
            double y_g = y_s * amplitude * radiusScaling;
            double z_g = z_s * amplitude * radiusScaling;

            glyph->GetPoints()->SetPoint(i, x_g, y_g, z_g);

            ///Compute color based on direction
            for(size_t j = 0; j < 3; j ++)
            {
                projectionType axis[3] = {0.0, 0.0, 0.0};
                const projectionType position[3] = {glyph->GetPoint(i)[0], glyph->GetPoint(i)[1], glyph->GetPoint(i)[2]};
                projectionType currentProjection[3] = {projections->GetTuple3(i)[0], projections->GetTuple3(i)[1], projections->GetTuple3(i)[2]}; //coordinate is vnl_vector_fixed

                axis[j] = 1.0;
                const projectionType projection = vtkMath::Dot(axis, position); //project to axis being done,
                currentProjection[j] = projection; //projection in each direction

                projections->SetTuple3(i, currentProjection[0], currentProjection[1], currentProjection[2]);
            }

            ///use point in input to displace sphere points accordingly
            x_g += m_x;
            y_g += m_y;
            z_g += m_z;
            glyph->GetPoints()->SetPoint(i, x_g, y_g, z_g);

            qApp->processEvents(); //keep UI responsive
        }

        ///Colour based on each of the gradient values in that direction
        vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
            scalars->SetNumberOfComponents(3);
            scalars->SetNumberOfTuples(glyph->GetNumberOfPoints());
            scalars->SetName("Colours");
            scalars->FillComponent(0, 0);
            scalars->FillComponent(1, 0);
            scalars->FillComponent(2, 0);

        for(vtkIdType i = 0; i < glyph->GetNumberOfPoints(); i ++)
        {
            coordinate currentProjection(projections->GetTuple3(i));
            coordinate currentProjSquared = element_product(currentProjection, currentProjection);
            projectionType maxProjection = currentProjSquared.max_value();
            currentProjSquared /= maxProjection;

            unsigned char colourOfPoint[3] = {0, 0, 0};
            colourOfPoint[0] = static_cast<unsigned char>( currentProjSquared[0]*255.0 );
            colourOfPoint[1] = static_cast<unsigned char>( currentProjSquared[1]*255.0 );
            colourOfPoint[2] = static_cast<unsigned char>( currentProjSquared[2]*255.0 );

            scalars->SetTupleValue(i, colourOfPoint);

            qApp->processEvents(); //keep UI responsive
        }
        //  cout << endl;
        glyph->GetPointData()->SetScalars(scalars);

    #if VTK_MAJOR_VERSION <= 5
        appendedPolyData->AddInput(glyph); //append to model
    #else
        appendedPolyData->AddInputData(glyph); //append to model
    #endif

        ++imageIterator;
    }
    cout << "Updating Polydata" << endl;
    appendedPolyData->Update();

    QPointer<milxQtModel> model = new milxQtModel;
        model->SetInput(appendedPolyData->GetOutput());
        model->setName("Diffusion Glyphs");
        model->generateModel();
        model->ImmediateModeRenderingOn(); //for large datasets

    emit done(-1);
    emit resultAvailable(model);
}

void milxQtDiffusionTensorImage::createActions()
{
    diffusionAct = new QAction(this);
        diffusionAct->setText(QApplication::translate("Model", "Diffusion Glyphs for Slice", 0, QApplication::UnicodeUTF8));
        diffusionAct->setShortcut(tr("Shift+Alt+d"));
    diffusion2Act = new QAction(this);
        diffusion2Act->setText(QApplication::translate("Model", "Diffusion Glyphs for Slice Brute Force", 0, QApplication::UnicodeUTF8));
        diffusion2Act->setShortcut(tr("Ctrl+Shift+d"));
}

void milxQtDiffusionTensorImage::createConnections()
{
    //Operations
    connect(diffusionAct, SIGNAL(triggered()), this, SLOT(diffusionGlyphs()));
    connect(diffusion2Act, SIGNAL(triggered()), this, SLOT(diffusionGlyphs2()));

    //milxQtImage::createConnections();
}

void milxQtDiffusionTensorImage::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    contextMenu = milxQtImage::basicContextMenu(); //!< Only exists for the duration of the context selection

    contextMenu->addSeparator()->setText(tr("Diffusion Tensor"));
    contextMenu->addAction(diffusionAct);
    contextMenu->addAction(diffusion2Act);
    contextMenu->addSeparator()->setText(tr("Extensions"));
    foreach(QAction *currAct, extActionsToAdd)
    {
        contextMenu->addAction(currAct);
    }
    contextMenu->addSeparator();
    ///Dont display extensions
    contextMenu->addMenu(milxQtRenderWindow::contourMenu);
    contextMenu->addMenu(milxQtRenderWindow::windowPropertiesMenu);
    contextMenu->addAction(milxQtRenderWindow::refreshAct);
    contextMenu->addAction(milxQtRenderWindow::resetAct);

    contextMenu->exec(currentEvent->globalPos());
}
