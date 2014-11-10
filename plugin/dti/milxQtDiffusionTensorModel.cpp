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
#include "milxQtDiffusionTensorModel.h"

#include <time.h>

#include <vtkMath.h>
#include <vtkLine.h>
#include <vtkSphereSource.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

milxQtDiffusionTensorModel::milxQtDiffusionTensorModel(QWidget *theParent) : milxQtModel(theParent)
{
    milxQtWindow::prefix = "DTI: ";

    createActions();

    createConnections();
}

milxQtDiffusionTensorModel::~milxQtDiffusionTensorModel()
{
    //dtor
}

void milxQtDiffusionTensorModel::colourByDirection()
{
    vtkSmartPointer<vtkPolyData> currentMesh = model.Result();

    typedef double projectionType;

    ///Determine colours based on axes directions
    vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
        scalars->SetNumberOfComponents(3);
        scalars->SetNumberOfTuples(currentMesh->GetNumberOfPoints());
        scalars->SetName("Fibre Colours");
        scalars->FillComponent(0, 0);
        scalars->FillComponent(1, 0);
        scalars->FillComponent(2, 0);
    vtkSmartPointer<vtkFloatArray> projections = vtkSmartPointer<vtkFloatArray>::New();
        projections->SetNumberOfComponents(3);
        projections->SetNumberOfTuples(currentMesh->GetNumberOfPoints());
        projections->SetName("Fibre Projections");
        projections->FillComponent(0, 0.0);
        projections->FillComponent(1, 0.0);
        projections->FillComponent(2, 0.0);

    emit working(-1);
    if(currentMesh->GetNumberOfLines() == 0)
    {
        printInfo("No lines in model. Computing colours for mesh.");

        ///Use the dot product in each axis
        for(size_t j = 0; j < 3; j ++)
        {
            projectionType axis[3] = {0.0, 0.0, 0.0};
            axis[j] = 1.0;
            cout << "Computing toward axis " << j << endl;

            ///Colour based on each of the gradient values in that direction
            for(vtkIdType k = 0; k < currentMesh->GetNumberOfPoints(); k ++)
            {
                coordinate currentProjection(projections->GetTuple3(k)), position(currentMesh->GetPoint(k));

                projectionType projection = vtkMath::Dot(axis, position.data_block()); //project to axis being done,
                currentProjection[j] = projection; //projection in each direction

                projections->SetTuple3(k, currentProjection[0], currentProjection[1], currentProjection[2]);
            }
        }

        ///Colour based on each of the gradient values in that direction
        for(vtkIdType k = 0; k < currentMesh->GetNumberOfPoints(); k ++)
        {
            coordinate currentProjection(projections->GetTuple3(k));
            coordinate currentProjSquared = element_product(currentProjection, currentProjection);
            projectionType maxProjection = currentProjSquared.max_value();
            currentProjSquared /= maxProjection;

            unsigned char colourOfPoint[3] = {0, 0, 0};
            colourOfPoint[0] = static_cast<unsigned char>( currentProjSquared[0]*255.0 );
            colourOfPoint[1] = static_cast<unsigned char>( currentProjSquared[1]*255.0 );
            colourOfPoint[2] = static_cast<unsigned char>( currentProjSquared[2]*255.0 );

            scalars->SetTupleValue(k, colourOfPoint);
        }

        currentMesh->GetPointData()->SetVectors(projections);
        currentMesh->GetPointData()->SetScalars(scalars);
    }
    else
    {
        printInfo("Re-colouring lines by axes directions");
        //Ensure lines done properly so can colour appropriately
        currentMesh->GetLines()->InitTraversal();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
        std::vector<vtkIdType> trackLengths;
        for(vtkIdType j = 0; j < currentMesh->GetNumberOfLines(); j ++)
        {
            currentMesh->GetLines()->GetNextCell(idList);
            trackLengths.push_back(idList->GetNumberOfIds());
            for(vtkIdType pointId = 0; pointId < idList->GetNumberOfIds(); pointId ++)
            {
    //            std::cout << idList->GetId(pointId) << " ";

                double position[3];
                currentMesh->GetPoint(idList->GetId(pointId), position);
                points->InsertNextPoint(position);
            }
        }

        //Re-stitch lines together
        vtkIdType step = 0;
        vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
        for(size_t j = 0; j < trackLengths.size(); j ++)
        {
            for(vtkIdType k = step; k < trackLengths[j]-1; k ++)
            {
                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    line->GetPointIds()->SetId(0, k);
                    line->GetPointIds()->SetId(1, k+1);
                lines->InsertNextCell(line);
            }
            qDebug() << trackLengths[j] << endl;
            step += trackLengths[j];
        }

        vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
            //Add the points to the dataset
            linesPolyData->SetPoints(points);
            //Add the lines to the dataset
            linesPolyData->SetLines(lines);
            linesPolyData->Modified();
        model.SetInput(linesPolyData);
    }

    emit done(-1);
    generateModel();
}

void milxQtDiffusionTensorModel::harmonics(QString ordersString)
{
    bool ok = false;
    if(ordersString.isEmpty())
    {
        ordersString = QInputDialog::getText(this, tr("Enter order magnitudes of harmonics separated by spaces"),
                                            tr("Orders: "), QLineEdit::Normal,
                                            "1 0 0 0 1 1 1 1 1", &ok);
    }
    if (!ok || ordersString.isEmpty())
        return;

    QStringList orders = ordersString.split(" ");
    std::vector<double> sHarmonicCoefficients;
    for(int i = 0; i < orders.size(); i++)
        sHarmonicCoefficients.push_back(orders[i].toDouble());

    int n_coefs = sHarmonicCoefficients.size();
    int l_max   = milxQtDiffusionTensorModel::LforN( n_coefs );

    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetThetaResolution( 64 );
        sphereSource->SetPhiResolution( 64 );
        sphereSource->SetRadius( 1.0 );
        sphereSource->Update();
    vtkSmartPointer<vtkPolyData> glyph = sphereSource->GetOutput();

    ///For every point in sphere mesh
    double m_x = 0, m_y = 0, m_z = 0;
    for(int i = 0; i < glyph->GetNumberOfPoints(); i++)
    {
        double point_sphere[3];
        glyph->GetPoint(i, point_sphere);
        double x_s = point_sphere[0];
        double y_s = point_sphere[1];
        double z_s = point_sphere[2];

        double x_g, y_g, z_g;

        ///Compute spherical harmonic amplitude
        double amplitude = computeAmplitude( sHarmonicCoefficients, x_s, y_s, z_s, l_max );

        if(amplitude < 0)
            amplitude = 0;

        ///use this to displace sphere points accordingly
        x_g = x_s * amplitude + m_x;
        y_g = y_s * amplitude + m_y;
        z_g = z_s * amplitude + m_z;

        glyph->GetPoints()->SetPoint(i, x_g, y_g, z_g);
    }

    model.SetInput(glyph);
    model.GenerateNormals(true);

    generateModel();
    milxQtRenderWindow::reset();
}

double milxQtDiffusionTensorModel::computeAmplitude(std::vector<double> SH, double x, double y, double z, int lmax)
{
    double az = atan2(y, x);
    double el = acos(z);
    double val = 0.0;
    for (int l = 0; l <= lmax; l+=2)
    {
        // val += SH[milxQtDiffusionTensorModel::index(l,0)] * boost::math::legendre_p<double>(l, 0, z);
        val += SH[milxQtDiffusionTensorModel::index(l,0)] * boost::math::spherical_harmonic_r<double>(l, 0, el,az);
    }

    for (int m = 1; m <= lmax; m++)
    {
//        float caz = cos(m*az);
//        float saz = sin(m*az);
        for (int l = 2*((m+1)/2); l <= lmax; l+=2)
        {
        // float buf = boost::math::legendre_p<double>(l, 0, z);
        // val += SH[index(l,m)]*buf*caz;
        // val += SH[index(l,-m)]*buf*saz;
            val += SH[milxQtDiffusionTensorModel::index(l,m)]*boost::math::spherical_harmonic_r(l,m,el,az);
            val += SH[milxQtDiffusionTensorModel::index(l,-m)]*boost::math::spherical_harmonic_i(l,m,el,az);
        }
    }

    return val;
}

void milxQtDiffusionTensorModel::createActions()
{
    colourDirectionAct = new QAction(this);
        colourDirectionAct->setText(QApplication::translate("Model", "Colour By Direction", 0, QApplication::UnicodeUTF8));
        colourDirectionAct->setShortcut(tr("Alt+c"));
    harmonicsAct = new QAction(this);
        harmonicsAct->setText(QApplication::translate("Model", "Show Spherical Harmonic ...", 0, QApplication::UnicodeUTF8));
        harmonicsAct->setShortcut(tr("Alt+s"));
}

void milxQtDiffusionTensorModel::createConnections()
{
    //Operations
    connect(colourDirectionAct, SIGNAL(triggered()), this, SLOT(colourByDirection()));
    connect(harmonicsAct, SIGNAL(triggered()), this, SLOT(harmonics()));

    //milxQtModel::createConnections();
}

void milxQtDiffusionTensorModel::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    contextMenu = milxQtModel::basicContextMenu(); //!< Only exists for the duration of the context selection

    contextMenu->addSeparator()->setText(tr("DiffusionTensor"));
    contextMenu->addAction(colourDirectionAct);
    contextMenu->addAction(harmonicsAct);
    contextMenu->addSeparator()->setText(tr("Extensions"));
    foreach(QAction *currAct, extActionsToAdd)
    {
        contextMenu->addAction(currAct);
    }
    contextMenu->addSeparator();
    ///Dont display extensions
    contextMenu->addAction(milxQtModel::scaleAct);
    contextMenu->addAction(milxQtRenderWindow::axesAct);
    contextMenu->addAction(milxQtRenderWindow::refreshAct);

    contextMenu->exec(currentEvent->globalPos());
}
