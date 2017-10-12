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
#include "milxQtPreferencesForm.h"

#include <QCheckBox>
#include <QLineEdit>
#include <QLabel>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QDesktopWidget>

milxQtPreferencesForm::milxQtPreferencesForm(milxQtMain *theParent) : QDialog(theParent)
{
    ui.setupUi(this);

    setWindowModality(Qt::ApplicationModal); //block user input
    setWindowTitle(tr("sMILX Preferences"));
	setFixedSize(this->size());

    MainWindow = theParent;

    setupPages();
    createConnections();
}

milxQtPreferencesForm::~milxQtPreferencesForm()
{
    //dtor
}

void milxQtPreferencesForm::setupPages()
{
    ///General
    QWidget *generalPage = new QWidget;
    QListWidgetItem *generalPageItem = new QListWidgetItem(ui.wdtOptions);
        generalPageItem->setText(tr("General"));
        generalPageItem->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    ///Background colour (check)
    backgroundCheckBox = new QCheckBox(tr("Prefer white background"));
        backgroundCheckBox->setChecked(MainWindow->isWhiteBackgroundPreferred());
    ///View (combo)
    ///Orientation standard (combo)
    ///Human glyph (check)
    humanCheckBox = new QCheckBox(tr("Show human orientation glyph"));
        humanCheckBox->setChecked(MainWindow->isHumanGlyphPreferred());
    ///Window size (text)
    QLabel *windowSizeLabel = new QLabel(tr("Preferred Window Size:"));
    QSize desktopSize = qApp->desktop()->availableGeometry().size();
    windowSizeEdit = new QSpinBox;
        windowSizeEdit->setMinimum(minWindowSize);
        windowSizeEdit->setMaximum(milx::Maximum<int>(desktopSize.width(), desktopSize.height()));
        windowSizeEdit->setValue(MainWindow->hasPreferredSubWindowSize());
    QHBoxLayout *windowSizeLayout = new QHBoxLayout;
        windowSizeLayout->addWidget(windowSizeLabel);
        windowSizeLayout->addWidget(windowSizeEdit);
    ///Number of processors (text)
    QLabel *processorsLabel = new QLabel(tr("Preferred Max. Processors:"));
    processorsEdit = new QSpinBox;
        processorsEdit->setMinimum(1);
        processorsEdit->setMaximum(milx::NumberOfProcessors());
        processorsEdit->setValue(MainWindow->hasMaximumProcessors());
    QHBoxLayout *processorsLayout = new QHBoxLayout;
        processorsLayout->addWidget(processorsLabel);
        processorsLayout->addWidget(processorsEdit);
    ///Magnification of Screenshots
    QLabel *magnifyLabel = new QLabel(tr("Screenshot Magnification Factor"));
    magnifyEdit = new QSpinBox;
        magnifyEdit->setMinimum(1);
        magnifyEdit->setMaximum(1024);
        magnifyEdit->setValue(MainWindow->hasScreenshotMagnifyFactor());
    QHBoxLayout *magnifyLayout = new QHBoxLayout;
        magnifyLayout->addWidget(magnifyLabel);
        magnifyLayout->addWidget(magnifyEdit);
    ///Timestamp in log? (check)
    timestampCheckBox = new QCheckBox(tr("Show timestamp in logs"));
        timestampCheckBox->setChecked(MainWindow->isTimestampsPreferred());
    //View Layout
    QVBoxLayout *viewLayout = new QVBoxLayout;
        viewLayout->addWidget(backgroundCheckBox);
        viewLayout->addWidget(humanCheckBox);
    QGroupBox *generalViewGroup = new QGroupBox(tr("View Options"));
        generalViewGroup->setLayout(viewLayout);
    //General layout
    QVBoxLayout *generalLayout = new QVBoxLayout;
        generalLayout->addLayout(windowSizeLayout);
        generalLayout->addLayout(processorsLayout);
        generalLayout->addLayout(magnifyLayout);
        generalLayout->addWidget(timestampCheckBox);
    QGroupBox *generalGroup = new QGroupBox(tr("General"));
        generalGroup->setLayout(generalLayout);
    QVBoxLayout *generalPageLayout = new QVBoxLayout;
        generalPageLayout->addWidget(generalViewGroup);
        generalPageLayout->addWidget(generalGroup);
        generalPage->setLayout(generalPageLayout);
    ui.wdtPages->insertWidget(0, generalPage);
    ui.wdtPages->setCurrentWidget(generalPage);
    ui.wdtOptions->setCurrentItem(generalPageItem);

    ///Imaging
    QWidget *imagingPage = new QWidget;
    QListWidgetItem *imagingPageItem = new QListWidgetItem(ui.wdtOptions);
        imagingPageItem->setText(tr("Images"));
        imagingPageItem->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    ///Colourmap (combo)
    ///Interpolation (check)
    interpolationCheckBox = new QCheckBox(tr("Interpolate images"));
        interpolationCheckBox->setChecked(MainWindow->isImageInterpolationPreferred());
    ///Apply Orientation (check)
    orientationCheckBox = new QCheckBox(tr("Apply orientation to images"));
        orientationCheckBox->setChecked(MainWindow->isOrientationPreferred());
    //Layout
    QVBoxLayout *imagingPageLayout = new QVBoxLayout;
        imagingPageLayout->addWidget(interpolationCheckBox);
        imagingPageLayout->addWidget(orientationCheckBox);
        imagingPage->setLayout(imagingPageLayout);
    ui.wdtPages->insertWidget(1, imagingPage);

    ///Models
    QWidget *modelPage = new QWidget;
    QListWidgetItem *modelPageItem = new QListWidgetItem(ui.wdtOptions);
        modelPageItem->setText(tr("Models"));
        modelPageItem->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    ///Colourmap (combo)
    ///Interpolation (check)
    interpolationModelCheckBox = new QCheckBox(tr("Interpolate (Phong Shading) models"));
        interpolationModelCheckBox->setChecked(MainWindow->isModelInterpolationPreferred());
    ///Scalar bar
    scalarBarCheckBox = new QCheckBox(tr("Always show scalar bar"));
        scalarBarCheckBox->setChecked(MainWindow->isScalarBarPreferred());
    //Layout
    QVBoxLayout *modelPageLayout = new QVBoxLayout;
        modelPageLayout->addWidget(interpolationModelCheckBox);
        modelPageLayout->addWidget(scalarBarCheckBox);
        modelPage->setLayout(modelPageLayout);
    ui.wdtPages->insertWidget(2, modelPage);

    ///Plugins
    QWidget *pluginsPage = new QWidget;
    QListWidgetItem *pluginsPageItem = new QListWidgetItem(ui.wdtOptions);
        pluginsPageItem->setText(tr("Plugins"));
        pluginsPageItem->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    ///List and Disable
    QListWidget *pluginsList = new QListWidget;
        pluginsList->setDisabled(true);
    foreach(QPointer<milxQtPluginInterface> plugin, MainWindow->getPlugins())
    {
        QListWidgetItem *pluginsItem = new QListWidgetItem(pluginsList);
            pluginsItem->setText(plugin->name());
            pluginsItem->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled );
    }
    //Layout
    QVBoxLayout *pluginsPageLayout = new QVBoxLayout;
        pluginsPageLayout->addWidget(pluginsList);
        pluginsPage->setLayout(pluginsPageLayout);

    ui.wdtPages->insertWidget(3, pluginsPage);
}

void milxQtPreferencesForm::changePage(QListWidgetItem *current, QListWidgetItem *previous)
{
    if (!current)
        current = previous;

//    cout << "Index of Page: " << ui.wdtOptions->row(current) << endl;
    ui.wdtPages->setCurrentIndex(ui.wdtOptions->row(current));
}

void milxQtPreferencesForm::accept()
{
    //General
    MainWindow->preferWhiteBackground(backgroundCheckBox->isChecked());
    MainWindow->preferHumanGlyph(humanCheckBox->isChecked());
    MainWindow->preferSubWindowSize(windowSizeEdit->value());
    MainWindow->preferMaximumProcessors(processorsEdit->value());
    MainWindow->preferScreenshotMagnifyFactor(magnifyEdit->value());
    MainWindow->preferTimestamps(timestampCheckBox->isChecked());
    
	//Imaging
    MainWindow->preferImageInterpolation(interpolationCheckBox->isChecked());
    MainWindow->preferOrientation(orientationCheckBox->isChecked());
    
	//Models
    MainWindow->preferModelInterpolation(interpolationModelCheckBox->isChecked());
    MainWindow->preferScalarBar(scalarBarCheckBox->isChecked());
    MainWindow->writeSettings();
    QDialog::accept();
}

void milxQtPreferencesForm::createConnections()
{
    connect(ui.wdtOptions, SIGNAL(currentItemChanged(QListWidgetItem*,QListWidgetItem*)),
                this, SLOT(changePage(QListWidgetItem*,QListWidgetItem*)));
}
