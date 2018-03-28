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
    ///General Page
	QWidget *generalPage = new QWidget;
	QListWidgetItem *generalPageItem = new QListWidgetItem(ui.wdtOptions);
		generalPageItem->setText(tr("General"));
		generalPageItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
	
	//Application options
	//Window size (text)
	QLabel *windowSizeLabel = new QLabel(tr("Preferred Window Size:"));
	QSize desktopSize = qApp->desktop()->availableGeometry().size();
		windowSizeEdit = new QSpinBox;
		windowSizeEdit->setMinimum(minWindowSize);
		windowSizeEdit->setMaximum(milx::Maximum<int>(desktopSize.width(), desktopSize.height()));
		windowSizeEdit->setValue(MainWindow->hasPreferredSubWindowSize());
	QHBoxLayout *windowSizeLayout = new QHBoxLayout;
		windowSizeLayout->addWidget(windowSizeLabel);
		windowSizeLayout->addWidget(windowSizeEdit);
	//Number of processors (text)
	QLabel *processorsLabel = new QLabel(tr("Preferred Max. Processors:"));
	processorsEdit = new QSpinBox;
		processorsEdit->setMinimum(1);
		processorsEdit->setMaximum(milx::NumberOfProcessors());
		processorsEdit->setValue(MainWindow->hasMaximumProcessors());
	QHBoxLayout *processorsLayout = new QHBoxLayout;
		processorsLayout->addWidget(processorsLabel);
		processorsLayout->addWidget(processorsEdit);
	//Magnification of Screenshots
	QLabel *magnifyLabel = new QLabel(tr("Screenshot Magnification Factor"));
	magnifyEdit = new QSpinBox;
		magnifyEdit->setMinimum(1);
		magnifyEdit->setMaximum(1024);
		magnifyEdit->setValue(MainWindow->hasScreenshotMagnifyFactor());
	QHBoxLayout *magnifyLayout = new QHBoxLayout;
		magnifyLayout->addWidget(magnifyLabel);
		magnifyLayout->addWidget(magnifyEdit);
	//Timestamp in log? (check)
	timestampCheckBox = new QCheckBox(tr("Show timestamp in logs"));
	timestampCheckBox->setChecked(MainWindow->isTimestampsPreferred());
	// Save settings option

	// Load settings option

	//Application options layout
	QVBoxLayout *generalLayout = new QVBoxLayout;
		generalLayout->addLayout(windowSizeLayout);
		generalLayout->addLayout(processorsLayout);
		generalLayout->addLayout(magnifyLayout);
		generalLayout->addWidget(timestampCheckBox);
	QGroupBox *applicationGroup = new QGroupBox(tr("Application"));
		applicationGroup->setLayout(generalLayout);

	//View options
    //Background colour (check)
    backgroundCheckBox = new QCheckBox(tr("Prefer white background"));
        backgroundCheckBox->setChecked(MainWindow->isWhiteBackgroundPreferred());
    //Human glyph (check)
    humanCheckBox = new QCheckBox(tr("Show human orientation glyph"));
        humanCheckBox->setChecked(MainWindow->isHumanGlyphPreferred());
	//View options layout
	QVBoxLayout *viewLayout = new QVBoxLayout;
		viewLayout->addWidget(backgroundCheckBox);
		viewLayout->addWidget(humanCheckBox);
	QGroupBox *generalViewGroup = new QGroupBox(tr("View Options"));
		generalViewGroup->setLayout(viewLayout);
       	
	//Imaging options
	//Interpolation (check)
	interpolationCheckBox = new QCheckBox(tr("Interpolate images"));
	interpolationCheckBox->setChecked(MainWindow->isImageInterpolationPreferred());
	//Apply Orientation (check)
	orientationCheckBox = new QCheckBox(tr("Apply orientation to images"));
	orientationCheckBox->setChecked(MainWindow->isOrientationPreferred());
	//TODO Add custom colour maps

	//Imaging options layout
	QVBoxLayout *imagingLayout = new QVBoxLayout;
		imagingLayout->addWidget(interpolationCheckBox);
		imagingLayout->addWidget(orientationCheckBox);
	QGroupBox *imagingGroup = new QGroupBox(tr("Imaging"));
		imagingGroup->setLayout(imagingLayout);
	
	//Model options
	//Interpolation (check)
	interpolationModelCheckBox = new QCheckBox(tr("Interpolate (Phong Shading) Models"));
	interpolationModelCheckBox->setChecked(MainWindow->isModelInterpolationPreferred());
	///Scalar bar (Check)
	scalarBarCheckBox = new QCheckBox(tr("Always show scalar bar"));
	scalarBarCheckBox->setChecked(MainWindow->isScalarBarPreferred());
	//layout
	QVBoxLayout *ModelLayout = new QVBoxLayout;
		ModelLayout->addWidget(interpolationModelCheckBox);
		ModelLayout->addWidget(scalarBarCheckBox);
	QGroupBox *modelGroup = new QGroupBox(tr("Models"));
		modelGroup->setLayout(ModelLayout);

	//General page layout
	QVBoxLayout *generalPageLayout = new QVBoxLayout;
		generalPageLayout->addWidget(generalViewGroup);
		generalPageLayout->addWidget(applicationGroup);
		generalPageLayout->addWidget(imagingGroup);
		generalPageLayout->addWidget(modelGroup);
		generalPageLayout->setAlignment(Qt::AlignTop);
	generalPage->setLayout(generalPageLayout);

	//Add General options page to preferences
    ui.wdtPages->insertWidget(0, generalPage);
    ui.wdtPages->setCurrentWidget(generalPage);
    ui.wdtOptions->setCurrentItem(generalPageItem);
	
    ///Plugins Page
    QWidget *pluginsPage = new QWidget;
	QVBoxLayout *pluginsPageLayout = new QVBoxLayout;
    QListWidgetItem *pluginsPageItem = new QListWidgetItem(ui.wdtOptions);
        pluginsPageItem->setText(tr("Plugins"));
        pluginsPageItem->setFlags( Qt::ItemIsSelectable | Qt::ItemIsEnabled );
	// Check if no plugins
	if (!MainWindow->getPlugins().size()) {
		noPluginMsg = new QLabel(tr("There are currently no plugins installed."));
		pluginsPageLayout->addWidget(noPluginMsg);
		pluginsPageLayout->setAlignment(Qt::AlignCenter);
	} else {
		///List and Disable
		QListWidget *pluginsList = new QListWidget;
		pluginsList->setDisabled(true);

		foreach(QPointer<milxQtPluginInterface> plugin, MainWindow->getPlugins())
		{
			QListWidgetItem *pluginsItem = new QListWidgetItem(pluginsList);
			pluginsItem->setText(plugin->name());
			pluginsItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
		}
		pluginsPageLayout->addWidget(pluginsList);
	}
    //Layout	
    pluginsPage->setLayout(pluginsPageLayout);
    ui.wdtPages->insertWidget(1, pluginsPage);
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
