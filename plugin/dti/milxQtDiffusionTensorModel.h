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
#ifndef milxQtDIFFUSIONTENSORMODEL_H
#define milxQtDIFFUSIONTENSORMODEL_H

#include "milxQtModel.h"

class MILXQT_PLUGIN_EXPORT milxQtDiffusionTensorModel : public milxQtModel
{
    Q_OBJECT

public:
    milxQtDiffusionTensorModel(QWidget *theParent = 0);
    virtual ~milxQtDiffusionTensorModel();

    static double computeAmplitude(std::vector<double> SH, double x, double y, double z, int lmax);
    static int NforL (int lmax) { return ((lmax+1)*(lmax+2)/2); }
    static int LforN (int N) { return (2*(((int) (sqrt((float) (1+8*N)))-3)/4)); }
    static int index (int l, int m) { return (l*(l+1)/2 + m); }

public slots:
    void colourByDirection();
    void harmonics(QString ordersString = "");

protected:
    QAction *colourDirectionAct; //!< colour by direction axes action
    QAction *harmonicsAct; //!< Display spherical harmonic

private:
    void createActions();
    void createConnections();
    void contextMenuEvent(QContextMenuEvent *currentEvent);

};

#endif // milxQtDIFFUSIONTENSORMODEL_H
