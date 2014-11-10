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
#ifndef milxQtDIFFUSIONTENSORIMAGE_H
#define milxQtDIFFUSIONTENSORIMAGE_H

#include "milxQtImage.h"
#include "milxQtModel.h"

class MILXQT_PLUGIN_EXPORT milxQtDiffusionTensorImage : public milxQtImage
{
    Q_OBJECT

public:
    milxQtDiffusionTensorImage(QWidget *theParent = 0);
    virtual ~milxQtDiffusionTensorImage();

public slots:
    void diffusionGlyphs(size_t subsampleFactor = 0, size_t resolution = 16);
    void diffusionGlyphs2(size_t subsampleFactor = 0, size_t resolution = 16);

signals:
    /*!
        \fn milxQtDiffusionTensorImage::resultAvailable(milxQtModel*)
        \brief Send signal that Resultant model is available for showing.
    */
    void resultAvailable(milxQtModel*);
    /*!
        \fn milxQtDiffusionTensorImage::resultAvailable(milxQtRenderWindow*)
        \brief Send signal that Resultant rendering is available for showing.
    */
    void resultAvailable(milxQtRenderWindow*);

protected:
    QAction *diffusionAct; //!< show diffusion glyphs
    QAction *diffusion2Act; //!< show diffusion glyphs

private:
    void createActions();
    void createConnections();
    void contextMenuEvent(QContextMenuEvent *currentEvent);
};

#endif // milxQtDIFFUSIONTENSORIMAGE_H
