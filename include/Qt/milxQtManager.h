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
#ifndef MILXQTMANAGER_H
#define MILXQTMANAGER_H

#include <QTabWidget>
#include <QTreeWidget>
//milxQt Specific
#include "milxQtAliases.h"

/**
    \class milxQtManager
    \brief A manager (tabbed) widget class for displaying information about data such as case ID etc.

    Plugins can utilise the class for case browsing etc.
*/
class MILXQT_EXPORT milxQtManager : public QTabWidget
{
    Q_OBJECT

public:
    /*!
        \fn milxQtManager::milxQtManager(QWidget *theParent = 0)
        \brief The standard constructor
    */
    milxQtManager(QWidget *theParent = 0);
    virtual ~milxQtManager();

public slots:
    /** \brief Creates a new tab in the manager with the title and headings provided
     *
     * \param tabTitle Title of the tab
     * \param headings Headings of the fields within the view
     * \return The index of the created tab
     */
    int newTab(QString tabTitle, QStringList headings);
    /** \brief Clear the current tab
     */
    void clearTab();
    /** \brief Clear the tab given by index
     *
     * \param tabIndex index for the tab to be cleared
     */
    void clearTab(int tabIndex);
    /** \brief Add all the entries in list to the view
     *
     * \param entries The entries for the current item or row to be add.
     * \param flags item properties, is it Qt::ItemIsEditable etc.
     */
    void addItem(QStringList entries, Qt::ItemFlags flags = Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    /** \brief Add all the entries in list to the view
     *
     * \param index tab index of the entry to be inserted
     * \param entries The entries for the current item or row to be add.
     * \param flags item properties, is it Qt::ItemIsEditable etc.
     */
    void addItem(int tabIndex, QStringList entries, Qt::ItemFlags flags = Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    /** \brief Add all the entries in list to the view, as well as the widget provided
     *
     * \param index tab index of the entry to be inserted
     * \param entries The entries for the current item or row to be add.
     * \param itemWidgetToAdd the widget to be added to the item
     * \param widgetColumn column where the widget is to go
     */
    void addItem(int tabIndex, QStringList entries, QWidget *itemWidgetToAdd, int widgetColumn);

protected:

private:

};

#endif // MILXQTMANAGER_H
