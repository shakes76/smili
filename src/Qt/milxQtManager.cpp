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
#include "milxQtManager.h"

milxQtManager::milxQtManager(QWidget *theParent) : QTabWidget(theParent)
{
    setTabsClosable(false);
    setObjectName("Manager");

//    newTab();
}

milxQtManager::~milxQtManager()
{
    //dtor
}

int milxQtManager::newTab(QString tabTitle, QStringList headings)
{
    QTreeWidget *treeWidget = new QTreeWidget(this);

    treeWidget->setColumnCount(headings.size());
    treeWidget->setHeaderLabels(headings);

    return addTab(treeWidget, tabTitle);
//    setCurrentWidget(treeWidget);

//    QList<QTreeWidgetItem *> items;
//    for (int i = 0; i < headings.size(); ++i)
//        items.append( new QTreeWidgetItem((QTreeWidget*)0, headings[i]) );
//    treeWidget->insertTopLevelItems(0, items);
}

void milxQtManager::clearTab()
{
    QTreeWidget *currentTree = qobject_cast<QTreeWidget *>( currentWidget() );

    if(!currentTree)
        return;

    currentTree->clear();
}

void milxQtManager::clearTab(int tabIndex)
{
    QTreeWidget *currentTree = qobject_cast<QTreeWidget *>( widget(tabIndex) );

    if(!currentTree)
        return;

    currentTree->clear();
}

void milxQtManager::addItem(QStringList entries, Qt::ItemFlags flags)
{
    QTreeWidget *currentTree = qobject_cast<QTreeWidget *>( currentWidget() );
    QTreeWidgetItem *newEntry = new QTreeWidgetItem(currentTree);

    for(int j = 0; j < entries.size(); j ++)
        newEntry->setText(j, entries[j]);
    newEntry->setFlags(flags);
    currentTree->itemBelow(newEntry);
}

void milxQtManager::addItem(int tabIndex, QStringList entries, Qt::ItemFlags flags)
{
    QTreeWidget *currentTree = qobject_cast<QTreeWidget *>( widget(tabIndex) );
    QTreeWidgetItem *newEntry = new QTreeWidgetItem(currentTree);

    for(int j = 0; j < entries.size(); j ++)
        newEntry->setText(j, entries[j]);
    newEntry->setFlags(flags);

    currentTree->itemBelow(newEntry);
}

void milxQtManager::addItem(int tabIndex, QStringList entries, QWidget *itemWidgetToAdd, int widgetColumn)
{
    QTreeWidget *currentTree = qobject_cast<QTreeWidget *>( widget(tabIndex) );
    QTreeWidgetItem *newEntry = new QTreeWidgetItem(currentTree);

    for(int j = 0; j < entries.size(); j ++)
        newEntry->setText(j, entries[j]);

    currentTree->setItemWidget(newEntry, widgetColumn, itemWidgetToAdd);
    currentTree->itemBelow(newEntry);
}
