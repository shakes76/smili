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
    setTabsClosable(true);
    setObjectName("Manager");

//    newTab();

    createActions();
    createConnections();
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

    int newIndex = addTab(treeWidget, tabTitle);
    setCurrentIndex(newIndex);
    return newIndex;
}

void milxQtManager::closeTab(int index)
{
    int newIndex = 0;

    if(count() > 1 || index > 0)
    {
        removeTab(index);

        if(index > 0 && currentIndex() == index)
            newIndex = currentIndex();
        else if(index > 0)
            newIndex = index-1; //Safe to decrement
        setCurrentIndex(newIndex); //else use zero
    }
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

void milxQtManager::exportTab(QString filename)
{
  QTreeWidget *currentTree = qobject_cast<QTreeWidget *>( currentWidget() );

  if(!currentTree)
    return;

  QPointer<QFileDialog> fileOpener = new QFileDialog;
  QSettings settings("Shekhar Chandra", "milxQt");
  QString path = settings.value("recentPath").toString();

  if(filename.isEmpty())
  {
      QFileDialog *fileOpener = new QFileDialog;
      filename = fileOpener->getSaveFileName(this,
                              tr("Select File to Save"),
                              path,
                              tr(openOtherExts.c_str()) ); //!< \todo Check and validate extensions support at Open in Main class
  }

  if(filename.isEmpty())
    return;

  QFile txtFile;
  txtFile.setFileName(filename);
  if(!txtFile.open(QIODevice::WriteOnly | QIODevice::Text))
  {
      milx::PrintError("Could not open text file for export.");
      return;
  }

  QTextStream outFile(&txtFile);
  QTreeWidgetItemIterator itemIterator(currentTree);
  while (*itemIterator)
  {
      QTreeWidgetItem *item = *itemIterator;
      for(int j = 0; j < item->columnCount(); j ++)
      {
          outFile << item->text(j);
          if(j < item->columnCount()-1)
             outFile << ", ";
      }
      outFile << endl;
      itemIterator ++;
  }
  txtFile.close();
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

void milxQtManager::addTreeItem(int tabIndex, QStringList topLevelName, QList<QStringList> entryList, Qt::ItemFlags flags)
{
    QTreeWidget *currentTree = qobject_cast<QTreeWidget *>(widget(tabIndex));
    //QList<QTreeWidgetItem *> items;

    //Add top level item first
    QTreeWidgetItem *topLevelEntry = new QTreeWidgetItem(currentTree, topLevelName);
    currentTree->addTopLevelItem(topLevelEntry);

    //Add leaves
    for (int j = 0; j < entryList.size(); j++)
    {
      QTreeWidgetItem *newEntry = new QTreeWidgetItem(topLevelEntry, entryList[j]);
      newEntry->setFlags(flags);
      //items.append(newEntry);
      currentTree->itemBelow(newEntry);
    }

    //currentTree->setCurrentItem(0, items);
}

void milxQtManager::createActions()
{
    actionExportTab = new QAction(this);
    actionExportTab->setText(QApplication::translate("Manager", "Export Contents ...", 0));
    actionExportTab->setShortcut(tr("Alt+e"));
    actionClearTab = new QAction(this);
    actionClearTab->setText(QApplication::translate("Manager", "Clear Tab", 0));
    actionClearTab->setShortcut(tr("Alt+t"));
    actionClear = new QAction(this);
    actionClear->setText(QApplication::translate("Manager", "Clear Manager", 0));
    actionClear->setShortcut(tr("Alt+c"));
}

void milxQtManager::createConnections()
{
    connect(this, SIGNAL(tabCloseRequested(int)), this, SLOT(closeTab(int)));

    connect(actionExportTab, SIGNAL(triggered()), this, SLOT(exportTab()));
    connect(actionClearTab, SIGNAL(triggered()), this, SLOT(clearTab()));
    connect(actionClear, SIGNAL(triggered()), this, SLOT(clear()));
}

void milxQtManager::contextMenuEvent(QContextMenuEvent *currentEvent)
{
    QMenu* contextMenu = new QMenu(this); //!< Only exists for the duration of the context selection

    contextMenu->addAction(actionExportTab);
    contextMenu->addAction(actionClearTab);
    contextMenu->addAction(actionClear);

    contextMenu->exec(currentEvent->globalPos());
}
