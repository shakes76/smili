# SMILI Release v1.0 Beta 2
**Highlights: Image Blending via a mixer, VTK 6.3 support, Image Cursors, DICOM Improvements, Image Autolevel, Bilateral Image filtering and surface flattening.**

Changelog:
Updated the apps that were out of date.
Fixed distance widget and others that were not working before.
Replaced crosshairs icon.
Fixed toolbar buttons not working for recent changes. Fixed human glyph preferences bug. Fixed instant display of crosshairs when enabled. Fixed window level reset in refresh for images.
Added refresh and cursors buttons for all windows in tab
Fixed VTK 5 compile errors. Added initial crosshair icon to be used later.
Added working linked and tracked cursors. Updated callback to invoke modified event for cursors. Added necessary API changes to access cursor components.
Added replacement cursor using ResliceCursor. Replaced Cursor3D with better looking one. Updated image class and main class to account for this. Added action for cursor.
Added copying of cursor position to other windows. Doesn't work but almost there I think.
Removed unused code for picking.
Added working version of an image cursor.
Added initial version of the cursor. Not working but code done.
Updated dev notes for recent changes.
Cleaned up auto level code for images. Removed Otsu version.
Completed working version of the auto level feature.
Properly reverted Otsu. Not great, needs fixing.
Added more code to do auto levelling using ITK. Still not working and reverted to Otsu levelling.
Added message for interpolator in orientation.
Merge branch 'topic-vtk6.3' into experimental
Fixed VTK 6 bug for orientation and output type. Added better histogram without zeros but doesn't make difference.
Added Otsu based level with autoLevel function. Augments previous approach with stats based approach.
Added notes and tweaked window leveling in autoLevel.
Completed working auto leveling and button. Levelling now works properly from inter-quartile ranges.
Added auto level button and better auto leveller. Moved info member to own section in image class. No code change. Moved histogram variable to class as needed throughout.
Fixed another VTK float type issue from 6.3. Need to check if change to float causes more issues.
Fixed render class for VTK 6.3 changes.

# SMILI Release v1.0 Beta 1
Fixed recent change to track view for multi-view. Replaced tracking with MMB instead of LMB.
Made the Setup Widgets member virtual for derived classes to redefine widgets.
Updated the widget default placement for images. 
Updated dev notes and todo list.
Added plane widget to interaction and annotation.
Added box widget for render windows.
Added bi directional and line widgets to render windows.
Added screenshot that shows updated multi-view and website.
Fixed recent changes of scalarstats breaking variance computation.
Improved scalar stats to ignore meshes with no scalars
Added more documentation for model app Updated Hausdorff script.
Fixed clipping of meshes in app and added more python scripts.
Fixed issues introduced by incorrect merging. This broke image orientation display and was confined to only the experimental branch.

# SMILI Release v1.0 Alpha
Added Save support for FTL PGMs image files. Reworked the patch scripts. Untested.
Added messages for view types boxes etc. Added refresh icon to reset action.
Added linking of view box in toolbar to window views. The view will now match those windows to view set immediately. 
Added working image contrast slider in toolbar. Works for all images windows per tab. Disabled dial, not sure if needed.
Fixed another issue with image size and load crash in the FTL plugin.
Fixed load crash for the FTL plugin.
Small fix to the about dialog to commodate different resolutions.
Fixed plots not working for model collection loaded SSMs.
Another fix to the script to fix output name.
Fixed Python Qt not being installed in packages properly. Fixed Python Qt paths settings.
Fix to script for executing custom module and Improved script case ID retrieval.

# SMILI Release v0.99 (Initial)

This is the first version of SMILI released to the public. A number of operating systems are supported including:
* Windows 8/7/Vista/XP 64/32-bit
* Mac OSX 10.6-10.8 64-bit
* Linux Mint 17 and Ubuntu LTS 64-bit
* OpenSUSE 13.1 64-bit

Notes:
* SSM versions of the installers have the Statistical Shape Model plugin built and ready for use.
* Due to EULA restrictions from Apple for virtual machines, Mac OSX releases will be made only when a Mac becomes available for builds.
* Windows 7 installer does not support Windows 8. These installers are marked explicitly as 'Windows7'
* 32-bit installers are not recommended and only provided for convenivence because you will likely run out of memory for even the most common tasks.

Known Issues:
* VTK 5.10.1 has a number of issues related to gamma of images, colour of axes etc. and is not recommended.