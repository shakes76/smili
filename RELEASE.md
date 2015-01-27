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

# SMILI Release v0.99 (Beta)

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