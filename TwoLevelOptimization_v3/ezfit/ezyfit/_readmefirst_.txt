
   EzyFit Toolbox for Matlab

   F. Moisy
   Lab. FAST - University Paris-Sud

   version 2.44, 24 May 2016
__________________________________________________________________________

See:  www.fast.u-psud.fr/ezyfit

Bug report:  moisy@fast.u-psud.fr
__________________________________________________________________________

The Ezyfit toolbox for Matlab enables you to perform simple curve fitting
of one-dimensional data using arbitrary fitting functions. It provides
command-line functions and a basic graphical user interface for interactive
selection of the data.
__________________________________________________________________________

Installation and system requirements

EzyFit needs Matlab 7.0 or higher. It has been tested under 7.0 to 8.6
(R2015b), but mainly under Windows. The command-line functions (e.g.
ezfit, showfit...) work equally well on all systems. However graphical
operations (e.g. selectfit, getslope...) may not be fully stable,
especially on non-Windows systems, and/or for most recent Matlab releases.

1. Download and unzip the EzyFit Toolbox in a directory somewhere in your
system. For instance, in a Windows installation, the directory
Documents/MATLAB/ezyfit may be a good location. Do NOT install
the toolbox in the Matlab directory itself (Program Files/Matlab directory
in Windows). If you upgrade from an older version, first empty the previous
directory.

2.Select 'Set Path' (available in the menu 'File' in Matlab 7, or in the
tab 'Home' in Matlab 8). In the dialog box, click on 'Add Folder' (NOT
'with subfolders') and select the 'ezyfit' directory. Click on 'Save' and
'Close'.

3. If you want to always have the Ezyfit menu in your figures, type
efmenu install. This will create or update your 'startup.m' file.
Careful: If you apply this step 4, all the figure files (.FIG) saved will
include the Ezyfit menu. See the function 'remove_efmenu_fig' to remove
the menu from the figure file.

Note: If you upgrade Matlab and you want to use your previous Ezyfit
installation, you just have to follow the steps 2-3.

__________________________________________________________________________

What's new?

Once installed, have a look to the Release Notes and the Known Bugs sections
in the help browser (type 'docezyfit').
__________________________________________________________________________

Copyrights

This toolbox is covered by the BSD License. See license.txt file.
__________________________________________________________________________

The End.
