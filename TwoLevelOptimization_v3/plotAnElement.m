function [] = plotAnElement(macro_meso_iteration,elementNumber)

folderNum=0;

  outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
  x = csvread(outname);
    p=plotResults;
       p. PlotArrayGenericWithBlueWhiteColors( x, sprintf('element %i and iteration %i',elementNumber,macro_meso_iteration));
       
       outname = sprintf('./mesoDesignforElement%i.png',elementNumber);
       print(outname,'-dpng');

    