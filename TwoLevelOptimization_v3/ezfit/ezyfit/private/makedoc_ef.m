%MAKEDOC_EF  Generate the HTML help files for the EzyFit toolbox
%   (internal use only - 19 oct 2006)

curdir = pwd;

cd ..

try
    makehtmldoc('*.m', 'upper', 'noad', 'quiet',...
        'title', '\f (Ezyfit Toolbox)',...
        'firstline', 'EzyFit Function Reference',...
        'lastline', '2005-2016 <a href="ezyfit.html">EzyFit Toolbox 2.44</a>');

    delete Contents.html;

    movefile *.html html;
    
    cd(curdir);
catch
    cd(curdir);
    error('makehtmldoc not found.'); 
end

clear curdir
