function [state,options,optchanged] = outPutFunction(options,state,flag)


csvwrite('populationVector.csv',state.Population );
csvwrite('popScoreVector.csv',state.Score  );
csvwrite('bestVector.csv',state.Best   );
csvwrite('SelectionVector.csv',state.Selection    );
csvwrite('generationNum.csv',state.Generation     );
sprintf('# geneations %i \n', state.Generation)


    optchanged = false;
