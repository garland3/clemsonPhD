function [averageDiffOfLastValues] = FindAvergeChangeOfLastValue(arrayOfValues, settings)
y2 = arrayOfValues;
t = settings.terminationAverageCount;
diff2 = y2-circshift(y2,[0,1]); % subtract the current one from the pervious on
diff2 = diff2(2:end-1);

% normalize and take ABS
% we want to normalize by max(abs) over all differences, so
diff2 = abs(diff2)/(max(abs(diff2)));

diff2 = diff2(end-t:end);
averagey2diff = sum(diff2)/t;
averageDiffOfLastValues = averagey2diff;