function state = gaCustomPlo(options, state, flag)
csvwrite('populationCurrent.csv',state.Population);
csvwrite('scoresCurrent.csv',state.Scores);

