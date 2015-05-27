 [X,Y] = deal(1:40);             % Create grid reference
       Z = peaks(40);                  % Create grid height
       stlwrite('test.stl',X,Y,Z,'mode','ascii')
