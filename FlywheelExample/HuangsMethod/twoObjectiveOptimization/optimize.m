function  [history,searchdir] = optimize(w1,w2,iterationNum)
LASTN = maxNumCompThreads(8)

w1 = str2double(w1)
w2 = str2double(w2)



% flywheel optimizatoin using fmincon
% do equation 5.15

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];

% x = numberInputs/2 control heights + numberInputs/2 vol fract heights
% x = [profile_ctrlpoints   vol_frac_ctrl_points]

numberInputs = 16;

half = numberInputs/2;


A = zeros(numberInputs,numberInputs);
b = zeros(numberInputs,1);
Aeq = [];
beq = [];
lb(1:half) = 0.02; % lower bounds for the profile control points
ub(1:half) = 0.1; % upper bounds for the profile control points

lb(half+1:numberInputs) = 0; % lower bounds for the vol frac control points
ub(half+1:numberInputs) = 1; % upper bounds for the vol frac  control points

x0(1:half) = (0.1-0.02)/2; % start in the middle of profile control points bounds
x0(half+1:numberInputs) = 0.5; % start in the middle of vol frac control points bounds

options =   optimoptions('fmincon','MaxFunEvals',1000,'TolCon',1,'TolFun',100,'TolX',0.01,'UseParallel','Always','Display','iter'); %,'OutputFcn',@outfun,'Display','iter','Algorithm','active-set');

%w1 = 0;
%w2 = 1;

 results = zeros(11,numberInputs);
%for i = 1:10
%i = evalNum
%w1 = i/10;
%w2 = 1-w1;
 message = sprintf('Optimize where w1 = %4.2f and w2 = %4.2f',w1,w2);
disp(message);

   [results,fval] = fmincon(@(x) Objective(x,w1,w2),x0,A,b,Aeq,beq,lb,ub,@nonlinearConstraints, options);
   filename= strcat('results',num2str(iterationNum), '.csv');
%end




csvwrite(filename,results)
results


% function stop = outfun(x,optimValues,state)
%      stop = false;
%  
%      switch state
%          case 'init'
%              hold on
%          case 'iter'
%          % Concatenate current point and objective function
%          % value with history. x must be a row vector.
%            history.fval = [history.fval; optimValues.fval];
%            history.x = [history.x; x];
%          % Concatenate current search direction with 
%          % searchdir.
%            searchdir = [searchdir;... 
%                         optimValues.searchdirection'];
%            plot(x(1),x(2),'o');
%          % Label points with iteration number and add title.
%          % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),... 
%                 num2str(optimValues.iteration));
%            title('Sequence of Points Computed by fmincon');
%          case 'done'
%              hold off
%          otherwise
%      end
%  end
end

