%% General info
 %--- FOUR STATES ---
 % S: Susceptible
 % I: Infected
 % R: Recovered
 % W: Vaccinated
 
 %--- TWO MODES ---
 % Mode 1: Conectivity factor is the same for every cell
 % Mode 2: Conectivity factor is diferent on west from east cells
 
 %Made on June 2020
%% Initialize variables regarding to the area, the population and experment's parameters 
x = 100;
y = 100;

%Every cell has the same population
S = single(ones(x,y));
I = single(zeros(x,y));
R = single(zeros(x,y));
% N = 50*ones(x,y);

%Population size varies among cells
N = ones(x,y);
N(:,1:ceil(y/2)) = 10;
N(:,ceil(y/2):y) = 100;

days = 100;

perdayS = zeros(1,days);
perdayI = zeros(1,days);
perdayR = zeros(1,days);

%Set the outbreak point
S(ceil(x/2),ceil(y/2)) = 0.7;
I(ceil(x/2),ceil(y/2)) = 0.3;

%Set the parameters for MODE 1
c = 1;
m = 0.5;
v = 0.6;
e = 0.2;
r = 0.03;

%Set the parameters for MODE 2
% c = ones(x,y);
% c(:,1:ceil(y/2)) = 0.3;
% c(:,ceil(y/2):y) = 1;
% 
% m = 0.5;
% v = 0.8;
% e = 0.2;
% r = 0.01;
%% Use SIR function to simulate the spreading of the virus
for day = 1:days
    
    % SIR model, more details inside function
    [S,I,R] = SIR(S,I,R,N,x,y,c,m,v,e,r);
    
    % Calculate the sum of S,I,R populations of each day
    sumS = 0;
    sumI = 0;
    sumR = 0;
    
    for i = 2:y-1
        sumS = sumS + sum(N(i,2:x-1).*S(i,2:x-1));
        sumI = sumI + sum(N(i,2:x-1).*I(i,2:x-1));
        sumR = sumR + sum(N(i,2:x-1).*R(i,2:x-1));
    end
    
    perdayS(day) = sumS;
    perdayI(day) = sumI;
    perdayR(day) = sumR;
   
    % Define the screen layout for plotting
    screen=zeros(x, y);

    [ blue , red, green ,black] = deal( 1, 2, 3, 4);

    screen(1,:) = 4;
    screen(x,:) = 4;
    screen(:,1) = 4;
    screen(:,y) = 4;

    for k = 2:(x-1)
          for l = 2:(y-1)
              
              if(S(k,l)==1)
                  screen(k,l)=blue;
              end
              
              if(I(k,l)>=0.4)
                  screen(k,l)=red;
              end
              
              if(R(k,l)>=0.6)
                  screen(k,l)=green;
              end
          end
    end
                  

    myColorMap = [ 0 0 1;
                   1 0 0
                   0 1 0
                   0 0 0];
 
    colormap(myColorMap);
    % On each day plot the heatmap of the S,I,R populations
    subplot(1,2,1)
    plot_matrix = imagesc(screen,[blue,black]); 
    drawnow

end

% Finaly, plot the overall change of S,I,R populations over time
subplot(1,2,2)
hold on
plot(perdayS)
hold on
plot(perdayI)
hold on
plot(perdayR)
hold off
legend('S', 'I', 'R');



