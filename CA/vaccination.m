%% General info
 %Four states
 % S: Susceptible
 % I: Infected
 % R: Recovered
 % W: Vaccinated
 
 %Made on June 2020
%% Initialize variables regarding to the area, the population and experment's parameters 
x = 50;
y = 50;

S = single(ones(x,y));
I = single(zeros(x,y));
R = single(zeros(x,y));
N = 500*ones(x,y);
pop = 500;

days = 50;
temp = single(zeros(x,y));

perdayS = zeros(1,days);
perdayI = zeros(1,days);
perdayR = zeros(1,days);

%Set initial conditions
S(ceil(x/2),ceil(y/2)) = 0.7;
I(ceil(x/2),ceil(y/2)) = 0.3;

%Set the experiment parameters
c = 1;
m = 0.9;
v = 0.9;
e = 0.2;
r = 0.01;
w = [0, 0.2, 0.3, 0.4, 0.7, 0.8];
perVacRate = zeros(length(w),days);
%% Use SIR function for every vaccination rate
for rate = 1:length(w)
    
    %For every experiment initialize the population
    S = single(ones(x,y));
    I = single(zeros(x,y));
    R = single(zeros(x,y));
    
    S(ceil(x/2),ceil(y/2)) = 0.7;
    I(ceil(x/2),ceil(y/2)) = 0.3;
   
    %Simulate vaccination
    S = S - w(rate);
    R = R + w(rate);
    
    for day = 1:days

        %SIR model, more details inside function
        [S,I,R] = SIR(S,I,R,N,x,y,c,m,v,e,r);

        % Calculate the sum of S,I,R populations of each day
        sumS = 0;
        sumI = 0;
        sumR = 0;

        for i = 2:y-1
            sumS = sumS + sum(pop*S(i,2:x-1));
            sumI = sumI + sum(pop*I(i,2:x-1));
            sumR = sumR + sum(pop*R(i,2:x-1));
        end
        %Store the number of the infected indiviluas for every day
        perdayI(day) = sumI;

    end
    %Store the perDayI vector of each vaccination rate
    perVacRate(rate,:) = perdayI;
end

plot(perVacRate(1,:),'LineWidth',2.0)
hold on
plot(perVacRate(2,:),'LineWidth',2.0)
hold on
plot(perVacRate(3,:),'LineWidth',2.0)
hold on
plot(perVacRate(4,:),'LineWidth',2.0)
hold on
plot(perVacRate(5,:),'LineWidth',2.0)
hold on
plot(perVacRate(6,:),'LineWidth',2.0)
hold off

legend('0%', '20%', '30%', '40%', '70%', '80%');


