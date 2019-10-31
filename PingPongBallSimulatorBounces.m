%{

The MIT License (MIT) 

Copyright (c) 2018 Andrew Wuller 

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal 
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions: 

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE. 



Ping Pong Ball Trajectory Simulator, with bouncing.
Created by: Andrew Wuller

Credit to Russell J. Phelan for showing me how to implement the Runge-Kutta
Method in this code. I implemented how he used his in-air trajectory
calculation method, with some modifications. However, the rest of the
script is entirely original.

Web Link for code: https://github.com/russphelan/matlab/tree/master/Baseball%20Simulation 


%}




% Z is up
% Thought about originally using [time, values] =
% ode45(@PingPongBallSimulatorBallMotion,tspan,position);
% for graphing it, but this way allows  effective graphing techniques.

clear all;
close all;  


% Input initial conditions
% Initial values are used for debugging purposes.

velocity = [5, 0, -2];
spin = [1000, 0, 100];
position = [0, 0, .2];
xMod = 20;
yMod = 10;


% If using above values for debugging purposes, then comment out the
% following 9 lines 

%{
velocity(1) = input('Initial velocity x (m/s)?\n');
velocity(2) = input('Initial velocity y (m/s)?\n');
velocity(3) = input('Initial velocity z (m/s)?\n');
spin(1) = input('Initial spin x (rad/s)?\n');
spin(2) = input('Initial spin y (rad/s)?\n');
spin(3) = input('Initial spin z (rad/s)?\n');
position(1) = input('Initial position x?\n');
position(2) = input('Initial position y?\n');
position(3) = input('Initial position z?\n');
xMod = input('Speed modifier x?\n');
yMod = input('Speed modifier y?\n');
%}



% Constants

timeStep = .0001;
timeMax = 3;
totalTime  = 0:timeStep:timeMax;
i  = length(totalTime);
functionStep = ceil(length(totalTime) / 200);

accelGravity  = 9.80665;
airDensity= .9846;

ballMass  = 0.145;
coefficientDrag = 0.29;

ballRadius  = 0.02;
crossSectionalArea  = pi*(ballRadius^2);

magnusCoefficient = .00065;

positionValues = zeros(i,3);
velocityValues = zeros(i,3);
velocityTotalValues = zeros(i,1);
spinValues = zeros(i,3);
spinTotalValues = zeros(i,1);
accelerationValues = zeros(i,3);
bounceLeft = 0;
bounceRight = 0;
winner = 'Sim TimeOut';
winnerLock = false;

restitutionCoefficient = -0.905;

% Runge-Kutta 4th-5th (or whatever it's called) Method
belowTableReversal = false;
bounceNumber = 0;

velocityTotal = (velocity(1) ^ 2 + velocity(2) ^ 2 + velocity(3) ^ 2) ^ .5;
spinTotal = (spin(1) ^ 2 + spin(2) ^ 2 + spin(3) ^ 2) ^ .5;

fprintf('\nAt Start:\nTime: 0\nVelocity: %s\t%s\nSpin: %s\t%s\nPosition: %s\n\n',num2str(velocity),num2str(velocityTotal),num2str(spin),num2str(spinTotal),num2str(position));


%Net Height: 0.915, 0.1525
for i=1:i
    
    % Compute in-air trajectory until in contact with table
   
    
    % Compute bounce velocities, then return to above section and continue
    % calculating trajectory
    
    % First solve for collision, then z dimension bounce
    
    
    if position(1) > 1.35 & position(1) < 1.39 & position(2) > -0.917 & position(2) < 0.917 & position(3) < 0.1545
        %velocity = 0.1*velocity;
        if position(1) > 1.37 & winnerLock == false
            winner = 'left: Hit net';
            winnerLock = true;
        elseif position(1) <= 1.37 & winnerLock == false
            winner = 'right: Hit net';
            winnerLock = true;
        end
        
    end
    
    
    if position(3) < 0 & belowTableReversal == false & position(2)>-0.7625 & position(2)<0.7625
        velocity(3) = restitutionCoefficient * velocity(3);
        bounceNumber = bounceNumber + 1;
        if position(2)>-0.7625 & position(2)<0.7625 & position(1)<2.74 & position(1)>0
            if position(1)<=1.37
                bounceLeft = bounceLeft + 1;
                bounceRight = 0;
                
                if velocity(1) < 0 & winnerLock == false & bounceLeft == 2
                    winner = 'right: Double bounce';
                    winnerLock = true;
                    %velocity = 0;
                    %spin = 0;
                    %position = [-999 -999 -999];
                end
                
            elseif position(1)>=1.37
                bounceRight = bounceRight + 1;
                bounceLeft = 0;
                
                if velocity(1) > 0 & winnerLock == false & bounceRight == 2
                    winner = 'left: Double Bounce';
                    winnerLock = true;
                    %velocity = 0;
                    %spin = 0;
                    %position = [-999 -999 -999]
                end
                
            end
            
           
        end
        
        
        for j=1:2
            % Spin equations, is necessary for loops so it doesn't overwrite
            % the third value.
            
            velocity2(j) = velocity(j) * (3 - 2 * restitutionCoefficient) / 5 + spin(j) * (2 * restitutionCoefficient * ballRadius + 2 * ballRadius) / 5;
            spin2(j) = spin(j) * (2 - 3 * restitutionCoefficient) / 5 + velocity(j) * (3 * restitutionCoefficient + 3) / 5;

            velocity(j) = velocity2(j);
            spin(j) = spin2(j);
            
        end
         
        velocity(3) = - restitutionCoefficient * velocity(3);
        
        belowTableReversal = true;
        
        position = position + velocity * timeStep;
        
        velocityTotal = (velocity(1) ^ 2 + velocity(2) ^ 2 + velocity(3) ^ 2) ^ .5;
        spinTotal = (spin(1) ^ 2 + spin(2) ^ 2 + spin(3) ^ 2) ^ .5;
        
        fprintf('\nOn Bounce %s:\nTime: %s\nVelocity: %s\t%s\nSpin: %s\t%s\nPosition: %s\n\n',num2str(bounceNumber),num2str(totalTime(i)),num2str(velocity),num2str(velocityTotal),num2str(spin),num2str(spinTotal),num2str(position));
    
    else
        % This section is easier since it can be directly put in.
        
        acceleration = accel(velocity,coefficientDrag,airDensity,crossSectionalArea,ballMass,accelGravity,magnusCoefficient,spin);

        tempVelocity = velocity + acceleration*timeStep/2;

        tempAcceleration = accel(tempVelocity,coefficientDrag,airDensity,crossSectionalArea,ballMass,accelGravity,magnusCoefficient,spin);

        velocity = velocity + tempAcceleration*timeStep;
        
        position = position + tempVelocity*timeStep;
        
        velocityTotal = (velocity(1) ^ 2 + velocity(2) ^ 2 + velocity(3) ^ 2) ^ .5;
        spinTotal = (spin(1) ^ 2 + spin(2) ^ 2 + spin(3) ^ 2) ^ .5;
    
    end
    

    
    if position(3) > 0
        belowTableReversal = false;
    end 
    
    if position(1)<0
        if position(2)>-0.7625 & position(2)<0.7625
            velocity(1) = -(velocity(1)*(randi([-xMod xMod],1,1)/100+1));
            velocity(2) = -(velocity(2)*(randi([-yMod yMod],1,1)/100+1));
            if velocity(3) < 0
                velocity(3) = -velocity(3);
            else
                velocity(3) = velocity(3)*1.5;
            end
            spin = spin.*randi([-200,200],1,3)./100;

        else

        end
        bounceLeft = 0;
        bounceRight = 0;
    elseif position(1)>2.74
        if position(2)>-0.7625 & position(2)<0.7625
            %This is Al. Check for clipping
            velocity(1) = -(velocity(1)*(randi([-20 20],1,1)/100+1));
            velocity(2) = -(velocity(2)*(randi([-10 10],1,1)/100+1));
            if velocity(3) < 0
                velocity(3) = -velocity(3);
            else
                velocity(3) = velocity(3)*1.5;
            end
            spin = spin.*randi([-200,200],1,3)./100;
        else

        end
        bounceLeft = 0;
        bounceRight = 0;
    end
    positionValues(i,:) = position;
    velocityValues(i,:) = velocity;
    velocityTotalValues(i,:) = velocityTotal;
    spinValues(i,:) = spin;
    spinTotalValues(i,:) = spinTotal;
    accelerationValues(i,:) = acceleration;
    
    if winnerLock == true
        break
    end
end







% Table height is 76 cm
% Net height is 15.25 cm, overhang is also 15.25 cm
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

netDimensions = [1.37, 0.915, 0; 1.37, 0.915, 0.1525; 1.37, -0.915, .1525; 1.37, -0.915, 0; 1.37, 0.915, 0];

tableDimensions = [0, 0.7625, 0; 2.74, 0.7625, 0; 2.74, -0.7625, 0; 0, -0.7625, 0; 0, 0.7625, 0];

plot3(tableDimensions(1:5,1), tableDimensions(1:5,2), tableDimensions(1:5,3));
hold on;
plot3(netDimensions(1:5,1), netDimensions(1:5,2), netDimensions(1:5,3));
hold on;
plot3(positionValues(1:functionStep:end,1), positionValues(1:functionStep:end,2), positionValues(1:functionStep:end,3));
xlim(axes1,[-.25 3]);
ylim(axes1,[-1.625 1.625]);
zlim(axes1,[-.25 3]);
xlabel('X');
ylabel('Y');
zlabel('Height');
view(axes1,[45.6 38]);

%winner
fprintf('Winner: %s\n',winner);





%{
% Querying system to find values
queryInput = input('Perform query (yes = y)?\n','s');

if queryInput == 'y'
    
    querying = 0;
    while querying == 0
      	
        whichTable = input('Select Table:\nPosition (p)\nVelocity (v)\nSpin (s)\nAcceleration (a)\nTime (t)\nQuit (q)\n','s');
        
        if whichTable == 'q'
            querying = 1;
            continue;
        end
        
        if whichTable ~= 't'
            if whichTable == 'p'
                myTable = positionValues;
                dontUse = 1;
            elseif whichTable == 'v'
                myTable = velocityValues;
                dontUse = 2;
            elseif whichTable == 's'
                myTable = spinValues;
                dontUse = 3;
            elseif whichTable == 'a'
                myTable = accelerationValues;
                dontUse = 4;
            else
                fprintf('Your input of %s is not a valid input; please try again.\n',whichTable)
            end
            
            queryColumn = input('Enter the column which you would like to search (x, y, z): ','s');
            if or(queryColumn == 'x', queryColumn == 1)
                queryColumn = 1;
            elseif or(queryColumn == 'y', queryColumn == 2)
                queryColumn = 2;
            elseif or(queryColumn == 'z', queryColumn == 3)
                queryColumn = 3;
            else
                fprintf('Your input of %s is not a valid input; please try again.\n',queryColumn)
            end
            
            inputTerm = input('Enter the number for which you wish to search: ');
            
            thresholdValue = input('Enter the threshold value: ');
            
            myTable = myTable(:,queryColumn);
            
            rowValues = find(abs(myTable - inputTerm) < thresholdValue);
            
            answerTable = zeros(length(rowValues),19);
            totalTime = totalTime.';
            q = 1;
            for q=1:length(rowValues)
                
                answerTable(q,:) = [totalTime(rowValues(q)),0,positionValues(rowValues(q),:),0,velocityValues(rowValues(q),:),velocityTotalValues(rowValues(q),:),0,spinValues(rowValues(q),:),spinTotalValues(rowValues(q),:),0,accelerationValues(rowValues(q),:)];
                
            end
            
        else
            inputTerm = input('Enter the time for which you wish to search: ');
            
            thresholdValue = input('Enter the threshold value: ');
            
            totalTime = totalTime.';
            rowValues = find(abs(totalTime - inputTerm) < thresholdValue);
            q = 1;
            for q=1:length(rowValues)
            
                answerTable(q,:) = [totalTime(rowValues(q)),0,positionValues(rowValues(q),:),0,velocityValues(rowValues(q),:),velocityTotalValues(rowValues(q),:),0,spinValues(rowValues(q),:),spinTotalValues(rowValues(q),:),0,accelerationValues(rowValues(q),:)];
            
            end
            
            
        end
        
        answerTable
        fprintf('Columns are time, position, velocity, spin, and acceleration, respectively. Answer is saved in variable answerTable.\n')
        queryInput = input('Conduct another query (yes = y)?\n','s');
        if queryInput ~= 'y'
            querying = 1;
        end
        
        
    end    
    
end
%}

function accel = accel(v,Cd,rho,A,m,g,Dm,w)
vel = [v(1) v(2) v(3)];
vmag  = norm(vel);
magnus = m*Dm*cross(w,v);
force = [-0.5*Cd*rho*A*vmag*v(1) + magnus(1), -0.5*Cd*rho*A*vmag*v(2) + magnus(2), -0.5*Cd*rho*A*vmag*v(3) - m*g + magnus(3)];
accel = force/m;
end