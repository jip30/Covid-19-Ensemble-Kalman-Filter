%% Obtaining the Data
clc
clear all
close all
% Reading the CSV Data. This gives the number of cases and deaths in the UK from
% 01-09-21 until 30-01-20. Column 4 = date, Column 5 =
% cumCasesByPublishDate, column 6 = cumDailyNsoDeathsByDeathDate
Data = readtable('Data/overview_2021-09-01.csv');

format long

%Ensemble Parameters
m = 1000; %number of ensemble members
tol = 0.001; %convergence criteria



% Formatting the data. Want the data from June 1st - July 31st 2021. This
% data is only measuring I and D and therefore will have a lower dimension
% that the state vector. 
Poscases = flip(Data.cumCasesByPublishDate(33:93));
Deaths = flip(Data.cumDailyNsoDeathsByDeathDate(33:93));


%% Parameters

dim = 4; %dimension of the problem (4 state variables, SIRD).

N= length(Poscases); %Number of time steps (days) (could be 60, idk yet)

%Parameters
beta = 3; %transmisson rate
gamma = 0.1; %recovery rate
delta = 0.01; %Death rate
Pop = 67e6; %population of the UK (or initial population, N = S+I+R+D)

param = [beta,gamma,delta,Pop];

%% Initialisation

%Initial data
I0 = Data.cumCasesByPublishDate(94);
D0 = Data.cumDailyNsoDeathsByDeathDate(94);
R0 = I0-D0;
S0 = Pop-I0-R0-D0;

xf = [S0; I0; R0; D0];


%Row = time step, Columns = SIRD. This is for plotting later
Xf = zeros(N,dim);
Xa = zeros(N,dim);
Y = zeros(N,dim);


%The operators M and H: M is the SIRD model, i.e. M= SIRDsolve.m
%The measurement operator. Only measuring 2/4 data points so
H = zeros(dim,dim);
H(2,2) = 1;
H(4,4) = 1;


%%


% %Noise: we have zero-mean Gaussian noise for both
q=0.4; %associated with model
r=0.1; %asscoiated with data

Q= (q^2)*eye(dim); %model error covariance
% R= (r^2)*eye(dim); %data error covariance


%Getting the corresponding covariance matrix 

Pf = Q; %initial Pf
XF = zeros(dim, N);

for i =1:N

    %Observation: getting the Data y
    j = linspace(-3,3,10000);
    idx1 = randi([1 length(j)],m,1);
    idx2 = randi([1 length(j)],m,1);
    idx3 = randi([1 length(j)],m,1);
    idx4 = randi([1 length(j)],m,1);
    
    %Rows = Ensemble, Columns = SIRD
    y = [0 ; Poscases(i); 0; Deaths(i)]; %this is the data, before the ensemble
    
    
    u = [normpdf(j(idx1),0,r^2); normpdf(j(idx2),0,r^2); normpdf(j(idx3),0,r^2); normpdf(j(idx4),0,r^2)];
    
    
    
    
    Yi = y*(ones(size(idx1))') + u; % this is the enesemble of data
    %Y(i,:) = y';
    
    %Computing the emperical error covariance matrix
    R = u*u';
    %R = zeros(4,4);
    
    
    %Analysis
    %Kalman Gain
    K = (Pf*(H'))/((H*Pf*H' +R));
    
    %Computing the Analysis Estimate for the entire ensemble
    xa = zeros(dim,m);
    for k =1:m
        xai = xf + K*(Yi(:,k)-H*xf);
        xa(:,k) = xai;
    end
%     xa = xf + K*(y-H*xf);
%     Xa(i,:) = xa';

    %Mean of the ensemble analyses
    xa_bar = zeros(4,1);
    for k = 1:m
        xa_bar = xa_bar+xa(:,k);
    end
    xa_bar = xa_bar/m;
    
    
    
    %Analysis Covariance
    %Pa = (eye(dim) - K*H)*Pf; 
    Pa = ((xa - xa_bar)*(xa - xa_bar)')/(m-1);
    
    
    %Forecasting
    %for the ensemble
    xf = zeros(dim,m);
    for k =1:m
        [~,xf_range] = SIRDsolve(xa(:,k), [i,i+1], param);
        xfi = (xf_range(end,:))';
        xf(:,k) = xfi;
    end
    
    
    %Mean of the ensemble forecasts
    xf_bar = zeros(4,1);
    for k = 1:m
        xf_bar = xf_bar+xa(:,k);
    end
    xf_bar = xf_bar/m;
    
    
    %Forecast Error Covariance 
    Pf = ((xf - xf_bar)*(xf - xf_bar)')/(m-1);
    
    xf = xf_bar;
    XF(:,i) = xf;
end


%%

% Plotting the results


Date = flip(Data.date(33:93));
date = datetime(Date); 

PopBar = Pop*ones(N,1);

figure(1);
plot(date,XF(1,:), 'LineWidth', 1)
hold on 
xtickformat('dd-MMM-yyyy')
xlabel('Date')
title('Susceptible')
axis tight
hold off

figure(2);
plot(date,XF(2,:), 'LineWidth', 1)
hold on 
plot(date,Poscases,'LineWidth', 1)
legend('Forecast', 'Data')
xtickformat('dd-MMM-yyyy')
xlabel('Date')
title('Infected')
axis tight
hold off

figure(3);
plot(date, XF(3,:), 'LineWidth', 1)
hold on 
xtickformat('dd-MMM-yyyy')
xlabel('Date')
title('Recovered')
axis tight
hold off

figure(4);
plot(date,XF(4,:), 'LineWidth', 1)
hold on 
plot(Deaths,'LineWidth', 1)
legend('Forecast','Data')
xtickformat('dd-MMM-yyyy')
xlabel('Date')
title('Dead')
axis tight
hold off


figure(5);
plot(date,XF(1,:), 'LineWidth', 1)
hold on 
plot(date,XF(2,:), 'LineWidth', 1)
plot(date,XF(3,:), 'LineWidth', 1)
plot(date,XF(4,:), 'LineWidth', 1)
plot(date,PopBar, 'LineWidth', 1)
legend('S','I','R','D','Pop')
xtickformat('dd-MMM-yyyy')
xlabel('Date')
title('Dead')

hold off