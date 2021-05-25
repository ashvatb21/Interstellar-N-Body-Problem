%% N-Body Problem where N = 6
clear all, close all, clc

% Nomenclature
% 1= Sun
% 2 = Mercury
% 3 = Venus
% 4 = Earth
% 5 = Mars
% 6 = Moon

% Orbital Velocities of Planets in km/s (Data from NASA Fact Sheet)
% Mercury
Mercury_Min = 38.86; Mercury_Max = 58.98; Mercury_Mean = 47.36;
% Venus
Venus_Min = 34.79; Venus_Max = 35.26; Venus_Mean = 35.02;
% Earth
Earth_Min = 29.29; Earth_Max = 30.29; Earth_Mean = 29.78;
% Mars
Mars_Min = 21.97; Mars_Max = 26.50; Mars_Mean = 24.07;
% Moon
Moon_Min = 0.970; Moon_Max = 1.082; Moon_Mean = 1.022;

% Gathering N-Body Information for the date of Jan 1, 2011
Julian = juliandate(2011,1,1);

% Position and Velocity Vectors
[r1, v1] = planetEphemeris(Julian,'Sun', 'Sun', '432t', 'km');
[r2, v2] = planetEphemeris(Julian,'Sun','Mercury','432t','km');
[r3, v3] = planetEphemeris(Julian,'Sun','Venus','432t','km');
[r4, v4] = planetEphemeris(Julian,'Sun','Earth','432t','km');
[r5, v5] = planetEphemeris(Julian,'Sun','Mars','432t','km');
[r6, v6] = planetEphemeris(Julian,'Sun','Moon','432t','km');

% Defining u0 (in meters)
u0 = 10^3 * [r1'; v1'; r2'; v2'; r3'; v3'; r4'; v4'; r5'; v5'; r6'; v6'];

% Masses (in kg)
m1 = 1.9891 * 10^30;
m2 = 3.285 * 10^23;
m3 = 4.867 * 10^24;
m4 = 5.97 * 10^24;
m5 = 7.35 * 10^22;
m6 = 7.34767 * 10^22;

% Defining the Parameters
G = 6.67*10^(-11); % Gravitational Constant
t0 = 0; % Inital Start Time
T = 60*60*24*687*1; % Final Time (1 Mars Revolution in Seconds)
dt = 60*60*24*0.002; % Time step (0.002 Days in Seconds)

% State-Space Form Vectors (For Reference)
%r1 = u(1:3); v1 = u(4:6);
%r2 = u(7:9); v2 = u(10:12);
%r3 = u(13:15); v3 = u(16:18);
%r4 = u(19:21); v4 = u(22:24);
%r5 = u(25:27); v5 = u(28;30);
%r6 = u(31:33); v6 = u(34:36);

% State-Space Model
f = @(u) [ u(4:6); (G*m2*(u(7:9)-u(1:3)))./(norm((u(7:9)-u(1:3))).^3) + (G*m3*(u(13:15)-u(1:3)))./(norm((u(13:15)-u(1:3))).^3) + (G*m4*(u(19:21)-u(1:3)))./(norm((u(19:21)-u(1:3))).^3) + (G*m5*(u(25:27)-u(1:3)))./(norm((u(25:27)-u(1:3))).^3) + (G*m6*(u(31:33)-u(1:3)))./(norm((u(31:33)-u(1:3))).^3); ...
        u(10:12); (G*m1*(u(1:3)- u(7:9)))./(norm((u(1:3)-u(7:9))).^3) + (G*m3*(u(13:15)-u(7:9)))./(norm((u(13:15)-u(7:9))).^3) + (G*m4*(u(19:21)-u(7:9)))./(norm((u(19:21)-u(7:9))).^3) + (G*m5*(u(25:27)-u(7:9)))./(norm((u(25:27)-u(7:9))).^3) + (G*m6*(u(31:33)-u(7:9)))./(norm((u(31:33)-u(7:9))).^3); ...
        u(16:18); (G*m1*(u(1:3)- u(13:15)))./(norm((u(1:3)-u(13:15))).^3) + (G*m2*(u(7:9)-u(13:15)))./(norm((u(7:9)-u(13:15))).^3) + (G*m4*(u(19:21)-u(13:15)))./(norm((u(19:21)-u(13:15))).^3) + (G*m5*(u(25:27)-u(13:15)))./(norm((u(25:27)-u(13:15))).^3) + (G*m6*(u(31:33)-u(13:15)))./(norm((u(31:33)-u(13:15))).^3); ...
        u(22:24); (G*m1*(u(1:3)- u(19:21)))./(norm((u(1:3)-u(19:21))).^3) + (G*m2*(u(7:9)-u(19:21)))./(norm((u(7:9)-u(19:21))).^3) + (G*m3*(u(13:15)-u(19:21)))./(norm((u(13:15)-u(19:21))).^3) + (G*m5*(u(25:27)-u(19:21)))./(norm((u(25:27)-u(19:21))).^3) + (G*m6*(u(31:33)-u(19:21)))./(norm((u(31:33)-u(19:21))).^3); ...
        u(28:30); (G*m1*(u(1:3)- u(25:27)))./(norm((u(1:3)-u(25:27))).^3) + (G*m2*(u(7:9)-u(25:27)))./(norm((u(7:9)-u(25:27))).^3) + (G*m3*(u(13:15)-u(25:27)))./(norm((u(13:15)-u(25:27))).^3) + (G*m4*(u(19:21)-u(25:27)))./(norm((u(19:21)-u(25:27))).^3) + (G*m6*(u(31:33)-u(25:27)))./(norm((u(31:33)-u(25:27))).^3); ...
        u(34:36); (G*m1*(u(1:3)- u(31:33)))./(norm((u(1:3)-u(31:33))).^3) + (G*m2*(u(7:9)-u(31:33)))./(norm((u(7:9)-u(31:33))).^3) + (G*m3*(u(13:15)-u(31:33)))./(norm((u(13:15)-u(31:33))).^3) + (G*m4*(u(19:21)-u(31:33)))./(norm((u(19:21)-u(31:33))).^3) + (G*m5*(u(25:27)-u(31:33)))./(norm((u(25:27)-u(31:33))).^3)] ;

%%%%%% AB2 Numerical Analysis %%%%%%%

tic

tk = 0; %starting time

%initialize iterates for AB2
u_AB2_k = u0;
u_AB2_km1 = u0;

%initialize vector that stores approximate soln at various times
u_AB2_approx = zeros( 36,T/dt );

%Initialize vector that stores the magnitude of velocities at T
u_AB2_vel = zeros(6,T/dt);

%advance to final time T
for i = 1 : T/dt
    
    if i < 2
        
        %advance with Heun's for 1st time step
        u_AB2_approx(:,i) = (u_AB2_k + 1/2*dt*(f(u_AB2_k) + f(u_AB2_k+dt*f(u_AB2_k))));
        
    else
        
        u_AB2_approx(:,i) = u_AB2_k + 1/2*dt*(-1*f(u_AB2_km1) + 3*f(u_AB2_k));
        
    end
    
    %update iterates
    tk = tk + dt;
    u_AB2_km1 = u_AB2_k;
    u_AB2_k = u_AB2_approx(:,i);
    
    %Calculating and storing the velocities of each body at T
    u_AB2_vel(1,i) = norm(norm(u_AB2_approx(4:6,i)));
    u_AB2_vel(2,i) = norm(norm(u_AB2_approx(10:12,i)));
    u_AB2_vel(3,i) = norm(norm(u_AB2_approx(16:18,i)));
    u_AB2_vel(4,i) = norm(norm(u_AB2_approx(22:24,i)));
    u_AB2_vel(5,i) = norm(norm(u_AB2_approx(28:30,i)));
    u_AB2_vel(6,i) = norm(norm(u_AB2_approx(34:36,i)));
        
end

t_AB2 = toc

% Gathering N-Body Information for 687 days in future from Jan 1, 2011 (T)
Julian2 = juliandate(2012,11,19);

% Position and Velocity Vectors at future position
[r1, v1] = planetEphemeris(Julian2,'Sun', 'Sun', '432t', 'km');
[r2, v2] = planetEphemeris(Julian2,'Sun','Mercury','432t','km');
[r3, v3] = planetEphemeris(Julian2,'Sun','Venus','432t','km');
[r4, v4] = planetEphemeris(Julian2,'Sun','Earth','432t','km');
[r5, v5] = planetEphemeris(Julian2,'Sun','Mars','432t','km');
[r6, v6] = planetEphemeris(Julian2,'Sun','Moon','432t','km');

% Calculating the Average, Minimum and Maximum Velocity of Each Planet
AB2_vel_avg = mean(u_AB2_vel,2);
AB2_vel_max = max(u_AB2_vel,[],2);
AB2_vel_min = min(u_AB2_vel,[],2);

% Assigning Variable Names
Planetary_Bodies = {'Mercury';'Venus';'Earth';'Mars';'Moon'};
Calculated_Position_AB2 = [u_AB2_approx(7,T/dt), u_AB2_approx(8,T/dt), u_AB2_approx(9,T/dt); u_AB2_approx(13,T/dt), u_AB2_approx(14,T/dt), u_AB2_approx(15,T/dt); u_AB2_approx(19,T/dt), u_AB2_approx(20,T/dt), u_AB2_approx(21,T/dt); u_AB2_approx(25,T/dt), u_AB2_approx(26,T/dt), u_AB2_approx(27,T/dt); u_AB2_approx(31,T/dt), u_AB2_approx(32,T/dt), u_AB2_approx(33,T/dt)];
JulianDate_Position_AB2 = 10^3 * [r2;r3;r4;r5;r6];
% Calculated_Position = [norm([u_AB2_approx(7,T/dt), u_AB2_approx(8,T/dt), u_AB2_approx(9,T/dt)]); norm([u_AB2_approx(13,T/dt), u_AB2_approx(14,T/dt), u_AB2_approx(15,T/dt)]); norm([u_AB2_approx(19,T/dt), u_AB2_approx(20,T/dt), u_AB2_approx(21,T/dt)]); norm([u_AB2_approx(25,T/dt), u_AB2_approx(26,T/dt), u_AB2_approx(27,T/dt)]); norm([u_AB2_approx(31,T/dt), u_AB2_approx(32,T/dt), u_AB2_approx(33,T/dt)])];
% JulianDate_Position = 10^3 * [norm(r2);norm(r3);norm(r4);norm(r5);norm(r6)];
Calculated_Final_V_AB2 = [norm([u_AB2_approx(10,T/dt), u_AB2_approx(11,T/dt), u_AB2_approx(12,T/dt)]); norm([u_AB2_approx(16,T/dt), u_AB2_approx(17,T/dt), u_AB2_approx(18,T/dt)]); norm([u_AB2_approx(22,T/dt), u_AB2_approx(23,T/dt), u_AB2_approx(24,T/dt)]); norm([u_AB2_approx(28,T/dt), u_AB2_approx(29,T/dt), u_AB2_approx(30,T/dt)]); norm([u_AB2_approx(34,T/dt), u_AB2_approx(35,T/dt), u_AB2_approx(36,T/dt)])];
JulianDate_Final_V_AB2 = 10^3 * [norm(v2);norm(v3);norm(v4);norm(v5);norm(v6)];
CalculatedMeanVelocity_AB2 = [AB2_vel_avg(2); AB2_vel_avg(3); AB2_vel_avg(4); AB2_vel_avg(5); AB2_vel_avg(6)];
KnownMeanVelocity_AB2 = 10^3 * [Mercury_Mean; Venus_Mean; Earth_Mean; Mars_Mean; Moon_Mean];

% Defining Tables
Table1 = table(Planetary_Bodies, Calculated_Position_AB2, JulianDate_Position_AB2);
Table2 = table(Planetary_Bodies, Calculated_Final_V_AB2, JulianDate_Final_V_AB2, CalculatedMeanVelocity_AB2, KnownMeanVelocity_AB2);

% Get the table in string form.
TString = evalc('disp(Table1)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(100)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Get the table in string form.
TString = evalc('disp(Table2)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(200)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Plotting the Orbits Estimated
figure(1)
scatter4(0,0,0,40,'filled')
grid on
hold on
plot3(u_AB2_approx(7,:), u_AB2_approx(8,:), u_AB2_approx(9,:),'LineWidth', 2)
scatter4(u_AB2_approx(7,1), u_AB2_approx(8,1), u_AB2_approx(9,1),40,'filled')
plot3(u_AB2_approx(13,:), u_AB2_approx(14,:), u_AB2_approx(15,:),'LineWidth', 2)
scatter4(u_AB2_approx(13,1), u_AB2_approx(14,1), u_AB2_approx(15,1),40,'filled')
plot3(u_AB2_approx(19,:), u_AB2_approx(20,:), u_AB2_approx(21,:),'LineWidth', 2)
scatter4(u_AB2_approx(19,1), u_AB2_approx(20,1), u_AB2_approx(21,1),40,'filled')
plot3(u_AB2_approx(25,:), u_AB2_approx(26,:), u_AB2_approx(27,:),'LineWidth', 2)
scatter4(u_AB2_approx(25,1), u_AB2_approx(26,1), u_AB2_approx(27,1),40,'filled')
plot3(u_AB2_approx(31,:), u_AB2_approx(32,:), u_AB2_approx(33,:),'LineWidth', 0.75)
scatter4(u_AB2_approx(31,1), u_AB2_approx(32,1), u_AB2_approx(33,1),10,'filled')
title('Inner Solar System (Adams–Bashforth Method)')
legend('Sun','Mercury Orbit','Mercury','Venus Orbit','Venus','Earth Orbit','Earth','Mars Orbit','Mars','Moon Orbit','Moon')


%%%%%% RK4 Numerical Analysis %%%%%%%

tk = 0; %starting time

tic

%initialize iterates for Huens
u_rk_k = u0;

%initialize vector that stores approximate soln at various times
u_rk_approx = zeros( 36,T/dt );

%Initialize vector that stores the magnitude of velocities at T
u_rk_vel = zeros(6,T/dt);

%advance to final time T
for i = 1 : T/dt
    
    %advance with RK4 for 1st time step
    y1 = f(u_rk_k);
    y2 = f(u_rk_k + 1/2*dt*y1);
    y3 = f(u_rk_k + (1/2)*dt*y2);
    y4 = f(u_rk_k + dt*y3);
    u_rk_approx(:,i) = u_rk_k + 1/6*dt*(y1 + 2*y2 + 2*y3 + y4);
    
    %update iterates
    tk = tk + dt; 
    u_rk_k = u_rk_approx(:,i);
    
    %Calculating and storing the velocities of each body at T
    u_rk_vel(1,i) = norm(norm(u_rk_approx(4:6,i)));
    u_rk_vel(2,i) = norm(norm(u_rk_approx(10:12,i)));
    u_rk_vel(3,i) = norm(norm(u_rk_approx(16:18,i)));
    u_rk_vel(4,i) = norm(norm(u_rk_approx(22:24,i)));
    u_rk_vel(5,i) = norm(norm(u_rk_approx(28:30,i)));
    u_rk_vel(6,i) = norm(norm(u_rk_approx(34:36,i)));   
    
end

t_r = toc

% Gathering N-Body Information for 687 days in future from Jan 1, 2011 (T)
Julian2 = juliandate(2012,11,19);

% Position and Velocity Vectors at future position
[r1, v1] = planetEphemeris(Julian2,'Sun', 'Sun', '432t', 'km');
[r2, v2] = planetEphemeris(Julian2,'Sun','Mercury','432t','km');
[r3, v3] = planetEphemeris(Julian2,'Sun','Venus','432t','km');
[r4, v4] = planetEphemeris(Julian2,'Sun','Earth','432t','km');
[r5, v5] = planetEphemeris(Julian2,'Sun','Mars','432t','km');
[r6, v6] = planetEphemeris(Julian2,'Sun','Moon','432t','km');

% Calculating the Average Velocity of Each Planet
rk_vel_avg = mean(u_rk_vel,2);
rk_vel_max = max(u_rk_vel,[],2);
rk_vel_min = min(u_rk_vel,[],2);

% Assigning Variable Names
Planetary_Bodies = {'Mercury';'Venus';'Earth';'Mars';'Moon'};
Calculated_Position_RK4 = [u_rk_approx(7,T/dt), u_rk_approx(8,T/dt), u_rk_approx(9,T/dt); u_rk_approx(13,T/dt), u_rk_approx(14,T/dt), u_rk_approx(15,T/dt); u_rk_approx(19,T/dt), u_rk_approx(20,T/dt), u_rk_approx(21,T/dt); u_rk_approx(25,T/dt), u_rk_approx(26,T/dt), u_rk_approx(27,T/dt); u_rk_approx(31,T/dt), u_rk_approx(32,T/dt), u_rk_approx(33,T/dt)];
JulianDate_Position_RK4 = 10^3 * [r2;r3;r4;r5;r6];
% Calculated_Position = [norm([u_rk_approx(7,T/dt), u_rk_approx(8,T/dt), u_rk_approx(9,T/dt)]); norm([u_rk_approx(13,T/dt), u_rk_approx(14,T/dt), u_rk_approx(15,T/dt)]); norm([u_rk_approx(19,T/dt), u_rk_approx(20,T/dt), u_rk_approx(21,T/dt)]); norm([u_rk_approx(25,T/dt), u_rk_approx(26,T/dt), u_rk_approx(27,T/dt)]); norm([u_rk_approx(31,T/dt), u_rk_approx(32,T/dt), u_rk_approx(33,T/dt)])];
% JulianDate_Position = 10^3 * [norm(r2);norm(r3);norm(r4);norm(r5);norm(r6)];
Calculated_Final_V_RK4 = [norm([u_rk_approx(10,T/dt), u_rk_approx(11,T/dt), u_rk_approx(12,T/dt)]); norm([u_rk_approx(16,T/dt), u_rk_approx(17,T/dt), u_rk_approx(18,T/dt)]); norm([u_rk_approx(22,T/dt), u_rk_approx(23,T/dt), u_rk_approx(24,T/dt)]); norm([u_rk_approx(28,T/dt), u_rk_approx(29,T/dt), u_rk_approx(30,T/dt)]); norm([u_rk_approx(34,T/dt), u_rk_approx(35,T/dt), u_rk_approx(36,T/dt)])];
JulianDate_Final_V_RK4 = 10^3 * [norm(v2);norm(v3);norm(v4);norm(v5);norm(v6)];
CalculatedMeanVelocity_RK4 = [rk_vel_avg(2); rk_vel_avg(3); rk_vel_avg(4); rk_vel_avg(5); rk_vel_avg(6)];
KnownMeanVelocity_RK4 = 10^3 * [Mercury_Mean; Venus_Mean; Earth_Mean; Mars_Mean; Moon_Mean];

% Defining Tables
Table1 = table(Planetary_Bodies, Calculated_Position_RK4, JulianDate_Position_RK4);
Table2 = table(Planetary_Bodies, Calculated_Final_V_RK4, JulianDate_Final_V_RK4, CalculatedMeanVelocity_RK4, KnownMeanVelocity_RK4);

% Get the table in string form.
TString = evalc('disp(Table1)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(300)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Get the table in string form.
TString = evalc('disp(Table2)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(400)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Plotting the Orbits Estimated
figure(2)
scatter4(0,0,0,40,'filled')
grid on
hold on
plot3(u_rk_approx(7,:), u_rk_approx(8,:), u_rk_approx(9,:),'LineWidth', 2)
scatter4(u_rk_approx(7,1), u_rk_approx(8,1), u_rk_approx(9,1),40,'filled')
plot3(u_rk_approx(13,:), u_rk_approx(14,:), u_rk_approx(15,:),'LineWidth', 2)
scatter4(u_rk_approx(13,1), u_rk_approx(14,1), u_rk_approx(15,1),40,'filled')
plot3(u_rk_approx(19,:), u_rk_approx(20,:), u_rk_approx(21,:),'LineWidth', 2)
scatter4(u_rk_approx(19,1), u_rk_approx(20,1), u_rk_approx(21,1),40,'filled')
plot3(u_rk_approx(25,:), u_rk_approx(26,:), u_rk_approx(27,:),'LineWidth', 2)
scatter4(u_rk_approx(25,1), u_rk_approx(26,1), u_rk_approx(27,1),40,'filled')
plot3(u_rk_approx(31,:), u_rk_approx(32,:), u_rk_approx(33,:),'LineWidth', 0.75)
scatter4(u_rk_approx(31,1), u_rk_approx(32,1), u_rk_approx(33,1),10,'filled')
title('Inner Solar System (Runge-Kutta 4 Method)')
legend('Sun','Mercury Orbit','Mercury','Venus Orbit','Venus','Earth Orbit','Earth','Mars Orbit','Mars','Moon Orbit','Moon')


%% N-Body Problem where N = 6 where Sun lost 25% of its mass
clear all, close all, clc

% Nomenclature
% 1 = Sun
% 2 = Jupiter
% 3 = Saturn
% 4 = Uranus
% 5 = Neptune
% 6 = Pluto

% Orbital Velocities of Planets in km/s (Data from NASA Fact Sheet)
% Jupiter
Jupiter_Min = 12.44; Jupiter_Max = 13.72; Jupiter_Mean = 13.06;
% Saturn
Saturn_Min = 9.09; Saturn_Max = 10.18; Saturn_Mean = 9.68;
% Uranus
Uranus_Min = 6.49; Uranus_Max = 7.11; Uranus_Mean = 6.80;
% Neptune
Neptune_Min = 5.37; Neptune_Max = 5.50; Neptune_Mean = 5.43;
% Pluto
Pluto_Min = 3.71; Pluto_Max = 6.10; Pluto_Mean = 4.67;

% Gathering N-Body Information for the date of Jan 1, 2011
Julian = juliandate(2011,1,1);

% Position and Velocity Vectors
[r1, v1] = planetEphemeris(Julian,'Sun', 'Sun', '432t', 'km');
[r2, v2] = planetEphemeris(Julian,'Sun','Jupiter','432t','km');
[r3, v3] = planetEphemeris(Julian,'Sun','Saturn','432t','km');
[r4, v4] = planetEphemeris(Julian,'Sun','Uranus','432t','km');
[r5, v5] = planetEphemeris(Julian,'Sun','Neptune','432t','km');
[r6, v6] = planetEphemeris(Julian,'Sun','Pluto','432t','km');


% Defining u0 (in meters)
u0 = 10^3 * [r1'; v1'; r2'; v2'; r3'; v3'; r4'; v4'; r5'; v5'; r6'; v6'];

% Masses (in kg)
m1 = 1.9891 * 10^30;
m2 = 1.898 * 10^27;
m3 = 5.683 * 10^26;
m4 = 8.681 * 10^25;
m5 = 1.024 * 10^26;
m6 = 1.309 * 10^22;


% Defining the Parameters
G = 6.67*10^(-11); % Gravitational Constant
t0 = 0; % Inital Start Time
T = 60*60*24*365*250; % Final Time (1 Pluto Years)
dt = 60*60*24*0.25; % Time step (0.01 Days in Seconds)

% State-Space Form Vectors (For Reference)
%r1 = u(1:3); v1 = u(4:6);
%r2 = u(7:9); v2 = u(10:12);
%r3 = u(13:15); v3 = u(16:18);
%r4 = u(19:21); v4 = u(22:24);
%r5 = u(25:27); v5 = u(28;30);
%r6 = u(31:33); v6 = u(34:36);

% State-Space Model
f = @(u) [ u(4:6); (G*m2*(u(7:9)-u(1:3)))./(norm((u(7:9)-u(1:3))).^3) + (G*m3*(u(13:15)-u(1:3)))./(norm((u(13:15)-u(1:3))).^3) + (G*m4*(u(19:21)-u(1:3)))./(norm((u(19:21)-u(1:3))).^3) + (G*m5*(u(25:27)-u(1:3)))./(norm((u(25:27)-u(1:3))).^3) + (G*m6*(u(31:33)-u(1:3)))./(norm((u(31:33)-u(1:3))).^3); ...
        u(10:12); (G*m1*(u(1:3)- u(7:9)))./(norm((u(1:3)-u(7:9))).^3) + (G*m3*(u(13:15)-u(7:9)))./(norm((u(13:15)-u(7:9))).^3) + (G*m4*(u(19:21)-u(7:9)))./(norm((u(19:21)-u(7:9))).^3) + (G*m5*(u(25:27)-u(7:9)))./(norm((u(25:27)-u(7:9))).^3) + (G*m6*(u(31:33)-u(7:9)))./(norm((u(31:33)-u(7:9))).^3); ...
        u(16:18); (G*m1*(u(1:3)- u(13:15)))./(norm((u(1:3)-u(13:15))).^3) + (G*m2*(u(7:9)-u(13:15)))./(norm((u(7:9)-u(13:15))).^3) + (G*m4*(u(19:21)-u(13:15)))./(norm((u(19:21)-u(13:15))).^3) + (G*m5*(u(25:27)-u(13:15)))./(norm((u(25:27)-u(13:15))).^3) + (G*m6*(u(31:33)-u(13:15)))./(norm((u(31:33)-u(13:15))).^3); ...
        u(22:24); (G*m1*(u(1:3)- u(19:21)))./(norm((u(1:3)-u(19:21))).^3) + (G*m2*(u(7:9)-u(19:21)))./(norm((u(7:9)-u(19:21))).^3) + (G*m3*(u(13:15)-u(19:21)))./(norm((u(13:15)-u(19:21))).^3) + (G*m5*(u(25:27)-u(19:21)))./(norm((u(25:27)-u(19:21))).^3) + (G*m6*(u(31:33)-u(19:21)))./(norm((u(31:33)-u(19:21))).^3); ...
        u(28:30); (G*m1*(u(1:3)- u(25:27)))./(norm((u(1:3)-u(25:27))).^3) + (G*m2*(u(7:9)-u(25:27)))./(norm((u(7:9)-u(25:27))).^3) + (G*m3*(u(13:15)-u(25:27)))./(norm((u(13:15)-u(25:27))).^3) + (G*m4*(u(19:21)-u(25:27)))./(norm((u(19:21)-u(25:27))).^3) + (G*m6*(u(31:33)-u(25:27)))./(norm((u(31:33)-u(25:27))).^3); ...
        u(34:36); (G*m1*(u(1:3)- u(31:33)))./(norm((u(1:3)-u(31:33))).^3) + (G*m2*(u(7:9)-u(31:33)))./(norm((u(7:9)-u(31:33))).^3) + (G*m3*(u(13:15)-u(31:33)))./(norm((u(13:15)-u(31:33))).^3) + (G*m4*(u(19:21)-u(31:33)))./(norm((u(19:21)-u(31:33))).^3) + (G*m5*(u(25:27)-u(31:33)))./(norm((u(25:27)-u(31:33))).^3)] ;

%%%%%% AB2 Numerical Analysis %%%%%%%

tic

tk = 0; %starting time

%initialize iterates for AB2
u_AB2_k = u0;
u_AB2_km1 = u0;

%initialize vector that stores approximate soln at various times
u_AB2_approx = zeros( 36,T/dt );

%Initialize vector that stores the magnitude of velocities at T
u_AB2_vel = zeros(6,T/dt);

%advance to final time T
for i = 1 : T/dt
    
    if i < 2
        
        %advance with Heun's for 1st time step
        u_AB2_approx(:,i) = (u_AB2_k + 1/2*dt*(f(u_AB2_k) + f(u_AB2_k+dt*f(u_AB2_k))));
        
    else
        
        u_AB2_approx(:,i) = u_AB2_k + 1/2*dt*(-1*f(u_AB2_km1) + 3*f(u_AB2_k));
        
    end
    
    %update iterates
    tk = tk + dt;
    u_AB2_km1 = u_AB2_k;
    u_AB2_k = u_AB2_approx(:,i);
    
    %Calculating and storing the velocities of each body at T
    u_AB2_vel(1,i) = norm(norm(u_AB2_approx(4:6,i)));
    u_AB2_vel(2,i) = norm(norm(u_AB2_approx(10:12,i)));
    u_AB2_vel(3,i) = norm(norm(u_AB2_approx(16:18,i)));
    u_AB2_vel(4,i) = norm(norm(u_AB2_approx(22:24,i)));
    u_AB2_vel(5,i) = norm(norm(u_AB2_approx(28:30,i)));
    u_AB2_vel(6,i) = norm(norm(u_AB2_approx(34:36,i)));
        
end

t_AB2 = toc

% Gathering N-Body Information for 249 years in future from Jan 1, 2011 (T)
Julian2 = juliandate(2260,12,31);

% Position and Velocity Vectors at future position
[r1, v1] = planetEphemeris(Julian2,'Sun', 'Sun', '432t', 'km');
[r2, v2] = planetEphemeris(Julian2,'Sun','Jupiter','432t','km');
[r3, v3] = planetEphemeris(Julian2,'Sun','Saturn','432t','km');
[r4, v4] = planetEphemeris(Julian2,'Sun','Uranus','432t','km');
[r5, v5] = planetEphemeris(Julian2,'Sun','Neptune','432t','km');
[r6, v6] = planetEphemeris(Julian2,'Sun','Pluto','432t','km');

% Calculating the Average, Minimum and Maximum Velocity of Each Planet
AB2_vel_avg = mean(u_AB2_vel,2);
AB2_vel_max = max(u_AB2_vel,[],2);
AB2_vel_min = min(u_AB2_vel,[],2);

% Assigning Variable Names
Planetary_Bodies = {'Jupiter';'Saturn';'Uranus';'Neptune';'Pluto'};
Calculated_Position_AB2 = [u_AB2_approx(7,T/dt), u_AB2_approx(8,T/dt), u_AB2_approx(9,T/dt); u_AB2_approx(13,T/dt), u_AB2_approx(14,T/dt), u_AB2_approx(15,T/dt); u_AB2_approx(19,T/dt), u_AB2_approx(20,T/dt), u_AB2_approx(21,T/dt); u_AB2_approx(25,T/dt), u_AB2_approx(26,T/dt), u_AB2_approx(27,T/dt); u_AB2_approx(31,T/dt), u_AB2_approx(32,T/dt), u_AB2_approx(33,T/dt)];
JulianDate_Position_AB2 = 10^3 * [r2;r3;r4;r5;r6];
% Calculated_Position = [norm([u_AB2_approx(7,T/dt), u_AB2_approx(8,T/dt), u_AB2_approx(9,T/dt)]); norm([u_AB2_approx(13,T/dt), u_AB2_approx(14,T/dt), u_AB2_approx(15,T/dt)]); norm([u_AB2_approx(19,T/dt), u_AB2_approx(20,T/dt), u_AB2_approx(21,T/dt)]); norm([u_AB2_approx(25,T/dt), u_AB2_approx(26,T/dt), u_AB2_approx(27,T/dt)]); norm([u_AB2_approx(31,T/dt), u_AB2_approx(32,T/dt), u_AB2_approx(33,T/dt)])];
% JulianDate_Position = 10^3 * [norm(r2);norm(r3);norm(r4);norm(r5);norm(r6)];
Calculated_Final_V_AB2 = [norm([u_AB2_approx(10,T/dt), u_AB2_approx(11,T/dt), u_AB2_approx(12,T/dt)]); norm([u_AB2_approx(16,T/dt), u_AB2_approx(17,T/dt), u_AB2_approx(18,T/dt)]); norm([u_AB2_approx(22,T/dt), u_AB2_approx(23,T/dt), u_AB2_approx(24,T/dt)]); norm([u_AB2_approx(28,T/dt), u_AB2_approx(29,T/dt), u_AB2_approx(30,T/dt)]); norm([u_AB2_approx(34,T/dt), u_AB2_approx(35,T/dt), u_AB2_approx(36,T/dt)])];
JulianDate_Final_V_AB2 = 10^3 * [norm(v2);norm(v3);norm(v4);norm(v5);norm(v6)];
CalculatedMeanVelocity_AB2 = [AB2_vel_avg(2); AB2_vel_avg(3); AB2_vel_avg(4); AB2_vel_avg(5); AB2_vel_avg(6)];
KnownMeanVelocity_AB2 = 10^3 * [Jupiter_Mean; Saturn_Mean; Uranus_Mean; Neptune_Mean; Pluto_Mean];

% Defining Tables
Table1 = table(Planetary_Bodies, Calculated_Position_AB2, JulianDate_Position_AB2);
Table2 = table(Planetary_Bodies, Calculated_Final_V_AB2, JulianDate_Final_V_AB2, CalculatedMeanVelocity_AB2, KnownMeanVelocity_AB2);

% Get the table in string form.
TString = evalc('disp(Table1)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(100)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Get the table in string form.
TString = evalc('disp(Table2)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(200)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Plotting the Orbits Estimated
figure(1)
scatter4(0,0,0,40,'filled')
grid on
hold on
plot3(u_AB2_approx(7,:), u_AB2_approx(8,:), u_AB2_approx(9,:),'LineWidth', 2)
scatter4(u_AB2_approx(7,1), u_AB2_approx(8,1), u_AB2_approx(9,1),40,'filled')
plot3(u_AB2_approx(13,:), u_AB2_approx(14,:), u_AB2_approx(15,:),'LineWidth', 2)
scatter4(u_AB2_approx(13,1), u_AB2_approx(14,1), u_AB2_approx(15,1),40,'filled')
plot3(u_AB2_approx(19,:), u_AB2_approx(20,:), u_AB2_approx(21,:),'LineWidth', 2)
scatter4(u_AB2_approx(19,1), u_AB2_approx(20,1), u_AB2_approx(21,1),40,'filled')
plot3(u_AB2_approx(25,:), u_AB2_approx(26,:), u_AB2_approx(27,:),'LineWidth', 2)
scatter4(u_AB2_approx(25,1), u_AB2_approx(26,1), u_AB2_approx(27,1),40,'filled')
plot3(u_AB2_approx(31,:), u_AB2_approx(32,:), u_AB2_approx(33,:),'LineWidth', 2)
scatter4(u_AB2_approx(31,1), u_AB2_approx(32,1), u_AB2_approx(33,1),40,'filled')
title('Outer Solar System (Adams–Bashforth Method)')
legend('Sun','Jupiter Orbit','Jupiter','Saturn Orbit','Saturn','Uranus Orbit','Uranus','Neptune Orbit','Neptune','Pluto Orbit','Pluto')



%%%%%% RK4 Numerical Analysis %%%%%%%

tk = 0; %starting time

tic

%initialize iterates for Huens
u_rk_k = u0;

%initialize vector that stores approximate soln at various times
u_rk_approx = zeros( 36,T/dt );

%Initialize vector that stores the magnitude of velocities at T
u_rk_vel = zeros(6,T/dt);

%advance to final time T
for i = 1 : T/dt
    
    %advance with RK4 for 1st time step
    y1 = f(u_rk_k);
    y2 = f(u_rk_k + 1/2*dt*y1);
    y3 = f(u_rk_k + (1/2)*dt*y2);
    y4 = f(u_rk_k + dt*y3);
    u_rk_approx(:,i) = u_rk_k + 1/6*dt*(y1 + 2*y2 + 2*y3 + y4);
    
    %update iterates
    tk = tk + dt; 
    u_rk_k = u_rk_approx(:,i);
    
    %Calculating and storing the velocities of each body at T
    u_rk_vel(1,i) = norm(norm(u_rk_approx(4:6,i)));
    u_rk_vel(2,i) = norm(norm(u_rk_approx(10:12,i)));
    u_rk_vel(3,i) = norm(norm(u_rk_approx(16:18,i)));
    u_rk_vel(4,i) = norm(norm(u_rk_approx(22:24,i)));
    u_rk_vel(5,i) = norm(norm(u_rk_approx(28:30,i)));
    u_rk_vel(6,i) = norm(norm(u_rk_approx(34:36,i)));   
    
end

t_r = toc

% Gathering N-Body Information for 249 years in future from Jan 1, 2011 (T)
Julian2 = juliandate(2260,12,31);

% Position and Velocity Vectors at future position
[r1, v1] = planetEphemeris(Julian2,'Sun', 'Sun', '432t', 'km');
[r2, v2] = planetEphemeris(Julian2,'Sun','Jupiter','432t','km');
[r3, v3] = planetEphemeris(Julian2,'Sun','Saturn','432t','km');
[r4, v4] = planetEphemeris(Julian2,'Sun','Uranus','432t','km');
[r5, v5] = planetEphemeris(Julian2,'Sun','Neptune','432t','km');
[r6, v6] = planetEphemeris(Julian2,'Sun','Pluto','432t','km');

% Calculating the Average Velocity of Each Planet
rk_vel_avg = mean(u_rk_vel,2);
rk_vel_max = max(u_rk_vel,[],2);
rk_vel_min = min(u_rk_vel,[],2);

% Assigning Variable Names
Planetary_Bodies = {'Jupiter';'Saturn';'Uranus';'Neptune';'Pluto'};
Calculated_Position_RK4 = [u_rk_approx(7,T/dt), u_rk_approx(8,T/dt), u_rk_approx(9,T/dt); u_rk_approx(13,T/dt), u_rk_approx(14,T/dt), u_rk_approx(15,T/dt); u_rk_approx(19,T/dt), u_rk_approx(20,T/dt), u_rk_approx(21,T/dt); u_rk_approx(25,T/dt), u_rk_approx(26,T/dt), u_rk_approx(27,T/dt); u_rk_approx(31,T/dt), u_rk_approx(32,T/dt), u_rk_approx(33,T/dt)];
JulianDate_Position_RK4 = 10^3 * [r2;r3;r4;r5;r6];
% Calculated_Position = [norm([u_rk_approx(7,T/dt), u_rk_approx(8,T/dt), u_rk_approx(9,T/dt)]); norm([u_rk_approx(13,T/dt), u_rk_approx(14,T/dt), u_rk_approx(15,T/dt)]); norm([u_rk_approx(19,T/dt), u_rk_approx(20,T/dt), u_rk_approx(21,T/dt)]); norm([u_rk_approx(25,T/dt), u_rk_approx(26,T/dt), u_rk_approx(27,T/dt)]); norm([u_rk_approx(31,T/dt), u_rk_approx(32,T/dt), u_rk_approx(33,T/dt)])];
% JulianDate_Position = 10^3 * [norm(r2);norm(r3);norm(r4);norm(r5);norm(r6)];
Calculated_Final_V_RK4 = [norm([u_rk_approx(10,T/dt), u_rk_approx(11,T/dt), u_rk_approx(12,T/dt)]); norm([u_rk_approx(16,T/dt), u_rk_approx(17,T/dt), u_rk_approx(18,T/dt)]); norm([u_rk_approx(22,T/dt), u_rk_approx(23,T/dt), u_rk_approx(24,T/dt)]); norm([u_rk_approx(28,T/dt), u_rk_approx(29,T/dt), u_rk_approx(30,T/dt)]); norm([u_rk_approx(34,T/dt), u_rk_approx(35,T/dt), u_rk_approx(36,T/dt)])];
JulianDate_Final_V_RK4 = 10^3 * [norm(v2);norm(v3);norm(v4);norm(v5);norm(v6)];
CalculatedMeanVelocity_RK4 = [rk_vel_avg(2); rk_vel_avg(3); rk_vel_avg(4); rk_vel_avg(5); rk_vel_avg(6)];
KnownMeanVelocity_RK4 = 10^3 * [Jupiter_Mean; Saturn_Mean; Uranus_Mean; Neptune_Mean; Pluto_Mean];

% Defining Tables
Table1 = table(Planetary_Bodies, Calculated_Position_RK4, JulianDate_Position_RK4);
Table2 = table(Planetary_Bodies, Calculated_Final_V_RK4, JulianDate_Final_V_RK4, CalculatedMeanVelocity_RK4, KnownMeanVelocity_RK4);

% Get the table in string form.
TString = evalc('disp(Table1)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(300)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Get the table in string form.
TString = evalc('disp(Table2)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure(400)
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',...
    FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

% Plotting the Orbits Estimated
figure(2)
scatter4(0,0,0,40,'filled')
grid on
hold on
plot3(u_rk_approx(7,:), u_rk_approx(8,:), u_rk_approx(9,:),'LineWidth', 2)
scatter4(u_rk_approx(7,1), u_rk_approx(8,1), u_rk_approx(9,1),40,'filled')
plot3(u_rk_approx(13,:), u_rk_approx(14,:), u_rk_approx(15,:),'LineWidth', 2)
scatter4(u_rk_approx(13,1), u_rk_approx(14,1), u_rk_approx(15,1),40,'filled')
plot3(u_rk_approx(19,:), u_rk_approx(20,:), u_rk_approx(21,:),'LineWidth', 2)
scatter4(u_rk_approx(19,1), u_rk_approx(20,1), u_rk_approx(21,1),40,'filled')
plot3(u_rk_approx(25,:), u_rk_approx(26,:), u_rk_approx(27,:),'LineWidth', 2)
scatter4(u_rk_approx(25,1), u_rk_approx(26,1), u_rk_approx(27,1),40,'filled')
plot3(u_rk_approx(31,:), u_rk_approx(32,:), u_rk_approx(33,:),'LineWidth', 2)
scatter4(u_rk_approx(31,1), u_rk_approx(32,1), u_rk_approx(33,1),40,'filled')
title('Outer Solar System (Runge-Kutta 4 Method)')
legend('Sun','Jupiter Orbit','Jupiter','Saturn Orbit','Saturn','Uranus Orbit','Uranus','Neptune Orbit','Neptune','Pluto Orbit','Pluto')
    

%% Convergence Test

% AB2
%--simulation params
% dtvect = [60*60*24*6.25, 60*60*24*1.25, 60*60*24*0.25, 60*60*24*0.05, 60*60*24*0.01, 60*60*24*0.002]; % For Inner Orbits
dtvect = 60*60*24 * [781.25, 156.25, 31.25, 6.25, 1.25, 0.25]; % For Outer Orbits
%--

%initialize vector that stores approximate soln at T for various dt
u_AB2_keep = zeros( 36,length( dtvect ) );
tic
%advance in time
for j = 1 : length(dtvect)
    
    %current dt
    dt = dtvect(j);
    
    tk = t0; %initialize time iterate
    
    %initialize iterates for various methods
    u_AB2_k = u0;
    u_AB2_km1 = u0;
    
    
    %run to final time T
    for jj = 1 : T/dt
        
        %AB2
        if jj < 2
            
            %advance with Heun's for 1st time step
            u_AB2_kp1 = (u_AB2_k + 1/2*dt*(f(u_AB2_k) + f(u_AB2_k + dt*f(u_AB2_k))));
            
        else
            u_AB2_kp1 = u_AB2_k + 1/2*dt*(-1*f(u_AB2_km1) + 3*f(u_AB2_k));
            
        end
        
        %update iterates
        u_AB2_km1 = u_AB2_k;
        u_AB2_k = u_AB2_kp1;
        
        %update time
        tk = tk + dt;
        
        
    end
    
    %store soln at T
    u_AB2_keep( :,j ) = u_AB2_kp1;
    
end

%compute difference between solution at smallest dt and the other dts

%initialize vector
u_AB2_diff = zeros( length( dtvect )-1,1 );

for j = 1 : length( dtvect )-1
    
    %AB2
    u_AB2_diff(j) = norm( u_AB2_keep(:,j) - u_AB2_keep(:,end) )/ ...
        norm( u_AB2_keep(:,end) );
    
end
toc
figure(50)
%AB2
loglog( dtvect(1:end-1), u_AB2_diff, '.-', 'markersize', 20, 'linewidth', 2 )

set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )

leg = legend( 'AB2' );
set( leg, 'fontsize', 12, 'interpreter', 'latex', 'location', 'southeast' )
xlabel('$\Delta t$', 'fontsize', 12, 'interpreter' , 'latex')
title('$\frac{||u_{\Delta t} - u_{2.5\times10^{-4}}||}{|| u_{2.5\times10^{-4}} ||}$', 'fontsize', 28,  'interpreter' , 'latex' )
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'ConvergencePlots1';
print( '-dpng', svnm, '-r300' )


% RK4
%--simulation params
% dtvect = [60*60*24*6.25, 60*60*24*1.25, 60*60*24*0.25, 60*60*24*0.05, 60*60*24*0.01, 60*60*24*0.002]; % For Inner Orbits
dtvect = 60*60*24 * [781.25, 156.25, 31.25, 6.25, 1.25, 0.25]; % For Outer Orbits
%--

%initialize vector that stores approximate soln at T for various dt
u_r_keep = zeros( 36,length( dtvect ) );
tic
%advance in time
for j = 1 : length(dtvect)
    
    %current dt
    dt = dtvect(j);
    
    tk = t0; %initialize time iterate
    
    %initialize iterates for various methods
    u_r_k = u0;
    
    
    %run to final time T
    for jj = 1 : T/dt
        
        %RK4
        y1 = f(u_r_k);
        y2 = f(u_r_k + 1/2*dt*y1);
        y3 = f(u_r_k + (1/2)*dt*y2);
        y4 = f(u_r_k + dt*y3);
        u_r_kp1 = u_r_k + 1/6*dt*(y1 + 2*y2 + 2*y3 + y4);
        
        %update iterates
        u_r_k = u_r_kp1;
        
        %update time
        tk = tk + dt;
        
        
    end
    
    %store soln at T
    u_r_keep( :,j ) = u_r_kp1;
    
end

%compute difference between solution at smallest dt and the other dts

%initialize vector
u_r_diff = zeros( length( dtvect )-1,1 );

for j = 1 : length( dtvect )-1
    
    %RK4
    u_r_diff(j) = norm( u_r_keep(:,j) - u_r_keep(:,end) )/ ...
        norm( u_r_keep(:,end) );
    
end
toc
figure(50), hold on
%RK4
loglog( dtvect(1:end-1), u_r_diff, '.-', 'markersize', 20, 'linewidth', 2 )

set( gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex' )

leg = legend( 'AB2', 'RK4' );
set( leg, 'fontsize', 16, 'interpreter', 'latex', 'location', 'southeast' )
xlabel('$\Delta t$', 'fontsize', 12, 'interpreter' , 'latex')
title('$\frac{||u_{\Delta t} - u_{2.5\times10^{-4}}||}{|| u_{2.5\times10^{-4}} ||}$', 'fontsize', 28,  'interpreter' , 'latex' )
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'ConvergencePlots2';
print( '-dpng', svnm, '-r300' )


% Euler
%--simulation params
% dtvect = [60*60*24*6.25, 60*60*24*1.25, 60*60*24*0.25, 60*60*24*0.05, 60*60*24*0.01, 60*60*24*0.002]; % For Inner Orbits
dtvect = 60*60*24 * [781.25, 156.25, 31.25, 6.25, 1.25, 0.25]; % For Outer Orbits
%--

%initialize vector that stores approximate soln at T for various dt
u_eu_approx = zeros( 36,length( dtvect ) );

%advance in time
for j = 1 : length(dtvect)
    
    %current dt
    dt = dtvect(j);
    
    tk = t0; %initialize time iterate
    
    %initialize iterates for various methods
    u_eu_k = u0;
    
    
    %run to final time T
    for jj = 1 : T/dt
        
        %Euler
        u_eu_kp1 = u_eu_k + dt*f(u_eu_k);
        
        %update iterates
        u_eu_k = u_eu_kp1;
        
        %update time
        tk = tk + dt;
        
        
    end
    
    %store soln at T
    u_eu_approx( :,j ) = u_eu_kp1;
    
end

%compute difference between solution at smallest dt and the other dts

%initialize vector
u_eu_diff = zeros( length( dtvect )-1,1 );

for j = 1 : length( dtvect )-1
    
    %Euler
    u_eu_diff(j) = norm( u_eu_approx(:,j) - u_eu_approx(:,end) )/ ...
        norm( u_eu_approx(:,end) );
    
end

figure(50), hold on
%Euler
loglog( dtvect(1:end-1), u_eu_diff, '.-', 'markersize', 20, 'linewidth', 2 )

set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )

leg = legend( 'AB2', 'RK4', 'Euler' );
set( leg, 'fontsize', 12, 'interpreter', 'latex', 'location', 'southeast' )
xlabel('$\Delta t$', 'fontsize', 12, 'interpreter' , 'latex')
title('$\frac{||u_{\Delta t} - u_{2.5\times10^{-4}}||}{|| u_{2.5\times10^{-4}} ||}$', 'fontsize', 28,  'interpreter' , 'latex' )
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'error_q2_Euler';
print( '-dpng', svnm, '-r200' )


% Heun's
%--simulation params
% dtvect = [60*60*24*6.25, 60*60*24*1.25, 60*60*24*0.25, 60*60*24*0.05, 60*60*24*0.01, 60*60*24*0.002]; % For Inner Orbits
dtvect = 60*60*24 * [781.25, 156.25, 31.25, 6.25, 1.25, 0.25]; % For Outer Orbits
%--

%initialize vector that stores approximate soln at T for various dt
u_h_keep = zeros( 36,length( dtvect ) );

%advance in time
for j = 1 : length(dtvect)
    
    %current dt
    dt = dtvect(j);
    
    tk = t0; %initialize time iterate
    
    %initialize iterates for various methods
    u_h_k = u0;
    
    
    %run to final time T
    for jj = 1 : T/dt
        
        %Huens
        u_h_kp1 = u_h_k + 1/2*dt*(f(u_h_k) + f(u_h_k + dt*f(u_h_k)));
        
        %update iterates
        u_h_k = u_h_kp1;
        
        %update time
        tk = tk + dt;
        
        
    end
    
    %store soln at T
    u_h_keep( :,j ) = u_h_kp1;
    
end

%compute difference between solution at smallest dt and the other dts

%initialize vector
u_h_diff = zeros( length( dtvect )-1,1 );

for j = 1 : length( dtvect )-1
    
    %Huens
    u_h_diff(j) = norm( u_h_keep(:,j) - u_h_keep(:,end) )/ ...
        norm( u_h_keep(:,end) );
    
end

figure(50), hold on
%Heun's
loglog( dtvect(1:end-1), u_h_diff, '.-', 'markersize', 20, 'linewidth', 2 )

set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )

leg = legend( 'AB2', 'RK4', 'Euler', 'Huens' );
set( leg, 'fontsize', 12, 'interpreter', 'latex', 'location', 'southeast' )
xlabel('$\Delta t$', 'fontsize', 12, 'interpreter' , 'latex')
title('$\frac{||u_{\Delta t} - u_{2.5\times10^{-4}}||}{|| u_{2.5\times10^{-4}} ||}$', 'fontsize', 28,  'interpreter' , 'latex' )
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'error_q2_Huens';
print( '-dpng', svnm, '-r200' )


%% Error Plots

% % Gathering N-Body Information for 687 days in future from Jan 1, 2011 (T)
% Julian2 = juliandate(2012,11,19);

% Gathering N-Body Information for 249 years in future from Jan 1, 2011 (T)
Julian2 = juliandate(2260,12,31);

% % Position and Velocity Vectors at future position
% [r1, v1] = planetEphemeris(Julian2,'Sun', 'Sun', '432t', 'km');
% [r2, v2] = planetEphemeris(Julian2,'Sun','Mercury','432t','km');
% [r3, v3] = planetEphemeris(Julian2,'Sun','Venus','432t','km');
% [r4, v4] = planetEphemeris(Julian2,'Sun','Earth','432t','km');
% [r5, v5] = planetEphemeris(Julian2,'Sun','Mars','432t','km');
% [r6, v6] = planetEphemeris(Julian2,'Sun','Moon','432t','km');


% Position and Velocity Vectors at future position
[r1, v1] = planetEphemeris(Julian2,'Sun', 'Sun', '432t', 'km');
[r2, v2] = planetEphemeris(Julian2,'Sun','Jupiter','432t','km');
[r3, v3] = planetEphemeris(Julian2,'Sun','Saturn','432t','km');
[r4, v4] = planetEphemeris(Julian2,'Sun','Uranus','432t','km');
[r5, v5] = planetEphemeris(Julian2,'Sun','Neptune','432t','km');
[r6, v6] = planetEphemeris(Julian2,'Sun','Pluto','432t','km');

[u_real] = 10^3 * [r1'; v1'; r2'; v2'; r3'; v3'; r4'; v4'; r5'; v5'; r6'; v6'];

% AB2
%--simulation params
% dtvect = [60*60*24*6.25, 60*60*24*1.25, 60*60*24*0.25, 60*60*24*0.05, 60*60*24*0.01, 60*60*24*0.002];
dtvect = 60*60*24 * [781.25, 156.25, 31.25, 6.25, 1.25, 0.25];
%--

%initialize vector that stores approximate soln at T for various dt
u_AB2_keep = zeros( 36,length( dtvect ) );
tic
%advance in time
for j = 1 : length(dtvect)
    
    %current dt
    dt = dtvect(j);
    
    tk = t0; %initialize time iterate
    
    %initialize iterates for various methods
    u_AB2_k = u0;
    u_AB2_km1 = u0;
    
    
    %run to final time T
    for jj = 1 : T/dt
        
        %AB2
        if jj < 2
            
            %advance with Heun's for 1st time step
            u_AB2_kp1 = (u_AB2_k + 1/2*dt*(f(u_AB2_k) + f(u_AB2_k + dt*f(u_AB2_k))));
            
        else
            u_AB2_kp1 = u_AB2_k + 1/2*dt*(-1*f(u_AB2_km1) + 3*f(u_AB2_k));
            
        end
        
        %update iterates
        u_AB2_km1 = u_AB2_k;
        u_AB2_k = u_AB2_kp1;
        
        %update time
        tk = tk + dt;
        
        
    end
    
    %store soln at T
    u_AB2_keep( :,j ) = u_AB2_kp1;
    
end

%compute difference between solution at smallest dt and the other dts

%initialize vector
u_AB2_error = zeros( length( dtvect ),1 );

for j = 1 : length( dtvect )
    
    %AB2
    u_AB2_error(j) = norm( u_AB2_keep(:,j) - u_real(:) )/ ...
        norm( u_real(:) );
    
end
toc
figure(50)
%AB2
loglog( dtvect(1:end), u_AB2_error, '.-', 'markersize', 20, 'linewidth', 2 )

set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )

leg = legend( 'AB2' );
set( leg, 'fontsize', 12, 'interpreter', 'latex', 'location', 'southeast' )
xlabel('$\Delta t$', 'fontsize', 12, 'interpreter' , 'latex')
ylabel('Error %', 'fontsize', 12, 'interpreter' , 'latex')
title('Error Plot', 'fontsize', 28,  'interpreter' , 'latex' )
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'ErrorPlots1';
print( '-dpng', svnm, '-r300' )



% RK4
%--simulation params
% dtvect = [60*60*24*6.25, 60*60*24*1.25, 60*60*24*0.25, 60*60*24*0.05, 60*60*24*0.01, 60*60*24*0.002];
dtvect = 60*60*24 * [781.25, 156.25, 31.25, 6.25, 1.25, 0.25];
%--

%initialize vector that stores approximate soln at T for various dt
u_r_keep = zeros( 36,length( dtvect ) );
tic
%advance in time
for j = 1 : length(dtvect)
    
    %current dt
    dt = dtvect(j);
    
    tk = t0; %initialize time iterate
    
    %initialize iterates for various methods
    u_r_k = u0;
    
    
    %run to final time T
    for jj = 1 : T/dt
        
        %RK4
        y1 = f(u_r_k);
        y2 = f(u_r_k + 1/2*dt*y1);
        y3 = f(u_r_k + (1/2)*dt*y2);
        y4 = f(u_r_k + dt*y3);
        u_r_kp1 = u_r_k + 1/6*dt*(y1 + 2*y2 + 2*y3 + y4);
        
        %update iterates
        u_r_k = u_r_kp1;
        
        %update time
        tk = tk + dt;
        
        
    end
    
    %store soln at T
    u_r_keep( :,j ) = u_r_kp1;
    
end

%compute the error

%initialize vector
u_r_error = zeros( length( dtvect ),1 );

for j = 1 : length( dtvect )
    
    %RK4
    u_r_error(j) = (norm( u_r_keep(:,j) - u_real(:) )/ ...
        norm( u_real(:) ));
    
end
toc
figure(50), hold on
%RK4
loglog( dtvect(1:end), u_r_error, '.-', 'markersize', 20, 'linewidth', 2 )

set( gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex' )

leg = legend( 'AB2', 'RK4' );
set( leg, 'fontsize', 16, 'interpreter', 'latex', 'location', 'southeast' )
xlabel('$\Delta t$', 'fontsize', 12, 'interpreter' , 'latex')
ylabel('Error %', 'fontsize', 12, 'interpreter' , 'latex')
title('Error Plot', 'fontsize', 28,  'interpreter' , 'latex' )
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'ErrorPlots2';
print( '-dpng', svnm, '-r300' )