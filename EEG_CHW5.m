clc
clear all
close all

%% Part 1
load('ElecPosXYZ') ;

%Forward Matrix
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;

% Plotting dipoles' locations
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:));
title("Location of dipoles");
xlabel("X [cm]");
ylabel("Y [cm]");
zlabel("Z [cm]");

% Saving gain matrix
save('GainMat.mat',"GainMat");

%% Part 2
% Finding electrodes' locations
electrode_name = string(size(ElecPos));
for i = 1:length(ElecPos)
    electrode_X(i) = ElecPos{i}.XYZ(1) * ModelParams.R(3);
    electrode_Y(i) = ElecPos{i}.XYZ(2) * ModelParams.R(3);
    electrode_Z(i) = ElecPos{i}.XYZ(3) * ModelParams.R(3);
    electrode_name(i) = ElecPos{i}.Name;
end

% Plotting dipoles' locations
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:));
hold on
title("Location of dipoles & electrodes");
xlabel("X [cm]");
ylabel("Y [cm]");
zlabel("Z [cm]");
hold on

% Plotting electrodes' locations
scatter3(electrode_X,electrode_Y,electrode_Z,'filled');
% Electrodes' names
text(electrode_X,electrode_Y,electrode_Z,electrode_name);

%% Part 3
% Random dipole
dip_r = sqrt(sum(abs(LocMat).^2,1)); % Dipoles' radius

surface_dips = find(dip_r == max(dip_r)); % On-surface dipoles
deep_dips = find(LocMat(1,:)>=-3 & LocMat(1,:)<=3 & LocMat(2,:)>=-3 & LocMat(2,:)<=3 & LocMat(3,:)>=2 & LocMat(3,:)<=4); %In the deep

%%%%%%%%%%%%%%%% Surface dipole:A=1   Deep dipole:A=2 %%%%%%%%%%%%%%%%
A = 2;
switch A
    case 1
        dip = surface_dips(randi(length(surface_dips))); % Random dipole on surface
    case 2
        dip = deep_dips(randi(length(deep_dips))); % Random dipole in the deep
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting the location of the dipole:
x = LocMat(1,dip);
y = LocMat(2,dip);
z = LocMat(3,dip);
R = sqrt(x^2 + y^2 + z^2);
delta_x = x/R;
delta_y = y/R;
delta_z = z/R;

% Plotting dipoles' locations
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:));
hold on
title("Location of dipoles & electrodes with a random dipole");
xlabel("X [cm]");
ylabel("Y [cm]");
zlabel("Z [cm]");
hold on

% Plotting electrodes' locations
scatter3(electrode_X,electrode_Y,electrode_Z,'filled');

% Electrodes' names
text(electrode_X,electrode_Y,electrode_Z,electrode_name);

% Adding the momentum vector:
quiver3(x,y,z,delta_x,delta_y,delta_z,'Green','LineWidth',3);

%% Part 4
load('Matlab\Interictal.mat');

% Random source for selected dipole
source_number = randi(size(Interictal,1));
source = Interictal(source_number,:);

% Making Q matrix
Q = zeros(3,size(source,2));
Q(1,:) = source*delta_x;
Q(2,:) = source*delta_y;
Q(3,:) = source*delta_z;

% Finding potentials on electrodes
M = GainMat(:,dip*3-2:dip*3)*Q;

% Plotting potentials
offset = max(max(abs(M)))/3 ;
disp_eeg(M,offset,250,electrode_name);
xlim('tight');
grid minor
title("Potentials from random dipole");

%% Part 5
% Defining peak thresholds
peaks_thresh = mean(M,2) + 2*std(M,[],2);

% Finding peaks
for i = 1:size(M,1)
    [peaks, locations] = findpeaks(M(i,:),'MinPeakHeight', peaks_thresh(i));
    peak_locations(i,1:size(locations,2)) = locations;
end

% Mean of spiky windows
for i = 1:size(peak_locations,1) %For each electrode
    for j = 1:size(peak_locations,2) %For each peak
        electrode_mean_potential(i) = mean(M(i,peak_locations(i,:)-3:peak_locations(i,:)+3));
    end
end

% Plotting dipoles' locations
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:));
hold on
title("Location of electrodes & their potentials");
xlabel("X [cm]");
ylabel("Y [cm]");
zlabel("Z [cm]");
hold on

% Normalized color for each electrode
color = normalize(electrode_mean_potential,2,"norm");

% Plotting electrodes' locations
scatter3(electrode_X,electrode_Y,electrode_Z,[],color,"filled");
colorbar

% Electrodes' names
text(electrode_X,electrode_Y,electrode_Z,electrode_name);

%% Part 6
figure;
Display_Potential_3D(ModelParams.R(3),electrode_mean_potential);
title("3D potential Display of Dipole number " + dip + " and source number " + source_number);
xlabel("X [cm]"); ylabel("Y [cm]"); zlabel("Z [cm]");

%% Part 7
% MNE
alpha = 1; % Tikhonov
N = size(M,1);

%Computing predicted momentum matrix
Q_predicted_MNE = GainMat.' * inv(GainMat*(GainMat.') + alpha*eye(N)) * M;

% wMNE
p = size(LocMat,2);

% defining Omega & W
omega = zeros(p);
for i = 1:p
    for j = 1:N
        omega(i,i) = omega(i,i) + GainMat(j,i)*GainMat(j,i).';
    end
end

W = kron(omega,eye(3));

%Computing predicted momentum matrix
Q_predicted_wMNE = inv(W.' * W) * GainMat.' * inv(GainMat * inv(W.' * W) * (GainMat.') + alpha*eye(N)) * M;

%% Part 8
% Calculating amplitudes of each dipole
squared_Q_predicted_MNE = Q_predicted_MNE.^2;
squared_Q_predicted_wMNE = Q_predicted_wMNE.^2;

% Arrays of momentum amplitudes
amolitudes_MNE = zeros(1,p); 
amolitudes_wMNE = zeros(1,p); 
for i = 1:p
    amolitudes_MNE(i) = sum(sum(squared_Q_predicted_MNE((3*i-2):(3*i),:)));
    amolitudes_wMNE(i) = sum(sum(squared_Q_predicted_wMNE((3*i-2):(3*i),:)));
end

% Predicting the dipoles (which maximizes the momentum amplitude)
dipode_predicted_MNE = find(amolitudes_MNE == max(amolitudes_MNE));
dipode_predicted_wMNE = find(amolitudes_wMNE == max(amolitudes_wMNE));

%Setting the locations of the predicted dipoles:
x_MNE = LocMat(1,dipode_predicted_MNE);
y_MNE = LocMat(2,dipode_predicted_MNE);
z_MNE = LocMat(3,dipode_predicted_MNE);
R_MNE = sqrt(x_MNE^2 + y_MNE^2 + z_MNE^2);
delta_x_MNE = x_MNE/R_MNE;
delta_y_MNE = y_MNE/R_MNE;
delta_z_MNE = z_MNE/R_MNE;

x_wMNE = LocMat(1,dipode_predicted_wMNE);
y_wMNE = LocMat(2,dipode_predicted_wMNE);
z_wMNE = LocMat(3,dipode_predicted_wMNE);
R_wMNE = sqrt(x_wMNE^2 + y_wMNE^2 + z_wMNE^2);
delta_x_wMNE = x_wMNE/R_wMNE;
delta_y_wMNE = y_wMNE/R_wMNE;
delta_z_wMNE = z_wMNE/R_wMNE;

% Plotting dipoles' locations
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:));
hold on
title("Location of Main and Predicted dipoles by MNE & wMNE");
xlabel("X [cm]");
ylabel("Y [cm]");
zlabel("Z [cm]");
hold on

% Plotting electrodes' locations
scatter3(electrode_X,electrode_Y,electrode_Z,'filled');

% Electrodes' names
text(electrode_X,electrode_Y,electrode_Z,electrode_name);

%Adding the momentum vectors (Main dipole)
quiver3(x,y,z,delta_x,delta_y,delta_z,'g','LineWidth',3);
%Adding the momentum vectors (Predicted dipoles)
quiver3(x_MNE,y_MNE,z_MNE,delta_x_MNE,delta_y_MNE,delta_z_MNE,'r','LineWidth',3);
quiver3(x_wMNE,y_wMNE,z_wMNE,delta_x_wMNE,delta_y_wMNE,delta_z_wMNE,'b','LineWidth',3);

%% Part 9
% Distance error
distance_error_MNE = sqrt((x-x_MNE)^2 + (y-y_MNE)^2 + (z-z_MNE)^2);
distance_error_wMNE = sqrt((x-x_wMNE)^2 + (y-y_wMNE)^2 + (z-z_wMNE)^2);

% Moment angle error
[fi, theta, radius] = cart2sph(x,y,z);
[fi_MNE, theta_MNE, radius_MNE] = cart2sph(x_MNE,y_MNE,z_MNE);
[fi_wMNE, theta_wMNE, radius_wMNE] = cart2sph(x_wMNE,y_wMNE,z_wMNE);

fi_error_MNE = fi_MNE - fi;
fi_error_wMNE = fi_wMNE - fi;
theta_error_MNE = theta_MNE - theta;
theta_error_wMNE = theta_wMNE - theta;

% Displaying errors:
disp("Dipole distance error for MNE = " + distance_error_MNE);
disp("Dipole distance error for wMNE = " + distance_error_wMNE);
disp("Dipole fi error for MNE = " + fi_error_MNE);
disp("Dipole fi error for wMNE = " + fi_error_wMNE);
disp("Dipole theta error for MNE = " + theta_error_MNE);
disp("Dipole theta error for wMNE = " + theta_error_wMNE);

%% Part 10
% For considering a deep dipole change A to 2

%% Part 13
% Random dipoles
dip_r = sqrt(sum(abs(LocMat).^2,1)); % Dipoles' radius
patch_dips = find(LocMat(1,:)>=1 & LocMat(1,:)<=4 & LocMat(2,:)>=1 & LocMat(2,:)<=4 & LocMat(3,:)>=3 & LocMat(3,:)<=4); %In the deep
patch_size = size(patch_dips,2);

% Setting the location of the dipole:
for i = 1:patch_size
    x_patch(i) = LocMat(1,patch_dips(i));
    y_patch(i) = LocMat(2,patch_dips(i));
    z_patch(i) = LocMat(3,patch_dips(i));
    R_patch(i) = sqrt(x_patch(i)^2 + y_patch(i)^2 + z_patch(i)^2);
    delta_x_patch(i) = x_patch(i)/R_patch(i);
    delta_y_patch(i) = y_patch(i)/R_patch(i);
    delta_z_patch(i) = z_patch(i)/R_patch(i);
end

% Plotting dipoles' locations
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:));
hold on
title("Location of patch dipoles & electrodes");
xlabel("X [cm]");
ylabel("Y [cm]");
zlabel("Z [cm]");
hold on

% Plotting electrodes' locations
scatter3(electrode_X,electrode_Y,electrode_Z,'filled');

% Electrodes' names
text(electrode_X,electrode_Y,electrode_Z,electrode_name);

% Adding the momentum vector:
for i = 1:patch_size
    quiver3(x_patch(i),y_patch(i),z_patch(i),delta_x_patch(i),delta_y_patch(i),delta_z_patch(i),'Green','LineWidth',3);
end

%% Part 14
%%%%%%%%%% section 1
% Random source for selected dipole
patch_source_numbers = randperm(size(Interictal,1),patch_size);

patch_sources = Interictal(patch_source_numbers,:);

% Making Q matrix
M_patch = zeros(21,10240);
Q_patch = zeros(3*patch_size,size(patch_sources,2));
for i = 1:patch_size
    Q_patch(3*i-2,:) = patch_sources(i,:)*delta_x_patch(i);
    Q_patch(3*i-1,:) = patch_sources(i,:)*delta_y_patch(i);
    Q_patch(3*i,:) = patch_sources(i,:)*delta_z_patch(i);

    % Finding potentials on electrodes
    M_patch = M_patch + GainMat(:,patch_dips(i)*3-2:patch_dips(i)*3)*Q_patch(3*i-2:3*i,:);
end

% Plotting potentials
offset = max(max(abs(M_patch)))/3 ;
disp_eeg(M_patch,offset,250,electrode_name);
xlim('tight');
grid minor
title("Potentials from patch dipoles");

%%%%%%%%%% section 2
figure;
Display_Potential_3D(ModelParams.R(3),electrode_mean_potential);
title("3D potential Display of Patch dipoles and source number");
xlabel("X [cm]"); ylabel("Y [cm]"); zlabel("Z [cm]");

%%%%%%%%%% section 3
% MNE
alpha = 1; % Tikhonov
N = size(M_patch,1); % Number of electrodes

%Computing predicted momentum matrix
Q_predicted_MNE_patch = GainMat.' * inv(GainMat*(GainMat.') + alpha*eye(N)) * M_patch;

% wMNE
p = size(LocMat,2);

% defining Omega & W
omega = zeros(p);
for i = 1:p
    for j = 1:N
        omega(i,i) = omega(i,i) + GainMat(j,i)*GainMat(j,i).';
    end
end

W = kron(omega,eye(3));

%Computing predicted momentum matrix
Q_predicted_wMNE_patch = inv(W.' * W) * GainMat.' * inv(GainMat * inv(W.' * W) * (GainMat.') + alpha*eye(N)) * M_patch;

%% Part 15
% Calculating amplitudes of each dipole
squared_Q_predicted_MNE_patch = Q_predicted_MNE_patch.^2;
squared_Q_predicted_wMNE_patch = Q_predicted_wMNE_patch.^2;

% Arrays of momentum amplitudes
amolitudes_MNE_patch = zeros(1,p); 
amolitudes_wMNE_patch = zeros(1,p); 
for i = 1:p
    amolitudes_MNE_patch(i) = sum(sum(squared_Q_predicted_MNE_patch((3*i-2):(3*i),:)));
    amolitudes_wMNE_patch(i) = sum(sum(squared_Q_predicted_wMNE_patch((3*i-2):(3*i),:)));
end

%% Part 16
% Making amplitudes vectors for ROC
amplitudes_main = zeros(1,p);
amplitudes_main(1,patch_dips) = 1;
amolitudes_MNE_patch_n = normalize(amolitudes_MNE_patch,2,"norm");
amolitudes_wMNE_patch_n = normalize(amolitudes_wMNE_patch,2,"norm");

% Labeling dipoles
n = 5001;
for i = 1:n
    i
    % Threshold defining
    thresh_roc = (1/(n-1))*(i-1);

    % Finding active dipoles based on threshold
    active_dipoles_MNE(i,:) = amolitudes_MNE_patch_n >= thresh_roc;
    active_dipoles_wMNE(i,:) = amolitudes_wMNE_patch_n >= thresh_roc;

    % TPR
    TPR_MNE(i) = sum(active_dipoles_MNE(i,:).*amplitudes_main(1,:) == 1)/patch_size;
    TPR_wMNE(i) = sum(active_dipoles_wMNE(i,:).*amplitudes_main(1,:) == 1)/patch_size;
    
    % FPR
    FPR_MNE(i) = (sum(active_dipoles_MNE(i,:)) - sum(active_dipoles_MNE(i,:).*amplitudes_main(1,:) == 1))/(p-patch_size);
    FPR_wMNE(i) = (sum(active_dipoles_wMNE(i,:)) - sum(active_dipoles_wMNE(i,:).*amplitudes_main(1,:) == 1))/(p-patch_size);
end

% ROCs
figure;
subplot(1,2,1);
plot(FPR_MNE,TPR_MNE,'LineWidth',2,'Color','b');
grid on
title("ROC of MNE");
xlabel("FPR");
ylabel("TPR");
xlim tight
ylim tight

subplot(1,2,2);
plot(FPR_wMNE,TPR_wMNE,'LineWidth',2,'Color','r');
grid on
title("ROC of wMNE");
xlabel("FPR");
ylabel("TPR");
xlim tight

%% (EXTRA) Part 11
%%%%%%%%%%%%% Section 1
% Loreta
alpha = 1; % Tikhonov
N = size(M,1);
p = size(LocMat,2);

% Defining Omega
omega = zeros(p);
for i = 1:p
    for j = 1:N
        omega(i,i) = omega(i,i) + GainMat(j,i)*GainMat(j,i).';
    end
end

% Defining B
d = 1;
oneP = ones(p,1);
a1 = zeros(p);
for i = 1:p
    for j = 1:p
        dist = sqrt( (LocMat(1,j)-LocMat(1,i))^2 + (LocMat(2,j)-LocMat(2,i))^2 + (LocMat(3,j)-LocMat(3,i))^2);
        if (dist == d) 
            a1(i,j) = 1/6;
        end
    end
end
a0 = inv(diag(a1*oneP)) * a1;
a = kron(a0,eye(3));
B = (6/d^2)*(a-eye(3*p));

% Defining W
new_omega = kron(omega,eye(3));
W = new_omega * (B.') * B * new_omega;

%Computing predicted momentum matrix
Q_predicted_LORETA = inv(W.' * W) * GainMat.' * inv(GainMat * inv(W.' * W) * (GainMat.') + alpha*eye(N)) * M;

%%%%%%%%%%%%% Section 2
% Calculating amplitudes of each dipole
squared_Q_predicted_LORETA = Q_predicted_LORETA.^2;

% Arrays of momentum amplitudes
amolitudes_LORETA = zeros(1,p); 
for i = 1:p
    amolitudes_LORETA(i) = sum(sum(squared_Q_predicted_LORETA((3*i-2):(3*i),:)));
end

% Predicting the dipoles (which maximizes the momentum amplitude)
dipode_predicted_LORETA = find(amolitudes_LORETA == max(amolitudes_LORETA));

%Setting the locations of the predicted dipoles:
x_LORETA = LocMat(1,dipode_predicted_LORETA);
y_LORETA = LocMat(2,dipode_predicted_LORETA);
z_LORETA = LocMat(3,dipode_predicted_LORETA);
R_LORETA = sqrt(x_LORETA^2 + y_LORETA^2 + z_LORETA^2);
delta_x_LORETA = x_LORETA/R_LORETA;
delta_y_LORETA = y_LORETA/R_LORETA;
delta_z_LORETA = z_LORETA/R_LORETA;

% Plotting dipoles' locations
figure;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:));
hold on
title("Location of Main and Predicted dipoles by LORETA");
xlabel("X [cm]");
ylabel("Y [cm]");
zlabel("Z [cm]");
hold on

% Plotting electrodes' locations
scatter3(electrode_X,electrode_Y,electrode_Z,'filled');

% Electrodes' names
text(electrode_X,electrode_Y,electrode_Z,electrode_name);

%Adding the momentum vectors (Main dipole)
quiver3(x,y,z,delta_x,delta_y,delta_z,'g','LineWidth',3);

%Adding the momentum vectors (Predicted dipoles)
quiver3(x_LORETA,y_LORETA,z_LORETA,delta_x_LORETA,delta_y_LORETA,delta_z_LORETA,'b','LineWidth',3);

%%%%%%%%%%%%% Section 3
% Distance error
distance_error_LORETA = sqrt((x-x_LORETA)^2 + (y-y_MNE)^2 + (z-z_LORETA)^2);

% Moment angle error
[fi, theta, radius] = cart2sph(x,y,z);
[fi_LORETA, theta_LORETA, radius_LORETA] = cart2sph(x_LORETA,y_LORETA,z_LORETA);

fi_error_LORETA = fi_LORETA - fi;
theta_error_LORETA = theta_LORETA - theta;

% Displaying errors:
disp("Dipole distance error for LORETA = " + distance_error_LORETA);
disp("Dipole fi error for LORETA = " + fi_error_LORETA);
disp("Dipole theta error for LORETA = " + theta_error_LORETA);

%%%%%%%%%%%%% Section 4
% By using A=2 in Part3 it will be studied

%% (EXTRA) Part 17
%%%%%%%%%%%%% Section 1
%Computing predicted momentum matrix
Q_predicted_LORETA_patch = inv(W.' * W) * GainMat.' * inv(GainMat * inv(W.' * W) * (GainMat.') + alpha*eye(N)) * M_patch;


%%%%%%%%%%%%% Section 2
% Calculating amplitudes of each dipole
squared_Q_predicted_LORETA_patch = Q_predicted_LORETA_patch.^2;

% Arrays of momentum amplitudes
amolitudes_LORETA_patch = zeros(1,p); 
for i = 1:p
    amolitudes_LORETA_patch(i) = sum(sum(squared_Q_predicted_LORETA_patch((3*i-2):(3*i),:)));
end

%%%%%%%%%%%%% Section 3
% Making amplitudes vectors for ROC
amplitudes_main = zeros(1,p);
amplitudes_main(1,patch_dips) = 1;
amolitudes_LORETA_patch_n = normalize(amolitudes_LORETA_patch,2,"norm");


% Labeling dipoles
n = 5001;
for i = 1:n
    i
    % Threshold defining
    thresh_roc = (1/(n-1))*(i-1);

    % Finding active dipoles based on threshold
    active_dipoles_LORETA(i,:) = amolitudes_LORETA_patch_n >= thresh_roc;

    % TPR
    TPR_LORETA(i) = sum(active_dipoles_LORETA(i,:).*amplitudes_main(1,:) == 1)/patch_size;
    
    % FPR
    FPR_LORETA(i) = (sum(active_dipoles_LORETA(i,:)) - sum(active_dipoles_LORETA(i,:).*amplitudes_main(1,:) == 1))/(p-patch_size);
end

% ROCs
figure;
plot(FPR_LORETA,TPR_LORETA,'LineWidth',2,'Color','b');
grid on
title("ROC of LORETA");
xlabel("FPR");
ylabel("TPR");
xlim tight
ylim tight

