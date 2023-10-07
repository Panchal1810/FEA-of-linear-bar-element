%% Intialization
close all;
clc;
clear all;

%% here we are following N for force and mm for displacement
% DATA FROM ABAQUS READ FROM EXCEL SHEETS
CCORD4=xlsread('bar_element_E4.xlsx',1,'A1:D6'); %Node number and corresponding cordinates are loaded
%NCA=xlsread('bar_element_E4.xlsx',1,'G1:I5');   % Element number and connected nodes with the elemnts are loaded

CCORD8=xlsread('bar_element_E8.xlsx',1,'A1:D10'); %Node number and corresponding cordinates are loaded
%NCA=xlsread('bar_element_E8.xlsx',1,'G1:I9');   % Element number and connected nodes with the elemnts are loaded

CCORD12=xlsread('bar_element_E12.xlsx',1,'A1:D14'); %Node number and corresponding cordinates are loaded
%NCA=xlsread('bar_element_E12.xlsx',1,'G1:I13');   % Element number and connected nodes with the elemnts are loaded

CCORD=xlsread('bar_element_E16.xlsx',1,'A1:D18'); %Node number and corresponding cordinates are loaded
NCA=xlsread('bar_element_E16.xlsx',1,'G1:I17');   % Element number and connected nodes with the elemnts are loaded


CCORD16 = CCORD;


L=0.2;
NNODES=length(CCORD);
NELEMENTS=length(NCA);       
DOFPN=1;          % Degree of freedom per node
nnp=DOFPN*NNODES; % Total degree of freedom
K=zeros(nnp,nnp); % Intializing stiffness matrix
F=zeros(nnp,1);   %Initializing force matrix
F_q = zeros(nnp,1); %Initializing distributed load
e = zeros(NELEMENTS,1);
sigma = zeros(NELEMENTS,1);

l=L/NELEMENTS;

A_0 = 0.1*0.02;
A_L = 0.01*0.02;
E = 80*1e+09;

 
%% Stifness Matrix Generation
for EN=1:NELEMENTS
     [B]=BCAL(l);           % Calculation of B matrix
    [ ke ] = Kel_bar(EN,A_0,A_L,L,l,CCORD,B,E);% Calculation of element stiffness matrix
    RN=[NCA(EN,2) NCA(EN,3)] ;%ROWN(EN,NCA);           % Global matrix location identification
    for i=1:2
        for  j=1:2
           K(RN(i),RN(j))=K(RN(i),RN(j))+ke(i,j);
        end
    end
end



%% Boundary Conditions
fixed_nodes = 1;     % Restricted degrees of freedom
free_nodes = setxor(1:nnp,fixed_nodes); % Free degrees of freedom
force_nodes =nnp;                       % Node at which load is applied
force_val =5000;
body_force = 10*1000;
F(force_nodes)= force_val;              % Assigning the load value at the particular node

for EN=1:NELEMENTS
   
 RN2 = [NCA(EN,2); NCA(EN,3)];

F_qe = [body_force*l/4*(-l/3+CCORD(EN,2)+CCORD(EN+1,2)); body_force*l/4*(l/3+CCORD(EN,2)+CCORD(EN+1,2))];

for i=1:2
F_q(RN2(i)) =  F_q(RN2(i))+F_qe(i);
end

end

%% Partitioning K and F matrix
Kpart = K(free_nodes,free_nodes);       
Fpart = F(free_nodes,1)+ F_q(free_nodes,1);

%% Solve the system of eqn
Up = Kpart\Fpart;
U=zeros(nnp,1);
U(free_nodes)=Up;


%% post-processing the data
% stress-strain calculation

for EN=1:NELEMENTS

e(EN) = B'*[U(EN,:); U(EN+1,:)] ;
sigma(EN) = E*e(EN);
end

%% plot convergence solution

U_ELE = [1.5347e-05 1.5981e-05 1.6140e-05 1.6202e-05];
for i = 1:3
    
ERROR(i) = ((U_ELE(i+1) - U_ELE(i))/U_ELE(i+1))*100;
end

N_ElE = [4 8 12 16];

U_ELE4 = [0
1.8295e-06
4.2686e-06
7.9263e-06
1.5347e-05];

U_ELE8 = [0
8.6075e-07
1.8368e-06
2.9630e-06
4.2932e-06
5.9187e-06
8.0110e-06
1.0957e-05
1.5981e-05];

U_ELE12 = [0
5.6272e-07
1.1727e-06
1.8382e-06
2.5703e-06
3.3834e-06
4.2979e-06
5.3429e-06
6.5624e-06
8.0277e-06
9.8663e-06
1.2340e-05
1.6140e-05];

U_ELE16 = [0
4.1799e-07
8.6152e-07
1.3338e-06
1.8387e-06
2.3810e-06
2.9665e-06
3.6028e-06
4.2996e-06
5.0696e-06
5.9302e-06
6.9061e-06
8.0337e-06
9.3700e-06
1.1012e-05
1.3146e-05
1.6202e-05];

figure(1);

plot(CCORD4(:,2),U_ELE4(:,1),'-or','MarkerSize',10);

hold on;

plot(CCORD8(:,2),U_ELE8(:,1),'-sk','MarkerSize',10);

hold on;


plot(CCORD12(:,2),U_ELE12(:,1),'-*g','MarkerSize',10);

hold on;

plot(CCORD16(:,2),U_ELE16(:,1),'-+b','MarkerSize',10);

%plot()

xlim([0 0.2]);

%% U profile
figure(2);

plot(CCORD(:,2),U(:,1),'-or','MarkerSize',10);


%% plotting the stress and strain with respect to x direction
figure(3);

for k=1:NELEMENTS

      plot([CCORD(k,2) CCORD(k+1,2)],e(k),'-sk');

       hold on;
end


figure(4);

for k=1:NELEMENTS

      plot([CCORD(k,2) CCORD(k+1,2)],sigma(k),'-or');

       hold on;
end



