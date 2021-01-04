
%-------------------------------------------------------------------%
% Software by -  Akash S. Shahade
%
% Title - Program for FE Analysis of Truss.
%-------------------------------------------------------------------%

clear;
clc; clf;
fprintf(' \nFEM ANALYSIS OF TRUSS \n\n');

%------------------------------------------------------%
% STEP 01 - PRE-PROCESSING
%------------------------------------------------------%

nel=2;

ele_nod=[1 2;2 3];                  % Element is connected with these nodes

nod_coor=[0 0;0.8 -0.6;0 -0.6];     % Node coordinates

num_ele=size(ele_nod,1);            % Element is connected with these nodes

ele_dof=[1 2 3 4;3 4 5 6];          % Degree of freedom (DOF) associated with nodes

num_nod=3;                          % Number of nodes

dof = 2;                            % DOF per node

displacement = zeros(dof*num_nod,1);    % Initialising Zero Matrix for Displacement
force = zeros(dof*num_nod,1);         % Initialising Zero Matrix for Force
stiffness = zeros(dof*num_nod);       % Initialising Zero Matrix for Stiffness


A(1)=5*10^(-5);                     % Area of element 01
A(2)=4*10^(-5);                     % Area of element 02
  

E(1)=200*10^9;                      % Young's Modulus of Element 01
E(2)=200*10^9;                      % Young's Modulus of Element 02

force(4)= -1000;                  % Load at node 2 in -ve Y direction


%------------------------------------------------------%
% Stiffness matrix calculation & ASSEMBLY
%------------------------------------------------------%

for e=1:num_ele                                                      % For i to number of elements
    
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...  % Calculate Length of element
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
  
 l=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);          % Calculate l = cos(theta)
 
 m=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);          % Calculate m = sin(theta)
 
 k=(A(e)*E(e)/L(e)*[l*l l*m -l*l -l*m;l*m m*m -l*m -m*m;...           % Calculate Stiffness Matrix
   -l*l -l*m l*l l*m; -l*m -m*m l*m m*m]);
   
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);

    for i=1:4
        for j=1:4
            
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);

        end
    end
end

fprintf('Global Stiffness Matrix: \n');
disp(stiffness);
fprintf('\n Global Load Vector: \n');
disp(force);
fprintf('\n------------------------------------------------\n');

%------------------------------------------------------%
% Boundary Conditions
%------------------------------------------------------%

fixed_node = [1 3];
stff = stiffness(:,3:4);
stiff_bc = stff(3:4,:);
force_bc = force(3:4,:);

%------------------------------------------------------%
% STEP 02 - SOLVE
%------------------------------------------------------%

dis_new = stiff_bc\force_bc ;

%------------------------------------------------------%
% STEP 03 - POST-PROCESSING
%------------------------------------------------------%

displacement = [0;0;dis_new;0;0];

% Stress in elements

for e=1:num_ele
    
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
  
 l=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 
 m=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
 
 stress(e)=(E(e)/L(e))*[-l -m l m]*displacement((ele_dof(e,:))');
 
end

% Reaction Forces

Reaction_Force = stiffness*displacement - force ;


%------------------------------------------------------%
% Plot 
%------------------------------------------------------%

for e=1:num_ele
    x=[nod_coor(ele_nod(e,1),1) nod_coor(ele_nod(e,2),1)];
    y=[nod_coor(ele_nod(e,1),2) nod_coor(ele_nod(e,2),2)];
plot(x,y,'b')
hold on
end

fprintf('\n Displacements: \n');
q={'q1x';'q1y';'q2x';'q2y';'q3x';'q3y'};
Displacement_mm = displacement*1000;
T = table(q,Displacement_mm);
disp(T);

fprintf('\n------------------------------------------------\n');

R = {'R1';'R2';'R3';'R4';'R5';'R6'};
P = table(R,Reaction_Force);
disp(P);

fprintf('\n------------------------------------------------\n');

Element = [1;2];
Stress_MPa = (stress.')/10^6;
B = table(Element,Stress_MPa);
disp(B);
fprintf('END OF PROGRAM.\n');
