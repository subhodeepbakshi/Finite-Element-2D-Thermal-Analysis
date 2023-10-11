clear all;
clc;
% Given Parameters
syms t4 t5 t6 t7 t8 t9 t10;
[kxx, kyy] = deal(45);
[t1, t2, t3] = deal(473); %% converted to kelvin
[t11, t12, t13, t14] = deal(298); %% converted to kelvin

% Nodes present at each element
ele = [1 2 5;
    2 3 6;
    2 6 5;
    1 5 4;
    3 7 6;
    6 7 10;
    6 10 9;
    5 6 9;
    5 9 8;
    4 5 8;
    8 12 11;
    8 9  12;
    9 13 12;
    9 10 13;
    10 14 13];
% x coordinate and y coordinate at each node
node = [0 0;
    0 1;
    0 2;
    1 0;
    1 0.5;
    1 1.5;
    1 2;
    2 0;
    2 1;
    2,2;
    3 0;
    3 0.5;
    3 1.5;
    3 2];

l = length(ele);
Kglobal = zeros(14, 14);
d = [t1; t2; t3; t4; t5; t6; t7; t8; t9; t10; t11; t12; t13; t14];

% Plotting the initial finite element mesh
figure(1);

for i = 1:l
    Na = ele(i, 1);
    Nb = ele(i, 2);
    Nc = ele(i, 3);

    % Plotting the element
    plot([node(Na, 1), node(Nb, 1), node(Nc, 1), node(Na, 1)], [node(Na, 2), node(Nb, 2), node(Nc, 2), node(Na, 2)], '-o');
    hold on;
end

title('Initial Finite Element Mesh'); % Corrected comment
xlabel('X-coordinate');
ylabel('Y-coordinate');

% Assigning nodal temperatures to the plot
for i = 1:size(node, 1)
    text(node(i, 1), node(i, 2), char(d(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

hold off;

for i = 1:l
    % for ith element to find corresponding node
    Na = ele(i, 1);
    Nb = ele(i, 2);
    Nc = ele(i, 3);

    % x coordinate for those nodes
    xa = node(Na, 1);
    xb = node(Nb, 1);
    xc = node(Nc, 1);
    % y coordinate for those nodes
    ya = node(Na, 2);
    yb = node(Nb, 2);
    yc = node(Nc, 2);

    b1 = (yb - yc);
    b2 = (yc - ya);
    b3 = (ya - yb);

    c1 = (xc - xb);
    c2 = (xa - xc);
    c3 = (xb - xa);

    % area of the element
    Ae = 1/2 * abs(det([1 xa ya; 1 xb yb; 1 xc yc]));
    % local stiffness matrix
    K = [kxx * b1^2 + kyy * c1^2 kxx * b1 * b2 + kyy * c1 * c2 kxx * b1 * b3 + kyy * c1 * c3;
        kxx * b1 * b2 + kyy * c1 * c2 kxx * b2^2 + kyy * c2^2 kxx * b2 * b3 + kyy * c2 * c3;
        kxx * b1 * b3 + kyy * c1 * c3 kxx * b2 * b3 + kyy * c2 * c3 kxx * b3^2 + kyy * c3^2] / (4 * Ae);

    % Populating the global matrix
    Kglobal(Na, Na) = Kglobal(Na, Na) + K(1, 1);
    Kglobal(Na, Nb) = Kglobal(Na, Nb) + K(1, 2);
    Kglobal(Na, Nc) = Kglobal(Na, Nc) + K(1, 3);

    Kglobal(Nb, Na) = Kglobal(Nb, Na) + K(2, 1);
    Kglobal(Nb, Nb) = Kglobal(Nb, Nb) + K(2, 2);
    Kglobal(Nb, Nc) = Kglobal(Nb, Nc) + K(2, 3);

    Kglobal(Nc, Na) = Kglobal(Nc, Na) + K(3, 1);
    Kglobal(Nc, Nb) = Kglobal(Nc, Nb) + K(3, 2);
    Kglobal(Nc, Nc) = Kglobal(Nc, Nc) + K(3, 3);

end

% Solving for nodal temperatures
R_Result = Kglobal * d;

% Forming 7 unknown equations from the matrix result
Eqn1 = R_Result(4) == 0;
Eqn2 = R_Result(5) == 0;
Eqn3 = R_Result(6) == 0;
Eqn4 = R_Result(7) == 0;
Eqn5 = R_Result(8) == 0;
Eqn6 = R_Result(9) == 0;
Eqn7 = R_Result(10) == 0;

sol = solve([Eqn1, Eqn2, Eqn3, Eqn4, Eqn5, Eqn6, Eqn7], [t4, t5, t6, t7, t8, t9, t10]);
t4 = double(sol.t4);
t5 = double(sol.t5);
t6 = double(sol.t6);
t7 = double(sol.t7);
t8 = double(sol.t8);
t9 = double(sol.t9);
t10 = double(sol.t10);

% Displaying the results
disp('-------------------------------------------------------------------');
disp('Nodal Temperatures are in kelvin:');
disp(t1);
disp(t2);
disp(t3);
disp(t4);
disp(t5);
disp(t6);
disp(t7);
disp(t8);
disp(t9);
disp(t10);
disp(t11);
disp(t12);
disp(t13);
disp(t14);

% Plotting the final finite element mesh
figure(2);

for i = 1:l
    Na = ele(i, 1);
    Nb = ele(i, 2);
    Nc = ele(i, 3);

    % Plotting the element
    plot([node(Na, 1), node(Nb, 1), node(Nc, 1), node(Na, 1)], [node(Na, 2), node(Nb, 2), node(Nc, 2), node(Na, 2)], '-o');
    hold on;
end

title('Final Finite Element Mesh'); % Corrected comment
xlabel('X-coordinate');
ylabel('Y-coordinate');

% Assigning nodal temperatures to the plot
for i = 1:size(node, 1)
    text(node(i, 1), node(i, 2), num2str(eval(['t', num2str(i)]), '%.2f'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

hold off;
