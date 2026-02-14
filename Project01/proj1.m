%../APPM2360/Project01/proj1.m
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

function exportGraph(name, figHandle)
    plotsFolder = 'Project01/plots';

    % Create folder if it doesn't exist
    if ~exist(plotsFolder, 'dir')
        mkdir(plotsFolder);
    end

    pdfFilename = fullfile(plotsFolder, [name, '.pdf']);
    exportgraphics(figHandle, pdfFilename, 'ContentType', 'vector');
end

function A = equation01(t, r, n, A0)
    % Discrete compounding
    A = A0 * (1 + r/n).^(n*t);
end

function A = equation02(t, r, A0)
    % Continuous compounding
    A = A0 * exp(r*t);
end

function dAdt = equation03(r, A, p, A0)
    % Differential equation 
    dAdt = r*A - 12*p;
end

%{
#1 
Examine the effect of continuous compounding on the value of a loan. Assuming that the interest rate is 3% (r = 0.03) and
the original loan is $750,000, compute the total cost of the loan after 5 years for loans compounded 1, 2, 4, and 12 times per
year, without any payments, using Equation (1). Use Equation (2) to compare these values to the value of the loan compounded
continuously. On the same graph, plot the value of the loan as a function of time compounded 4 times a year and 12 times a year
as well as the value of the loan when the interest is compounded continuously for 0 ≤ t ≤ 30 years.
%}

ns = [1 2 4 12]; % number of compounds per year
comp_n = zeros(size(ns));

% Total cost after 5 years
for i = 1:length(ns)
    comp_n(i) = equation01(5, 0.03, ns(i), 750000);
end

fprintf('Total cost after 5 years:\n');
for i = 1:length(ns)
    fprintf('  n=%d: $%.2f\n', ns(i), comp_n(i));
end

continuous_value = equation02(5, 0.03, 750000);
fprintf('  Continuous: $%.2f\n\n', continuous_value);

% Plot for 0 to including 30 years
t = 0:0.1:30;

f = figure;

plot(t, equation01(t, 0.03, 4, 750000), 'Color', '#CFB87C', 'LineWidth', 2);
hold on;
plot(t, equation01(t, 0.03, 12, 750000), 'Color', '#A2A4A3', 'LineWidth', 2);
plot(t, equation02(t, 0.03, 750000), 'Color', '#0A3758', 'LineWidth', 2);

xlabel('Time (years)');
ylabel('Loan Value ($)');
title('Loan Value: Compounded 4x/year, 12x/year, and Continuously');

legend('n = 4', 'n = 12', 'Continuous', 'Location', 'northwest');

grid on;
hold off;

exportGraph('3.1.1', f);


%{
Total cost after 5 years:
  n=1: $869455.56
  n=2: $870405.62
  n=4: $870888.11
  n=12: $871212.59
  Continuous: $871375.68

%}

%{
#2
Next, gain a broad understanding of the behavior of the loan value by determining whether there any equilibrium solutions to Eq.
(3). If so, what are they, and what is their stability? What do these equilibria represent in real-world terms?
%}


%{
#3
Determine the exact behavior of the loan in your friends’ situation by solving (3) with A(0) = A0 and r and p arbitrary. Be sure
to show your work so that your friends are confident that you have the correct solution.
%}


%{
#4
The size of the monthly payment p that your friends are willing to make plays a large role in deciding the type of loan they should
choose. Use the solution to (3) to find the correct p to pay off a 10-year fixed rate mortgage with rate of 3% and initial debt of
$750,000. Do the same for a 30-year fixed rate mortgage with an interest rate of 5%. Hint: you want to find p such that A(tl) = 0,
were tl is the duration (years) of your mortgage. Find this analytically, not numerically using a root finding routine.
%}
A0 = 750000;

% 10-year at 3%
r1 = 0.03;
t1 = 10;
p1 = (A0*r1) / (12*(1 - exp(-r1*t1)));

% 30-year at 5%
r2 = 0.05;
t2 = 30;
p2 = (A0*r2) / (12*(1 - exp(-r2*t2)));

fprintf('10-year @ 3%% monthly payment p = $%.2f\n', p1);
fprintf('30-year @ 5%% monthly payment p = $%.2f\n', p2);

% Optional: verify analytically that A(t_l)=0 using the closed form
A_t1 = (12*p1)/r1 + (A0 - (12*p1)/r1)*exp(r1*t1);
A_t2 = (12*p2)/r2 + (A0 - (12*p2)/r2)*exp(r2*t2);

fprintf('Check A(t1)= %.6f (should be ~0)\n', A_t1);
fprintf('Check A(t2)= %.6f (should be ~0)\n', A_t2);

%{
#5
While having a low monthly payment is nice, you should warn your friends that there is quite literally a price to pay for this
convenience. We can determine the total amount paid by summing each monthly payment over the duration of the loan. How
much interest is paid in the 30-year fixed rate mortgage? The 10-year?
%}

%{
#6
Buyers often choose to pay as much of the cost as they can up front (make a down payment) so that they don’t have to borrow quite
so much. Might this option be worth it for your friends? How much money would they save in each case if they paid $100,000
down on the house, i.e., the original loan amount was $650,000? Use the interest rates and loan periods from part (4).
%}

%{
#7
What are the advantages and disadvantages of taking out a 30-year fixed rate mortgage as opposed to a 10-year mortgage?
%}








%Consider a mortgage for $750,000 with a constant interest rate of 5% (r = 0.05) and a monthly payment p = $4000.

    % 1. Implement Euler’s method for Eq. (3) with step size h = 0.5. Run the method until the mortgage is paid off and determine when it is paid off. Note: in reality, the mortgage is paid off when its value is zero. However, due to errors in the computations (both discretization and roundoff), it is likely that Euler’s method will not produce an exact value of 0 for the mortgage value for any time. To account for this, consider the mortgage to be paid off when its value first becomes negative.
    

    % 2. Plot the numerical solution A(t) and the true solution to Eq. (3) with the parameters given here on the same graph and compare the two. 
    

    %3. Repeat the previous item for a step size h = 0.01 and comment on the difference

A0 = 750000;
r  = 0.05;
p  = 4000;

% (1) Euler with h = 0.5 until payoff  
h = 0.5;
t = 0;
A = A0;
k = 1;
maxYears = 200; 

while A(k) >= 0 && t(k) < maxYears
    t(k+1) = t(k) + h;
    A(k+1) = A(k) + h*(r*A(k) - 12*p);
    k = k + 1;
end

tpay_h05 = t(end);

% (2) Plot Euler and true solution together h=0.5
A_true = (12*p)/r + (A0 - (12*p)/r).*exp(r*t);

f = figure;
plot(t, A, 'LineWidth', 2, 'Color', '#565A5C'); hold on;
plot(t, A_true, 'LineWidth', 2, 'Color', '#CFB87C');
yline(0, 'k--');
xlabel('t (years)'); ylabel('A(t) ($)');
title(sprintf('Fixed rate: Euler vs True (h = %.2f)', h));
legend('Euler', 'True', 'Location', 'best');
grid on; hold off;
exportGraph('3.2.2', f);

% (3) Repeat for h=0.01  
h = 0.01;
t2 = 0;
A2 = A0;
k = 1;

while A2(k) >= 0 && t2(k) < maxYears
    t2(k+1) = t2(k) + h;
    A2(k+1) = A2(k) + h*(r*A2(k) - 12*p);
    k = k + 1;
end

tpay_h001 = t2(end);

fprintf('Fixed rate (r=%.2f, p=$%.0f): Euler h=%.2f payoff ~ %.2f years\n\n', r, p, h, tpay_h001);

A_true2 = (12*p)/r + (A0 - (12*p)/r).*exp(r*t2);

f = figure;
plot(t2, A2, 'LineWidth', 2, 'Color', '#565A5C'); hold on;
plot(t2, A_true2, 'LineWidth', 2, 'Color', '#CFB87C');
yline(0, 'k--');
xlabel('t (years)'); ylabel('A(t) ($)');
title(sprintf('Fixed rate: Euler vs True (h = %.2f)', h));
legend('Euler', 'True', 'Location', 'best');
grid on; hold off;
exportGraph('3.2.3', f);


%Now we turn to the adjustable rate mortgages. Suppose that for the same $750,000 mortgage a bank offers an adjustable rate mortgage, which starts with an initial lower fixed rate of 3% (r = 0.03) for the first 5 years and is tied to credit markets after that. Let’s assume that after the first 5 years the rate increases as r(t) = 0.03 + 0.015√t − 5, so
%Use Euler’s method with h = 0.01 to answer the following.

    % 1. Suppose your friends pay $4000 per month. How long will it take them to pay off the mortgage?


    % 2. What about if they pay $4500 per month?


    % 3. How much interest is paid in each case?


    % 4. Plot the numerical solution A(t) for both scenarios on the same graph. How does the variable interest rate affect the graph, compared to the fixed rate? How do the different payment sizes affect the graph?


% 3.2.2 Adjustable rate mortgage (Euler h=0.01)
% r(t)=0.03 for t<=5
% r(t)=0.03 + 0.015*sqrt(t-5) for t>5

A0 = 750000;
h  = 0.01;

% Scenario 1: p = 4000   
p1 = 4000;

tA = 0;
AA = A0;
interestA = 0;

k = 1;
while AA(k) >= 0 && tA(k) < maxYears
    if tA(k) <= 5
        rk = 0.03;
    else
        rk = 0.03 + 0.015*sqrt(tA(k) - 5);
    end

    interestA = interestA + h*(rk*AA(k));           % interest accumulated
    tA(k+1) = tA(k) + h;
    AA(k+1) = AA(k) + h*(rk*AA(k) - 12*p1);         % Euler step
    k = k + 1;
end

tpay_p4000 = tA(end);

% Scenario 2: p = 4500   
p2 = 4500;

tB = 0;
AB = A0;
interestB = 0;

k = 1;
while AB(k) >= 0 && tB(k) < maxYears
    if tB(k) <= 5
        rk = 0.03;
    else
        rk = 0.03 + 0.015*sqrt(tB(k) - 5); 
    end

    interestB = interestB + h*(rk*AB(k));
    tB(k+1) = tB(k) + h;
    AB(k+1) = AB(k) + h*(rk*AB(k) - 12*p2);
    k = k + 1;
end

tpay_p4500 = tB(end);

fprintf('ARM (Euler h=%.2f):\n', h);
fprintf('  p=$%d payoff ~ %.2f years, interest paid ~ $%.2f\n', p1, tpay_p4000, interestA);
fprintf('  p=$%d payoff ~ %.2f years, interest paid ~ $%.2f\n', p2, tpay_p4500, interestB);


f = figure;
plot(tA, AA, 'LineWidth', 2,'Color', '#565A5C'); 
hold on;
plot(tB, AB, 'LineWidth', 2, 'Color', '#CFB87C');
yline(0, 'k--');
xlabel('t (years)'); ylabel('A(t) ($)');
title('Adjustable rate mortgage: Euler h=0.01');
legend(sprintf('p=$%d (payoff ~ %.2f yrs)', p1, tpay_p4000), ...
       sprintf('p=$%d (payoff ~ %.2f yrs)', p2, tpay_p4500), ...
       'Location', 'best');
grid on; hold off;
exportGraph('3.2.4', f);