%../APPM2360/Project01/proj1.m

s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

function ex = exportGraph(name, figure)
    plotsFolder = 'Project01/plots';
    pdfFilename = fullfile(plotsFolder, [name, '.pdf']);
    exportgraphics(figure, pdfFilename, 'ContentType', 'vector');
end

function A = equation01(t, r, n, A0)
    A = A0 * (1 + r/n).^(n*t);
end

function A = equation02(t, r, A0)
    A = A0*exp(r*t);
end

function dAdt = equation03(r, A, p, A0)
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

ns = [1 2 4 12]; %number compounds per year
comp_n = zeros(size(ns));

% total cost after 5 years
for i = 1:length(ns)
    comp_n(i) = equation01(5, 0.03, ns(i), 750000);
end

fprintf('Tot cost after 5 years:\n');
for i = 1:length(ns)
    fprintf('  n=%d: $%.2f\n', ns(i), comp_n(i));
end
fprintf('  Continuous: $%.2f\n\n', equation02(5, 0.03, 750000));

% Plot for 0 <= t <= 30 years
t = 0:0.1:30; % time in years

f = figure;
plot(t, equation01(t, 0.03, 4, 750000), 'b-', 'LineWidth', 2);
hold on;
plot(t, equation01(t, 0.03, 12, 750000), 'g-', 'LineWidth', 2);
plot(t, equation02(t, 0.03, 750000), 'r-', 'LineWidth', 2);

xlabel('Time (years)');
ylabel('Loan Value ($)');
title('Loan Value: Compounded 4x/year, 12x/year, and Continuously');
legend('n=4', 'n=12', 'Continuous');
grid on;

exportGraph('3.1.1', f);


%{
Tot cost after 5 years:
  n=1: $869455.56
  n=2: $870405.62
  n=4: $870888.11
  n=12: $871212.59
  Continuous: $871375.68%}

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