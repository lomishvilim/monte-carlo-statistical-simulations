%%
% Lottery Match Distribution project (theory + simulations)

%this project:
% 1) Computes the theoretical distribution of the number of matches X
% between a fixed lottery ticket and the winning 6/49 numbers.
% 2) Computes a binomial approximation for comparison.
% 3) Simulates many lottery drawings using randperm (without replacement).
% 4) Estimates a second player's matches and the joint distribution (X1, X2).
% 5) Uses CLT to study the sampling distribution of P_hat = P(X >= 3).
% 6) Generates figures and CSV tables for the final report.



clear; clc; close all;
%%
% lottery parameters
N = 49; % population size (numbers 1..49)
K = 6; % number of winning numbers in the population
n = 6; % numbers drawn each time

% tickets for two players, no repeats
ticket1 = sort([3 11 17 25 36 42]);
ticket2 = sort([5 11 19 25 33 42]);

% simulation parameters
M = 1000; % draws per block
B = 200; % number of blocks (total trials = B * M)
Ntrials = B * M; % total number of simulated lottery drawings = 200,000
threshold_ge = 3; % we'll focus on event {X >= 3}

% reproducibility - initialize random functions' seed
rng(12345);


%%
% 1. Theoretical distribution of X using hypergeometric

k_vals = 0:n; % possible number of matches 0..6
num_k  = numel(k_vals);
pmf_hyper = zeros(1, num_k);

%hypergeometric pmf P(X = k)
for i = 1:num_k
   k = k_vals(i);
   pmf_hyper(i) = hypergeo_pmf(k, N, K, n); %helper function, defined below
end

% Cumulative Distribution Function CDF from pmf values
cdf_hyper = cumsum(pmf_hyper);

% theoretical mean and variance from the distribution (for checking)
mean_hyper = sum(k_vals .* pmf_hyper);
var_hyper  = sum((k_vals - mean_hyper).^2 .* pmf_hyper);

% Hypergeometric theory: event probabilities
P_eq0  = pmf_hyper(k_vals == 0); % P(X = 0)
P_eq1  = pmf_hyper(k_vals == 1); % P(X = 1)
P_eq2  = pmf_hyper(k_vals == 2); % P(X = 2)
P_eq3  = pmf_hyper(k_vals == 3); % P(X = 3)
P_ge3  = sum(pmf_hyper(k_vals >= 3)); % P(X >= 3)
P_le2  = sum(pmf_hyper(k_vals <= 2)); % P(X <= 2) = 1 - P_ge3



%%
% 2. Binomial distribution - for comparison

p_success = K / N; % probability of matching a given winning number
pmf_binom = zeros(1, num_k);

for i = 1:num_k
   k = k_vals(i);
   pmf_binom(i) = binomial_pmf(k, n, p_success); %helper function, defined below
end



%%
% 3. Simulation of lottery drawings

%Simulate Ntrials amount of independent drawings of 6 numbers from 1 to 49 
% (without replacement) using randperm)

X1 = zeros(Ntrials, 1); %matches for ticket1
X2 = zeros(Ntrials, 1); %matches for ticket2
X = zeros(Ntrials, 1); %same as X1, for single ticket case!

for t = 1:Ntrials
   % random winning numbers (with no replacement)
   winning = randperm(N, n);
  
   % count matches for each ticket
   matches1 = sum(ismember(ticket1, winning));
   matches2 = sum(ismember(ticket2, winning));
  
   X1(t) = matches1;
   X2(t) = matches2;
   X(t)  = matches1; %for single ticket distribution
end



%%
% 4. Empirical distribution and basic sample statistics


% histogram-based pmf estimate for X 
edges = -0.5 : 1 : (n + 0.5); % bin edges centered on integers
counts_X = histcounts(X, edges);
pmf_sim  = counts_X / Ntrials;

% sample mean and variance for X
mean_sim = mean(X);
var_sim  = var(X, 1);

% compare theoretical vs simulated mean/variance
fprintf('Single-ticket X: number of matches\n');

fprintf('Theoretical mean E[X] = %.6f\n', mean_hyper);
fprintf('Simulated  mean (sample) = %.6f\n', mean_sim);
fprintf('Theoretical variance Var[X] = %.6f\n', var_hyper);
fprintf('Simulated  variance (sample) = %.6f\n\n', var_sim);

% empirical event probabilities from simulation
P_sim_eq0 = pmf_sim(k_vals == 0);
P_sim_eq1 = pmf_sim(k_vals == 1);
P_sim_eq2 = pmf_sim(k_vals == 2);
P_sim_eq3 = pmf_sim(k_vals == 3);
P_sim_ge3 = sum(pmf_sim(k_vals >= 3));
P_sim_le2 = sum(pmf_sim(k_vals <= 2));

fprintf('P(X = 0): theory = %.6f, simulation = %.6f\n', P_eq0, P_sim_eq0);
fprintf('P(X = 1): theory = %.6f, simulation = %.6f\n', P_eq1, P_sim_eq1);
fprintf('P(X = 2): theory = %.6f, simulation = %.6f\n', P_eq2, P_sim_eq2);
fprintf('P(X = 3): theory = %.6f, simulation = %.6f\n', P_eq3, P_sim_eq3);
fprintf('P(X >= 3): theory = %.8f, simulation = %.8f\n\n', P_ge3, P_sim_ge3);




%%
% 5. Joint distribution of X1 and X2 (two players' tickets), covariance, correlation

% Joint frequency table for X1, X2 in 0…6
% but add 1 to use as indices 1…7
joint_counts = accumarray([X1 + 1, X2 + 1], 1, [n + 1, n + 1]);
joint_probs  = joint_counts / Ntrials;


%sample means, variances, covariance, correlation
mean_X1 = mean(X1);
mean_X2 = mean(X2);
var_X1  = var(X1, 1);
var_X2  = var(X2, 1);
cov_X1X2 = mean((X1 - mean_X1) .* (X2 - mean_X2));
corr_X1X2 = cov_X1X2 / sqrt(var_X1 * var_X2);


% probability both players have at least 3 matches
P_both_ge3_sim = mean( (X1 >= threshold_ge) & (X2 >= threshold_ge) );
fprintf('Two-player analysis\n');
fprintf('Mean(X1) = %.6f, Mean(X2) = %.6f\n', mean_X1, mean_X2);
fprintf('Var(X1)  = %.6f, Var(X2)  = %.6f\n', var_X1, var_X2);
fprintf('Cov(X1, X2) = %.6f\n', cov_X1X2);
fprintf('Corr(X1, X2) = %.6f\n', corr_X1X2);
fprintf('P(X1 >= 3 and X2 >= 3) (simulated) = %.10f\n\n', P_both_ge3_sim);


%%
% 6. CLT analysis for P_hat = P(X >= 3)


% Break the Ntrials observations into B blocks of size M (Ntrials = B*M)
% For each block b,   P_hat_b = (# of trials with X >= 3) / M

if B * M ~= Ntrials % check for no errors
   error('B * M must equal Ntrials.');
end


indicator_ge3 = (X >= threshold_ge);
P_hat_blocks  = zeros(B, 1);

for b = 1:B
   idx_start = (b - 1) * M + 1;
   idx_end   = b * M;
   block_data = indicator_ge3(idx_start : idx_end);
   P_hat_blocks(b) = mean(block_data);
end


% theoretical CLT mean and std
p_star = P_ge3;
sigma_clt = sqrt(p_star * (1 - p_star) / M);


% empirical mean and std of P_hat across blocks
mean_P_hat = mean(P_hat_blocks);
std_P_hat  = std(P_hat_blocks, 1); % population std across B blocks


fprintf('CLT analysis for P_hat = P(X >= %d)\n', threshold_ge);
fprintf('Theoretical p* = P(X >= %d) = %.8f\n', threshold_ge, p_star);
fprintf('Mean(P_hat) across blocks (sim) = %.8f\n', mean_P_hat);
fprintf('Theoretical std(P_hat) = %.8f\n', sigma_clt);
fprintf('Empirical  std(P_hat) (sim) = %.8f\n\n', std_P_hat);




%%
% 7. Tables (MATLAB tables and csv for the report)
%needed gpt's help for this part


% Table 1: comparison of hypergeometric theory, binomial approximation,
% and simulated pmf for X (0..6).

Sim_minus_Theory = pmf_sim(:) - pmf_hyper(:);
Binom_minus_Hyper = pmf_binom(:) - pmf_hyper(:);
T_pmf = table(k_vals(:), pmf_hyper(:), pmf_binom(:), pmf_sim(:), ...
             Binom_minus_Hyper, Sim_minus_Theory, ...
   'VariableNames', {'k', 'P_hyper', 'P_binom', 'P_sim', ...
                     'P_binom_minus_P_hyper', 'P_sim_minus_P_hyper'});
disp('Table: pmf comparison (hypergeometric vs binomial vs simulation)');
disp(T_pmf);
writetable(T_pmf, 'pmf_comparison_table.csv');


% Table 2: summary of CLT statistics
T_clt_summary = table(p_star, mean_P_hat, sigma_clt, std_P_hat, ...
   'VariableNames', {'p_star_theory', 'mean_P_hat_sim', ...
                     'sigma_clt_theory', 'std_P_hat_sim'});
disp('Table: CLT summary');
disp(T_clt_summary);
writetable(T_clt_summary, 'clt_summary_table.csv');


% 3. Joint probability matrix written to csv 
% Rows: X1 = 0..6, Columns: X2 = 0..6
writematrix(joint_probs, 'joint_distribution_probs.csv');
% Also write the block-wise P_hat values
T_clt_blocks = table((1:B).', P_hat_blocks, ...
   'VariableNames', {'BlockIndex', 'P_hat'});
writetable(T_clt_blocks, 'clt_blocks_table.csv');


%%
% 8. Figures: pmf comparison, binomial vs hypergeometric, joint heatmap,
% CLT histogram with normal overlay


% Theoretical vs simulated pmf of X
figure;
bar(k_vals, [pmf_hyper(:), pmf_sim(:)], 'grouped');
xlabel('Number of matches k');
ylabel('Probability');
legend('Hypergeometric theory', 'Simulation', 'Location', 'best');
title('Lottery 6/49: P(X = k) – Theory vs Simulation');
grid on;
saveas(gcf, 'pmf_theory_vs_simulation.png');


% Hypergeometric vs binomial approximation
figure;
bar(k_vals, [pmf_hyper(:), pmf_binom(:)], 'grouped');
xlabel('Number of matches k');
ylabel('Probability');
legend('Hypergeometric theory', 'Binomial approximation', 'Location', 'best');
title('Hypergeometric vs Binomial Approximation (n = 6, p = 6/49)');
grid on;
saveas(gcf, 'pmf_hyper_vs_binomial.png');


% Joint distribution heatmap for (X1, X2)
figure;
imagesc(joint_probs);
colorbar;
axis equal tight;
set(gca, 'XTick', 1:(n+1), 'XTickLabel', 0:n, ...
        'YTick', 1:(n+1), 'YTickLabel', 0:n);
xlabel('X_2 (matches for ticket 2)');
ylabel('X_1 (matches for ticket 1)');
title('Joint Distribution of (X_1, X_2) – Simulated Probabilities');
saveas(gcf, 'joint_distribution_heatmap.png');


% CLT: Histogram of P_hat with normal approximation overlay
figure;
h = histogram(P_hat_blocks, 'Normalization', 'pdf');
hold on;
x_min = min(P_hat_blocks);
x_max = max(P_hat_blocks);
x_grid = linspace(x_min, x_max, 300);
normal_pdf = (1 / (sqrt(2*pi)*sigma_clt)) * ...
            exp(-(x_grid - p_star).^2 / (2 * sigma_clt^2));
plot(x_grid, normal_pdf, 'LineWidth', 2);
xlabel('P\_hat (proportion of trials with X \ge 3 in each block)');
ylabel('Density');
legend('Empirical distribution of P\_hat', 'Normal approximation', ...
      'Location', 'best');
title(sprintf('CLT: Distribution of P\\_hat (Block size M = %d)', M));
grid on;
hold off;
saveas(gcf, 'clt_phat_histogram.png');


%%
% Helper functions


function p = hypergeo_pmf(k, N, K, n) 
% Hypergeometric pmf P(X = k).
% N,K,n = Population, success states/winning numbers, draws
% X is number of successes in the sample

%   P(X = k) = [C(K, k) * C(N-K, n-k)] / C(N, n)
   if k < 0 || k > n || k > K || (n - k) > (N - K)
       p = 0;
       return;
   end
  
   p = nchoosek(K, k) * nchoosek(N - K, n - k) / nchoosek(N, n);
end



function p = binomial_pmf(k, n, p_success)
%binomial pmf P(Y = k) for Y = Binomial(n, p_success).

%   P(Y = k) = C(n, k) * p^k * (1-p)^(n-k)
   if k < 0 || k > n
       p = 0;
       return;
   end
  
   p = nchoosek(n, k) * (p_success^k) * ((1 - p_success)^(n - k));
end
