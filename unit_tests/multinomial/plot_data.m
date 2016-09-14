close all
clear all

%% test 2D multinomial CDF

load('mn_2D_cdf.dat') 

input = read_input_parameters('mn_2D_parameters.in'); 
NN = input.N; 
pp = input.p; 
xx = mn_2D_cdf(:,1);
bn_cdf = binocdf(xx,NN,pp); 

figure
plot(mn_2D_cdf(:,1),mn_2D_cdf(:,2),'reds-') 
hold on
plot(xx,bn_cdf,'greens-') 
legend('2D multinomial (C++)','binomial (matlab)') 

%% test 3D multinomial

% read parameters
input = read_input_parameters('mn_3D_parameters.in'); 
NN = input.N; 
pp = [input.p0 input.p1 input.p2]; 

% load samples 
load('mn_3D_samples.dat') 
x1x2 = mn_3D_samples(:,1:2);
number_samples = size(x1x2,1);

% count samples on a 2D mesh
edges1 = -0.5:(NN+0.5);
edges2 = -0.5:(NN+0.5);
[count, ~, mid, ~] = histcn(x1x2, edges1, edges2);
probability_empirical = count/number_samples;
x1 = mid{1};
x2 = mid{2};

% plot histogram of samples
figure
bar3(probability_empirical)
zlabel('Probability Mass')
title('Trinomial Distribution from 3D multinomial (C++)')

% create matlab's 3D multinomial pdf
probability_exact = zeros(length(x1),length(x2));
for ii = 1:length(x1) 
    for jj = 1:length(x2) 
        xx = [x1(ii), x2(jj), NN - (x1(ii)+x2(jj))]; 
        probability_exact(ii,jj) = mnpdf(xx,pp);
    end
end
        
% plot bar graph of matlab's pdf
figure
bar3(probability_exact)
zlabel('Probability Mass')
title('Trinomial Distribution from mnpdf (matlab)')

% 3D multinomial error 
error = abs(probability_empirical - probability_exact)./probability_exact;
corrected_error = error; 
prob_cutoff = 1e-2*max(probability_exact(:));
corrected_error(probability_exact < prob_cutoff) = 0;

% plot error
figure
set(gca,'fontsize',20)
bar3(corrected_error) 
title('corrected % error; 3D multinomial')


%% test 4D multinomial

% read parameters
input = read_input_parameters('mn_4D_parameters.in'); 
NN = input.N; 
pp = [input.p0 input.p1 input.p2 input.p3]; 

% load samples 
load('mn_4D_samples.dat') 
x1x2x3 = mn_4D_samples(:,1:3);
number_samples = size(x1x2x3,1);

% compare predicted with empirical means
disp('compare predicted with empirical moments in 4D multinomial')
mean_exact = pp*NN;
disp(['exact means = ' num2str(mean_exact)]) 
mean_empirical = mean(mn_4D_samples);
disp(['empirical means = ' num2str(mean_empirical)]) 

% compare predicted with empirical variances
variance_exact = pp.*(1-pp)*NN;
disp(['exact variances = ' num2str(variance_exact)]) 
variance_empirical = var(mn_4D_samples);
disp(['empirical variances = ' num2str(variance_empirical)]) 

% count samples on a 3D mesh
edges1 = -0.5:(NN+0.5);
edges2 = -0.5:(NN+0.5);
edges3 = -0.5:(NN+0.5);
[count, ~, mid, ~] = histcn(x1x2x3, edges1, edges2, edges3);
probability_empirical = count/number_samples;
x1 = mid{1};
x2 = mid{2};
x3 = mid{3};

% fixed slice in dim 3
fixed_slice_3 = int64(mean_exact(3)); 

% plot histogram of samples with fixed dimension
figure
set(gca,'fontsize',20)
bar3(probability_empirical(:,:,fixed_slice_3))
zlabel('Probability Mass')
title('4D multinomial (C++); one dimension fixed')

% create matlab's 4D multinomial pdf on a 3D mesh
probability_exact = zeros(length(x1),length(x2),length(x3));
for ii = 1:length(x1)
    for jj = 1:length(x2)
        for kk = 1:length(x3)
            xx = [x1(ii), x2(jj), x3(kk), NN - (x1(ii)+x2(jj)+x3(kk))];
            probability_exact(ii,jj,kk) = mnpdf(xx,pp);
        end
    end
end
    
% plot bar graph of matlab's pdf with fixed dimension
figure
set(gca,'fontsize',20)
bar3(probability_exact(:,:,fixed_slice_3))
zlabel('Probability Mass')
title('4D multinomial Distribution from mnpdf (matlab); one dimension fixed')

% calculate 4D multinomial error 
error = abs(probability_empirical - probability_exact)./probability_exact;
corrected_error = error; 
prob_cutoff = 1e-2*max(probability_exact(:));
corrected_error(probability_exact < prob_cutoff) = 0;

disp(['4D multinomial max corrected % error = ' num2str(max(corrected_error(:)))]) 

