function [Otsu_result]=Function_otsu_algorithm(histogram,number_of_class)
% The Otsu's algorithm is well explained in J. Villanova et al. JPS 2013
% The same notation is used in this program, so please refer to the article for a detailed parameter definition

%% CONTINUM HISTOGRAM
x_min=1;
x_max=max(histogram(:,1));
new_histogram=zeros(x_max,2);
new_histogram(:,1)=1:1:x_max;

%unique_histogram=unique(histogram(:,1));
%n_=length(unique_histogram);
[n_,~] = size(histogram);
for k=1:1:n_
    kk=histogram(k,1);
    new_histogram(kk,2)=histogram(k,2);
end
histogram=new_histogram;
   
%% MEAN VALUE AND TOTAL VARIANCE OF THE DISTRIBUTION

I = histogram(:,1); % Axe x of the histogram
ni = histogram(:,2); % Axe y of the histogram
N = sum(ni);
pi = ni/N; % Probalibty of I, Sum(pi)=1
mut = sum(I.*pi); % Mean value
sigmat2 = sum ( ((I-mut).^2).*pi ); % Total variance of the distribution


%% OTSU ALGORITHM

% Number of different value
n = length(histogram(:,1));
% Indice assigned for each value
indice = (1:1:n)';

% Left threshold of the class
k_left=zeros(number_of_class,1);
% Right threshold of the class
k_right=zeros(number_of_class,1);
% Permanent value (always true, whatever is the numbe of class)
k_left(1)=1; k_right(end)=n;

%k_left(1)=min(I); k_right(end)=max(I);
% k_left(i)<=Value<=k_right(i) belongs to the class i

% For each class, there is a loop to choose the threshold value, and the bound of the loop is defined by the precedent loop.
% It's recursive!
% A solution is to compute all the different permutations, and then
% calculate the ratio for each of the permutations
all_permutations = nchoosek( indice(2:end-1) ,number_of_class-1);
% Note that Matlab has already sorted the permutations in both row and
% columns, with an increasing order

% Number of permutation
n_permutation = size(all_permutations,1);

% Add k_left(1) and k_right(end)
all_permutations = [ones(n_permutation,1)*k_left(1) all_permutations ones(n_permutation,1)*k_right(end)];

% Number of possible case
n_case = size(all_permutations,1);

% Store all the sigmab2/sigmat2
ratio_=zeros(n_case,1);

% Iteration on all the different cases
for current_case = 1:1:n_case
    % Current permutation
    current_permutation = all_permutations(current_case,:);
    % Set the threshold
    k_right(1)=current_permutation(2);
    for current_class=2:1:number_of_class
        k_left(current_class)=k_right(current_class-1)+1;
        k_right(current_class)=current_permutation(current_class+1);
    end
    % Check error
    if min(k_right-k_left>=0)==0
        disp 'Error on the Otsu threshold: left>right!'
    end
    % Calculate the sigmab2
    sum_class=0;
    for class_=1:1:number_of_class
        l_= k_left(class_); r_ = k_right(class_); % Threshold
        muc_ = sum( I(l_:r_).*pi(l_:r_) ) / sum(pi(l_:r_));
        diff_mu = (muc_-mut)^2;
        sum_class=sum_class+ (sum(pi(l_:r_))*diff_mu);
    end
    ratio_(current_case,1)=sum_class/sigmat2;
    
end

% Find the case that has maximised the ratio
index = find(ratio_==max(ratio_));

% Save the result for this number of class
Otsu_result.number_of_class = number_of_class;
Otsu_result.sigmab2_sigmat2=ratio_(index(1));
Otsu_result.threshold=all_permutations(index(1),:);
Otsu_result.values=I(all_permutations(index(1),:));

Otsu_result.allratio=ratio_;
Otsu_result.allpermuation=all_permutations;



end


