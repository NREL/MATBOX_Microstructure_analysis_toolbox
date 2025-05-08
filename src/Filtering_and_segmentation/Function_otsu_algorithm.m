function [Otsu_result]=Function_otsu_algorithm(histogram,number_of_class)
% The Otsu's algorithm is well explained in J. Villanova et al. JPS 2013
% The same notation is used in this program, so please refer to the article for a detailed parameter definition

noninteger = false;
if sum(histogram(:,1)==round(histogram(:,1)))~=length(histogram(:,1))
    histogram(:,1) = round(histogram(:,1)*1000);
    noninteger = true;
end

%% CONTINUM HISTOGRAM
x_min=min(histogram(:,1));
x_max=max(histogram(:,1));
N = x_max-x_min+1;
new_histogram=zeros(N,2);
new_histogram(:,1)=x_min:1:x_max;

%unique_histogram=unique(histogram(:,1));
%n_=length(unique_histogram);
[n_,~] = size(histogram);
a = (N-1)/(x_max-x_min);
b = N-a*x_max;
for k=1:1:n_
    %kk = round(histogram(k,1));
    x = round(histogram(k,1));
    kk = round(a*x + b);
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
sigmaw2 = zeros(n_case,1); % Within variance per class
thresholds_bounds = zeros(number_of_class,2);

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
        sigmaw2(current_case,class_) = sum( pi(l_:r_) .* (I(l_:r_)-muc_).^2);
    end
    ratio_(current_case,1)=sum_class/sigmat2;  
end

% Find the case that has maximised the ratio
index = find(ratio_==max(ratio_));

all_permutations = all_permutations+x_min-1;

% Save the result for this number of class
Otsu_result.number_of_class = number_of_class;
Otsu_result.sigmab2_sigmat2=ratio_(index(1));
Otsu_result.sigmaw2=sigmaw2(index(1),:);
Otsu_result.threshold=all_permutations(index(1),:);


thresholds_bounds(1,1) = -1;
thresholds_bounds(1,2) = Otsu_result.threshold(2);
for k=2:number_of_class
    thresholds_bounds(k,1) = thresholds_bounds(k-1,2)+1;
    thresholds_bounds(k,2) = Otsu_result.threshold(k+1);
end

Otsu_result.thresholdbounds=thresholds_bounds;

%Otsu_result.values=I(all_permutations(index(1),:));

Otsu_result.allratio=ratio_;
Otsu_result.allpermuation=all_permutations;

if noninteger
    Otsu_result.threshold = Otsu_result.threshold/1000;
    Otsu_result.allpermuation = Otsu_result.allpermuation/1000;
end

end


