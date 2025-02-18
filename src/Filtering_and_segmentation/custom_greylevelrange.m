function [GL] = custom_greylevelrange(GL,Custom_greyscale)
max_=max(max(max(GL)));
min_=min(min(min(GL)));
initial_delta=max_-min_;
final_delta=Custom_greyscale-1;
domain_size=size(GL);
GL=round( ((double((GL-min_)) ./ double(initial_delta)) .* final_delta)+1 );
end