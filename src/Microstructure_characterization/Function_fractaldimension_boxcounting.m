function [N,box_lengths,fractal_dimension,fractal_dimension_convergence] = Function_fractaldimension_boxcounting(binary_array,p)
% Method from: https://doi.org/10.1016/j.rineng.2020.100106

Domain_size = size(binary_array);
number_dimension = length(Domain_size);

% Resolution
box_minlength=1;
box_maxlength=min(p.boxmaxlength, floor(min(Domain_size(1:number_dimension))/2));
box_length_step = 1;

% Number of box, and length
n_tmp=(box_maxlength-box_minlength)/box_length_step+1;
box_lengths = unique(round(linspace(box_minlength,box_maxlength,n_tmp)));
n_box_length = length(box_lengths);

N = zeros(n_box_length,1);
% Loop over box length
for kbox=1:1:n_box_length
    current_box_length = box_lengths(kbox);
    for kdim = 1:1:number_dimension
        boxbounds.dim(kdim).as =  1:current_box_length:Domain_size(kdim)-(current_box_length-1);
        boxbounds.dim(kdim).bs =  boxbounds.dim(kdim).as+current_box_length-1;
        boxbounds.dim(kdim).n = length(boxbounds.dim(kdim).as);
    end
    if number_dimension==2
        boxbounds.dim(3).as =  1;
        boxbounds.dim(3).bs =  1;
        boxbounds.dim(3).n = 1;
    end
    for kx=1:1:boxbounds.dim(1).n
        xa = boxbounds.dim(1).as(kx); xb = boxbounds.dim(1).bs(kx);
        for ky=1:1:boxbounds.dim(2).n
            ya = boxbounds.dim(2).as(ky); yb = boxbounds.dim(2).bs(ky);
            for kz=1:1:boxbounds.dim(3).n
                za = boxbounds.dim(3).as(kz); zb = boxbounds.dim(3).bs(kz);
                Current_Box = binary_array(xa:xb,ya:yb,za:zb);
                if sum(sum(sum(Current_Box))) > 0
                    N(current_box_length,1) = N(current_box_length,1) + 1;
                end
            end
        end
    end
end

fractal_dimension_convergence = zeros(n_box_length-1,3);
for k=2:1:n_box_length
    x=log(1./box_lengths(1:k));
    y=log(N(1:k));
    [p_,S]=polyfit(x,y,1);
    fractal_dimension_convergence(k-1,1) = k;
    fractal_dimension_convergence(k-1,2) = p_(1);
    fractal_dimension_convergence(k-1,3) = S.normr;
end

fractal_dimension = zeros(4,1);
fractal_dimension(1,1) = 2;
fractal_dimension(2,1) = fractal_dimension_convergence(1,2);
fractal_dimension(3,1) = p.topology_dimension;
fractal_dimension(4,1) = abs(fractal_dimension(3,1) - fractal_dimension(2,1));

if p.plot
    opt.save.savefig = false;
    opt.format.autoclosefig = false;
    p.todo = 1;
    infovol.phasecolor(1,:) = [0    0.4470    0.7410];
    opt.format.linewidth = 2;
    infovol.phasename(1,1) = {'Binary phase'};
    opt.format.fontname = 'Times new roman';
    opt.format.sgtitlefontsize = 16;
    opt.format.titlefontsize = 14;
    opt.format.axefontsize = 12;
    opt.format.legendfontsize = 12;
    function_fractalfig(' ','n/a', 'n/a',box_lengths,N,fractal_dimension_convergence,p.topology_dimension,1,p,opt,infovol)
end

end