function [normal] = function_fitplane3dpoints(points,method,display)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if strcmp(method,'Orthogonalregression_PCA')
 
    %http://mathforum.org/library/drmath/view/63765.html
    
    % see https://www.mathworks.com/help/stats/examples/fitting-an-orthogonal-regression-using-principal-components-analysis.html
    [coeff,score,roots] = pca(points); % Principal Component Analysis
    normal = coeff(:,3); % Normalized normal
    
    if display==true
        
        plot3(points(:,1),points(:,2),points(:,3),'bo');
        grid on;
        maxlim = max(abs(points(:)))*1.1;
        axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
        axis square
        view(-9,12);
        
        meanX=mean(points,1);
        
        hold on;
        [xgrid,ygrid] = meshgrid(linspace(min(points(:,1)),max(points(:,1)),5), ...
            linspace(min(points(:,2)),max(points(:,2)),5));
        zgrid = (1/normal(3)) .* (meanX*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
        h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);
        hold off;
        
    end
    
end

% Minimize z distance
% https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
% http://www.ilikebigbits.com/blog/2015/3/2/plane-from-points
% 
% Const = ones(size(all_x)); % Vector of ones for constant term
% Coefficients_plane = double([all_x all_y Const])\double(all_z); % Find the coefficients
% % Equation plane
% % a*x + b+y + c*z = d
% % with c=1
% % z = -a*x - b*y + d
% 
% % a*x + b*y + z = d
% % z=C(1)*x + C(2)*y + C(3)
% 
% a = -Coefficients_plane(1)
% b = -Coefficients_plane(2)
% c = 1
% d = Coefficients_plane(3)