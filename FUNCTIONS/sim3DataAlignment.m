%% align data by sm3 method
%   Author:         Xiaochen Qiu from Beihang Univ.
%   Reference:      "Least-Squares Estimation of Transformation Parameters 
%                   between Two Point Patterns"
%   Announcement:   Not debug free. Feel free to do any modifications and
%                   use it at will
%-Inputs:
% @X: 3*n data to be aligned with @Y
% @Y: 3*n data to be aligned with @X
% @mod: sensors setup mode, deciding the alignment method
%-Outputs:
% @R: rotation
% @t: transition
% @s: scale
function [R, t, s] = sim3DataAlignment (X, Y, mod)
% Assertion
if ( size(X,1)~=3 || size(Y,1)~=3 )
    error('TERMINATED BY FUNC Sm3DataAlignment: INPUT IS NOT 3D VECTOR');
end
if ( size(X,2)~=size(Y,2) )
    error('TERMINATED BY FUNC Sm3DataAlignment: INPUTS SIZE UNMATCHED');
end
if nargin==2
    mod = 'mono';
end
if ~strcmp(mod,'mono') && ~strcmp(mod,'stereo') && ~strcmp(mod,'vio')
    error('SENSORS SETUP IS NOT EXISTED, SHOULD BE ''mono'' (DEFAULT), ''stereo'' or ''vio''');
end
% calculate some necesssary quantity
n = size(X,2);
mu_X = sum(X,2)/n;  % mean of X
mu_Y = sum(Y,2)/n;  % mean of Y
sigma2_X = sum(sum(((X-repmat(mu_X,1,size(X,2))).^2),1))/n;  % std of X
% sigma2_Y = sum(sum(((Y-repmat(mu_Y,1,size(Y,2))).^2),1))/n;  % std of Y
covXY = zeros(3,3);
for i = 1:n
    cov = (Y(:,i)-mu_Y) * (X(:,i)-mu_X)';
    covXY = covXY + cov;
end
if strcmp(mod,'vio')
    a = covXY(2,1) - covXY(1,2);
    b = covXY(1,1) + covXY(2,2);
    c = sqrt(a*a+b*b);
    sin_ = a/c;
    cos_ = b/c;
    R = [cos_,-sin_,0; sin_,cos_,0; 0,0,1];
    s = 1;
else
    covXY = covXY/n;    % covariance of X and Y
    % SVD of covXY
    [U,D,V] = svd(covXY);
    % calculate @R, @s and @t
    r = rank(covXY);
    S = eye(3);
    if det(covXY)<0
        S(3,3) = -1;
    end
    R = U*S*V';
    if strcmp(mod,'mono')
        s = trace(D*S)/sigma2_X;
    else
        s = 1;
    end
    if r<2
        fprintf('Solution is not unique\n');
    end
end
t = mu_Y-s*R*mu_X;
end