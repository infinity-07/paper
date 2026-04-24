function reconstructed_value = wenoExp3(uavem, uave, uavep, gp, lambda)
    % WENO (Weighted Essentially Non-oscillatory) 重构
    %
    % 输入参数:
    %   v1, v2, v3, v4, v5: 重构所需的五个连续点的值
    %
    % 输出参数:
    %   reconstructed_value: WENO重构值

    
    % WENO方法参数
    EPSILON = 1e-6;  % 避免除零的小参数
    
    % 计算光滑指示器 (beta系数)
    % beta值越小表示相应区域越光滑
    beta0 = (uavem-uave).^2;
    beta1 = (uave-uavep).^2;

    h=1;
    
switch gp
    case 0
        stencil_val0 = (h*lambda*exp(h*lambda)*(uave - uavem))/(exp(h*lambda) - 1)^2 - (uave - uavem*exp(h*lambda))/(exp(h*lambda) - 1);
        stencil_val1 = (uavep - uave*exp(h*lambda) + uave*exp(2*h*lambda) - uavep*exp(h*lambda) - h*lambda*uave + h*lambda*uavep)/(exp(h*lambda) - 1)^2;

        LINEAR_WEIGHT0 =  -(exp(-h*lambda)*(h*lambda - exp(2*h*lambda) + h*lambda*exp(2*h*lambda) + 1))/(2*sinh(h*lambda)*(h*lambda - exp(h*lambda) + 1)); 
        LINEAR_WEIGHT1 = (exp(h*lambda)*(2*h*lambda*exp(h*lambda) - exp(2*h*lambda) + 1))/((exp(h*lambda) - 1)*(exp(h*lambda) + 1)*(h*lambda - exp(h*lambda) + 1));  
    case 1
        stencil_val0 = (h*lambda*exp(-(h*lambda*(5^(1/2) - 15))/10)*(uave - uavem))/(exp(h*lambda) - 1)^2 - (uave - uavem*exp(h*lambda))/(exp(h*lambda) - 1);
        stencil_val1 = -(uave*exp((h*lambda*(5^(1/2) + 10))/10) - uave*exp((h*lambda*(5^(1/2) + 20))/10) + uavep*exp((h*lambda*(5^(1/2) + 10))/10) - uavep*exp((5^(1/2)*h*lambda)/10) + h*lambda*uave*exp((h*lambda)/2) - h*lambda*uavep*exp((h*lambda)/2))/(exp((h*lambda*(5^(1/2) + 20))/10) - 2*exp((h*lambda*(5^(1/2) + 10))/10) + exp((5^(1/2)*h*lambda)/10));

        LINEAR_WEIGHT0 = -(exp((h*lambda*(5^(1/2) - 15))/10)*((exp(lambda*(h + xa)) + h*lambda*exp((lambda*(15*h + 10*xa - 5^(1/2)*h))/10))/(2*exp(lambda*(h + xa)) - exp(lambda*xa) - 2*exp(lambda*(3*h + xa)) + exp(lambda*(4*h + xa))) - (exp(-(lambda*(5*h + 20*xa - 5^(1/2)*h))/10)*(exp((lambda*(35*h + 20*xa - 5^(1/2)*h))/10) - h*lambda*exp(lambda*(3*h + 2*xa))))/((exp(h*lambda) - 1)^3*(exp(h*lambda) + 1)))*(exp(h*lambda) - 1)^2)/(exp((h*lambda*(5^(1/2) - 5))/10) - exp((h*lambda*(5^(1/2) + 5))/10) + h*lambda);
        LINEAR_WEIGHT1 = ((exp(lambda*(h + xa))/(2*exp(lambda*(h + xa)) - exp(lambda*xa) - 2*exp(lambda*(3*h + xa)) + exp(lambda*(4*h + xa))) + (exp(-(lambda*(5*h + 20*xa - 5^(1/2)*h))/10)*(h*lambda*exp((lambda*(15*h + 10*xa - 5^(1/2)*h))/5) - exp((lambda*(35*h + 20*xa - 5^(1/2)*h))/10) + h*lambda*exp(2*lambda*(h + xa))))/((exp(h*lambda) - 1)^3*(exp(h*lambda) + 1)))*(exp((h*lambda*(5^(1/2) + 20))/10) - 2*exp((h*lambda*(5^(1/2) + 10))/10) + exp((5^(1/2)*h*lambda)/10)))/(exp((5^(1/2)*h*lambda)/10) - exp((h*lambda*(5^(1/2) + 10))/10) + h*lambda*exp((h*lambda)/2));
    case 2
        stencil_val0 = (h*lambda*exp((h*lambda*(5^(1/2) + 15))/10)*(uave - uavem))/(exp(h*lambda) - 1)^2 - (uave - uavem*exp(h*lambda))/(exp(h*lambda) - 1);
        stencil_val1 = (uavep - uave*exp(h*lambda) + uave*exp(2*h*lambda) - uavep*exp(h*lambda) - h*lambda*uave*exp((h*lambda*(5^(1/2) + 5))/10) + h*lambda*uavep*exp((h*lambda*(5^(1/2) + 5))/10))/(exp(h*lambda) - 1)^2;

        LINEAR_WEIGHT0 = -(((exp(h*lambda) + h*lambda*exp((h*lambda*(5^(1/2) + 15))/10))/((exp(h*lambda) - 1)^3*(exp(h*lambda) + 1)) - (exp(-(lambda*(5*h + 20*xa + 5^(1/2)*h))/10)*(exp((lambda*(35*h + 20*xa + 5^(1/2)*h))/10) - h*lambda*exp(lambda*(3*h + 2*xa))))/((exp(h*lambda) - 1)^3*(exp(h*lambda) + 1)))*(exp(h*lambda) - 1)^2)/(exp(h*lambda) - exp(2*h*lambda) + h*lambda*exp((h*lambda*(5^(1/2) + 15))/10));
        LINEAR_WEIGHT1 = (exp(-(lambda*(5*h + 20*xa + 5^(1/2)*h))/10)*(exp((lambda*(15*h + 20*xa + 5^(1/2)*h))/10) - exp((lambda*(35*h + 20*xa + 5^(1/2)*h))/10) + h*lambda*exp((lambda*(15*h + 10*xa + 5^(1/2)*h))/5) + h*lambda*exp(2*lambda*(h + xa))))/((exp(h*lambda) - 1)*(exp(h*lambda) + 1)*(h*lambda*exp((h*lambda*(5^(1/2) + 5))/10) - exp(h*lambda) + 1));
    case 3
        stencil_val0 = (h*lambda*exp(2*h*lambda)*(uave - uavem))/(exp(h*lambda) - 1)^2 - (uave - uavem*exp(h*lambda))/(exp(h*lambda) - 1);    % 左边模板
        stencil_val1 = (uavep - uave*exp(h*lambda) + uave*exp(2*h*lambda) - uavep*exp(h*lambda) - h*lambda*uave*exp(h*lambda) + h*lambda*uavep*exp(h*lambda))/(exp(h*lambda) - 1)^2;      % 中间模板

        LINEAR_WEIGHT0 = -(exp(-h*lambda)*(2*h*lambda*exp(h*lambda) - exp(2*h*lambda) + 1))/(2*sinh(h*lambda)*(h*lambda*exp(h*lambda) - exp(h*lambda) + 1)); 
        LINEAR_WEIGHT1 = (exp(h*lambda)*(h*lambda - exp(2*h*lambda) + h*lambda*exp(2*h*lambda) + 1))/((exp(h*lambda) - 1)*(exp(h*lambda) + 1)*(h*lambda*exp(h*lambda) - exp(h*lambda) + 1));  
end

    % 计算非线性权重
    alpha0 = LINEAR_WEIGHT0./((EPSILON + beta0).^2);
    alpha1 = LINEAR_WEIGHT1./((EPSILON + beta1).^2);
   
    % 权重归一化
    weight_sum = alpha0 + alpha1;
    weight0 = alpha0./weight_sum;
    weight1 = alpha1./weight_sum;

    % 计算最终的重构值
    reconstructed_value = weight0.*stencil_val0 + ...
                         weight1.*stencil_val1;
end