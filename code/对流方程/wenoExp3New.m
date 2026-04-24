function [reconstructed_value, lambda] = wenoExp3New(uavemm, uavem, uave, uavep, uavepp)
    % WENO (Weighted Essentially Non-oscillatory) 重构
    %
    % 输入参数:
    %   v1, v2, v3, v4, v5: 重构所需的五个连续点的值
    %
    % 输出参数:
    %   reconstructed_value: WENO重构值

    lambda = adaptive_lambda(uavemm, uavem, uave, uavep, uavepp);
    h = 1;

    % WENO方法参数
    EPSILON = 1e-6;  % 避免除零的小参数
    
    % 计算每个模板上的重构值（三点模板）
    stencil_val0 = (uave - uave*exp((h*lambda)) - uavem*exp((h*lambda)) + uavem*exp((2*h*lambda)) + h*lambda*uave*exp((2*h*lambda)) - h*lambda*uavem*exp((2*h*lambda)))/(exp(h*lambda) - 1)^2;    % 左边模板
    stencil_val1 = (uavep - uave*exp((h*lambda)) + uave*exp((2*h*lambda)) - uavep*exp((h*lambda)) - h*lambda*uave*exp((h*lambda)) + h*lambda*uavep*exp((h*lambda)))/(exp(h*lambda) - 1)^2;      % 中间模板

    
    % 计算光滑指示器 (beta系数)
    % beta值越小表示相应区域越光滑
    beta0 = (uavem-uave).^2;
    beta1 = (uave-uavep).^2;

    % 线性权重（理想权重）
    LINEAR_WEIGHT0 = (1831*h^6*lambda^6)/6123600 - (1549*h^5*lambda^5)/1020600 - (197*h^4*lambda^4)/68040 + (13*h^3*lambda^3)/810 + (7*h^2*lambda^2)/270 - (2*h*lambda)/9 + 1/3;  % d0 = 1/10
    LINEAR_WEIGHT1 = - (1831*h^6*lambda^6)/6123600 + (1549*h^5*lambda^5)/1020600 + (197*h^4*lambda^4)/68040 - (13*h^3*lambda^3)/810 - (7*h^2*lambda^2)/270 + (2*h*lambda)/9 + 2/3;  % d1 = 6/10

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