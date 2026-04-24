function reconstructed_value = wenoPol3(uavem, uave, uavep, gp)
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
    
switch gp
    case 0
        stencil_val0 = uave/2 + uavem/2;
        stencil_val1 = (3*uave)/2 - uavep/2;

        LINEAR_WEIGHT0 = 2/3; 
        LINEAR_WEIGHT1 = 1/3;  
    case 1
        stencil_val0 = uave - (5^(1/2)*uave)/10 + (5^(1/2)*uavem)/10;
        stencil_val1 = uave + (5^(1/2)*uave)/10 - (5^(1/2)*uavep)/10;

        LINEAR_WEIGHT0 = 1/2 - 5^(1/2)/30;
        LINEAR_WEIGHT1 = 1/2 + 5^(1/2)/30;
    case 2
        stencil_val0 = uave + (5^(1/2)*uave)/10 - (5^(1/2)*uavem)/10;
        stencil_val1 = uave - (5^(1/2)*uave)/10 + (5^(1/2)*uavep)/10;

        LINEAR_WEIGHT0 = 1/2 + 5^(1/2)/30;
        LINEAR_WEIGHT1 = 1/2 - 5^(1/2)/30;
    case 3
        stencil_val0 = (3*uave)/2 - uavem/2;    % 左边模板
        stencil_val1 = uave/2 + uavep/2;      % 中间模板

        LINEAR_WEIGHT0 = 1/3; 
        LINEAR_WEIGHT1 = 2/3;  
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