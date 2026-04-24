function numFlux = lfFlux(flux, alpha, uL, uR)
    % 计算 Lax-Friedrichs 数值通量
    % flux  : 物理通量函数句柄
    % alpha : 最大特征速度 (标量)
    % uL    : 左侧重构值
    % uR    : 右侧重构值
    % numFlux : 输出的数值通量

    numFlux(:) = 0.5 * (flux(uL) + flux(uR)) - 0.5 * alpha * (uR - uL);

end