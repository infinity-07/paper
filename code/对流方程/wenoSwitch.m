function [value, lambda] = wenoSwitch(v1, v2, v3, v4, v5)

global config 
lambda = 0;
switch config
    case 1
    value = weno(v1, v2, v3, v4, v5);

    case 2
    value = wenoPol3(v2, v3, v4);

    case 3
    value = wenoExp3(v2, v3, v4);

    case 5
    [value, lambda] = wenoExp3New(v1, v2, v3, v4, v5);

    case 6
    [value, lambda] = wenoExpThinc3(v1, v2, v3, v4, v5);
end

end