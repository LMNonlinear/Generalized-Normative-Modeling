function [isallow,freemem]=test_gpu(bitSize,count)
if nargin<1 || isempty(bitSize)
    bitSize=4e9;
end
if nargin<2
    count=1;
end


if gpuDeviceCount>0    
    try
        gpuInfo=gpuDevice(count);
    catch
        fprintf('gpuDevice() has error, skipped.\n');
        gpuInfo=gpuDevice(1);
    end    
    freemem=gpuInfo.AvailableMemory;
    if freemem>bitSize
        isallow=true;
    else
        isallow=false;
    end
else
    isallow=false;
end

end
