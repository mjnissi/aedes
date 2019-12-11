function w_fid=linebroad(fid,at,lb)
%line braodening for csi data set

dims=size(fid);
[np,loc]=max(dims);

if length(dims)==2
  dims(3)=1;
end

sw = np/at;
t1 = ((1:np) .* 1/sw);

if (loc == 1)   
    for idx = 1:dims(2)
        for idx2 = 1:dims(3)
            fid(:,idx,idx2) = squeeze(fid(:,idx,idx2)) .* exp(-t1 .* pi * lb)';
        end    
    end
end 

if (loc == 2)
    for idx = 1:dims(1)
        for idx2 = 1:dims(3)
            fid(:,idx,idx2) = fid(:,idx,idx2) .* exp(-t1 .* pi * lb)';
        end    
    end
end 

if (loc == 3)
    for idx = 1:dims(1)
        for idx2 = 1:dims(2)
            fid(:,idx,idx2) = fid(:,idx,idx2) .* exp(-t1 .* pi * lb)';
        end    
    end
end 

w_fid=fid;
