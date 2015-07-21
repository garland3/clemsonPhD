%% ---- VJC's VERSION OF "padarray" ----
function padded_array = padarray_vjc(array)
% The call padarray_vjc(array) is equivalent to padarray(array,[1,1],'replicate')
padded_array = zeros(size(array)+2);
padded_array(2:end-1,2:end-1) = array;
padded_array(2:end-1,      1) = array(1:end,  1);
padded_array(2:end-1,    end) = array(1:end,end);
padded_array(      1,  1:end) = padded_array(    2,1:end);
padded_array(    end,  1:end) = padded_array(end-1,1:end);
