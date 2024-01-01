function quantized = quantization(input, resolution, quantized_max_num)
                quantized = ( double(resolution - 1) / quantized_max_num ) * input;
                round_value = round(quantized);
                quantized = double(quantized_max_num * round_value) / double(resolution - 1);
end