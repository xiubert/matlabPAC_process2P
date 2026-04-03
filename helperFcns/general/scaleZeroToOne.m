function scaledInput = scaleZeroToOne(input)

    scaledInput = input-min(input(:));
    scaledInput = scaledInput / max(input(:));
    
end