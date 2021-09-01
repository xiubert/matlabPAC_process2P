function RGBout = redGreenMerge(redImg,greenImg)
    if size(redImg)~=size(greenImg)
        error('Image sizes must be the same')
    end
    RGBout = cat(3,scaleZeroToOne(redImg),scaleZeroToOne(greenImg),nan(size(redImg)));