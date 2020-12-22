
classdef Synapse < handle
    
    properties
        angle = 0; % double, the rotation angle
        areaBassoon; % pi * std_x * std_y   -> std = standard deviation
        areaBassoonDescription = 'pi * stdev(xLocs) * stdev(yLocs)';
        areaHomer; % pi * std_x * std_y     -> std = standard deviation
        areaHomerDescription = 'pi * stdev(xLocs) * stdev(yLocs)';
        areaProteins; % pi * std_x * std_y  -> std = standard deviation
        areaProteinsDescription = 'pi * stdev(xLocs) * stdev(yLocs)';
        areaUnits; % string
        bassoon; % 2D array, localizations
        bassoonRotated; % 2D array, rotated localizations
        bassoonRotatedResult; % the gaussian fit result to the histogram
        bassoonRotatedResultDescription = '1D Gaussian fit [bg; sigma; area; mean ...] to the rotated bassoon histogram';
        bassoonRotatedX; % the x values used to fit the histogram to a Gaussian
        bassoonRotatedXDescription = 'The x centers used to fit the rotated bassoon histogram to a 1D Gaussian function';
        bassoonRotatedXFit; % the x values from the fit to the histogram
        bassoonRotatedXFitDescription = 'The x centers from the line of best fit of the rotated bassoon histogram'
        bassoonRotatedY; % the y values used to fit the histogram to a Gaussian
        bassoonRotatedYDescription = 'The y counts used to fit the rotated bassoon histogram to a 1D Gaussian function';
        bassoonRotatedYFit; % the y values from the fit to the histogram
        bassoonRotatedYFitDescription = 'The y counts from the line of best fit of the rotated bassoon histogram'
        distances; % the distances between bassoon and homer and each protein and homer
        homer; % 2D array, localizations
        homerRotated; % 2D array, rotated localizations
        homerRotatedResult; % the gaussian fit result to the histogram
        homerRotatedResultDescription = '1D Gaussian fit [bg; sigma; area; mean ...] to the rotated homer histogram';
        homerRotatedX; % the x values used to fit the histogram to a Gaussian
        homerRotatedXDescription = 'The x centers used to fit the rotated homer histogram to a 1D Gaussian function';
        homerRotatedXFit; % the x values from the fit to the histogram
        homerRotatedXFitDescription = 'The x centers from the line of best fit of the rotated homer histogram'
        homerRotatedY; % the y values used to fit the histogram to a Gaussian
        homerRotatedYDescription = 'The y counts used to fit the rotated homer histogram to a 1D Gaussian function';
        homerRotatedYFit; % the y values from the fit to the histogram
        homerRotatedYFitDescription = 'The y counts from the line of best fit of the rotated homer histogram'
        origin; % double, the origin for all rotations
        originDescription; % string, text describing how the origin is calculated
        proteins; % cell array, of 2D array of localizations
        proteinsRotated; % cell array, of 2D array of rotated localizations        
        proteinsRotatedResult; % the gaussian fit result to the histogram
        proteinsRotatedResultDescription = '1D Gaussian fit [bg; sigma; area; mean ...] to each rotated proteins histogram';
        proteinsRotatedX; % the x values used to fit the histogram to a Gaussian
        proteinsRotatedXDescription = 'The x centers used to fit each rotated proteins histogram to a 1D Gaussian function';
        proteinsRotatedXFit; % the x values from the fit to the histogram
        proteinsRotatedXFitDescription = 'The x centers from the line of best fit of each rotated proteins histogram'
        proteinsRotatedY; % the y values used to fit the histogram to a Gaussian
        proteinsRotatedYDescription = 'The y counts used to fit each rotated proteins histogram to a 1D Gaussian function';
        proteinsRotatedYFit; % the y values from the fit to the histogram
        proteinsRotatedYFitDescription = 'The y counts from the line of best fit of each rotated proteins histogram'
        sigmaBassoon; % = max(std_x, std_y)  -> std = standard deviation
        sigmaBassoonDescription; % string
        sigmaHomer; % = max(std_x, std_y)    -> std = standard deviation
        sigmaHomerDescription; % string 
        sigmaXProteins; % the standard deviation in the X dimension
        sigmaXProteinsDescription = 'The standard deviation of the xLocs for each protein, in nm';
        sigmaYProteins; % the standard deviation in the Y dimension
        sigmaYProteinsDescription = 'The standard deviation of the yLocs for each protein, in nm';
        widthBassoon; % Gaussian fit full width (i.e., 2 * sqrt(2*log(2)) * sigma * nm_per_pixel), in nm
        widthBassoonDescription = 'FWHM of the Gaussian fit -> 2 * sqrt(2*log(2)) * sigma_value_of_the_Gaussian_with_the_largest_peak * nm_per_pixel)';
        widthHomer; % Gaussian fit full width (i.e., 2 * sqrt(2*log(2)) * sigma * nm_per_pixel), in nm
        widthHomerDescription = 'FWHM of the Gaussian fit -> 2 * sqrt(2*log(2)) * sigma_value_of_the_Gaussian_with_the_largest_peak * nm_per_pixel)';
    end
            
end