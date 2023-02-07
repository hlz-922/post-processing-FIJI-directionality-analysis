# Post-processing-FIJI-directionality-analysis
This is to document the python codes for post processing of saved histogram data by FIJI directionality analysis.

The current postprocessing functions include:
1. Build the 3D direction distribution based on the 2D direction distribution on individual slices
2. Normalise the 3D direction distribution by setting the maximum intensity as 1
3. Linearly scale the 3D direction distribution by setting the maximum intensity as 1 and the minimum intensity as 0
4. Create a pseudo dataset based on the a histogram data
5. Create customised box plots for comparison
