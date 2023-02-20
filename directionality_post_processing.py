import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

'''
This module allows user to do post-processings based on the results of Directionality analysis in FIJI.

Functions include read_sample_roi, organise, select_plot, find_max_shift, pseudo_dataset, create_pseudo_boxplot. Further details can be found in the documentation of the specific function.

'''

# 1. Extract all samples and all rois (produced using the same conidtion) to a list of pandas dataframe
# list layer 1: sample (no.i)
# list layer 2: roi (no.j)
# list[i][j] stores the pandas dataframe of scaffold no.i, roi no.j

def read_sample_roi(sample_code:str, sample_no:int, roi_no:int):
    
    '''
    This function allows the user to extract directionality data from many sample folders that individually contain roi data
    
    ---input---
    sample_code: the code for the sample
    sample_no: the number of samples to be read
    roi_no: the number of rois to be read
    
    ---output---
    A list of dataframes
    List layer 1: sample no (index starts from 0)
    List layer 2: roi no (index starts from 0)

    '''
   
    scaffold_list = [0]*(sample_no)

    for i in range(sample_no):
        roi_list = [0]*(roi_no)
        for j in range(roi_no):
            roi_list[j] = pd.read_csv(f'{sample_code}_{i+1}/r{j+1}.csv',encoding= 'unicode_escape')
        scaffold_list[i] = roi_list
    
    return scaffold_list





# 2. Organise a dataframe
# Raw dataframe contains the individual 2D distribution histogram results for every slice in the 3D distribution
# For each slice: (1) Raw histogram (2) Fitted Gaussian distribution (discarded for future calculations)

def organise(df):
    
    '''
    This function smash 2D distributions of multiple slices to 3D distribution
    
    ---input
    df: the raw dataframe
    
    ---ouput
    Add 6 new columns in the previous dataframe
    'total': add the frequency of all individual slices together 
    'total_fit': add the frequency of all individual slices (Gaussian fit) together 
    'total_norm': normalise column 'total' by setting the maximum count in the column to 1
    'total_fit_norm': normalise column 'total_fit' by setting the maximum count in the column to 1
    'total_frequency': the frequency values for column 'total'
    'total_linearscale': linearly scale column 'total' so that the maximum count in the column is 1 and the minimum is 0
    
    '''
    
    column_no = df.shape[1]
    
    # initialise the columns 'total' and 'total_fit'
    # 'total': the sum of intensity for all 2D slices at specific angle (bin midpoint)
    # 'total_fit': the sum of intensity for all 2D slices based on the Gaussian fit
    df['total']=np.zeros(df.shape[0])
    df['total_fit']=np.zeros(df.shape[0])
    
    # sum by looping through every two columns (default save by FIJI)
    for i in np.arange(1,column_no,2):
        df['total'] = df['total']+df.iloc[:,i]
    for i in np.arange(2,column_no,2):
        df['total_fit'] = df['total_fit']+df.iloc[:,i]
        
    # normalise column 'total' and 'total_fit' by setting the maximum being 1
    df['total_norm'] = df['total']/max(df['total'])
    df['total_fit_norm'] = df['total_fit']/max(df['total_fit'])
    
    # normalise column 'total' by setting the overal sum is 1
    df['total_frequency'] = df['total']/sum(df['total'])
    
    # linearly scale column 'total' by setting the maximum being 1 and the minimum being 0
    a = 1 / (max(df['total']) - min(df['total']))
    b = - a * min(df['total'])
    df['total_linearscale'] = df['total'].apply(lambda x: a*x+b)
    
    return df
    
    
    


# 3. Choose which scaffold to plot and which rois to plot
# The default setting is to plot all 5 rois but can be specified
# The default colour combination for the histogram plot is set but can be customised if needed

def select_plot(df_list:list, plot_kind:str, width=1.1, roi_list=[1,2,3,4,5], 
                colour_list=['lightcoral', 'gray', 'cornflowerblue', 'gray', 'lightcoral'], 
                rotate=90):
    
    '''
    This function allows the user to choose certain roi(s) to plot for one sample

    ---input
    df_list: the list of dataframes, which includes all rois from one sample
    plot_kind: 'norm' or 'linearscale', more details can be found in dpp.clean function
    *width: the bin width | default: 1.1, which performs well for the original dataset with a bin resolution of 1/deg
    *roi_list: the roi(s) that will be shown on the plot (index starts from 1) | default: [1,2,3,4,5]
    *colour_list: the alternating colour for different roi(s) | default: ['lightcoral', 'gray', 'cornflowerblue', 'gray', 'lightcoral']
    *rotate: rotatation degree of the reference of frame (anti-clockwise) | default: 90
    
    ---output
    A bar plot disgused to be a histogram plot
    
    '''
    
    label_list = ['B1', 'B2', 'B3', 'B4', 'B5']
    for i in roi_list:
        # Extract direction (bin centre) from the roi (index 0 column)
        x = df_list[i-1].iloc[:,0]
        # Plot histogram plot by using bar plot
        # Convert the direction from 0~180/deg to -90~90/deg
        plt.bar(np.multiply(np.subtract(x,rotate),-1),df_list[i-1][f'total_{plot_kind}'],label=label_list[i-1],
        width=width,color=colour_list[i-1],linewidth=0)
    
    plt.xlabel('Pore orientation (°)')
    plt.ylabel('Normalised intensity (%)')

    plt.xlim([-90,90])
    # plt.ylim([20,110])
    plt.yticks([0,1])
    plt.legend(loc='upper left')
    plt.show()

    
    
    

# 4. Find the maximum orientation change using the normalised or the median
# The default method is to use the maximum peak in the histogram
# The optional method is to use the pseudo dataset created from the histogram to find the median

def find_max_shift(df_list:list, method='peak', size=5000):
    '''
    This function returns the maximum directionality change in the datasets (dataframes in list format)

    ---input
    df_list: the list of dataframes, which includes all rois from one sample
    *method: the method used to calculate the maximum directionalty change. (Options: 'peak' and 'median') | default: 'peak'
    ->The 'peak' method is to use the peak positions in the histogram plot. 
    ->The 'median' method is to use the median in the distribution based on the pseudo dataset.
    *size: the size of the pseudo dataset | default: 5000 (more details in pseudo_dataset function documentation)
    
    ---output
    Returns a number that is the maximum directionality change
    
    '''
    index_list = []
    
    if method == 'peak':
        for df in df_list:
            index_list.append(df['total_norm'].idxmax())
        return max(index_list) - min(index_list)
    
    elif method == 'median':
        for df in df_list:
            dataset, _ = pseudo_dataset(df, size)
            index_list.append(np.median(dataset))
        return max(index_list) - min(index_list)
    
    else:
        print('Wrong useage of the function. :)')
    




# 5. Create pseudo dataset from histogram data (for one dataframe)
# The resolution for the pseudodataset is also returned

def pseudo_dataset(df, size=5000):
    
    '''
    This function creates a pseudo dataset from the recorded histogram data from FIJI directionality analysis.
    
    ---input--- 
    df: the dataframe containing histogram data from FIJI directionality analysis
    *size: the size of the pseudo dataset | default: 5000, which produces an identical histogram to the original directionality analysis with a bin resolution of 1/deg. The ideal size is recommended to tune based on specific datasets.
    
    ---output---
    pseudo_dataset: the created pseudo dataset
    res: the bin resolution for the original dataset for reference
    
    '''
    
    pseudo_dataset = []
    for i in range(df.shape[0]):
        bin_centre = df['Direction (°)'][i]
        N = int(df['total_frequency'][i] * size)
        # Convert the bin centre from 0~180/deg to -90~90/deg
        pseudo_dataset.extend( [(90 - bin_centre)] * N )
    # The 'resolution' of the pseudodataset is defined as the the distance between bin centres
    res = df['Direction (°)'][1] - df['Direction (°)'][0]
    return pseudo_dataset, res




        
# 6. Create pseudo dataset from histogram data for boxplots
# The default colour_list is set but can be customised if needed
# The default size for the dataset was scaled to 5000. The ideal size is recommended to tune based on specific dataset. With a resolution of 

def create_pseudo_boxplot(df_list:list, size=5000, colour_list=['pink', 'lightgray', 'lightblue', 'lightgray', 'pink'],
                         labels = ['Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5']):
    
    '''
    This function create boxplots based on pseudo datasets (more details in pseudo_dataset function documentation
    
    ---input---
    df_list: the list of dataframes, which includes all rois from one sample
    *size: the size of the pseudo dataset | default: 5000 (more details in pseudo_dataset function documentation)
    *colour_list: the alternating colour for different roi(s) | default: ['pink', 'lightgray', 'lightblue', 'lightgray', 'pink']
    *labels: the labels for the roi(s) | default: ['Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5']
    
    ---ouput---
    The boxplots for different rois
    
    '''
    
    pseudo_data_list = []
    for df in df_list:
        pseudo_ds, _ = pseudo_dataset(df, size)
        pseudo_data_list.append(pseudo_ds)
    
    # styling the median and mean
    medianprops = dict(linestyle='-.', linewidth=1.5, color='black')
    meanpointprops = dict(marker='D', markeredgecolor='black', markersize=10, markerfacecolor='firebrick')
    
    bplot = plt.boxplot(pseudo_data_list, patch_artist=True, labels=labels, showfliers=False, showmeans=True, medianprops=medianprops, meanline=False, meanprops=meanpointprops)
        
    # Fill the box with colours
    for patch, color in zip(bplot['boxes'], colour_list):
        patch.set_facecolor(color)
 
    plt.ylim([-100,100])
    plt.ylabel('Pore orientation (°)')
    plt.show()
    

