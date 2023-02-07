import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 1. Extract all scaffolds and all rois (produced using the same conidtion) to a list of pandas dataframe
# list layer 1: scaffold (no.i)
# list layer 2: roi (no.j)
# list[i][j] stores the panda dataframe of scaffold no.i, roi no.j

def read_scaffold_roi(scaffold_code, scaffold_no, roi_no):
   
    scaffold_list = [0]*(scaffold_no)

    for i in range(scaffold_no):
        roi_list = [0]*(roi_no)
        for j in range(roi_no):
            roi_list[j] = pd.read_csv(f'{scaffold_code}_{i+1}/r{j+1}.csv',encoding= 'unicode_escape')
        scaffold_list[i] = roi_list
    
    return scaffold_list



# 2. Organise a dataframe
# Raw dataframe contains the individual 2D distribution histogram results for every slice in the 3D distribution
# For each slice: (1) Raw histogram (2) Fitted Gaussian distribution (discarded for future calculations)
def total_normal_linearscale(df):
    
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
    df['total_norm2'] = df['total']/sum(df['total'])
    
    # linearly scale column 'total' by setting the maximum being 1 and the minimum being 0
    a = 1 / (max(df['total']) - min(df['total']))
    b = - a * min(df['total'])
    df['total_linearscale'] = df['total'].apply(lambda x: a*x+b)
    
    return df
    


# Choose which scaffold to plot and which rois to plot
# The default colour_list is set but can be customised if needed
def select_plot(df_list:list, number_list:list,plot_kind:str, colour_list=['lightcoral', 'gray', 'cornflowerblue', 'gray', 'lightcoral']):
    
    label_list = ['B1', 'B2', 'B3', 'B4', 'B5']
    for i in number_list:
        plt.bar(np.multiply(np.subtract(df_list[i-1].iloc[:,0],90),-1),df_list[i-1][f'total_{plot_kind}'],label=label_list[i-1],
        width=1.1,color=colour_list[i-1],linewidth=0)
    
    plt.xlabel('Pore orientation (°)')
    plt.ylabel('Normalised intensity (%)')

    plt.xlim([-90,90])
    # plt.ylim([20,110])
    plt.yticks([0,1])
    plt.legend(loc='upper left')
    plt.show()


# Find the maximum orientation change using the normalised
def find_max_shift(df_list):
    index_list = []
    for df in df_list:
        index_list.append(df['total_norm'].idxmax())
    return max(index_list) - min(index_list)


# Create pseudo dataset from histogram data for boxplots
# The default colour_list is set but can be customised if needed
def create_pseudo_boxplot(df_list:list, colour_list=['pink', 'lightgray', 'lightblue', 'lightgray', 'pink']):
    pseudo_data_list = []
    for df in df_list:
        pseudo_dataset = []
        for i in range(df.shape[0]):
            bin_centre = df['Direction (°)'][i]
            N = int(df['total_norm2'][i]*1000)
            pseudo_dataset.extend( [bin_centre] * N )
        pseudo_data_list.append(pseudo_dataset)
    bplot = plt.boxplot(pseudo_data_list, patch_artist=True, labels=['Base 1', 'Base 2', 'Base 3', 'Base 4', 'Base 5'])
    
    # fill with colours
    for patch, color in zip(bplot['boxes'], colour_list):
        patch.set_facecolor(color)
 
    plt.ylim([-10,200])
    plt.ylabel('Pore orientation (°)')
    plt.show()
    

